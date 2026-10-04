// ============================================================================
//  dsr_io.hpp  -  Small POSIX helpers: pipes to child processes, TCP, timing
// ============================================================================
#pragma once

#include <arpa/inet.h>
#include <fcntl.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

#include <poll.h>
#include <atomic>
#include <condition_variable>
#include <deque>
#include <memory>
#include <mutex>
#include <thread>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace dsr::io {

// Write the complete buffer (handles partial writes and EINTR).
inline bool write_all(int fd, const void* buf, size_t n) {
  const char* p = static_cast<const char*>(buf);
  while (n) {
    ssize_t w = ::write(fd, p, n);
    if (w < 0) { if (errno == EINTR) continue; return false; }
    p += w; n -= static_cast<size_t>(w);
  }
  return true;
}

// Blocking read of at most n bytes; returns 0 on EOF, <0 on error.
inline ssize_t read_some(int fd, void* buf, size_t n) {
  for (;;) {
    ssize_t r = ::read(fd, buf, n);
    if (r < 0 && errno == EINTR) continue;
    return r;
  }
}

// Read exactly n bytes unless EOF; returns bytes read.
inline size_t read_full(int fd, void* buf, size_t n) {
  char* p = static_cast<char*>(buf);
  size_t got = 0;
  while (got < n) {
    ssize_t r = read_some(fd, p + got, n - got);
    if (r <= 0) break;
    got += static_cast<size_t>(r);
  }
  return got;
}

// Start a program (no shell, argv[0] is searched in PATH) whose stdout we read.
inline int spawn_reader(const std::vector<std::string>& argv, pid_t* pid) {
  int p[2];
  if (pipe(p) < 0) return -1;
  pid_t c = fork();
  if (c < 0) return -1;
  if (c == 0) {
    dup2(p[1], 1);
    close(p[0]); close(p[1]);
    int dn = open("/dev/null", O_RDONLY);
    if (dn >= 0) { dup2(dn, 0); close(dn); }
    std::vector<char*> a;
    for (auto& s : argv) a.push_back(const_cast<char*>(s.c_str()));
    a.push_back(nullptr);
    execvp(a[0], a.data());
    _exit(127);
  }
  close(p[1]);
  *pid = c;
  return p[0];
}

// Start "sh -c cmd" whose stdin we write to (used for the audio player).
inline int spawn_writer(const std::string& cmd, pid_t* pid) {
  int p[2];
  if (pipe(p) < 0) return -1;
  pid_t c = fork();
  if (c < 0) return -1;
  if (c == 0) {
    dup2(p[0], 0);
    close(p[0]); close(p[1]);
    execl("/bin/sh", "sh", "-c", cmd.c_str(), static_cast<char*>(nullptr));
    _exit(127);
  }
  close(p[0]);
  *pid = c;
  return p[1];
}

// TCP server on 127.0.0.1:port; blocks until one client connected.
inline int tcp_listen_accept(int port) {
  int s = socket(AF_INET, SOCK_STREAM, 0);
  int one = 1;
  setsockopt(s, SOL_SOCKET, SO_REUSEADDR, &one, sizeof one);
  sockaddr_in a{};
  a.sin_family = AF_INET;
  a.sin_port = htons(static_cast<uint16_t>(port));
  inet_pton(AF_INET, "127.0.0.1", &a.sin_addr);
  if (bind(s, reinterpret_cast<sockaddr*>(&a), sizeof a) < 0 || listen(s, 1) < 0) {
    perror("bind/listen");
    return -1;
  }
  std::fprintf(stderr, "waiting for a client on 127.0.0.1:%d ...\n", port);
  int c = accept(s, nullptr, nullptr);
  close(s);
  int sndbuf = 1 << 22;
  setsockopt(c, SOL_SOCKET, SO_SNDBUF, &sndbuf, sizeof sndbuf);
  return c;
}

// TCP client; retries for `wait_s` seconds (the encoder may not be up yet).
inline int tcp_connect(const std::string& host, int port, int wait_s = 30) {
  for (int t = 0; t <= wait_s; ++t) {
    int s = socket(AF_INET, SOCK_STREAM, 0);
    sockaddr_in a{};
    a.sin_family = AF_INET;
    a.sin_port = htons(static_cast<uint16_t>(port));
    inet_pton(AF_INET, host.c_str(), &a.sin_addr);
    if (connect(s, reinterpret_cast<sockaddr*>(&a), sizeof a) == 0) return s;
    close(s);
    sleep(1);
  }
  return -1;
}

// ----------------------------------------------------------------------------
// TcpBroadcaster: a TCP server that stays open for the whole lifetime of the
// encoder.  Any number of clients may connect and disconnect at any time; each
// one gets the sample stream from the moment it joins.  Every client has its own
// sender thread and byte queue, so a slow or stalled client can never block the
// real-time encoder: if its queue exceeds `max_queue_bytes` it is disconnected.
// ----------------------------------------------------------------------------
class TcpBroadcaster {
 public:
  using Chunk = std::shared_ptr<std::vector<uint8_t>>;

  ~TcpBroadcaster() { stop(); }

  bool start(int port, size_t max_queue_bytes) {
    max_q_ = max_queue_bytes;
    lfd_ = socket(AF_INET, SOCK_STREAM, 0);
    int one = 1;
    setsockopt(lfd_, SOL_SOCKET, SO_REUSEADDR, &one, sizeof one);
    sockaddr_in a{};
    a.sin_family = AF_INET;
    a.sin_port = htons(static_cast<uint16_t>(port));
    inet_pton(AF_INET, "127.0.0.1", &a.sin_addr);
    if (bind(lfd_, reinterpret_cast<sockaddr*>(&a), sizeof a) < 0 || listen(lfd_, 8) < 0) {
      perror("bind/listen");
      return false;
    }
    acceptor_ = std::thread([this] { accept_loop(); });
    return true;
  }

  // Queue one chunk of the stream for every connected client (never blocks).
  void send(const void* data, size_t n) {
    Chunk c;
    std::lock_guard<std::mutex> l(m_);
    reap();
    if (clients_.empty()) return;
    c = std::make_shared<std::vector<uint8_t>>(static_cast<const uint8_t*>(data), static_cast<const uint8_t*>(data) + n);
    for (auto& cl : clients_) {
      std::lock_guard<std::mutex> ql(cl->m);
      if (cl->bytes > max_q_) { cl->dead = true; shutdown(cl->fd, SHUT_RDWR); std::fprintf(stderr, "client too slow - disconnected\n"); continue; }
      cl->q.push_back(c);
      cl->bytes += n;
      cl->cv.notify_one();
    }
  }

  size_t clients() { std::lock_guard<std::mutex> l(m_); return clients_.size(); }

  void stop() {
    if (stop_.exchange(true)) return;
    if (lfd_ >= 0) { shutdown(lfd_, SHUT_RDWR); close(lfd_); }
    if (acceptor_.joinable()) acceptor_.join();
    std::lock_guard<std::mutex> l(m_);
    for (auto& cl : clients_) cl->dead = true;
    reap(true);
  }

 private:
  struct Client {
    int fd = -1;
    std::thread th;
    std::mutex m;
    std::condition_variable cv;
    std::deque<Chunk> q;
    size_t bytes = 0;
    bool dead = false, finish = false;
  };

  void accept_loop() {
    while (!stop_) {
      pollfd p{lfd_, POLLIN, 0};
      if (poll(&p, 1, 200) <= 0) continue;
      sockaddr_in ca{};
      socklen_t cl = sizeof ca;
      int fd = accept(lfd_, reinterpret_cast<sockaddr*>(&ca), &cl);
      if (fd < 0) continue;
      int sndbuf = 1 << 22;
      setsockopt(fd, SOL_SOCKET, SO_SNDBUF, &sndbuf, sizeof sndbuf);
      auto c = std::make_unique<Client>();
      c->fd = fd;
      Client* raw = c.get();
      c->th = std::thread([this, raw] { sender(raw); });
      std::lock_guard<std::mutex> l(m_);
      clients_.push_back(std::move(c));
      std::fprintf(stderr, "client connected (%zu now)\n", clients_.size());
    }
  }

  // One thread per client: blocking writes, isolated from the encoder.
  void sender(Client* c) {
    for (;;) {
      Chunk ch;
      {
        std::unique_lock<std::mutex> l(c->m);
        c->cv.wait(l, [&] { return !c->q.empty() || c->dead || c->finish; });
        if (c->dead || (c->q.empty() && c->finish)) break;
        ch = c->q.front();
        c->q.pop_front();
        c->bytes -= ch->size();
      }
      if (!write_all(c->fd, ch->data(), ch->size())) break;
    }
    std::lock_guard<std::mutex> l(c->m);
    c->dead = true;
  }

  // Remove finished clients (caller holds m_).
  void reap(bool all = false) {
    for (size_t i = 0; i < clients_.size();) {
      Client& c = *clients_[i];
      bool dead;
      { std::lock_guard<std::mutex> l(c.m); dead = c.dead; }
      if (dead || all) {
        { std::lock_guard<std::mutex> l(c.m); c.dead = true; c.cv.notify_one(); }
        shutdown(c.fd, SHUT_RDWR);
        if (c.th.joinable()) c.th.join();
        close(c.fd);
        clients_.erase(clients_.begin() + i);
        if (!all) std::fprintf(stderr, "client disconnected (%zu left)\n", clients_.size());
      } else ++i;
    }
  }

  int lfd_ = -1;
  size_t max_q_ = 0;
  std::atomic<bool> stop_{false};
  std::thread acceptor_;
  std::mutex m_;
  std::vector<std::unique_ptr<Client>> clients_;
};

inline int64_t now_ns() {
  timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return int64_t(ts.tv_sec) * 1000000000 + ts.tv_nsec;
}
inline void sleep_until_ns(int64_t t) {
  timespec ts{time_t(t / 1000000000), long(t % 1000000000)};
  clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, &ts, nullptr);
}

}  // namespace dsr::io
