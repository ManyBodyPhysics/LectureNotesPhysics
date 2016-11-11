#ifndef TIMER_H
#define TIMER_H

class timer {
  public:
    timer() : isStarted(false) {}
    void start() {
      clock_gettime(CLOCK_MONOTONIC, &started);
      isStarted = true;
    }

    double get_elapsed() {
      if (notStarted() ) { return 0.0; }
      struct timespec stopped;
      clock_gettime(CLOCK_MONOTONIC, &stopped);
      isStarted = false;

      double elapsed = (stopped.tv_sec - started.tv_sec);
      elapsed += (stopped.tv_nsec - started.tv_nsec)/1000000000.0;

      return elapsed;
    }
  private:
    bool notStarted() { return !isStarted; }
    bool isStarted;
    struct timespec started;
};
#endif
