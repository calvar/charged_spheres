#ifndef __TIMER_HPP
#define __TIMER_HPP

#include <ctime>
#include <sys/time.h>

/**
 * Class to keep track of time
 */
class Timer {
 private:
  timeval m_tim;
  long m_tcpu;
  
 public:
  Timer();
  ~Timer() {}
  
   void start();
   void stop(double&, long&);
};


#endif
