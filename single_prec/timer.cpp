#include "timer.hpp"


Timer::Timer() {
  gettimeofday(&m_tim, NULL); 
  m_tcpu = clock();
}

void Timer::start() {
  gettimeofday(&m_tim, NULL); 
  m_tcpu = clock();
}

void Timer::stop(float& telap, long& twork) { 
  float t1 = m_tim.tv_sec * 1e+6 + m_tim.tv_usec;
  long start = m_tcpu;
  gettimeofday(&m_tim, NULL);
  float t2 = m_tim.tv_sec * 1e+6 + m_tim.tv_usec;
  m_tcpu = clock();
  long end = m_tcpu;
  telap = t2 - t1;
  twork = end - start;
}
