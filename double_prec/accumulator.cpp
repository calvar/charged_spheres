#include "accumulator.hpp"


accumulator::accumulator() : m_K(0.), m_U(0.), m_P(0.), m_Usr(0.), m_Psr(0.),
		  m_Uc(0.), m_Pc(0.) {}

void accumulator::set_K(double K) { m_K = K; }

double accumulator::get_K() { return m_K; }

double accumulator::get_K() const{ return m_K; }

void accumulator::set_U(double U) { m_U = U; }

double accumulator::get_U() { return m_U; }

double accumulator::get_U() const{ return m_U; }

void accumulator::set_P(double P) { m_P = P; }

double accumulator::get_P() { return m_P; }

double accumulator::get_P() const{ return m_P; }

void accumulator::set_Uss(double Usr) { m_Usr = Usr; }

double accumulator::get_Usr() { return m_Usr; }

double accumulator::get_Usr() const{ return m_Usr; }

void accumulator::set_Psr(double Psr) { m_Psr = Psr; }

double accumulator::get_Psr() { return m_Psr; }

double accumulator::get_Psr() const{ return m_Psr; }

void accumulator::set_UchR(double Uc) { m_Uc = Uc; }

double accumulator::get_Uc() { return m_Uc; }

double accumulator::get_Uc() const{ return m_Uc; }

void accumulator::set_Pc(double Pc) { m_Pc = Pc; }

double accumulator::get_Pc() { return m_Pc; }
 
double accumulator::get_Pc() const{ return m_Pc; }

void accumulator::set_zero() {
  m_K = 0.; m_U = 0.; m_P = 0.;
  m_Usr = 0.; m_Psr = 0.; m_Uc = 0.; m_Pc = 0.; 
}

void accumulator::operator=(const accumulator& acc) {
  m_K = acc.m_K;
  m_U = acc.m_U;
  m_P = acc.m_P;
  m_Usr = acc.m_Usr;
  m_Psr = acc.m_Psr;
  m_Uc = acc.m_Uc;
  m_Pc = acc.m_Pc;
}

void accumulator::operator+=(const accumulator& acc) {
  m_K += acc.m_K;
  m_U += acc.m_U;
  m_P += acc.m_P;
  m_Usr += acc.m_Usr;
  m_Psr += acc.m_Psr;
  m_Uc += acc.m_Uc;
  m_Pc += acc.m_Pc;
}
  
const accumulator accumulator::operator*(const accumulator& acc) {
  accumulator result;
  result.m_K = m_K * acc.m_K;
  result.m_U = m_U * acc.m_U;
  result.m_P = m_P * acc.m_P;
  result.m_Usr = m_Usr * acc.m_Usr;
  result.m_Psr = m_Psr * acc.m_Psr;
  result.m_Uc = m_Uc * acc.m_Uc;
  result.m_Pc = m_Pc * acc.m_Pc;
  return result;
}

const accumulator accumulator::operator/(long a) {
  accumulator result;
  result.m_K = m_K / a;
  result.m_U = m_U / a;
  result.m_P = m_P / a;
  result.m_Uss = m_Usr / a;
  result.m_Pss = m_Psr / a;
  result.m_UchR = m_Uc / a;
  result.m_PchR = m_Pc / a;
  return result;
}
