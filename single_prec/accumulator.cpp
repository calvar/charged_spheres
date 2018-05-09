#include "accumulator.hpp"


accumulator::accumulator() : m_K(0.), m_U(0.), m_P(0.), m_Usr(0.), m_Psr(0.),
		  m_Uc(0.), m_Pc(0.) {}

void accumulator::set_K(float K) { m_K = K; }

float accumulator::get_K() { return m_K; }

float accumulator::get_K() const{ return m_K; }

void accumulator::set_U(float U) { m_U = U; }

float accumulator::get_U() { return m_U; }

float accumulator::get_U() const{ return m_U; }

void accumulator::set_P(float P) { m_P = P; }

float accumulator::get_P() { return m_P; }

float accumulator::get_P() const{ return m_P; }

void accumulator::set_Uss(float Usr) { m_Usr = Usr; }

float accumulator::get_Usr() { return m_Usr; }

float accumulator::get_Usr() const{ return m_Usr; }

void accumulator::set_Psr(float Psr) { m_Psr = Psr; }

float accumulator::get_Psr() { return m_Psr; }

float accumulator::get_Psr() const{ return m_Psr; }

void accumulator::set_UchR(float Uc) { m_Uc = Uc; }

float accumulator::get_Uc() { return m_Uc; }

float accumulator::get_Uc() const{ return m_Uc; }

void accumulator::set_Pc(float Pc) { m_Pc = Pc; }

float accumulator::get_Pc() { return m_Pc; }
 
float accumulator::get_Pc() const{ return m_Pc; }

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
