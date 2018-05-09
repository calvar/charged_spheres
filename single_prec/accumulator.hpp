#ifndef __ACCUMULATOR_HPP
#define __ACCUMULATOR_HPP


/**
 *Class to accumulate the values of the variables
 */
class accumulator {
private:
  float m_K, m_U, m_P;
  float m_Usr, m_Psr;
  float m_Uc, m_Pc;
  
public:
  accumulator();
  ~accumulator() {}
  
  void set_K(float K);
  float get_K();
  float get_K() const;
  void set_U(float U);
  float get_U();
  float get_U() const;
  void set_P(float P);
  float get_P();
  float get_P() const;
  void set_Usr(float Usr);
  float get_Usr();
  float get_Usr() const;
  void set_Psr(float Psr);
  float get_Psr();
  float get_Psr() const;
  void set_Uc(float Uc);
  float get_Uc();
  float get_Uc() const;
  void set_Pc(float Pc);
  float get_Pc();
  float get_Pc() const;

  void set_zero();
  void operator=(const accumulator&);
  void operator+=(const accumulator&);
  const accumulator operator*(const accumulator&);
  const accumulator operator/(long);
};




#endif
