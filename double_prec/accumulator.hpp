#ifndef __ACCUMULATOR_HPP
#define __ACCUMULATOR_HPP


/**
 *Class to accumulate the values of the variables
 */
class accumulator {
private:
  double m_K, m_U, m_P;
  double m_Usr, m_Psr;
  double m_Uc, m_Pc;
  
public:
  accumulator();
  ~accumulator() {}
  
  void set_K(double K);
  double get_K();
  double get_K() const;
  void set_U(double U);
  double get_U();
  double get_U() const;
  void set_P(double P);
  double get_P();
  double get_P() const;
  void set_Usr(double Usr);
  double get_Usr();
  double get_Usr() const;
  void set_Psr(double Psr);
  double get_Psr();
  double get_Psr() const;
  void set_Uc(double Uc);
  double get_Uc();
  double get_Uc() const;
  void set_Pc(double Pc);
  double get_Pc();
  double get_Pc() const;

  void set_zero();
  void operator=(const accumulator&);
  void operator+=(const accumulator&);
  const accumulator operator*(const accumulator&);
  const accumulator operator/(long);
};




#endif
