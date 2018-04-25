#include "particles.hpp"

/**
 * Default constructor
 */
Particles::Particles() {
  m_N = 0;
  m_R = 1.;
  m_M = 1.;
  m_Q = 0.;

  m_pos = new double[3];
  m_mom = new double[3];
  m_force = new double[3];
}

/**
 * Constructor that accepts the number of particles
 * and their parameters
 */
Particles::Particles(int N, double R, double M, double Q) {
  m_N = N;
  m_R = R;
  m_M = M;
  m_Q = Q;

  m_pos = new double[3*N];
  m_mom = new double[3*N];
  m_force = new double[3*N];
}

/**
 * Destructor
 */
Particles::~Particles() {
  delete[] m_pos;
  delete[] m_mom;
  delete[] m_force;
}

/**
 * Getter methods
 */
int Particles::getNoOfParticles() {
  return m_N;
}

int Particles::getNoOfParticles() const {
  return m_N;
}

double Particles::getRadius() {
  return m_R;
}

double Particles::getRadius() const {
  return m_R;
}

double Particles::getMass() {
  return m_M;
}

double Particles::getMass() const {
  return m_M;
}

double Particles::getCharge() {
  return m_Q;
}

double Particles::getCharge() const {
  return m_Q;
}

/**
 * Get position, velocity or force for the dim dimension 
 * of the ith particle
 */
double Particles::getPosition(int i, int dim) {
  //check boundaries
  if(i < 0 || i > m_N-1) printf("Particle index out of bounds.\n");
  if(dim < 0 || dim > 2) printf("Dimension out of bounds.\n");

  //return component "dim" of the position of particle "i"
  return m_pos[i*3+dim];
}

double Particles::getPosition(int i, int dim) const {
  //check boundaries
  if(i < 0 || i > m_N-1) printf("Particle index out of bounds.\n");
  if(dim < 0 || dim > 2) printf("Dimension out of bounds.\n");

  //return component "dim" of the position of particle "i"
  return m_pos[i*3+dim];
}

double Particles::getMomentum(int i, int dim) {
  //check boundaries
  if(i < 0 || i > m_N-1) printf("Particle index out of bounds.\n");
  if(dim < 0 || dim > 2) printf("Dimension out of bounds.\n");

  //return component "dim" of the momentum of particle "i"
  return m_mom[i*3+dim];
}

double Particles::getMomentum(int i, int dim) const {
  //check boundaries
  if(i < 0 || i > m_N-1) printf("Particle index out of bounds.\n");
  if(dim < 0 || dim > 2) printf("Dimension out of bounds.\n");

  //return component "dim" of the momentum of particle "i"
  return m_mom[i*3+dim];
}

double Particles::getForce(int i, int dim) {
  //check boundaries
  if(i < 0 || i > m_N-1) printf("Particle index out of bounds.\n");
  if(dim < 0 || dim > 2) printf("Dimension out of bounds.\n");

  //return component "dim" of the force of particle "i"
  return m_force[i*3+dim];
}

double Particles::getForce(int i, int dim) const {
  //check boundaries
  if(i < 0 || i > m_N-1) printf("Particle index out of bounds.\n");
  if(dim < 0 || dim > 2) printf("Dimension out of bounds.\n");

  //return component "dim" of the force of particle "i"
  return m_force[i*3+dim];
}


/**
 * Get the pointer to the positions, momenta or forces of
 * all the particles (to pass to the device).
 */
double* Particles::getPositions() {
  return m_pos;
}

double* Particles::getPositions() const {
  return m_pos;
}

double* Particles::getMomenta() {
  return m_mom;
}

double* Particles::getMomenta() const {
  return m_mom;
}

double* Particles::getForces() {
  return m_force;
}

double* Particles::getForces() const {
  return m_force;
}


/**
 * Setter methods
 */
void Particles::setNoOfParticles(int N) {
  m_N = N;
}

void Particles::setRadius(double R) {
  m_R = R;
}

void Particles::setMass(double M) {
  m_M = M;
}

void Particles::setCharge(double Q) {
  m_Q = Q;
}

/**
 * Set the position, velocity or force for the dim dimension 
 * of the ith particle
 */
void Particles::setPosition(int i, int dim, double val) {
  //check boundaries
  if(i < 0 || i > m_N-1) printf("Particle index out of bounds.\n");
  if(dim < 0 || dim > 2) printf("Dimension out of bounds.\n");

  //set component "dim" of the position of particle "i" to "val"
  m_pos[i*3+dim] = val;
}

void Particles::setMomentum(int i, int dim, double val) {
  //check boundaries
  if(i < 0 || i > m_N-1) printf("Particle index out of bounds.\n");
  if(dim < 0 || dim > 2) printf("Dimension out of bounds.\n");

  //set component "dim" of the momentum of particle "i" to "val"
  m_mom[i*3+dim] = val;
}

void Particles::setForce(int i, int dim, double val) {
  //check boundaries
  if(i < 0 || i > m_N-1) printf("Particle index out of bounds.\n");
  if(dim < 0 || dim > 2) printf("Dimension out of bounds.\n");

  //set component "dim" of the force of particle "i" to "val"
  m_force[i*3+dim] = val;
}

/**
 * Set the pointer to the forces of all the particles 
 * (to copy from the device).
 */
void Particles::setForces(double* F, int size) {
  //check size
  if(size != 3*m_N) printf("Wrong array size at setForces\n");
  
  for(int i = 0; i < 3*m_N; ++i)
    m_force[i] = F[i];
}

/**
 * Add contribution to the total force on a particle 
 */
void Particles::addToForce(int i, int dim, double val) {
  //check boundaries
  if(i < 0 || i > m_N-1) printf("Particle index out of bounds.\n");
  if(dim < 0 || dim > 2) printf("Dimension out of bounds.\n");

  //Add force contribution "val" to component "dim" of particle "i"
  m_force[i*3+dim] += val;
}

/**
 * Set forces to zero 
 */
void Particles::setForcesToZero() {
  for(int i = 0; i < 3*m_N; ++i)
    m_force[i] = 0.;
}
