#include "particles.hpp"

//Particles******************************************************
Particles::Particles(): Nkinds(1) {
  Npart = new int[1];
  M = new double[1];
  R = new double[1];
  Q = new double[1];

  pos = new double[3];
  mom = new double[3];
  force = new double[3];
}

Particles::Particles(int kinds, const vector<int>& N, const vector<double>& m, const vector<double>& r, const vector<double>& q): Nkinds(kinds) {
  Npart = new int[Nkinds];
  M = new double[Nkinds];
  R = new double[Nkinds];  
  Q = new double[Nkinds];
  
  int tot = 0;
  for(int i = 0; i < Nkinds; ++i){
    Npart[i] = N[i];
    M[i] = m[i];  
    R[i] = r[i];
    Q[i] = q[i];
    tot += N[i];
  }

  pos = new double[3*tot];
  mom = new double[3*tot];
  force = new double[3*tot];
}

Particles::~Particles() {
  delete[] Npart;
  delete[] M;
  delete[] R;
  delete[] Q;
  
  delete[] pos;
  delete[] mom;
  delete[] force;
}

int Particles::get_Nkinds() {
  return Nkinds;
}

int Particles::get_N(int i) {
  return Npart[i];
}

int* Particles::get_Npart() {
  return Npart;
}

double Particles::get_mass(int i) {
  return M[i];
}

double* Particles::get_M() {
  return M;
}

double Particles::get_rad(int i) {
  return R[i];
}

double* Particles::get_R() {
  return R;
}

double Particles::get_charge(int i) {
  return Q[i];
}

double* Particles::get_Q() {
  return Q;
}

int Particles::get_Nkinds() const {
  return Nkinds;
}

int Particles::get_N(int i) const {
  return Npart[i];
}

int* Particles::get_Npart() const {
  return Npart;
}

double Particles::get_mass(int i) const {
  return M[i];
}

double* Particles::get_M() const{
  return M;
}

double Particles::get_rad(int i) const {
  return R[i];
}

double* Particles::get_R() const {
  return R;
}

double Particles::get_charge(int i) const {
  return Q[i];
}

double* Particles::get_Q() const {
  return Q;
}

void Particles::set_charge(int i, double q) {
  Q[i] = q;
}


//Get coordinate "a" from particle "i" of species "n"; 
double Particles::get_pos(int n, int i, int a) {
  //check boundaries
  if(i > Npart[n]-1) printf("Particle index out of bounds.\n");
  if(a > 2) printf("Component out of bounds.\n");

  //get the index of the particle
  int part = 0;
  for(int l = 0; l < n; ++l)
    part += Npart[l];
  part += i;

  //get component "a" of the position
  return pos[part*3+a];
}

//Get coordinate "a" from particle "i" of species "n"; 
double Particles::get_pos(int n, int i, int a) const {
  //check boundaries
  if(i > Npart[n]-1) printf("Particle index out of bounds.\n");
  if(a > 2) printf("Component out of bounds.\n");

  //get the index of the particle
  int part = 0;
  for(int l = 0; l < n; ++l)
    part += Npart[l];
  part += i;

  //get component "a" of the position
  return pos[part*3+a];
}

//Set coordinate "a" from particle "i" of species "n"
void Particles::set_pos(int n, int i, int a, double val){
  //check boundaries
  if(i > Npart[n]-1) printf("Particle index out of bounds.\n");
  if(a > 2) printf("Component out of bounds.\n");
  
  //get the index of the particle
  int part = 0;
  for(int l = 0; l < n; ++l)
    part += Npart[l];
  part += i;

  pos[part*3+a] = val;
}

double* Particles::get_X(){
  return pos;
}

double* Particles::get_X() const {
  return pos;
}

//Get momentum component "a" from particle "i" of species "n"; 
double Particles::get_mom(int n, int i, int a) {
  //check boundaries
  if(i > Npart[n]-1) printf("Particle index out of bounds.\n");
  if(a > 2) printf("Component out of bounds.\n");

  //get the index of the particle
  int part = 0;
  for(int l = 0; l < n; ++l)
    part += Npart[l];
  part += i;

  //get component "a" of the position
  return mom[part*3+a];
}

//Get momentum component "a" from particle "i" of species "n"; 
double Particles::get_mom(int n, int i, int a) const {
  //check boundaries
  if(i > Npart[n]-1) printf("Particle index out of bounds.\n");
  if(a > 2) printf("Component out of bounds.\n");

  //get the index of the particle
  int part = 0;
  for(int l = 0; l < n; ++l)
    part += Npart[l];
  part += i;

  //get component "a" of the position
  return mom[part*3+a];
}

//Set momentum component "a" from particle "i" of species "n"
void Particles::set_mom(int n, int i, int a, double val){
  //check boundaries
  if(i > Npart[n]-1) printf("Particle index out of bounds.\n");
  if(a > 2) printf("Component out of bounds.\n");
  
  //get the index of the particle
  int part = 0;
  for(int l = 0; l < n; ++l)
    part += Npart[l];
  part += i;

  mom[part*3+a] = val;
}

double* Particles::get_P(){
  return mom;
}

double* Particles::get_P() const {
  return mom;
}

//Get force component "a" on particle "i" of species "n"; 
double Particles::get_F(int n, int i, int a) {
  //check boundaries
  if(i > Npart[n]-1) printf("Particle index out of bounds.\n");
  if(a > 2) printf("Component out of bounds.\n");

  //get the index of the particle
  int part = 0;
  for(int l = 0; l < n; ++l)
    part += Npart[l];
  part += i;

  //get component "a" of the position
  return force[part*3+a];
}

//Get force component "a" on particle "i" of species "n"; 
double Particles::get_F(int n, int i, int a) const{
  //check boundaries
  if(i > Npart[n]-1) printf("Particle index out of bounds.\n");
  if(a > 2) printf("Component out of bounds.\n");

  //get the index of the particle
  int part = 0;
  for(int l = 0; l < n; ++l)
    part += Npart[l];
  part += i;

  //get component "a" of the position
  return force[part*3+a];
}

//Get force vector 
double* Particles::get_F() {
  return force;
}

//Get force vector 
double* Particles::get_F() const{
  return force;
}

//Set forces to zero
void Particles::set_zero_F(){
  int tot = 0;
  for(int i = 0; i < Nkinds; ++i)
    tot += Npart[i];
  
  for(int i = 0; i < 3*tot; ++i)
    force[i] = 0.;
}

//Add to force component "a" on particle "i" of species "n"
void Particles::add_F(int n, int i, int a, double val){
  //check boundaries
  if(i > Npart[n]-1) printf("Particle index out of bounds.\n");
  if(a > 2) printf("Component out of bounds.\n");
  
  //get the index of the particle
  int part = 0;
  for(int l = 0; l < n; ++l)
    part += Npart[l];
  part += i;

  #pragma omp atomic
  force[part*3+a] += val;
}

//Set the force vector
void Particles::set_F(double *Fin, int size){
  int Tpart = 0;
  for(int n =0; n < Nkinds; ++n)
    Tpart += Npart[n];
  if(size != 3*Tpart) printf("Wrong array size in set_F!\n");

  for(int i = 0; i < size; ++i)
    force[i] = Fin[i];
}

//Get the kind and number of the particle
void Particles::get_kind(int id, int& n, int& i){
  int Nto=0, Ntn=0;
  for(int nn = 0; nn < Nkinds; ++nn){
    Nto = Ntn;
    Ntn += Npart[nn];
    if(id < Ntn){
      n = nn;
      i = id - Nto;
      break;
    }
  }
}

//Get the kind and number of the particle
void Particles::get_kind(int id, int& n, int& i) const{
  int Nto=0, Ntn=0;
  for(int nn = 0; nn < Nkinds; ++nn){
    Nto = Ntn;
    Ntn += Npart[nn];
    if(id < Ntn){
      n = nn;
      i = id - Nto;
      break;
    }
  }
}
