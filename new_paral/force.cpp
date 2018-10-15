#include "force.hpp"

//Compute the terms including the erfc function and store in a table for
//certain values.
void tabul(double *Tab, double alpha){
  double extL = 200 * 1.001;
  double alsq = alpha * alpha;
  double spi = sqrt(M_PI);
  for(int i = 1; i < TAB_SIZE+1; i++){
    double r =  0.5 * extL * i / TAB_SIZE;
    double rsq = r * r;
    double sr1 = 1.0 / r;
    double ERFC = erfc(alpha * r);
    double aer = 2.0 * alpha * exp(- alsq * rsq) / spi;
    Tab[i-1] = ERFC * sr1;
    Tab[TAB_SIZE+(i-1)] = Tab[i-1] + aer;
  }
}

bool force_R(Particles& part, const double L[3], const double *Tab, double& U_sr, double& W_sr, double& U_c, double& W_c, double Dt){
  //get total number of particles
  int Tpart = 0;
  for(int i = 0; i < part.get_Nkinds(); ++i)
    Tpart += part.get_N(i);
  
  //set forces to 0
  part.set_zero_F();
  //short range
  double Usr = 0.;
  double Wsr = 0.;
  double srt = pow(2., 1./6);
  //Coulomb
  double Uc = 0.;
  double Wc = 0.;
  double Lh = 0.5 * (200*1.001) / TAB_SIZE;
  
  //Loop over all the particles****************************************************
  bool abort = false; //If there is an overlap abort the summation
  for(int i = 0; i < Tpart; ++i){
    if(!abort){
      double ri[3], rj[3], f[3];
      double vi[3], vj[3];
      double qi, qj;
      double massi, massj;
      
      int m, mi;
      part.get_kind(i, m, mi);
      massi = part.get_mass(m);
      for(int a = 0; a < 3; ++a){
	ri[a] = part.get_pos(m, mi, a);
	vi[a] = part.get_mom(m, mi, a);
	vi[a] /= massi;
      }
      qi = part.get_charge(m);
      
      //even out the charge among threads*****************************
      int mx = static_cast<int>(ceil(static_cast<float>(Tpart-1)/2));
      if(fmod(static_cast<float>(Tpart),2) == 0. && i >= Tpart/2)
	mx = static_cast<int>(floor(static_cast<float>(Tpart-1)/2));
      
      int j = i+1 - Tpart*static_cast<int>(floor((i+1)/Tpart + 0.5));
      int cnt = 0;
      //**************************************************************
      while(cnt < mx){
	int n, nj;
	part.get_kind(j, n, nj);
	massj = part.get_mass(n);
	for(int a = 0; a < 3; ++a){
	  rj[a] = part.get_pos(n, nj, a);
	  vj[a] = part.get_mom(n, nj, a);
	  vj[a] /= massj;
	}
	qj = part.get_charge(n);
	
	double rij[3];
	for(int a = 0; a < 3; ++a){
	  rij[a] = ri[a] - rj[a];
	  rij[a] -= L[a] * floor(rij[a] / L[a] + 0.5);
	}
	
	//cutoff distance for the real part (cube)
	if(fabs(rij[0])<L[0]/2 && fabs(rij[1])<L[1]/2 && fabs(rij[2])<L[2]/2){
	  double RIJSQ;
	  dot(rij, rij, RIJSQ);
	  double RIJ = sqrt(RIJSQ);
	  
	  double sigma = part.get_rad(m) + part.get_rad(n);
	  double sig_sq = sigma * sigma;
	  double rsh = srt * sigma;
	  
	  if(RIJ < MX_OVRLP * sigma){ //check for overlaps
	    abort = true;
	  }
	  
	  //Lennard-Jones shifted potential, virial and force----------------------
	  double SR2s = sig_sq / RIJSQ;
	  double SR6s = SR2s * SR2s * SR2s;
	  double pot = 4 * SR6s * (SR6s - 1.);
	  double UIJ = 0.;
	  double WIJ = 0.;
	  if(RIJ < rsh){
	    UIJ = pot + 1.;//shifted potential (epsilon = 1)
	    WIJ = 6 * (pot + 4*SR6s*SR6s);
	  }
	  Usr += UIJ;
	  Wsr += WIJ;
	  
	  //short range force
	  double FIJsr = WIJ / RIJSQ;
	  for(int a = 0; a < 3; ++a){
	    f[a] = FIJsr * rij[a];
	    part.add_F(m, mi, a, f[a]);
	  }
	  for(int a = 0; a < 3; ++a){
	    f[a] *= -1;
	    part.add_F(n, nj, a, f[a]);
	  }
	  //----------------------------------------------------------------------
	  
	  //Coulomb potential, virial and force-----------------------------------
	  double qq = qi * qj;
	  if(fabs(qq) > 1.e-12){
	    //make a linear interpolation from the values on the table
	    int ri = static_cast<int>(floor(RIJ / Lh));
	    double R_L0 = RIJ - Lh * ri;
	    double R_L1 = RIJ - Lh * (ri+1);
	    double A = (Tab[ri] * R_L0 - Tab[ri-1] * R_L1) / Lh;
	    double B = (Tab[ri+TAB_SIZE] * R_L0 - Tab[(ri-1)+TAB_SIZE] * R_L1) / Lh;

	    //potential and virial
	    UIJ = qq * A;
	    WIJ = qq * B;
	    Uc += UIJ;
	    Wc += WIJ;
	    
	    //Coulomb force
	    for(int a = 0; a < 3; ++a){
	      f[a] = WIJ * rij[a] / RIJSQ;
	      part.add_F(m, mi, a, f[a]);
	    }
	    for(int a = 0; a < 3; ++a){
	      f[a] *= -1;
	      part.add_F(n, nj, a, f[a]);
	    }
	  }
	  //----------------------------------------------------------------------
	}
	//update j
	j += 1 - Tpart*static_cast<int>(floor((j+1)/Tpart + 0.5));
	cnt++;
      }
      //**************************************************************
    }
  }
  if(abort)
    return false;
  
  U_sr = Usr;
  W_sr = Wsr;
  U_c = Uc;
  W_c = Wc;
  
  return true;
}


//Compute the wave vectors
void wave_v(double *Kvec, const double L[3], double alpha, int K_max){
  double c = 4. * alpha * alpha;
  double P2[3] = {2.*M_PI/L[0], 2.*M_PI/L[1], 2.*M_PI/L[2]};
  int sz = K_max+1;
  int sz2 = sz*sz;
  for(int nz = 0; nz <= K_max; ++nz){
    double Kz = P2[2] * nz;
    for(int ny = 0; ny <= K_max; ++ny){
      double Ky = P2[1] * ny;
      for(int nx = 0; nx <= K_max; ++nx){
	double Kx = P2[0] * nx;
	double Ksq = Kx * Kx + Ky * Ky + Kz * Kz;
	if(Ksq != 0.){
	  Kvec[nz*sz2 + ny*sz + nx] = 4. * M_PI * exp(- Ksq / c) / Ksq; 
	}else{
	  Kvec[nz*sz2 + ny*sz + nx] = 0.;
	}
      }
    }
  }
}

void force_K(Particles& part, const double L[3], int K_max, double alpha, const double *Kvec, double& U, double& W){
  int Tpart = 0;
  for(int i = 0; i < part.get_Nkinds(); ++i)
    Tpart += part.get_N(i);
  
  double V = L[0] * L[1] * L[2];
  double P2[3] = {2.*M_PI/L[0], 2.*M_PI/L[1], 2.*M_PI/L[2]};

  K_max += 1;
  int K_max2 = K_max * K_max;
  int K_max3 = K_max2 * K_max;
  double Uc = 0.;
  double Wc = 0.;

  for(int kn = 0; kn < K_max3; ++kn){
    //Integer vectors, nx, ny, nz
    int nz = static_cast<int>(floor(static_cast<float>(kn)/K_max2));
    float l = kn - K_max2 * nz;
    int ny = static_cast<int>(floor(l/K_max));
    int nx = static_cast<int>(l) - K_max * ny;
    //K vectors
    double Kz = P2[2] * nz;
    double Ky = P2[1] * ny;
    double Kx = P2[0] * nx;
    double Ksq = Kx * Kx + Ky * Ky + Kz * Kz;
    
    if(Ksq != 0){
      //
      int vali = nz*K_max2+ny*K_max+nx;
      if(kn != vali){
	printf("kn != nz*K_max2+ny*K_max+nx ! ");
	printf("%d != %d\n",kn,vali);
      }
      //
      double KK = Kvec[kn];
      double val = 2. * KK / V;  //mult by 2 for symmetry
      
      double K[4][3];
      K[0][0] = Kx; K[0][1] = Ky; K[0][2] = Kz;
      K[1][0] = -Kx; K[1][1] = Ky; K[1][2] = Kz;
      K[2][0] = Kx; K[2][1] = -Ky; K[2][2] = Kz;
      K[3][0] = -Kx; K[3][1] = -Ky; K[3][2] = Kz;
      
      //sum array
      complex<double> cc(0., 0.);
      complex<double> sum[8];
      for(int b = 0; b < 8; ++b)
	sum[b] = cc;
      
      double self = 0.; //extract the self energy corresponding to each term (Kmax finite)
      
      //individual values array for each particle
      vector<vector<complex<double> > >  M(4, vector<complex<double> >(Tpart, cc));
      
      double ri[3], qi;
      
      //generate sum of k vectors************************************
      for(int i = 0; i < Tpart; ++i){
	int n, ni;
	part.get_kind(i, n, ni);
	for(int a = 0; a < 3; ++a)
	  ri[a] = part.get_pos(n, ni, a);
	qi = part.get_charge(n);
	
	if(qi != 0.){
	  double RIK[4];
	  complex<double> z[4];
	  for(int b = 0; b < 4; ++b){
	    dot(ri, K[b], RIK[b]);
	    z[b] = polar(1., RIK[b]);
	  }
	  
	  for(int b = 0; b < 4; ++b){
	    M[b][i] = conj(z[b]);      
	    sum[b] += qi * M[b][i];	
	    sum[b+4] += qi * RIK[b] * z[b];
	  }
	  self += qi * qi;
	}
      }
      //***********************************************************
      
      //compute energy - virial--------------------------------
      double Upar = 0, Wpar = 0;
      for(int b = 0; b < 4; ++b){
	Upar += val * 0.5 * (norm(sum[b]) - self);
	Wpar -= val * imag(sum[b+4] * sum[b]);
      }
      if((nx == 0 && ny == 0) || (nx == 0 && nz == 0) || (ny == 0 && nz == 0)){
	Upar /= 4.;
	Wpar /= 4.;
      }else if((nz == 0 && nx != 0 && ny != 0) || (ny == 0 && nz != 0 && nx != 0) || (nx == 0 && ny != 0 && nz != 0)){
	Upar /= 2.;
	Wpar /= 2.;
      }
      Uc += Upar;       
      Wc += Wpar;
	
      //force (needs a cicle)**************************************
      double f[3];
      for(int i = 0; i < Tpart; ++i){
	int n, ni;
	part.get_kind(i, n, ni);
	for(int a = 0; a < 3; ++a)
	  ri[a] = part.get_pos(n, ni, a);
	qi = part.get_charge(n);
	
	if(qi != 0.){
	  double fact = 0.;
	  //as force is not scalar, different values of "a" are not always equivalent.
	  for(int a = 0; a < 3; ++a)
	    f[a] = 0.;
	  for(int b = 0; b < 4; ++b){
	    fact = -val * imag(qi * M[b][i] * conj(sum[b])); //self force term = 0
	    for(int a = 0; a < 3; ++a)
	      f[a] += fact * K[b][a];
	  }
	  if((nx == 0 && ny == 0) || (nx == 0 && nz == 0) || (ny == 0 && nz == 0)){
	    for(int a = 0; a < 3; ++a)
	      f[a] /= 4.;
	  }else if((nz == 0 && nx != 0 && ny != 0) || (ny == 0 && nz != 0 && nx != 0) || (nx == 0 && ny != 0 && nz != 0)){
	    for(int a = 0; a < 3; ++a)
	      f[a] /= 2.;
	    
	  }
	  for(int a = 0; a < 3; ++a)
	    part.add_F(n, ni, a, f[a]);  
	}
      }
      //*********************************************************
      
    }
  }
  U = Uc;
  W = Wc;
}


void dot(const double a[3], const double b[3], double& ans){
  ans = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
