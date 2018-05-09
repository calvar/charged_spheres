#include "force.hpp"

//Compute the total short range + Coloumb force on each particle
bool force_R(Particles& part, float L, float& Usr, float& Wsr, float& Uc, float& Wc){
  int N = part.getNoOfParticles();
  
  //set forces to 0
  part.setForcesToZero();

  float sigma = 2*part.getRadius();
  float sig_sq = sigma * sigma;
  float srt = pow(2., 1./6);
  float rsh = srt * sigma;
  
  //short range
  Usr = 0.;
  Wsr = 0.;
  //Coulomb
  Uc = 0.;
  Wc = 0.;

  float ri[3], rj[3], f[3];
  float q = part.getCharge();
  
  for(int i = 0; i < N; ++i){
    
    for(int a = 0; a < 3; ++a)
      ri[a] = part.getPosition(i, a);
    float qi;
    (i%2 == 0) ?  qi = q : qi = -q;

    //even out the load among threads
    int mx = static_cast<int>(ceil(static_cast<float>(N-1)/2));
    if(N%2 == 0. && i >= N/2)
       mx = static_cast<int>(floor(static_cast<float>(N-1)/2));
       
    int j = i+1 - N*static_cast<int>(floor((i+1)/N + 0.5));
    int cnt = 0;
    while(cnt < mx){
      
      for(int a = 0; a < 3; ++a)
	rj[a] = part.getPosition(j, a);
      float qj;
      (i%2 == 0) ?  qj = q : qj = -q;
      
      float rij[3];
      for(int a = 0; a < 3; ++a){
	rij[a] = ri[a] - rj[a];
	rij[a] -= L * floor(rij[a] / L + 0.5);
      }
      float RIJSQ = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
      float RIJ = sqrt(RIJSQ);

      if(RIJ < MX_OVRLP * sigma){ //check for overlaps
	printf("Overlap between particles %d and %d.\n",i,j);
	return false;
      }

      //Lennard-Jones shifted potential, virial and force----------------------
      float SR2s = sig_sq / RIJSQ;
      float SR6s = SR2s * SR2s * SR2s;
      float pot = 4 * SR6s * (SR6s - 1.);
      float UIJ = 0.;
      float WIJ = 0.;
      if(RIJ < rsh){
	UIJ = pot + 1.;//shifted potential (epsilon = 1)
	WIJ = 6 * (pot + 4*SR6s*SR6s);
      }
      Usr += UIJ;
      Wsr += WIJ;

      //short range force
      float FIJsr = WIJ / RIJSQ;
      for(int a = 0; a < 3; ++a){
	f[a] = FIJsr * rij[a];
	part.addToForce(i, a, f[a]);
      }
      for(int a = 0; a < 3; ++a){
	f[a] *= -1;
	part.addToForce(j, a, f[a]);
      }
      //----------------------------------------------------------------------

      float qq = qi * qj;
      //Coulomb potential, virial and force-----------------------------------
      if(fabs(qq) > 1.e-12){
	UIJ = qq / RIJ;
	WIJ = qq / RIJ;
	Uc += UIJ;
	Wc += WIJ;
	
	//Coulomb force
	float FIJq = WIJ / RIJSQ;
	for(int a = 0; a < 3; ++a){
	  f[a] = FIJq * rij[a];
	  part.addToForce(i, a, f[a]);
	}
	for(int a = 0; a < 3; ++a){
	  f[a] *= -1;
	  part.addToForce(j, a, f[a]);
	}
      }
      //----------------------------------------------------------------------
      j += 1 - N*static_cast<int>(floor((j+1)/N + 0.5));
      cnt++;
    }
  }
  
  return true;
}
