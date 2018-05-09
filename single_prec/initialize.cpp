
#include "initialize.hpp"

//Test for overlaps
bool overlap_test(const Particles& part, float L, float frac){
  int N = part.getNoOfParticles();
  float rad = part.getRadius();
  
  float sigma = 2*rad;
  float srt = pow(2., 1./6);
  float sig_sq = sigma * sigma;
  float rsh = srt * sigma;
  
  float ri[3], rj[3];
  
  for(int i = 0; i < N-1; ++i){
    for(int a = 0; a < 3; ++a)
      ri[a] = part.getPosition(i, a);

    //even out the load among threads
    int mx = static_cast<int>(ceil(static_cast<float>(N-1)/2));
    if(N%2 == 0. && i >= N/2)
       mx = static_cast<int>(floor(static_cast<float>(N-1)/2));
       
    int j = i+1 - N*static_cast<int>(floor((i+1)/N + 0.5));
    int cnt = 0;
    while(cnt < mx){
      for(int a = 0; a < 3; ++a)
	rj[a] = part.getPosition(j, a);
      
      float rij[3];
      for(int a = 0; a < 3; ++a){
	rij[a] = ri[a] - rj[a];
	//Periodic boudary conditions
	rij[a] -= L * floor(rij[a] / L + 0.5);
      }
      float RIJSQ = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
      float RIJ = sqrt(RIJSQ);

      if(RIJ < frac * sigma){ //check for overlaps
	printf("Overlap between particles %d and %d.\n",i,j);
	return true;
      }
      j += 1 - N*static_cast<int>(floor((j+1)/N + 0.5));
      cnt++;
    }
  }
  return false;
}


//Initial positions
int ini_pos(Particles& part, float L, float dis){
  int N = part.getNoOfParticles();
  float rad = part.getRadius();
  
  float r[3] = {0., 0., 0.};
  float sq = sqrt(3.) / 2.;
  float sqt = 1./ (2 * sqrt(3.));
  float sqtt = sqrt(6.) / 3;

  float ssigma = 2*rad+dis;
  
  int num1 = static_cast<int>(max(floor((L-1.e-2)/(ssigma)), 1.));
  int num2 = static_cast<int>(max(floor((L-1.e-2)/(ssigma*sq)), 1.));
  int num3 = static_cast<int>(max(floor((L-1.e-2)/(ssigma*sqtt)), 1.));

  int cont = 0;
  while(cont < N){
    bool escape = false;
    for(int nz = 0; nz < num3; ++nz){
      for(int ny = 0; ny < num2; ++ny){
	for(int nx = 0; nx < num1; ++nx){
	  // //
	  // cout<<cont<<endl;
	  // //
	  if(cont < N){
	    
	    if(fmod(static_cast<float>(nz),2.) == 0){
	      if(fmod(static_cast<float>(ny),2.) == 0){
		r[0] = (ssigma) * static_cast<float>(nx);
	      }else{
		r[0] = (ssigma) * (static_cast<float>(nx) + 0.5);
	      }
	      r[1] = (ssigma) * static_cast<float>(ny) * sq;
	    }else{
	      if(fmod(static_cast<float>(ny),2.) != 0){
		r[0] = (ssigma) * static_cast<float>(nx);
	      }else{
		r[0] = (ssigma) * (static_cast<float>(nx) + 0.5);
	      }
	      r[1] = (ssigma) * (static_cast<float>(ny) * sq + sqt);
	    }
	    r[2] = (ssigma) * static_cast<float>(nz) * sqtt;
	      
	    for(int a = 0; a < 3; ++a)
	      r[a] -= L * floor(r[a] / L + 0.5);
	
	    for(int a = 0; a < 3; ++a)
	      part.setPosition(cont, a, r[a]);
	    
	    cont++;
	    
	  }else{
	    escape = true;
	  }
	  
	  // //
	  // cout<<"esc "<<escape<<endl;
	  // //
	  
	  if(escape)
	    break;
	}
	if(escape)
	  break;
      }
      if(escape)
	break;
    }
  }
 
  //print file
  ofstream InPos("iniPos.dat");
  if(! InPos){
    InPos.close();
    return -1;
  }
  InPos << "0 "; //charge
  InPos << endl;
  InPos << L << endl;
  
  for(int i = 0; i < N; i++){
    for(int a = 0; a < 3; ++a)
      r[a] = part.getPosition(i, a);
    
    InPos << r[0] << " " << r[1] << " " << r[2];
    InPos << " 0 " << "0 " << "0 " << endl; 
  }
  
  InPos.close();

  //test for overlaps
  if(overlap_test(part, L, 1.)){
    return 1;
  }
  
  return 0;
}



void ini_mom(Particles& part, float temp, const gsl_rng* ran){
  int N = part.getNoOfParticles();
  
  //for translation velocities************************************************
  vector<vector<float> > rnd(3,vector<float>(1));
  for(int a = 0; a < 3; ++a)
    rnd[a].resize(N+1);
  
  for(int i = 0; i < N; ++i){
    for(int a = 0; a < 3; ++a)
      rnd[a][i] = gsl_rng_uniform(ran);
  }
  
  float mtemp = temp * part.getMass();
  for(int i = 0; i < N; i += 2){
    if((fmod(static_cast<float>(N),2) != 0.) && (i == N-1)){
      for(int a = 0; a < 3; a++){
	float rr = gsl_rng_uniform(ran);
	gauss(mtemp, rr, rnd[a][i]);
      }
    }else{
      for(int a = 0; a < 3; a++){
	gauss(mtemp, rnd[a][i], rnd[a][i+1]);
      }
    }
  }

  float P[3] = {0., 0., 0.};
  //float pp[3][Ntot+1];
  vector<vector<float> > pp(3, vector<float>(N+1));

  //total momentum on each direction
  for(int i = 0; i < N; ++i){
    for(int a = 0; a < 3; a++){
      pp[a][i] = rnd[a][i];
      P[a] += pp[a][i];
    }
  }

  float pr = 1.;
  
//   //
//   cout<<"begin linear X..\n";
//   //
  
 //  // X part/////////////////////////////////////////////////
  int II = 0;/*, nn=0, ii=0;*/
  long scont = 0;
  while(fabs(P[0]) > pr){ //NOTE: in fabsP > x, if the value of x is very 
                           // small the program will have an error. 
                           //(what is small depends on Ntot) 
                           // for N=256 x=0.1 is fine.
    while(P[0] > pr){
      if(pp[0][II] > 0){ 
	pp[0][II] = -pp[0][II];
	P[0] += 2 * pp[0][II];
      }
      II++;
      if(II == N - 1) II = 0;
    }
    while(P[0] < -pr){
      if(pp[0][II] < 0){ 
	pp[0][II] = -pp[0][II];
	P[0] += 2 * pp[0][II];
      }
      II++;
      if(II == N - 1) II = 0;
    }
    if(scont > 1e6) break;
    scont++;
  }
  P[0] -= pp[0][N-1];
  pp[0][N-1] = - P[0];
  P[0] = 0.;
  
  for(int i = 0; i < N; ++i){
    rnd[0][i] = pp[0][i];
    P[0] += rnd[0][i];
  }
  
  if(fabs(P[0]) > 1e-12){
    printf("Px!=0 !!\n");
  }


  //   //
//   cout<<"begin linear Y..\n";
//   //
  
 //  // Y part/////////////////////////////////////////////////
  II = 0;/*, nn=0, ii=0;*/
  scont = 0;
  while(fabs(P[1]) > pr){ //NOTE: in fabsP > x, if the value of x is very 
                           // small the program will have an error. 
                           //(what is small depends on Ntot) 
                           // for N=256 x=0.1 is fine.
    while(P[1] > pr){
      if(pp[1][II] > 0){ 
	pp[1][II] = -pp[1][II];
	P[1] += 2 * pp[1][II];
      }
      II++;
      if(II == N - 1) II = 0;
    }
    while(P[1] < -pr){
      if(pp[1][II] < 0){ 
	pp[1][II] = -pp[1][II];
	P[1] += 2 * pp[1][II];
      }
      II++;
      if(II == N - 1) II = 0;
    }
    if(scont > 1e6) break;
    scont++;
  }
  P[1] -= pp[1][N-1];
  pp[1][N-1] = - P[1];
  P[1] = 0.;
  
  for(int i = 0; i < N; ++i){
    rnd[1][i] = pp[1][i];
    P[1] += rnd[1][i];
  }
  
  if(fabs(P[1]) > 1e-12){
    printf("Py!=0 !!\n");
  }


    //   //
//   cout<<"begin linear Y..\n";
//   //
  
 //  // Y part/////////////////////////////////////////////////
  II = 0;/*, nn=0, ii=0;*/
  scont = 0;
  while(fabs(P[2]) > pr){ //NOTE: in fabsP > x, if the value of x is very 
                           // small the program will have an error. 
                           //(what is small depends on Ntot) 
                           // for N=256 x=0.1 is fine.
    while(P[2] > pr){
      if(pp[2][II] > 0){ 
	pp[2][II] = -pp[2][II];
	P[2] += 2 * pp[2][II];
      }
      II++;
      if(II == N - 1) II = 0;
    }
    while(P[2] < -pr){
      if(pp[2][II] < 0){ 
	pp[2][II] = -pp[2][II];
	P[2] += 2 * pp[2][II];
      }
      II++;
      if(II == N - 1) II = 0;
    }
    if(scont > 1e6) break;
    scont++;
  }
  P[2] -= pp[2][N-1];
  pp[2][N-1] = - P[2];
  P[2] = 0.;
  
  for(int i = 0; i < N; ++i){
    rnd[2][i] = pp[2][i];
    P[2] += rnd[2][i];
  }
  
  if(fabs(P[2]) > 1e-12){
    printf("Pz!=0 !!\n");
  }

  
  //set momenta
  for(int i = 0; i < N; i++){
    for(int a = 0; a < 3; a++)
      part.setMomentum(i, a, rnd[a][i]);
  }
}


//Box-Muller algorithm
void gauss(float variance, float& rnd1, float& rnd2){
  float rndA, rndB;
  rndA = sqrt(- 2 * log(rnd1)) * cos(2 * M_PI * rnd2);
  rndB = sqrt(- 2 * log(rnd1)) * sin(2 * M_PI * rnd2);
  float sqrtVar = sqrt(variance);
  rnd1 = rndA * sqrtVar;
  rnd2 = rndB * sqrtVar;
}


