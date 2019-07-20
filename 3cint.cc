#include <cstdio>
#include <stdlib.h>
#include <iostream>

#include "3CInt/angpow_exceptions.h"
#include "3CInt/angpow_numbers.h"
#include "3CInt/angpow_func.h"
#include "3CInt/angpow_tools.h"
#include "3CInt/angpow_chebyshevInt.h"
#include "3CInt/walltimer.h"


namespace Angpow {

struct PARAM {
  int Lmax_for_xmin;
  r_8 jl_xmin_cut;
  int Lmax;
  r_8 kMin;
  r_8 kMax;
  r_8 R1;
  r_8 R2;
  int chebyshev_order_1;
  int chebyshev_order_2;
  int n_sub_intervals;
} para ;


class FuncType1: public ClassFunc1D {

public:
  FuncType1(int ell, r_8 R): ell_(ell), R_(R) {}
  inline virtual r_8 operator()(r_8 x) const {
    return cos(x*R_ -ell_*M_PI*0.5 - M_PI*0.25);
  }
  virtual ~FuncType1() {}
private:
  int ell_;
  r_8 R_; 

};//ProdJBess





void test0() {
  r_8 kMin = para.kmin;
  r_8 kMax = para.kmax;
  int Lmax = para.Lmax; //ell<Lmax
  int nSubInterv = para.n_sub_intervals;


  std::vector<r_8> R(2);
  R[0] = para.R1;
  R[1] = para.R2;
  //Total Nber of Radius
  int nbreOfRadius = R.size();

  r_8 Rcur = std::accumulate(R.begin(), R.end(), 0.)/((r_8)R.size());

  r_8 kscale = 1./Rcur;

  
  //Start here for a specific ell
  int ell = 20;
  std::cout << "ell = " << ell << std::endl;

  //Function to be integrated
  FuncType1* f1 = new FuncType1(ell,R[0]);
  FuncType1* f2 = new FuncType1(ell,R[1]);

  //final k-integral bounds
  std::vector<r_8> klp(nSubInterv);
  r_8 dK = kMax-kMin;
  for(int i=0; i<=nSubInterv; i++){
    klp[i] = kMin + dK * i/((r_8)nSubInterv);
  }


  std::cout << "Dump klp values:\n";
  for(int i=0;i<klp.size();i++){
    std::cout << klp[i] << " "; 
  }
  std::cout << std::endl;

  
  printf("ell=%d, Rcur=%f, Nintervales=%d, Nradius=%d\n",
	 ell,Rcur,nSubInterv,nbreOfRadius);



  // Chebyshev machinery 
  int iOrd1 = para.chebyshev_order_1;
  int iOrd2 = para.chebyshev_order_2;
  ChebyshevIntNonSym* cheInt =  new ChebyshevIntNonSym(iOrd1, iOrd2);


  //Integration
  r_8 integral = 0.;

  for(int p = 1; p<nSubInterv; p++){ //init at p=1

    //get the bounds
    r_8 lowBound = klp[p-1];
    r_8 uppBound = klp[p];

    std::cout << "current interval: [" << lowBound << ", " << uppBound << "]" 
	      << std::endl;
    
    if(lowBound > uppBound)
      throw AngpowError("KIntegrator::Compute uppBound < lowBound Fatal");
    
    //get the Chebyshev coeffs in the final space for each fonctions
    std::vector<r_8> ChebTrans1;
    cheInt->ChebyshevTransform(f1,0,ChebTrans1,lowBound,uppBound);
    std::vector<r_8> ChebTrans2;
    cheInt->ChebyshevTransform(f2,0,ChebTrans2,lowBound,uppBound);
    
    //compute local integral
    integral += 
      cheInt->ComputeIntegral(ChebTrans1,ChebTrans2,lowBound,uppBound);
  }//p-loop 

  std::cout << "Approx. Integ = " << integral << std::endl;

  //truth
  r_8 true_int = 0.;
  if(R[0] != R[1]){
    r_8 Rdiff = R[0]-R[1];
    r_8 Rsum  = R[0]+R[1];
    true_int = (Rdiff*(cos(ell*M_PI-kmin*Rsum) - cos(ell*M_PI-kMax*Rsum)) 
		- Rsum*( sin(kmin*Rdiff)-sin(kMax*Rdiff)) )/(2.*Rdiff*Rsum);
  }
  std::cout << "True Integ = " << true_int 
	    << " diff= " << true_int - integral
	    << std::endl;



  //-------
  // clean
  //-------
  // Bessel stuff
  delete jlPtr;
  delete jliniPtr;
  // Che
  delete cheInt;
  // func
  delete f1;
  delete f2;

    

  std::cout << "End test0......" << std::endl;

}

}//namespace

//----------------------------------------------
//               Main 
//----------------------------------------------

int main(int narg, char *arg[]) {
  
  using namespace Angpow;

  //  unsigned int maxmemsize = getMemorySize()/1e6;
  //std::cout << "Max Memory size: " <<  maxmemsize << " MBytes" <<  std::endl;

  //The cosmological distance tool
   
  int test=0;
  int Lmax = 1000;
  r_8 R1 = 3350.; //Mpc z=1.0
  r_8 R2 = 3592.; //Mpc z=1.1
  r_8 kMin = 0.;
  r_8 kMax = 1.0; //Mpc^(-1)

  int chebyshev_order_1 = 9;
  int chebyshev_order_2 =  chebyshev_order_1;
  int n_sub_intervals = 10;

  int ka=1;
  while (ka<narg) {
    if (strcmp(arg[ka],"-h")==0) {
      return 0;
    }
    else if (strcmp(arg[ka],"-lmax")==0) {
      Lmax=atoi(arg[ka+1]);
      ka+=2;
    }
    else if (strcmp(arg[ka],"-kmin")==0) {
      kMin=atof(arg[ka+1]);
      ka+=2;
    }
    else if (strcmp(arg[ka],"-kmax")==0) {
      kMax=atof(arg[ka+1]);
      ka+=2;
    }
    else if (strcmp(arg[ka],"-r1")==0) {
      R1 =  atof(arg[ka+1]);
      ka+=2;      
    }    
    else if (strcmp(arg[ka],"-r2")==0) {
      R2 =  atof(arg[ka+1]);
      ka+=2;      
    }    
    else if (strcmp(arg[ka],"-che1")==0) {
      chebyshev_order_1 =  atof(arg[ka+1]);
      ka+=2;      
    }    
    else if (strcmp(arg[ka],"-che2")==0) {
      chebyshev_order_2 =  atof(arg[ka+1]);
      ka+=2;      
    }    
    else if (strcmp(arg[ka],"-nroot")==0) {
      n_bessel_roots_per_interval  = atoi(arg[ka+1]);
      ka+=2;      
    }    
    else ka++;
  }//eo while


  para.Lmax = Lmax;
  para.kMin = kMin;
  para.kMax = kMax;
  para.R1 = R1;
  para.R2 = R2;

  para.chebyshev_order_1  = chebyshev_order_1;
  para.chebyshev_order_2  = chebyshev_order_2;
  para.n_sub_intervals = n_sub_intervals;
    

  std::cout << "Configuration parameters are set to: " << std::endl;

  int rc=0;
  try {
    
    switch(test) {
    case 0:
      test0();
      break;
    default:
      throw AngpowError("EROOR Test type unkwown");
    }//end of switch


    std::cout << "---/ Fin bloc try ---- " << std::endl;
  }
    
  catch (AngpowError & e) {
    std::cerr << " besssht.cc: Catched Exception (AngpowError)" << (std::string)typeid(e).name() 
	 << " - Msg= " << e.what() << std::endl;
    rc = 99;
  }
  catch (std::exception & e) {
    std::cerr << " besssht.cc: Catched std::xception "  
	 << " - what()= " << e.what() << std::endl;
    rc = 98;
  }
  catch (...) {
    std::cerr << " besssht.cc: some other exception (...) was caught ! " << std::endl;
    rc = 97;
  }
  std::cout << " ---- Programme besssht.cc -  FIN  (Rc=" << rc << ") --- " << std::endl;
  return rc;
}//main