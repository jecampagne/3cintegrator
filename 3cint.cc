#include <cstdio>
#include <stdlib.h>
#include <iostream>

#include "3CInt/angpow_exceptions.h"
//#include "3CInt/angpow_utils.h" //getMemorySize
#include "3CInt/angpow_numbers.h"
#include "3CInt/angpow_func.h"
#include "3CInt/angpow_bessel.h"
#include "3CInt/angpow_tools.h"
#include "3CInt/angpow_chebyshevInt.h"
#include "3CInt/walltimer.h"


namespace Angpow {

struct PARAM {
  int Lmax_for_xmin;
  r_8 jl_xmin_cut;
  int Lmax;
  r_8 kmax;
  r_8 R1;
  r_8 R2;
  int chebyshev_order_1;
  int chebyshev_order_2;
  int n_bessel_roots_per_interval;
} para ;

class FuncType2: public ClassFunc1D {

public:
  FuncType2(BesselJImp* jl, int ell, r_8 R) {
    jl_ = new JBess1(jl,ell,R);
  }
  inline virtual r_8 operator()(r_8 x) const {
    return (*jl_)(x) * 1000. * sqrt(x) * exp(-x*3.0);
  }
  virtual ~FuncType2() {}
private:
  JBess1* jl_;      //owner

};

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
  r_8 kMax = para.kmax;
  int Lmax = para.Lmax; //ell<Lmax
  int nRootPerInt = para.n_bessel_roots_per_interval;


  r_8 jl_xmin_cut   = para.jl_xmin_cut;
  int Lmax_for_xmin = para.Lmax_for_xmin;  
  if(Lmax_for_xmin<Lmax){
    printf("ProcessTest0: Warning: auto-reset Lmax_for_xmin > Lmax\n");
    Lmax_for_xmin = Lmax+1;
  }
  //Initialize some stuff for the bessel functions
  BesselJInit* jliniPtr = new BesselJInit(Lmax_for_xmin, jl_xmin_cut);
  BesselJImp*  jlPtr    = new BesselJImp(jliniPtr);



  //Bessel roots 
  std::vector<r_8> R(2);
  R[0] = para.R1;
  R[1] = para.R2;
  //Total Nber of Radius
  int nbreOfRadius = R.size();

  r_8 Rcur = std::accumulate(R.begin(), R.end(), 0.)/((r_8)R.size());

  r_8 kscale = 1./Rcur;
  int Pmaxmax = (kMax*Rcur)/M_PI;

  BesselRoot broots(Lmax,Pmaxmax,nRootPerInt);
  int NbessRoots = broots.NRootsL(); //last bessel root is NbessRoots*nRootPerInt_
    
  int maxRoots = NbessRoots;
  
  //Start here for a specific ell
  int ell = 20;
  std::cout << "ell = " << ell << std::endl;

  //Function to be integrated
  FuncType1* f1 = new FuncType1(ell,R[0]);
  FuncType1* f2 = new FuncType1(ell,R[1]);

  std::vector<r_8> qlp; broots.GetVecRoots(qlp, ell);	


  //get the last value = kMax
  
  r_8 kLast = qlp.back()*kscale;
  std::cout << "(1) kLast= " << kLast << std::endl;
  while(kLast>kMax && !qlp.empty()){
    qlp.pop_back();
    kLast = qlp.back()/Rcur;
  }
  if(kLast<kMax){
    qlp.push_back(Rcur*kMax);
  }
  kLast = qlp.back()*kscale;
  std::cout << "(2) kLast= " << kLast << std::endl;
  
  //final k-integral bounds
  std::vector<r_8> klp;
  klp.push_back(jlPtr->Xmin(ell));
  klp.insert(klp.end(),qlp.begin(),qlp.end());
  std::transform(klp.begin(),klp.end(),klp.begin(),std::bind1st(std::multiplies<r_8>(),kscale));


  std::cout << "Dump klp values:\n";
  for(int i=0;i<klp.size();i++){
    std::cout << klp[i] << " "; 
  }
  std::cout << std::endl;

  
  maxRoots = klp.size();
  
  printf("ell=%d, Rcur=%f, Nintervales=%d, Nradius=%d\n",
	 ell,Rcur,maxRoots-1,nbreOfRadius);



  // Chebyshev machinery 
  int iOrd1 = para.chebyshev_order_1;
  int iOrd2 = para.chebyshev_order_2;
  ChebyshevIntNonSym* cheInt =  new ChebyshevIntNonSym(iOrd1, iOrd2);


  //Integration
  r_8 integral = 0.;
  //Split the [kmin, kmax] intervalle into sub-intervales
  for(int p = 1; p<maxRoots; p++){ //init at p=1

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

  std::cout << "Integ = " << integral << std::endl;


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
  r_8 kmax = 1.0; //Mpc^(-1)

  int chebyshev_order_1 = 9;
  int chebyshev_order_2 =  chebyshev_order_1;
  int n_bessel_roots_per_interval = 100;

  int ka=1;
  while (ka<narg) {
    if (strcmp(arg[ka],"-h")==0) {
      return 0;
    }
    else if (strcmp(arg[ka],"-lmax")==0) {
      Lmax=atoi(arg[ka+1]);
      ka+=2;
    }
    else if (strcmp(arg[ka],"-kmax")==0) {
      kmax=atof(arg[ka+1]);
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


  //Bessel parameters
  para.Lmax_for_xmin = 5000;
  para.jl_xmin_cut   = 5e-10;

  para.Lmax = Lmax;
  para.kmax = kmax;
  para.R1 = R1;
  para.R2 = R2;

  para.chebyshev_order_1  = chebyshev_order_1;
  para.chebyshev_order_2  = chebyshev_order_2;
  para.n_bessel_roots_per_interval = n_bessel_roots_per_interval;
    

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
