#include "3CInt/angpow_chebyshevInt.h"
#include <iostream>
#include <sstream>

#ifdef PROFILING
#include "3CInt/walltimer.h" //profiling 
#endif

namespace Angpow {

void ChebyshevIntBase::Init(int nOrdProd) {
  
  nOrdProd_ = nOrdProd; //may be not necessary to store it in fact JEC 14/2/17

  vecDCTInv_.resize(nOrdProd_,0.);
  planInv_ = new FFTPlanning(nOrdProd_,vecDCTInv_);
  
  //Clenshow-Curtis single quadrature weight
  wCC_.resize(nOrdProd_,0);
  planCC_ = new FFTPlanning(nOrdProd_, wCC_);

  ClenshawCurtisWeightsFast();

}//ChebyshevIntBase::Init

  
r_8 ChebyshevIntBase::ComputeIntegral(std::vector<r_8>& v1, std::vector<r_8>& v2, r_8 lowBnd, r_8 uppBnd){

  int nelem = v1.size();

  if(v1.size() != v2.size()){
    std::stringstream ss1; ss1<<v1.size(); 
    std::stringstream ss2; ss2<<v2.size();
    std::string msg = "ChebyshevInt::ComputeIntegral  size error (1): " + ss1.str() + ", "+ss2.str();
    throw AngpowError(msg);
  }

  if(nelem != (int)wCC_.size()){
    std::stringstream ss1; ss1<<nelem; 
    std::stringstream ss2; ss2<<wCC_.size();
    std::string msg = "ChebyshevInt::ComputeIntegral  size error (2): " + ss1.str() + ", "+ss2.str();
    throw AngpowError(msg);
  }
 

  std::vector<r_8>invCoefProd(v1);
  for (size_t i=0;i<invCoefProd.size();i++) invCoefProd[i]*=v2[i];

  
  //Compute integral

  r_8 integral = inner_product(invCoefProd.begin(),invCoefProd.end(),wCC_.begin(),0.);

  integral *= (uppBnd - lowBnd)*0.5; //the 2 division comme from CC quadrature computation throw FFT
  
  return integral;
  
}//ChebyshevIntBase::ComputeIntegral



void ChebyshevIntSym::Init() {

  //plan for single function
  nOrdFunc_=pow(2.,(r_8)ordFunc_)+1;
  vecDCTFunc_.resize(nOrdFunc_,0.);
  planFunc_ = new FFTPlanning(nOrdFunc_,vecDCTFunc_);
  
 //Unique Product plan
  nOrdProd_= 2*pow(2.,(r_8)ordFunc_) +1 ;

  ChebyshevIntBase::Init(nOrdProd_);

}//ChebyshevIntSym::Init




void ChebyshevIntNonSym::Init() {

  //plan for 1st function

  nOrdFunc1_=pow(2.,(r_8)ordFunc1_)+1;
  vecDCTFunc1_.resize(nOrdFunc1_,0.);
  planFunc1_ = new FFTPlanning(nOrdFunc1_,vecDCTFunc1_);
  
  //plan for 2nd function 
  nOrdFunc2_ = pow(2.,(r_8)ordFunc2_)+1;
  vecDCTFunc2_.resize(nOrdFunc2_,0.);
  planFunc2_ = new FFTPlanning(nOrdFunc2_,vecDCTFunc2_);
  
 //Unique Product plan
  nOrdProd_= pow(2.,(r_8)ordFunc1_)+ pow(2.,(r_8)ordFunc2_) +1;

  ChebyshevIntBase::Init(nOrdProd_);

}//ChebyshevIntSym::Init









}//namespace



