#include <cstdio>
#include <stdlib.h>
#include <iostream>

#include "3CInt/angpow_exceptions.h"
#include "3CInt/angpow_utils.h" //getMemorySize
#include "3CInt/angpow_numbers.h"
#include "3CInt/angpow_func.h"
#include "3CInt/angpow_bessel.h"
#include "3CInt/angpow_tools.h"
#include "3CInt/angpow_chebyshevInt.h"
#include "3CInt/walltimer.h"

class ProdJBess: public ClassFunc1D {

public:
  ProdJBess(BesselJImp* jl, int el, r_8 R1, r_8 R2):jl_(jl) {
    j1_ = new JBess1(jl,el,R1);
    j2_ = new JBess1(jl,el,R2);    
  }
  virtual ~ProdJBess() {}
private:
  BesselJImp* jl_; //not owner
  JBess1* j1_;      //owner
  JBess1* j2_;      //owner

};//ProdJBess


void test1() {

  r_8 jl_xmin_cut   = para.jl_xmin_cut;
  int Lmax_for_xmin = para.Lmax_for_xmin;  
  if(Lmax_for_xmin<Lmax){
    printf("Process: Warning: auto-reset Lmax_for_xmin > Lmax\n");
    Lmax_for_xmin = Lmax+1;
  }
  //Initialize some stuff for the bessel functions
  BesselJInit* jliniPtr = new BesselJInit(Lmax_for_xmin, jl_xmin_cut);
  BesselJImp*  jlPtr    = new BesselJImp(jliniPtr);


}
