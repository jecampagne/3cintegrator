#ifndef ANGPOW_BESSEL_SEEN
#define ANGPOW_BESSEL_SEEN
/*
 *  This file is part of Angpow.
 *
 *  Angpow is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Angpow is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libsharp; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  Angpow is being developed at the Linear Accelerateur Laboratory (LAL)
 *  91898 ORSAY CEDEX - FRANCE
 *  main author: J.E Campagne
 *    co-authors: J. Neveu, S. Plaszczynski
 *
 * fucntion bessel_j(el,x,&val) is adapted from CLASSgal code (http://cosmology.unige.ch/content/classgal)
 * which in turn uses:
 * Asymptotic approximations 8.11.5, 8.12.5, and 8.42.7 from
 * G.N.Watson, A Treatise on the Theory of Bessel Functions,
 * 2nd Edition (Cambridge University Press, 1944).
 * Higher terms in expansion for x near l given by
 * Airey in Phil. Mag. 31, 520 (1916).

 */

//STD
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <iterator>
#include <math.h>
#include <algorithm>




#include "angpow_numbers.h"
#include "angpow_func.h"
#include "angpow_exceptions.h"

#define _GAMMA1_ 2.6789385347077476336556
#define _GAMMA2_ 1.3541179394264004169452

/*
POLICY = 2 => CLASSgal bessel_j function
       = 1 => BOOST + policy reduite
       = 0 => BOOST call avec r_4 arguments
*/
#define POLICY 2


//Boost: see how to cut this dependancy (see zeros computation)
#include "boost/math/special_functions/bessel.hpp"



namespace Angpow {






class BesselRoot {
   public:
  //Ctor
  BesselRoot(int lmax, int nmax, int step): lmax_(lmax),nmax_(lmax), step_(step) { 
    Init(lmax,nmax,step); 
  }

  //Dtor
  virtual ~BesselRoot() {}

  //Initializor
  void Init(int lmax, int nmax, int step) {
    
    lmax_= lmax;  
    nmax_ = nmax;
    step_ = step;
    nroots_ = nmax/step+2;

    //Set proper dimension of qln zeros matrix
    qln_.resize(lmax_*nroots_);

    //compute the Bessel roots
    CmpRoots();

  }//end Init


  //Helper
  r_8 operator()(int l, int n) const { 
    return qln_[l*nroots_+n];
  } //qln(l,n)

  void GetVecRoots(std::vector<r_8>& roots, int l) const {
    roots.resize(nroots_);
    std::copy(qln_.begin()+l*nroots_,qln_.begin()+(l+1)*nroots_,roots.begin());
  }

  int NRootsL() const { return nroots_;}
  int NRoots()  const { return qln_.size();}

  //Find Bessel root (qln) values such that j_l(qln) = 0 for l=O,...,Lmax and nroots from 1, nmax
  void CmpRoots();
 
 protected:

  int lmax_;  //!< l=0,...,lmax-1  order of Spherical Bessel indice
  int nmax_;  //!< n=0,...,nmax-1  index of the roots of Spherical Bessel func
  int step_;
  int nroots_; 
  
  std::vector<r_8> qln_;  //!< roots: Care qln(l,n) l=0,...,lmax-1; n=0,...,nroots_

private:
  
  /*
    Explicit prohibit Copy Ctor
  */
  BesselRoot(const BesselRoot& );
  
  /*
    Explicit prohibit  Assignment Operator
  */
  BesselRoot& operator=(const BesselRoot&);

};//BesselRoot


// j_l(x) hub
inline void bessel_j(int l, double x, double* jl) {
    
  double nu,nu2,beta,beta2;
  double x2,sx,sx2,cx;
  double cotb,cot3b,cot6b,secb,sec2b;
  double trigarg,expterm,fl;
  double l3,cosb;
       
  fl = (double)l;
    
  x2 = x*x;
    
  /************* Use closed form for l<7 **********/
    
  if (l < 7) {
      
    sx=sin(x);
    cx=cos(x);
      
    if(l == 0) {
      if (x > 0.1) *jl=sx/x;
      else *jl=1.-x2/6.*(1.-x2/20.);
      return;
    }
      
    if(l == 1) {
      if (x > 0.2) *jl=(sx/x -cx)/x;
      else *jl=x/3.*(1.-x2/10.*(1.-x2/28.));
      return;
    }

    if (l == 2) {
      if (x > 0.3) *jl=(-3.*cx/x-sx*(1.-3./x2))/x;
      else *jl=x2/15.*(1.-x2/14.*(1.-x2/36.));
      return;
    }

    if (l == 3) {
      if (x > 0.4) *jl=(cx*(1.-15./x2)-sx*(6.-15./x2)/x)/x;
      else *jl=x*x2/105.*(1.-x2/18.*(1.-x2/44.));
      return;
    }

    if (l == 4) {
      if (x > 0.6) *jl=(sx*(1.-45./x2+105./x2/x2) +cx*(10.-105./x2)/x)/x;
      else *jl=x2*x2/945.*(1.-x2/22.*(1.-x2/52.));
      return;
    }
    
    if (l == 5) {
      if (x > 1) *jl=(sx*(15.-420./x2+945./x2/x2)/x -cx*(1.0-105./x2+945./x2/x2))/x;
      else *jl=x2*x2*x/10395.*(1.-x2/26.*(1.-x2/60.));
      return;
    }

    if (l == 6) {
      if (x > 1) *jl=(sx*(-1.+(210.-(4725.-10395./x2)/x2)/x2)+
		      cx*(-21.+(1260.-10395./x2)/x2)/x)/x;
      else *jl=x2*x2*x2/135135.*(1.-x2/30.*(1.-x2/68.));
      return;
    }

  }

  else {

    if (x <= 1.e-40) {
      *jl=0.0;
      return;
    }

    nu= fl + 0.5;
    nu2=nu*nu;

    if ((x2/fl) < 0.5) {
      *jl=exp(fl*log(x/nu/2.)+nu*(1-log(2.))-(1.-(1.-3.5/nu2)/nu2/30.)/12./nu)
	/nu*(1.-x2/(4.*nu+4.)*(1.-x2/(8.*nu+16.)*(1.-x2/(12.*nu+36.))));
      return;
    }

    if ((fl*fl/x) < 0.5) {

      beta = x - M_PI/2.*(fl+1.);
      *jl = (cos(beta)*(1.-(nu2-0.25)*(nu2-2.25)/8./x2*(1.-(nu2-6.25)*(nu2-12.25)/48./x2))
	     -sin(beta)*(nu2-0.25)/2./x* (1.-(nu2-2.25)*(nu2-6.25)/24./x2*(1.-(nu2-12.25)*(nu2-20.25)/80./x2)) )/x;
      
      return;

    }

    l3 = pow(nu,0.325);

    if (x < nu-1.31*l3) {

      cosb=nu/x;
      sx=sqrt(nu2-x2);
      cotb=nu/sx;
      secb=x/nu;
      beta=log(cosb+sx/x);
      cot3b=cotb*cotb*cotb;
      cot6b=cot3b*cot3b;
      sec2b=secb*secb;
      expterm=((2.+3.*sec2b)*cot3b/24.
	       - ((4.+sec2b)*sec2b*cot6b/16.
		  + ((16.-(1512.+(3654.+375.*sec2b)*sec2b)*sec2b)*cot3b/5760.
		     + (32.+(288.+(232.+13.*sec2b)*sec2b)*sec2b)*sec2b*cot6b/128./nu)*cot6b/nu)/nu)/nu;
      *jl=sqrt(cotb*cosb)/(2.*nu)*exp(-nu*beta+nu/cotb-expterm);

      return;

    }

    if (x > nu+1.48*l3) {

      cosb=nu/x;
      sx=sqrt(x2-nu2);
      cotb=nu/sx;
      secb=x/nu;
      beta=acos(cosb);
      cot3b=cotb*cotb*cotb;
      cot6b=cot3b*cot3b;
      sec2b=secb*secb;
      trigarg=nu/cotb-nu*beta-M_PI/4.
	-((2.+3.*sec2b)*cot3b/24.
	  +(16.-(1512.+(3654.+375.*sec2b)*sec2b)*sec2b)*cot3b*cot6b/5760./nu2)/nu;
      expterm=((4.+sec2b)*sec2b*cot6b/16.
	       -(32.+(288.+(232.+13.*sec2b)*sec2b)*sec2b)*sec2b*cot6b*cot6b/128./nu2)/nu2;
      *jl=sqrt(cotb*cosb)/nu*exp(-expterm)*cos(trigarg);

      return;
    }
    
    /* last possible case */

    beta=x-nu;
    beta2=beta*beta;
    sx=6./x;
    sx2=sx*sx;
    secb=pow(sx,1./3.);
    sec2b=secb*secb;
    *jl=(_GAMMA1_*secb + beta*_GAMMA2_*sec2b
	 -(beta2/18.-1./45.)*beta*sx*secb*_GAMMA1_
	 -((beta2-1.)*beta2/36.+1./420.)*sx*sec2b*_GAMMA2_
	 +(((beta2/1620.-7./3240.)*beta2+1./648.)*beta2-1./8100.)*sx2*secb*_GAMMA1_
	 +(((beta2/4536.-1./810.)*beta2+19./11340.)*beta2-13./28350.)*beta*sx2*sec2b*_GAMMA2_
	 -((((beta2/349920.-1./29160.)*beta2+71./583200.)*beta2-121./874800.)*
	   beta2+7939./224532000.)*beta*sx2*sx*secb*_GAMMA1_)*sqrt(sx)/12./sqrt(M_PI);

    return;

  }
    
  std::stringstream ss1; ss1 <<l; 
  std::stringstream ss2; ss2 <<x; 
  std::string msg = "l or x out of bounds: ";
  msg += ss1.str() + ", " + ss2.str();
  throw AngpowError(msg);
  //    std::cout << "l or x out of bounds: " << l << ", " << x << std::endl;

  //    return -1;

}//bessel_j

class BesselJInit {
public:
  BesselJInit(int Lmax=10000, r_8 cut=5e-10):Lmax_(Lmax), cut_(cut) {
    Init(Lmax, cut);
  }//Ctor

  void Init(int Lmax, r_8 cut) {
    vec_.resize(Lmax);
    
    for(int el=0;el<Lmax;el++){
      if(el==0){//j_0(x)=sin(x)/x
	vec_[el]=0;
	continue;
      }
      r_8 xmin=0;
      r_8 xmax=el+0.5;
      r_8 fmax;
      bessel_j(el,xmax,&fmax);

      if(fmax < cut)
	throw AngpowError("BesselJImp::Xmin cut not appropriate. FATAL");
    
      r_8 eps=0.5e-10;
      r_8 xmiddle=(xmax+xmin)*0.5;
      r_8 fmiddle;
      while((xmax-xmin)/(xmax+xmin)> eps){
	xmiddle = (xmax+xmin)*0.5;
	bessel_j(el,xmax,&fmiddle);
	if(fmiddle>cut)
	  xmax = xmiddle;
	else
	  xmin = xmiddle;
      }//while

      vec_[el]=xmiddle;

    }//el-loop
  }//Init

public:
  std::vector<r_8> vec_;
  int Lmax_;
  r_8 cut_;
};


class BesselJImp {
public:

  BesselJImp(BesselJInit* jlinit): jlinit_(jlinit) {}
  
  
  //return x such that j_l(x)=cut in [0,l+1/2] 
  //  static r_8 Xmin(int el) { return xmin_[el]; }
  r_8 Xmin(int el) { return jlinit_->vec_[el]; }


  //Heart of the class  (avoid virtual)   
  inline r_8 operator()(int el, r_8 x) const{

    //    if(el>=(int)xmin_.size())
    if(el>=(int)jlinit_->vec_.size()) {
      throw AngpowError("BesselJImp used with el > xmin size");
    }
    //    if(x<xmin_[el])return 0;
    if(x<jlinit_->vec_[el])return 0;

    r_8 val=0;
#if POLICY == 2
    bessel_j(el,x,&val);
#elif POLICY == 1
    typedef boost::math::policies::policy<boost::math::policies::digits10<5>, boost::math::policies::promote_float<false>, boost::math::policies::promote_double<false>> pol;
    val =  boost::math::sph_bessel(el, (r_4)(x), pol()); //JEc 19mai16 flaot arg/return
#else
    val = boost::math::sph_bessel(el,x);
#endif
    return val;
  }

 protected:
  //  static std::vector<r_8> xmin_;
  //  int Lmax_;
  BesselJInit* jlinit_; //not owner

}; //BesselJImp



/*! 
  Class for Spherical Bessel sampling
  Sqrt(2/Pi) j_l(kR)
*/
class JBess1 : public ClassFunc1D {
 public:
  JBess1(BesselJImp* jl, int el, r_8 R, r_8 norm=1.0): jl_(jl),el_(el),R_(R) {
    norm_ = sqrt(2./M_PI) * norm;
  }
  inline virtual r_8 operator()(r_8 x) const {
    r_8 val = (*jl_)(el_,x*R_);
    return norm_*val;
  }
  inline r_8 Norm() const { return norm_;}

  virtual ~JBess1() {}
private:
  BesselJImp* jl_; //not owner
  int el_;
  r_8 R_;
  r_8 norm_;
};




#undef POLICY



}//end namespace

#undef _GAMMA1_
#undef _GAMMA2_


#endif //ANGPOW_BESSEL_SEEN
