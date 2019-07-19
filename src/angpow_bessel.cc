#include "Angpow/angpow_bessel.h"
#include "Angpow/angpow_exceptions.h"
#include "Angpow/angpow_parameters.h"
#include <boost/exception/diagnostic_information.hpp> 
#include <boost/exception_ptr.hpp> 
#include <boost/math/policies/policy.hpp> //test
#include <boost/mpl/bool.hpp> //test devra etre enleve


namespace Angpow {



void BesselRoot::CmpRoots() {
  typedef boost::math::policies::policy<boost::math::policies::digits10<5> > pol;

    try{
      for (int lcur =0; lcur<lmax_; lcur++){
	//j(l,x) ~ J(l+1/2,x)/sqrt(x). 
	r_4 order = (r_4)(lcur + 0.5);
	std::vector<r_8> zeros(nroots_);
	int keepFirstRoots = 0; //was 2 JEC 9/10/16
	for(int i=0; i<nroots_; i++) {
	  int rIdx = (i<keepFirstRoots) 
	    ? i 
	    : keepFirstRoots-1+(i-keepFirstRoots+1)*step_;
	  // zeros[i] = boost::math::cyl_bessel_j_zero(order,rIdx+1,pol()); //1-based index of zero
	  //JEC 20mai16: use the approximated roots given by boost inside the
	  //boost::math::cyl_bessel_j_zero. It means do not use newton_raphson_iterate process
	  zeros[i] = boost::math::detail::bessel_zero::cyl_bessel_j_zero_detail::initial_guess(order, rIdx+1, pol());
	}
	std::copy(zeros.begin(),zeros.end(),qln_.begin() + lcur*nroots_);
      }
    }
    catch (boost::exception& ex) {
      // error handling
      std::cerr << boost::diagnostic_information(ex);
    }

}// BesselRoot::CmpRoots

  
}//namespace 

