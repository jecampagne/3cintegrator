
import Angpow

class FuncType0(Angpow.ClassFunc1D):
  dummy = 2
  def __init__(self):
    dummy = 2
    
#public Angpow.ClassFunc1D
#  inline virtual r_8 operator()(r_8 x) const {
#    return x * (x-1.)*(x-1.);
#  }
#private:
#  r_8 p_;
#  r_8 scale_; 
#};

def test0():
  f0 = FuncType0()
  
  iOrd0 = 2
  
  farr = Angpow.std_vector_CheFunc()
  #che_func_0 = Angpow.CheFunc(f0,iOrd0)
  #farr.push_back(Angpow.CheFunc(f0,iOrd0))
  
#  CheAlgo cheAlgo(farr);
  print('End test0......');

test0()
