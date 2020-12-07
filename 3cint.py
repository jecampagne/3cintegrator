
import Angpow

class FuncType0(Angpow.Py_ClassFunc1D):
#  dummy = 2
  def __init__(self):
    return
  def get_value(x):    
    return x * (x-1.)*(x-1.)
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
