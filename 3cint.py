
import Angpow

#class FuncType0(Angpow.Py_ClassFunc1D):
class FuncType0(Angpow.ClassFunc1D):
  def get_value(a_x):    
    return a_x * (a_x-1.)*(a_x-1.)

import math
  
class FuncType1(Angpow.ClassFunc1D):
  m_ell = 0
  m_R = 0
  #def __init__(self,a_ell,a_R):
  #  #Angpow.ClassFunc1D()
  #  self.m_ell = a_ell
  #  self.m_R = a_R
  def get_value(self,a_x):
    return math.cos(a_x*self.m_R - self.m_ell*math.pi*0.5 - math.pi*0.25)

def test0():
  ell = 20
  R1 = 2000.
  R2 = 2200.
  chebyshev_order_1 = 8
  chebyshev_order_2 =  chebyshev_order_1

  #f1 = FuncType1(ell,R1)
  f1 = FuncType1()
  f1.m_ell = ell
  f1.m_R = R1
  
  #f2 = FuncType1(ell,R2)
  f2 = FuncType1()
  f2.m_ell = ell
  f2.m_R = R2

  f0 = FuncType0()
  
  #print(f1.get_value(1))
  
  iOrd0 = 2
  iOrd1 = chebyshev_order_1
  iOrd2 = chebyshev_order_2
  
  farr = Angpow.std_vector_CheFunc()
  farr.push_back(Angpow.CheFunc(f1, iOrd1))
  farr.push_back(Angpow.CheFunc(f2, iOrd2))
  farr.push_back(Angpow.CheFunc(f0, iOrd0))
  
#  CheAlgo cheAlgo(farr);
  print('End test0......')

test0()
