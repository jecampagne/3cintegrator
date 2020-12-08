
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
  n_sub_intervals = 5

  # k-integral bounds
  kMin = 0.
  kMax = 1.0 #Mpc^(-1)
  klp = Angpow.std_vector_double(n_sub_intervals+1)
  dK = kMax-kMin
  for i in range(0,n_sub_intervals+1):
    klp[i] = kMin + dK * i/n_sub_intervals
    print "klp[",i,"]",klp[i]
  print "ell=",ell," , Nintervales=",n_sub_intervals

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
  print "debug : add CheFunc : 000"
  che_fun_0 = Angpow.CheFunc(f1, iOrd1)
  farr.push_back(che_fun_0)

  print "debug : add CheFunc : 001"
  che_fun_1 = Angpow.CheFunc(f2, iOrd2)
  farr.push_back(che_fun_1)

  print "debug : add CheFunc : 002"
  che_fun_2 = Angpow.CheFunc(f0, iOrd0)
  farr.push_back(che_fun_2)
  print "debug : add CheFunc : 003"
  
  # Initialisation of the Clenshow-Curtis quadrature
  cheAlgo = Angpow.CheAlgo(farr)

  # Integration
  integral = 0.

  for p in range(1,n_sub_intervals+1):
    # get the bounds
    lowBound = klp[p-1]
    uppBound = klp[p]
    
    if lowBound > uppBound:
      print 'KIntegrator::Compute uppBound < lowBound Fatal'
      return
    
    # Loop on each function to compute their  Foward Chebyshev coefficents
    print 'debug : test0 : 000-0 ',farr.size()
    for i in range(0,farr.size()):
      print 'debug : test0 : 000-1 ',i
      farr[i].ChebyshevTransform(lowBound, uppBound)
    print 'debug : test0 : 001'

    # Compute the sampling of all the functions in the final space dimension
    cheAlgo.InverseChebyshevTransform()

    # Compute the integral thanks to CC quadrature and the function sampling 
    integral += (uppBound - lowBound) * cheAlgo.ComputeIntegralUnscaled()
    #i_unscaled = cheAlgo.ComputeIntegralUnscaled()
    #print i_unscaled
    #integral += i_unscaled
    
  print "Approx. Integ = ",integral

  print 'End test0......'

test0()
