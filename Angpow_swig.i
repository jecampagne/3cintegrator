%module(directors="1") Angpow

class wallTimer {
public:
  wallTimer();
  void start();
  void stop();
};

//////////////////////////////
/// angpow_func.h : //////////
//////////////////////////////

namespace Angpow {

class ClassFunc1D {
public:  
  ClassFunc1D();
  virtual ~ClassFunc1D();
  //public:
  //  virtual double operator()(double x) const = 0;
};

%feature("director") Py_ClassFunc1D;
class Py_ClassFunc1D : public ClassFunc1D {
public:  
  Py_ClassFunc1D();
  virtual ~Py_ClassFunc1D();
public:
  virtual double get_value(double x) const;
};
 
//////////////////////////////
/// 3cint_chefunc.h : ////////
//////////////////////////////

class CheFunc {
public:
  CheFunc(ClassFunc1D* func, int ordFunc);
  virtual ~CheFunc();
public:  
  int orderFunc();
  void ChebyshevCoeffFFT();
//void ChebyshevTransform(Angpow::r_8 a, Angpow::r_8 b);
  void ChebyshevTransform(double a, double b);
};

//////////////////////////////
/// 3cint_chealgo.h : ////////
//////////////////////////////

class CheAlgo {
public:
  CheAlgo(std::vector<CheFunc*> farr);
  virtual ~CheAlgo();
public:
  void ClenshawCurtisWeightsFast();
  void InverseChebyshevCoeffFFT();
  void InverseChebyshevTransform();
  //Angpow::r_8 ComputeIntegralUnscaled();
  double ComputeIntegralUnscaled();
};

//////////////////////////////
//////////////////////////////
//////////////////////////////

}//end namespace Andpow

%include std_vector.i
%template(std_vector_CheFunc)  std::vector<Angpow::CheFunc*>;
%template(std_vector_double)  std::vector<double>;
