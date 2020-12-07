%module Angpow

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

  //class GenericFunction {
  //public:
  //  GenericFunction();
  //  virtual ~GenericFunction();
  //};

class ClassFunc1D {
public:  
  ClassFunc1D();
  virtual ~ClassFunc1D();
  //public:
  //  virtual double operator()(double x) const = 0;
};

//class Py_ClassFunc1D : public ClassFunc1D {
//public:  
//  Py_ClassFunc1D();
//  virtual ~Py_ClassFunc1D();
//public:
//  virtual double get_value(double x) const = 0;
//};
 
//class Function1D
//: public ClassFunc1D
//{
//public:
//  Function1D(double (*g)(double));
//  virtual double operator()(double x) const;
//};

//class ClassFunc2D : public GenericFunction {
//public:
//  virtual double operator()(double x, double y) const =0;
//};

//class Function2D
//: public ClassFunc2D
//{
//public:
//  Function2D(double (*g)(double, double));
//  virtual double operator()(double x, double y) const;
//};

//////////////////////////////
/// 3cint_chefunc.h : ////////
//////////////////////////////

class CheFunc {
 public:
  CheFunc(ClassFunc1D* func, int ordFunc);
  virtual ~CheFunc();
  //getter
  int orderFunc();
  //std::vector<r_8>& vec();
  //std::vector<r_8>& vecFinal();

  //void updateVec(std::vector<r_8>& vec);

  // Single value evaluation
  virtual Angpow::r_8 operator()(Angpow::r_8 x) const;

  // Vectorized evaluations
  //virtual void getValues(const std::vector<double>& vin,std::vector<double>& vout) const;

  //void ChebyshevSampling(r_8 a, r_8 b);
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
  void ClenshawCurtisWeightsFast();
  void InverseChebyshevCoeffFFT();
  void InverseChebyshevTransform();
  Angpow::r_8 ComputeIntegralUnscaled();
};

//////////////////////////////
//////////////////////////////
//////////////////////////////

}//end namespace Andpow

%include std_vector.i
%template(std_vector_CheFunc)  std::vector<Angpow::CheFunc*>;
%template(std_vector_double)  std::vector<double>;
