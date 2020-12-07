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

class GenericFunction {
public:
  GenericFunction();
  virtual ~GenericFunction();
};

class ClassFunc1D : public GenericFunction {
public:
  virtual double operator()(double x) const = 0;
  virtual void getValues(const std::vector<double>& vin, std::vector<double>& vout) const;
};

class Function1D
//: public ClassFunc1D
{
public:
  Function1D(double (*g)(double));
  virtual double operator()(double x) const;
};

//class ClassFunc2D : public GenericFunction {
//public:
//  virtual double operator()(double x, double y) const =0;
//};

class Function2D
//: public ClassFunc2D
{
public:
  Function2D(double (*g)(double, double));
  virtual double operator()(double x, double y) const;
};

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
  virtual r_8 operator()(r_8 x) const;

  // Vectorized evaluations
  //virtual void getValues(const std::vector<double>& vin,std::vector<double>& vout) const;

  //void ChebyshevSampling(r_8 a, r_8 b);
  void ChebyshevCoeffFFT();
  void ChebyshevTransform(r_8 a, r_8 b);
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
  r_8 ComputeIntegralUnscaled();
};

//////////////////////////////
//////////////////////////////
//////////////////////////////

}//end namespace Andpow
