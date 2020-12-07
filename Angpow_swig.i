%module Angpow

class wallTimer {
public:
  wallTimer();
  void start();
  void stop();
};

namespace Angpow {

class GenericFunction {
public:
  GenericFunction();
  virtual ~GenericFunction();
};

//class ClassFunc1D : public GenericFunction {
//public:
//  virtual double operator()(double x) const = 0;
//  virtual void getValues(const std::vector<double>& vin, std::vector<double>& vout) const;
//};

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

}//end namespace Andpow
