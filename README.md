Code extracted from AngPow to perform the k-integration.

The purpose of the exercice is to perform the following integration 

```math
I = \int_{k_{Min}}^{k_{kMax}} f_\ell(k,R_1) f_\ell(k, R_2) dk
```
with

```math
f_\ell(k,R) = \cos( kR - \ell \pi/2 - \pi/4) 
```

The user has to provide the function parameters, 
the  integral bounds kMin and kMax, 
the number of k-sub intervales, the chebyshev orders.

Then, on each k-sub interval, the chebyshev coefficients of 
the two functions are computed separately 
and then multiplied to compute the integral on the k-sub interval

The "cos" product function leads to exact computation that can be used
as accuracy benchmark.

Test: after cloning the code, and ```FFTW``` installed

    > make
    > ./bin/3cint

Output

    Configuration parameters are set to: 
    ell=20, Nintervales=5
    Approx. Integ = -0.0019508
    True Integ = -0.0019508 diff= 3.67111e-15
    End test0......
    ---/ Fin bloc try ---- 
    ---- Programme 3cint.cc -  FIN  (Rc=0) --- 




