# butcher
Calculates the Butcher Tableau for a given order.

Calculates A, B, and C for a given order and method. Works up to order = 17.

This class calculates the transformation matrix, T, and its respective inverse. (However small errors may exist in the method used)

It can also calculate the interpolator coefficients, P, for a cubic polynomial spline.

# Example

```
order = 5
X = butcher(order, 15)
A, B, C = X.radau() 
Ainv = X.inv(A)        
T, TI = X.Tmat(Ainv)  
P = X.P(C)

X.printm(A)
X.printm(Ainv)
X.printm(T)
X.printm(TI)
X.printm(P)

```
