# Butcher Tableau

##### Calculates the Butcher Tableau for a given order.

##### Calculates A, B, and C for a given order and method. 

##### This class calculates the transformation matrix, T, and its respective inverse. (However small errors may exist in the method used)

##### It can also calculate the interpolator coefficients, P, for a cubic polynomial spline.

# Example



```python
order = 5
X = butcher(order, 15)
A, B, C = X.radau() 
Ainv = X.inv(A)        
T, TI = X.Tmat(Ainv)  
P = X.P(C)
              
```




