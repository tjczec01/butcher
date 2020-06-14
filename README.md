# Butcher Tableau

##### Calculates the [Butcher Tableau](https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Radau_IIA_methods) for a given order.

##### Calculates A, B, and C for a given order and method. 

##### This class calculates the transformation matrix, T, and its respective inverse. (However small errors may exist in the method used)

##### It can also calculate the interpolator coefficients, P, for a cubic polynomial spline.

###### [butchertableau.py](https://github.com/tjczec01/butcher/blob/master/butchertableau.py)

## [Example](https://github.com/tjczec01/butcher/blob/master/butchertableauf.ipynb) 

## Code

```python
order = 5
X = butcher(order, 15)
A, B, C = X.radau() 
Ainv = X.inv(A)        
T, TI = X.Tmat(Ainv)  
P = X.P(C)
              
```




