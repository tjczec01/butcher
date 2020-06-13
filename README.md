# Butcher Tableau

##### Calculates the Butcher Tableau for a given order.

##### Calculates A, B, and C for a given order and method. Works up to order = 17.

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

display(Latex(r"$\{}{} = {}{} \\ {} \\ {}{}$".format('matrix', '{A}', '{', A[0], A[1], A[2], '}')))
display(Latex("${} \ = \ {}[ {}, \ {}, \ {} ]{}$".format( '{B}', '{', B[0], B[1], B[2], '}'))) 
display(Latex("${} \ = \ {}[ {}, \ {}, \ {} ]{}$".format( '{C}', '{', C[0], C[1], C[2], '}'))) 
display(Latex(r"$ \{}{}^{} = {}{} \\ {} \\ {}{}$".format('matrix', '{A}', '{-1}', '{', Ainv[0], Ainv[1], Ainv[2], '}')))
display(Latex(r"$\{}{} = {}{} \\ {} \\ {}{}$".format('matrix', '{T}', '{', T[0], T[1], T[2], '}')))
display(Latex(r"$ \{}{}^{} = {}{} \\ {} \\ {}{}$".format('matrix', '{T}', '{-1}', '{', TI[0], TI[1], TI[2], '}')))
display(Latex(r"$\{}{} = {}{} \\ {} \\ {}{}$".format('matrix', '{P}', '{', P[0], P[1], P[2], '}')))
              
```








# Results


\begin{align}
\nabla \cdot \vec{\mathbf{E}} & = 4 \pi \rho \\
\nabla \times \vec{\mathbf{E}}\, +\, \frac1c\, \frac{\partial\vec{\mathbf{B}}}{\partial t} & = \vec{\mathbf{0}} \\
\nabla \cdot \vec{\mathbf{B}} & = 0
\end{align}

[equation][A = \matrix{[0.19681547722366058,  -0.06553542585019845,  0.023770974348220134], [0.39442431473908746,  0.29207341166522777,  -0.04154875212599763], [0.376403062700466,  0.512485826188424,  0.111111111111110]}]

```math
$\matrix{B} = \ {[0.376403062700466, \ 0.512485826188424, \ 0.111111111111110]}$


$\matrix{C} = \ {[0.15505102572168228, \ 0.6449489742783174, \ 1.0]}$


$ \matrix{A}^{-1} = \ {[3.22474487139158, \ 1.16784008469041, \ -0.253197264742181] \\ [-3.56784008469039, \ 0.775255128608406, \ 1.05319726474218] \\ [5.53197264742194, \ -7.53197264742187, \ 5.00000000000001]}$


$\matrix{T} = \ {[0.0944387624889745, \ -0.141255295020953, \ 0.0300291941051473] \\ [0.250213122965334, \ 0.204129352293798, \ -0.38294211275726] \\ [1.0, \ 1.0, \ 0.0]}$


$ \matrix{T}^{-1} = \ {[4.178718591551935, \ 0.32768282076106514, \ 0.5233764454994487] \\ [-4.178718591551935, \ -0.3276828207610649, \ 0.476623554500551] \\ [0.5028726349458223, \ -2.571926949855616, \ 0.5960392048282263]}$


$\matrix{P} = \ {[10.048809399827414, \ -25.62959144707665, \ 15.580782047249254] \\ [-1.38214273316075, \ 10.29625811374331, \ -8.914115380582556] \\ [0.3333333333333328, \ -2.6666666666666616, \ 3.333333333333328]}$

```

