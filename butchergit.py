# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 23:21:54 2020

@author: tjcze
"""

import numpy as np
import scipy as sc
import sympy as sp
import pprint as pp
from scipy import linalg
pprint = pp.pprint

flatten = lambda l: [item for sublist in l for item in sublist] 

def factorial(n):
     if n == 0:
            return 1
     else:
            nn = n
            ne = n - 1
            while ne >= 1:
                   nn = nn*ne
                   ne -= 1
            return nn

def Gamma(n):
       if n <= 0:
              return None
       else:
              return factorial(n-1)
class butcher:
       
       def __init__(self, order):
              self.s = int(order)
       
       def Px(self, s, evalf=True):
                     xp = sp.symbols("x")
                     if s == 0:
                            return 1.0
                     else:
                            term1 = 1.0/((2.0**s)*factorial(s))
                            term2 = sp.expand((xp**2 - 1)**s)
                            term3 = sp.expand(sp.diff(term2,xp,s))
                            term4 = sp.Mul(term1,term3)
                            return term4
       # Shifed legendre
       def PxS(self, s, evalf=True):
                     xps = sp.symbols("x")
                     if s == 0:
                            return 1.0
                     else:
                            term1 = 1.0/factorial(s)
                            term2 = sp.expand((xps**2 - xps)**s)
                            term3 = sp.expand(sp.diff(term2,xps,s))
                            term4 = sp.Mul(term1,term3)
                            return term4
                     
       def legendreP(self, s, listall=False, expand=True):
              listc = [i for i in range(s+1)]
              terms = []
              for i in listc:
                     if expand is True:
                            termm = sp.expand(self.Px(i))
                     elif expand is False:
                            termm = self.Px(i, False)
                     terms.append(termm)
              if listall is True:
                     return terms
              elif listall is False:
                     return terms[-1]
       # Shifed legendre      
       def legendrePS(self, s, listall=False, expand=True):
              listc = [i for i in range(s+1)]
              terms = []
              for i in listc:
                     if expand is True:
                            termm = sp.expand(self.PxS(i))
                     elif expand is False:
                            termm = self.PxS(i, False)
                     terms.append(termm)
              if listall is True:
                     return terms
              elif listall is False:
                     return terms[-1]    
       
       # Shifed legendre 
       def CiLPS(self):
              x = sp.symbols('x')
              eq1 = self.legendrePS(self.s,False,True)
              eq2 = self.legendrePS(self.s - 1,False,True)
              eq3 = sp.Add(eq1,sp.Mul(sp.Integer(-1),eq2))
              listeq = sp.nroots(sp.Poly(eq3,x),15)
              return listeq[:]
              
       def CiLP(self):
              x = sp.symbols('x')
              eq1 = self.legendrePS(self.s,False,True)
              eq2 = self.legendrePS(self.s-1,False,True)
              eq3 = sp.Add(eq1,sp.Mul(sp.Integer(-1),eq2))
              listeq = sp.nroots(sp.Poly(eq3,x),15)
              return listeq
       
       def Bi(self, C):
              bsyms = [sp.symbols("b{}".format(i+1)) for i in range(len(C))]
              eqs = []
              for i in range(len(C)):
                     eqi = [sp.Mul(ij**(i),j) for ij,j in zip(C,bsyms)]
                     eqif = sp.Add(sum(eqi),sp.Mul(sp.Integer(-1),1.0/(i+1)))
                     eqs.append(eqif)
              bs = sp.solve(eqs,bsyms)
              return list(bs.values())
       
       def Ai(self, B, C):
              n = self.s
              Vsyms = sp.Matrix(sp.symarray('a',(n-1,n)))
              Asyms = []
              Asymsl = []
              eqsa = []
              symsd = []
              for sy in range(n-1):
                     for sk in range(n):
                            sd = sp.symbols("a_{}_{}".format(sy,sk))
                            symsd.append(sd)
              for i in range(n):
                     asm = []
                     for j in range(n):
                            s1 = sp.symbols("a_{}_{}".format(i,j))
                            asm.append(s1)
                            Asymsl.append(s1)
                     Asyms.append(asm[:])
                     asm.clear()
              rterms = []
              lterms = []
              for l in range(0,n-1,1):
                     iv = 0
                     for i in range(1,n+1,1):
                            termeq = []
                            termsr = []
                            rv = (C[l]**(i))/(i)
                            termsr.append(rv)
                            rterms.append(rv)
                            for j in range(1,n+1,1):
                                   term1 = Asyms[l][j-1]
                                   term2 = C[j-1]**(i-1)
                                   term3 = sp.Mul(term1,term2)
                                   termeq.append(term3)
                            term4 = sp.sympify(sum(termeq))
                            lterms.append(term4)
                            term5 = sp.Add(term4,sp.Mul(sp.Integer(-1),rv))
                            eqsa.append(term5)
                            termeq.clear()
                     iv += 3
              vs = flatten(sp.matrix2numpy(Vsyms).tolist())
              Aslin = sp.linsolve(eqsa, vs)
              linans = list(Aslin.args[:][0])
              lin2 = np.array(linans, dtype=float).reshape((n-1,n))
              lin3 = lin2.tolist()
              lin3.append(B[:])
              return lin3
       
       def Bfunc(self, s, Bi,Ci, k):
              sk = int(2*k + 1)
              leftlb = []
              rightlb = []
              for j in range(1,sk+1,1):
                     leftllb = []
                     rightb = j**-1
                     rightlb.append(rightb)
                     for i in range(0,s+1,1):
                            Bval = Bi[i]
                            Cval = Ci[i]**(j-1)
                            Fvalb = Bval*Cval
                            leftllb.append(Fvalb)
                     finallb = sum(leftllb)
                     leftlb.append(finallb)
                     leftllb.clear()
              return [leftlb, rightlb]
              
       def Cfunc(self, s, Aij,Ci, k):
              leftlc = []
              rightlc = []
              for l in range(k):
                     for i in range(s):
                            leftllc = []
                            rightc = (Ci[i]**(l+1))/(l+1)
                            rightlc.append(rightc)
                            for j in range(0,s,1):
                                   Avals = Aij[i][j]
                                   Cval = Ci[j]**l
                                   Fvalc = Avals*Cval
                                   leftllc.append(Fvalc)
                            finallc = sum(leftllc)
                            leftlc.append(finallc)
                            leftllc.clear()
              return [leftlc, rightlc]
                            
       def Dfunc(self, s, Aij,Bi,Ci,k):
              leftld = []
              rightld = []
              for l in range(k):
                     for j in range(s):
                            rightd = (Bi[j]*(1.0 - (Ci[j]**(l+1))))/(l+1)
                            rightld.append(rightd)
                            leftlld = []
                            for i in range(s):
                                   Aval = Aij[i][j]
                                   Cs = Ci[i]**(l) 
                                   Bv = Bi[i]
                                   Vall = Aval*Cs*Bv
                                   leftlld.append(Vall)
                            finalld = sum(leftlld)
                            leftld.append(finalld)
                            leftlld.clear()
              return [leftld, rightld]
       
       def Efunc(self, s, Aij,Bi,Ci):
              leftf = []
              rightf = []
              for m in range(s):
                     rightee = []
                     lefte = []
                     for n in range(s):
                            leftle = []
                            rightle = []
                            mm = m + 1
                            nn = n + 1
                            righte = 1.0/((mm+nn)*nn)
                            rightle.append(righte)
                            for i in range(s):
                                   for j in range(s):
                                          Bval = Bi[i]
                                          Cval = Ci[i]**m
                                          Aval = Aij[i][j]
                                          C2val = Ci[j]**n
                                          leftv = Bval*Aval*Cval*C2val
                                          leftle.append(leftv)
                            finalle = sum(leftle)
                            lefte.append(finalle)
                            leftle.clear()
                            rightee.append(righte)
                     rightf.append(rightee[:])
                     rightee.clear()
                     leftf.append(lefte[:])
                     lefte.clear()
              return [leftf, rightf]
                                          
       def inv(self, A):
              Ainv = sc.linalg.inv(np.array(A, dtype=float))
              return Ainv
       
       def evals(self, Am):
              A = sc.linalg.eigvals(Am)
              return A.tolist()
       
       def evects(self, Am):
              A = sp.simplify(sp.Matrix(Am))
              Eigs = A.eigenvects()
              return Eigs
       
       def diag(self, Am):
              A = sp.Matrix(Am)
              (P, D) = A.diagonalize()
              return P, D
       
       
X = butcher(17)
C = X.CiLPS()
B = X.Bi(C)
A = X.Ai(B, C)
pprint(A)





