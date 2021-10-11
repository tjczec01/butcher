# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 23:21:54 2020

@author: tjcze
"""

import time
import numpy as np
import scipy as sc
import sympy as sp
import pprint as pp
from sympy import I
from scipy import linalg
from scipy import optimize
from scipy.optimize import fsolve
from sympy.matrices import *
from numpy.linalg import *
pprint = pp.pprint

start = time.time() #Real time when the program starts to run

flatten = lambda l: [item for sublist in l for item in sublist] 
λ = sp.Symbol('λ')

# s = Stage Order
def sym2num(x):
              try:
                     v = float(sp.re(x))*1.0 + float(sp.im(x))*1j
                     if v.imag ==0:
                            return v.real
                     else:
                            return complex(v.real, v.imag)
              except:
                     return float(x)
              
       

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
       
def isodd(num):
       if num % 2 == 0:
           return False # Even 
       else:
           return True # Odd       

def stage(p):
        s = sp.symbols("s")
        odd = isodd(p)
        if odd is True:
              Eq = 2*s - 1 - int(p)
        elif odd is False:
              Eq = 2*s - int(p)
        S = int(sp.solve(Eq,s)[0])
        return S


class butcher:
       
       def stage(self, p):
              ss = sp.symbols("s")
              odd = isodd(p)
              if odd is True:
                     Eq = 2*ss - 1 - int(p)
              elif odd is False:
                     Eq = 2*ss - int(p)
              S = int(sp.solve(Eq,ss)[0])
              return S
       
       def __init__(self, order):
              self.order = int(order)
              self.s = self.stage(order)
              
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
       def CiLPS(self, s):
              x = sp.symbols('x')
              eq1 = self.legendrePS(self.s,False,True)
              eq2 = self.legendrePS(self.s - 1,False,True)
              eq3 = sp.Mul(sp.Integer(-1),eq2)
              eq4 = sp.Add(eq1,eq3)
              P = sp.Poly(eq4,x)
              Pc = P.all_coeffs()
              roots = list(np.roots(Pc))[::-1]
              listeq = sp.nroots(sp.Poly(eq4,x),15)
              return roots, listeq
              
       def CiLP(self, s):
              x = sp.symbols('x')
              eq1 = self.legendreP(self.s,False,True)
              eq2 = self.legendreP(self.s-1,False,True)
              eq3 = sp.Mul(sp.Integer(-1),eq2)
              eq4 = sp.Add(eq1,eq3)
              P = sp.Poly(eq4,x)
              Pc = P.all_coeffs()
              roots = list(np.roots(Pc))[::-1]
              listeq = sp.nroots(sp.Poly(eq4,x),15)
              return roots, listeq
       
       def Bi(self, s, C):
              bsyms = [sp.symbols("b{}".format(i+1)) for i in range(len(C))]
              eqs = []
              for i in range(len(C)):
                     eqi = [sp.Mul(ij**(i),j) for ij,j in zip(C,bsyms)]
                     eqif = sp.Add(sum(eqi),sp.Mul(sp.Integer(-1),1.0/(i+1)))
                     eqs.append(eqif)
              bs = sp.solve(eqs,bsyms)
              return list(bs.values())
       
       def Ai(self, s, B, C):
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
                            eqsa.append([term5])
                            termeq.clear()
                     iv += 3
              vs = flatten(sp.matrix2numpy(Vsyms).tolist())
              Aslin = sp.linsolve(sp.Matrix(eqsa), vs)
              linans = list(Aslin.args[:][0])
              lin2 = np.array(linans, dtype=float).reshape((n-1,n))
              lin3 = lin2.tolist()
              lin3.append(B[:])
              return lin3
       
       def radau(self):
              Cs = self.CiLPS(self.s)[0]
              Cs[-1] = 1.0
              Bs = self.Bi(self.s, Cs)
              As = self.Ai(self.s, Bs, Cs)
              return [As, Bs, Cs]
       
       def gauss(self):
              Cs = self.CiLP(self.s)[0]
              Bs = self.Bi(self.s, Cs)
              As = self.Ai(self.s, Bs, Cs)
              return [As, Bs, Cs]
       
       def Bfunc(self, s, k, Bi,Ci):
              leftlb = []
              rightlb = []
              for j in range(k):
                     jj = int(j + 1)
                     leftllb = []
                     rightb = round(float(1.0/jj),15)
                     rightlb.append(rightb)
                     for i in range(s):
                            Bval = Bi[i]
                            Cval = Ci[i]**j
                            Fvalb = float(Bval*Cval)
                            leftllb.append(Fvalb)
                     finallb = round(sum(leftllb),15)
                     leftlb.append(finallb)
                     leftllb.clear()
              return [leftlb, rightlb]
              
       def Cfunc(self, s, k, Aij,Ci):
              leftlcf = []
              rightlcf = []
              for l in range(k):
                     leftlc = []
                     rightlc = []
                     ll = l + 1
                     for i in range(s):
                            leftllc = []
                            Cvalct = Ci[i]**ll
                            Cvali = Cvalct/ll
                            rightc = Cvali
                            rightlc.append(rightc)
                            for j in range(s):
                                   Avalsf = Aij[i][j]
                                   Cvalj = Ci[j]**l
                                   Fvalc = Avalsf*Cvalj
                                   leftllc.append(Fvalc)
                            finallc = sum(leftllc)
                            leftlc.append(finallc)
                            leftllc.clear()
                     leftlcf.append(leftlc[:])
                     leftlc.clear()
                     rightlcf.append(rightlc[:])
                     rightlc.clear()
              return [leftlcf, rightlcf]
                            
       def Dfunc(self, s, k, A,Bi,Ci):
              leftlf = []
              rightlf = []
              for l in range(k):
                     ll = l + 1
                     leftld = []
                     rightld = []
                     for j in range(s):
                            Bvald = Bi[j]
                            Cvald = Ci[j]
                            top = Bvald*(1.0 - Cvald**ll)
                            rightd = top/ll
                            rightld.append(rightd)
                            leftlld = []
                            for i in range(s):
                                   Aij = A[i][j]
                                   Cs = Ci[i]**(ll-1) 
                                   Bv = Bi[i]
                                   Vall = Aij*Cs*Bv
                                   leftlld.append(Vall)
                            finalld = round(sum(leftlld), 15)
                            leftld.append(finalld)
                            leftlld.clear()
                     leftlf.append(leftld[:])
                     leftld.clear()
                     rightlf.append(rightld[:])
                     rightld.clear()
              return [leftlf, rightlf]
       
       def Efunc(self, s, k, l, A,Bi,Ci):
              leftf = []
              rightf = []
              for m in range(k):
                     rightee = []
                     lefte = []
                     for n in range(l):
                            leftle = []
                            rightle = []
                            mm = m + 1
                            nn = n + 1
                            bottom = (mm + nn) * nn
                            righte = round(1.0/bottom, 15)
                            rightle.append(righte)
                            for i in range(s):
                                   Cival = Ci[i]**m
                                   for j in range(s):
                                          Bval = Bi[i]
                                          
                                          Aij = A[i][j]
                                          Cjval = Ci[j]**n
                                          leftv = round(Bval*Cival*Aij*Cjval, 15)
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
              A = sp.Matrix(Am)
              Eigs = A.eigenvects()
              return Eigs
       
       def diag(self, Am):
              A = sp.Matrix(Am)
              (P, D) = A.diagonalize()
              return P, D
       
       def jordan(self, Am, calc=True):
              A = sp.Matrix(Am)
              if calc is True:
                     Pa, Ja = A.jordan_form(calc_transform=calc)
                     return Pa, Ja
              if calc is False:
                     Jb = A.jordan_form(calc_transform=calc)
                     return Jb
          
       def char(self, Am):
              A = sp.Matrix(Am)
              M, N = np.array(Am).shape
              λ = sp.Symbol('λ')
              II = sp.eye(M)
              Ad = A - II*λ
              Deta = Ad.det()
              return Deta 
       
       def sroots(self, Am, tol=None, nn=None, char=None):
              #returns sympy format, python format
              λ = sp.Symbol('λ')
              if tol == None:
                     tol = 10**-15
              if nn == None:
                     nv = 15
              elif nn != None:
                     nv = int(nn)
              if char == None:
                     A = sp.Matrix(Am)
                     M, N = np.array(Am).shape
                     II = sp.eye(M)
                     Ad = A - II*λ
                     Deta = Ad.det()
              elif char != None:
                     Deta = char
              Detsimp = sp.nfloat(sp.nsimplify(Deta, tolerance=tol, full=True), n=nv)
              # rlist = list(sp.solve(Detsimp, λ))
              rlist = list(sp.solve(Detsimp, λ,  **{'set': True ,'particular': True}))
              rootsd = [sp.nsimplify(i, tolerance=tol, full=True) for i in rlist]
              sympl = []
              numpr = []
              numpi = []
              for i, j in enumerate(rootsd):
                     sympl.append(sp.simplify(sp.nfloat(rootsd[i], n=nv), rational=False))
                     vals = j.as_real_imag()
                     vali = sp.nfloat(vals[0], n=nv)
                     valj = sp.nfloat(vals[1], n=nv)
                     numpr.append(vali)
                     numpi.append(valj)
              reals = []
              iss = []
              simps = []
              for i, j in zip(numpr, numpi):
                    reale = i
                    compe = j
                    reals.append(reale)
                    if compe != 0.0:
                           iss.append(compe)
                           simps.append(complex(i,j))
                    else:
                           simps.append(reale)  
                           
              return rlist, simps 
       
       def alpha(self, eigl):
              eigs = list(eigl)
              lambds = []
              alps = []
              betas = []
              for i in eigs:
                     rr = i.real
                     ii = i.imag
                     ii2 = -ii
                     if rr not in alps and ii != 0.0:
                            alps.append(rr)
                     if ii not in betas and ii2 not in betas and ii != 0.0:
                            betas.append(ii)
                     if ii == 0.0:
                            lambds.append(rr)
                     
              return lambds, alps, betas
       
       def block(self, ls, als, bs):
              matrices = []
              for i in range(len(als)):
                     matrixi = sp.Matrix([[als[i], -bs[i]], [bs[i], als[i]]])
                     
                     matrices.append(matrixi)
              B = sp.BlockDiagMatrix(sp.Matrix([ls]),*matrices)
              return B
       
       def Tmat(self, Am, eis):
              A = sp.Matrix(Am)
              M, N = np.array(Am).shape
              Xi = sp.symarray('x',M)
              listd = [sp.symbols('x_{}'.format(i)) for i in range(M)]
              llist = [sp.symbols('x_{}'.format(i)) for i in range(M-1)]
              T = []
              realss = []
              comps = []
              imagine = bool(False)
              count = 0
              for i in eis:
                     II = sp.eye(M)
                     Ad = A - i*II
                     AA = sp.matrix2numpy(Ad*sp.Matrix(Xi))
                     AAf = flatten(AA)
                     ss = sp.nonlinsolve(AAf, *llist)
                     # print(ss)
                     Xvec = list(sp.simplify(ss.args[0].subs({listd[-1] : -1.0})))
                     for iss in Xvec:
                            indexx = Xvec.index(iss)
                            Xvec[indexx] = sp.simplify(iss)
                     XXvec = []
                     for ii,jb in enumerate(Xvec):
                            vall = jb*1.0
                            Xvec[ii] = vall 
                            vali = sp.re(vall)
                            valj = sp.im(vall)
                            realss.append(vali)
                            comps.append(valj)
                            if valj != 0.0:
                                   imagine = bool(True)
                     if imagine == True:
                            count += 1
                     realss.insert(len(realss), 1)  
                     comps.insert(len(comps), 0)
                     if count % 2 == 0 and imagine == False:
                            T.append(realss[:])
                            realss.clear()
                            comps.clear()
                     elif count % 2 != 0 and imagine == True:
                             T.append(realss[:])
                             realss.clear()
                             comps.clear()
                     elif count % 2 == 0 and imagine == True:
                            T.append(comps[:])
                            realss.clear()
                            comps.clear()
                     Xvec.clear()
                     XXvec.clear()
              Tfixed = np.array(T, dtype=float).T
              Tinv = sc.linalg.inv(Tfixed)
              return Tfixed, Tinv
              
       
E = [-10.0488094 ,   1.38214273,  -0.33333333]       
# order = int(input("Enter order of desired butcher tableau --> "))
order = 17
ni = 15 # Decimal places for nfloat
toll = 10E-15 # Tolerence for nsimplify
X = butcher(order)
A, B, C = X.radau()
An = np.array(A, dtype=float)
Ainv = X.inv(A)
Aninv = np.array(Ainv)
Ans = sp.Matrix(Ainv)
Eigs = sc.linalg.eigvals(Aninv)
ljd = list(Eigs)
eigenvals = list(np.sort(np.array(Eigs)))
eigenvals.insert(0, eigenvals.pop(-1))
ljd2 = eigenvals.copy()
Tf, TIf = X.Tmat(Ainv, eigenvals)
T = Tf.tolist()
TI = TIf.tolist()
print(T)
print(TI)
