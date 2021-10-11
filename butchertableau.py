# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 23:21:54 2020

@author: tjcze01@gmail.com
"""
import numpy as np
import scipy as sc
import sympy as sp
from sympy import I
from scipy.interpolate import CubicSpline

class butcher:
       
       def __init__(self, order, decs):
              
              # s = Stage Order
              def stage(p):
                     ss = sp.symbols("s")
                     odd = self.isodd(p)
                     if odd is True:
                            Eq = 2*ss - 1 - int(p)
                     elif odd is False:
                            Eq = 2*ss - int(p)
                     S = int(sp.solve(Eq,ss)[0])
                     return S
              
              self.order = int(order)
              self.s = stage(order)
              self.decs = int(decs)
              self.nd = int(self.decs)
       
       λ = sp.Symbol('λ')
       
       def flatten(self, lm):
              flatten = lambda l: [item for sublist in l for item in sublist] 
              return flatten(lm)
       
       def factorial(self, n):
            if n == 0:
                   return 1
            else:
                   nn = n
                   ne = n - 1
                   while ne >= 1:
                          nn = nn*ne
                          ne -= 1
                   return nn
       
       def Gamma(self, n):
              if n <= 0:
                     return None
              else:
                     return self.factorial(n-1)
              
       def isodd(self, num):
              if num % 2 == 0:
                  return False # Even 
              else:
                  return True # Odd  
           
       def Tt(self, M):
           """
           Returns a transpose of a matrix.
               :param M: The matrix to be transposed
               :return: The transpose of the given matrix
           """
           # Section 1: if a 1D array, convert to a 2D array = matrix
           if not isinstance(M[0], list):
               M = [M]
       
           # Section 2: Get dimensions
           rows = len(M)
           cols = len(M[0])
       
           # Section 3: MT is zeros matrix with transposed dimensions
           MT = self.zeros_matrix(cols, rows)
       
           # Section 4: Copy values from M to it's transpose MT
           for i in range(rows):
               for j in range(cols):
                   MT[j][i] = M[i][j]
       
           return MT
              
       def Px(self, s, evalf=True):
                     xp = sp.symbols("x")
                     if s == 0:
                            return 1.0
                     else:
                            term1 = 1.0/((2.0**s)*self.factorial(s))
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
                            term1 = 1.0/self.factorial(s)
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
              listeq = sp.nroots(sp.Poly(eq4,x), self.nd)
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
              listeq = sp.nroots(sp.Poly(eq4,x), self.nd)
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
              vs = self.flatten(sp.matrix2numpy(Vsyms).tolist())
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
       
       def radauref(self, A, B, C):
              Cn = [1.0 - C[i] for i in range(len(C))]
              AN = [[j*0.0 for j in range(len(A[0]))] for i in range(len(A))]
              for i in range(len(A)):
                     for j in range(len(A[0][:])):
                            aij = B[j] - A[i][j]
                            AN[j][i] = aij
              return AN, B, Cn
                            
       def gauss(self):
              Cs = self.CiLP(self.s)[0]
              Bs = self.Bi(self.s, Cs)
              As = self.Ai(self.s, Bs, Cs)
              return [As, Bs, Cs]
       
       def Bfunc(self, s, k, Bi, Ci):
              leftlb = []
              rightlb = []
              for j in range(k):
                     jj = int(j + 1)
                     leftllb = []
                     rightb = round(float(1.0/jj),15)
                     rightlb.append(rightb)
                     for i in range(s):
                            Bval = Bi[i]
                            Cval = Ci[i]**jj
                            Fvalb = float(Bval*Cval)
                            leftllb.append(Fvalb)
                     finallb = round(sum(leftllb),15)
                     leftlb.append(finallb)
                     leftllb.clear()
              return [leftlb, rightlb]
              
       def Cfunc(self, s, k, Aij, Ci):
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
                                   Cvalj = Ci[j]**ll
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
                            
       def Dfunc(self, s, k, A, Bi, Ci):
              leftlf = []
              rightlf = []
              for l in range(k):
                     ll = l + 1
                     leftld = []
                     rightld = []
                     for jj in range(s):
                             Cvald = Ci[jj]**ll
                             Bvald = Bi[jj]
                             fval = (Bvald*(1.0 - Cvald))/ll
                             rightld.append(fval)
                     del rightld[-1]
                     rightlf.append(sum(rightld[:]))
                     rightld.clear()
                     for j in range(s):
                            Cvald = Ci[j]**ll
                            Bvald = Bi[j]
                            top = Bvald*(1 - Cvald)
                            rightd = top/ll
                            leftlld = []
                            for i in range(s):
                                   Aij = A[i][j]
                                   Cs = Ci[i]**l
                                   Bv = Bi[i]
                                   Vall = Aij*Cs*Bv
                                   leftlld.append(Vall)
                            finalld = round(sum(leftlld), 15)
                            leftld.append(finalld)
                            leftlld.clear()
                     leftlf.append(sum(leftld[:]))
                     leftld.clear()
              
              return leftlf, rightlf
       
       def Dfuncb(self, s, A, B, C):
              leftl = []
              rightl = []
              for j in range(s):
                     x = sp.symbols('x')
                     sj = j + 1
                     llist = []
                     Bjval = B[j]
                     Cjval = C[j]
                     eq1a = self.legendreP(self.s,False,True)
                     eq2a = self.legendreP(self.s - 1,False,True)
                     eq3a = sp.Mul(sp.Integer(-1),eq2a)
                     eq4a = sp.Add(eq1a,eq3a)
                     fff = sp.integrate(eq4a, (x, 1.0 - Cjval, 1.0))
                     Pa = sp.Poly(eq4a,x)
                     Paa = Pa.all_coeffs()
                     PPa = np.polynomial.legendre.legval(Cjval, Paa)
                     Pfunc = list(np.polynomial.laguerre.lagfromroots(C))
                     Pintb = np.polynomial.laguerre.Laguerre(Pfunc)
                     Pint = list(Pintb.coef)
                     Pfb = Pintb.integ()
                     Pf = list(Pfb.coef)
                     top = np.polynomial.laguerre.lagval(1.0, Pf)
                     bot = np.polynomial.laguerre.lagval(Cjval, Pf)
                     Fjval = fff*Bjval
                     rightl.append(Fjval)
                     for i in range(s):
                            x = sp.symbols('x')
                            ss = i + 1
                            Aval = A[i][j]
                            Bval = B[i]
                            Cval = C[i]
                            eq1 = self.legendreP(self.s,False,True)
                            eq2 = self.legendreP(self.s - 1,False,True)
                            eq3 = sp.Mul(sp.Integer(-1),eq2)
                            eq4 = sp.Add(eq1,eq3)
                            P = sp.Poly(eq4,x)
                            Pc = P.all_coeffs()
                            PP = np.polynomial.legendre.legval(1.0 - Cval, Pc)
                            Pvalc = list(np.polynomial.laguerre.lagfromroots(C))
                            Pval = np.polynomial.laguerre.lagval(1.0 - Cval, Pvalc)
                            Fval = PP*(Bjval - Aval)*Bval
                            llist.append(Fval)
                     leftl.append(sum(llist))
                     llist.clear()
              return leftl, rightl
       
       def Efunc(self, s, k, l, A, Bi, Ci):
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
                                          
       def invs(self, A):
              Ainv = sc.linalg.inv(A)
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
       
       def eigs(self, Mm):
              rs = []
              cs = []
              ff = []
              Eigs = list(sc.linalg.eigvals(np.array(Mm, dtype=float)))
              
              for ei in Eigs:
                     ri = ei.real
                     ci = ei.imag
                     if ci == 0.0 or ci == 0:
                            ff.append(ri)
                     elif len(rs) > 0:
                            ril = rs[-1]
                            if ril != ri:
                                   rs.append(ri)
                                   cs.append(ci)
                            else:
                                   rs.append(ri)
                                   cs.append(-1.0*ci)
                     else:
                            rs.append(ri)
                            cs.append(ci)
              
              if len(rs) > 0:
                     for ris, cis in zip(rs, cs):
                            ind = rs.index(ris)
                            r = ris
                            iz = cis
                            z = complex(r, iz)
                            zc = z.conjugate()
                            ff.append(z)
              return ff
       
       def Tmat(self, Am):
              A = sp.Matrix(Am)
              eis = self.eigs(Am)
              M, N = A.shape
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
                     AAf = self.flatten(AA)
                     ss = sp.nonlinsolve(AAf, llist)
                     Xvec = list(ss.args[0].subs({listd[-1] : 1.00})) # One is used by default but may be wrong in certain situations
                     for iss in Xvec:
                            indexx = Xvec.index(iss)
                            Xvec[indexx] = sp.simplify(iss)
                     XXvec = []
                     for ii,jb in enumerate(Xvec):
                            vall = jb
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
              T_matrix = self.Tt(T)
              for i in range(len(T_matrix)):
                     for j in range(len(T_matrix[0])):
                            ijval = float("{:.25f}".format(T_matrix[i][j]))
                            T_matrix[i][j] = ijval
              TI_matrix = self.inv(T_matrix)
              return T_matrix , TI_matrix
       
       def zeros_matrix(self, rows, cols):
           """
           Creates a matrix filled with zeros.
               :param rows: the number of rows the matrix should have
               :param cols: the number of columns the matrix should have
               :return: list of lists that form the matrix
           """
           M = []
           while len(M) < rows:
               M.append([])
               while len(M[-1]) < cols:
                   M[-1].append(0.0)

           return M


       def identity_matrix(self, n):
           """
           Creates and returns an identity matrix.
               :param n: the square size of the matrix
               :return: a square identity matrix
           """
           IdM = self.zeros_matrix(n, n)
           for i in range(n):
               IdM[i][i] = 1.0
       
           return IdM    
          
       def copy_matrix(self, M):
           """
           Creates and returns a copy of a matrix.
               :param M: The matrix to be copied
               :return: A copy of the given matrix
           """
           # Section 1: Get matrix dimensions
           rows = len(M)
           cols = len(M[0])
       
           # Section 2: Create a new matrix of zeros
           MC = self.zeros_matrix(rows, cols)
       
           # Section 3: Copy values of M into the copy
           for i in range(rows):
               for j in range(cols):
                   MC[i][j] = M[i][j]
       
           return MC
       
       def check_matrix_equality(self, Am, Bm, tol=None):
           """
           Checks the equality of two matrices.
               :param A: The first matrix
               :param B: The second matrix
               :param tol: The decimal place tolerance of the check
               :return: The boolean result of the equality check
           """
           # Section 1: First ensure matrices have same dimensions
           if len(Am) != len(Bm) or len(Am[0]) != len(Bm[0]):
               return False
       
           # Section 2: Check element by element equality
           #            use tolerance if given
           for i in range(len(Am)):
               for j in range(len(Am[0])):
                   if tol is None:
                       if Am[i][j] != Bm[i][j]:
                           return False
                   else:
                       if round(Am[i][j], tol) != round(Bm[i][j], tol):
                           return False
       
           return True
       
       def check_squareness(self, Am):
           """
           Makes sure that a matrix is square
               :param A: The matrix to be checked.
           """
           if len(Am) != len(Am[0]):
               raise ArithmeticError("Matrix must be square to inverse.")
               
       def matrix_multiply(self, Am, Bm):
           """
           Returns the product of the matrix A * B
               :param A: The first matrix - ORDER MATTERS!
               :param B: The second matrix
               :return: The product of the two matrices
           """
           # Section 1: Ensure A & B dimensions are correct for multiplication
           rowsA = len(Am)
           colsA = len(Am[0])
           rowsB = len(Bm)
           colsB = len(Bm[0])
           if colsA != rowsB:
               raise ArithmeticError(
                   'Number of A columns must equal number of B rows.')
       
           # Section 2: Store matrix multiplication in a new matrix
           C = self.zeros_matrix(rowsA, colsB)
           for i in range(rowsA):
               for j in range(colsB):
                   total = 0
                   for ii in range(colsA):
                       total += Am[i][ii] * Bm[ii][j]
                   C[i][j] = total
       
           return C     
         
       def detr(self, Am, total=0):
           """
           Find determinant of a square matrix using full recursion
               :param A: the matrix to find the determinant for
               :param total=0: safely establish a total at each recursion level
               :returns: the running total for the levels of recursion
           """
           # Section 1: store indices in list for flexible row referencing
           indices = list(range(len(Am)))
       
           # Section 2: when at 2x2 submatrices recursive calls end
           if len(Am) == 2 and len(Am[0]) == 2:
               val = Am[0][0] * Am[1][1] - Am[1][0] * Am[0][1]
               return val
       
           # Section 3: define submatrix for focus column and call this function
           for fc in indices:  # for each focus column, find the submatrix ...
               As = self.copy_matrix(Am)  # make a copy, and ...
               As = As[1:]  # ... remove the first row
               height = len(As)
       
               for i in range(height):  # for each remaining row of submatrix ...
                   As[i] = As[i][0:fc] + As[i][fc+1:]  # zero focus column elements
       
               sign = (-1) ** (fc % 2)  # alternate signs for submatrix multiplier
               sub_det = self.detr(As)  # pass submatrix recursively
               total += sign * Am[0][fc] * sub_det  # total all returns from recursion
       
           return total
       
       def detf(self, Am):
           # Section 1: Establish n parameter and copy A
           n = len(Am)
           AM = self.copy_matrix(Am)
        
           # Section 2: Row ops on A to get in upper triangle form
           for fd in range(n): # A) fd stands for focus diagonal
               for i in range(fd+1,n): # B) only use rows below fd row
                   if AM[fd][fd] == 0: # C) if diagonal is zero ...
                       AM[fd][fd] == 1.0e-18 # change to ~zero
                   # D) cr stands for "current row"
                   crScaler = AM[i][fd] / AM[fd][fd] 
                   # E) cr - crScaler * fdRow, one element at a time
                   for j in range(n): 
                       AM[i][j] = AM[i][j] - crScaler * AM[fd][j]
            
           # Section 3: Once AM is in upper triangle form ...
           product = 1.0
           for i in range(n):
               # ... product of diagonals is determinant
               product *= AM[i][i] 
        
           return product
    
       def check_non_singular(self, Am):
           """
           Ensure matrix is NOT singular
               :param A: The matrix under consideration
               :return: determinant of A - nonzero is positive boolean
                         otherwise, raise ArithmeticError
           """
           det = self.detf(Am)
           if det != 0:
               return det
           else:
               raise ArithmeticError("Singular Matrix!")
               
       def invert_mAmtrix(self, Am, tol=None):
           """
           Returns the inverse of the pAmssed in mAmtrix.
               :pAmrAmm Am: The mAmtrix to be inversed
        
               :return: The inverse of the mAmtrix Am
           """
           # Section 1: MAmke sure Am cAmn be inverted.
           self.check_squareness(Am)
           self.check_non_singular(Am)
        
           # Section 2: MAmke copies of Am & I, AmM & IM, to use for row ops
           n = len(Am)
           AmM = self.copy_matrix(Am)
           Id = self.identity_matrix(n)
           IM = self.copy_matrix(I)
        
           # Section 3: Perform row operAmtions
           indices = list(range(n)) # to Amllow flexible row referencing ***
           for fd in range(n): # fd stAmnds for focus diAmgonAml
               fdScAmler = 1.0 / AmM[fd][fd]
               # FIRST: scAmle fd row with fd inverse. 
               for j in range(n): # Use j to indicAmte column looping.
                   AmM[fd][j] *= fdScAmler
                   IM[fd][j] *= fdScAmler
               # SECOND: operAmte on Amll rows except fd row Ams follows:
               for i in indices[0:fd] + indices[fd+1:]: 
                   # *** skip row with fd in it.
                   crScAmler = AmM[i][fd] # cr stAmnds for "current row".
                   for j in range(n): 
                       # cr - crScAmler * fdRow, but one element Amt Am time.
                       AmM[i][j] = AmM[i][j] - crScAmler * AmM[fd][j]
                       IM[i][j] = IM[i][j] - crScAmler * IM[fd][j]
        
           # Section 4: MAmke sure IM is Amn inverse of Am with specified tolerAmnce
           if self.check_matrix_equality(Id, self.matrix_multiply(Am,IM), tol):
               return IM
           else:
                  # return IM
                   raise ArithmeticError("MAmtrix inverse out of tolerAmnce.")
       
       def inv(self, Am):
           """
           Returns the inverse of the pAmssed in mAmtrix.
               :pAmrAmm Am: The mAmtrix to be inversed
        
               :return: The inverse of the mAmtrix Am
           """
           # Section 1: MAmke sure Am cAmn be inverted.
           self.check_squareness(Am)
           self.check_non_singular(Am)
        
           # Section 2: MAmke copies of Am & I, AmM & IM, to use for row ops
           n = len(Am)
           AmM = self.copy_matrix(Am)
           I = self.identity_matrix(n)
           IM = self.copy_matrix(I)
        
           # Section 3: Perform row operAmtions
           indices = list(range(n)) # to Amllow flexible row referencing ***
           for fd in range(n): # fd stAmnds for focus diAmgonAml
               fdScAmler = 1.0 / AmM[fd][fd]
               # FIRST: scAmle fd row with fd inverse. 
               for j in range(n): # Use j to indicAmte column looping.
                   AmM[fd][j] *= fdScAmler
                   IM[fd][j] *= fdScAmler
               # SECOND: operAmte on Amll rows except fd row Ams follows:
               for i in indices[0:fd] + indices[fd+1:]: 
                   # *** skip row with fd in it.
                   crScAmler = AmM[i][fd] # cr stAmnds for "current row".
                   for j in range(n): 
                       # cr - crScAmler * fdRow, but one element Amt Am time.
                       AmM[i][j] = AmM[i][j] - crScAmler * AmM[fd][j]
                       IM[i][j] = IM[i][j] - crScAmler * IM[fd][j]
        
           return IM
    
       def printm(self, Mm, decimals=25):
           """
           Print a matrix one row at a time
               :param M: The matrix to be printed
           """
           for row in Mm:
               print([round(x, decimals) + 0 for x in row])
               
       def P(self, Cm):
              Ps = []
              Cmb = [i for i in Cm]
              c0 = Cmb[0]
              if c0 == 0.0 or c0 == 0:
                     pass
              else:
                     Cmb.insert(0, 0.0)
              for i in range(len(Cmb)-1):
                     ys = [0.0 for i in range(len(Cmb))]
                     ys[i+1] = 1.0
                     CC = CubicSpline(Cmb, ys)
                     coeffs = CC.c
                     coeffs2 = (coeffs.T)[0]
                     coeffs3 = list(coeffs2[::-1])
                     del coeffs3[0]
                     Ps.append(coeffs3)
              return Ps
       
       def dot(self, v1, v2):
            return sum([x*y for x, y in zip(v1, v2)])
       
       # def dotd(self, v1, v2, pr):
       #        vv = sum([x*y for x, y in zip(v1, v2)])
       #        aa = '{}:.{}f{}'.format('{', pr,'}')
       #        aa1 = '{}'.format(aa)
       #        aa2 = str(aa1)
       #        aa3 = str(aa2.format(vv))
       #        aa4 = De(aa3)
       #        return aa4
               
       
__all__ = ["butcher"]                            
              
 
