# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 17:00:34 2020

@author: Travis Czechorski tjczec01@gmail.com
"""
import numpy as np
import pprint as pp
import sympy as sp
import scipy.linalg as la
import scipy as sc
pp = pp.pprint
EPS = np.finfo(float).eps
print("{:.52f}".format(EPS))
As = sp.Matrix([[3, 2, 3],
     [4, 6, 6],
     [7, 4, 9]])
Bs = sp.Matrix([[5, 5], [6, 7], [9, 9]])
AM = [[3, 2, 3],      [4, 6, 6],      [7, 4, 9]]
B = [[5, 5], [6, 7], [9, 9]]
AN = np.array(AM)
# print(np.array(A).T)
print(len("0000000000000002220446049250313080847263336181640625"))
flatten = lambda l: [item for sublist in l for item in sublist] 

def permutations(iterable, r=None):
    # permutations('ABCD', 2) --> AB AC AD BA BC BD CA CB CD DA DB DC
    # permutations(range(3)) --> 012 021 102 120 201 210
    pool = tuple(iterable)
    n = len(pool)
    r = n if r is None else r
    if r > n:
        return
    indices = list(range(n))
    cycles = list(range(n, n-r, -1))
    yield tuple(pool[i] for i in indices[:r])
    while n:
        for i in reversed(range(r)):
            cycles[i] -= 1
            if cycles[i] == 0:
                indices[i:] = indices[i+1:] + indices[i:i+1]
                cycles[i] = n - i
            else:
                j = cycles[i]
                indices[i], indices[-j] = indices[-j], indices[i]
                yield tuple(pool[i] for i in indices[:r])
                break
        else:
            return
     
def rotate_left(x, y):
    """
    Left rotates a list x by the number of steps specified
    in y.
    Examples
    ========
    >>> from sympy.utilities.iterables import rotate_left
    >>> a = [0, 1, 2]
    >>> rotate_left(a, 1)
    [1, 2, 0]
    """
    if len(x) == 0:
        return []
    y = y % len(x)
    return x[y:] + x[:y]


def rotate_right(x, y):
    """
    Right rotates a list x by the number of steps specified
    in y.
    Examples
    ========
    >>> from sympy.utilities.iterables import rotate_right
    >>> a = [0, 1, 2]
    >>> rotate_right(a, 1)
    [2, 0, 1]
    """
    if len(x) == 0:
        return []
    y = len(x) - y % len(x)
    return x[y:] + x[:y]  

def minlex(seq, directed=True, is_set=False, small=None):
    """
    Return a tuple where the smallest element appears first; if
    ``directed`` is True (default) then the order is preserved, otherwise
    the sequence will be reversed if that gives a smaller ordering.
    If every element appears only once then is_set can be set to True
    for more efficient processing.
    If the smallest element is known at the time of calling, it can be
    passed and the calculation of the smallest element will be omitted.
    Examples
    ========
    >>> from sympy.combinatorics.polyhedron import minlex
    >>> minlex((1, 2, 0))
    (0, 1, 2)
    >>> minlex((1, 0, 2))
    (0, 2, 1)
    >>> minlex((1, 0, 2), directed=False)
    (0, 1, 2)
    >>> minlex('11010011000', directed=True)
    '00011010011'
    >>> minlex('11010011000', directed=False)
    '00011001011'
    """
    is_str = type(seq)
    seq = list(seq)
    if small is None:
       small = min(seq)
    if is_set:
        i = seq.index(small)
        if not directed:
            n = len(seq)
            p = (i + 1) % n
            m = (i - 1) % n
            if seq[p] > seq[m]:
                seq = list(reversed(seq))
                i = n - i - 1
        if i:
            seq = rotate_left(seq, i)
        best = seq
    else:
        count = seq.count(small)
        if count == 1 and directed:
            best = rotate_left(seq, seq.index(small))
        else:
            # if not directed, and not a set, we can't just
            # pass this off to minlex with is_set True since
            # peeking at the neighbor may not be sufficient to
            # make the decision so we continue...
            best = seq
            for i in range(count):
                seq = rotate_left(seq, seq.index(small, count != 1))
                if seq < best:
                    best = seq
                # it's cheaper to rotate now rather than search
                # again for these in reversed order so we test
                # the reverse now
                if not directed:
                    seq = rotate_left(seq, 1)
                    seq = list(reversed(seq))
                    if seq < best:
                        best = seq
                    seq = list(reversed(seq))
                    seq = rotate_right(seq, 1)
    # common return
    if is_str == str:
        return ''.join(best)
    return tuple(best)       

def cyclic_form(Pl):
        """
        This is used to convert to the cyclic notation
        from the canonical notation. Singletons are omitted.
        Examples
        ========
        >>> from sympy.combinatorics.permutations import Permutation
        >>> Permutation.print_cyclic = False
        >>> p = Permutation([0, 3, 1, 2])
        >>> p.cyclic_form
        [[1, 3, 2]]
        >>> Permutation([1, 0, 2, 4, 3, 5]).cyclic_form
        [[0, 1], [3, 4]]
        See Also
        ========
        array_form, full_cyclic_form
        """
        pt = type(Pl)
        if pt == tuple:
               Pl = list(Pl)
               # return list(Pl)
        elif pt == str:
               raise Exception('Given Value must be either tuple or list')
        elif pt == int:
               raise Exception('Given Value must be either tuple or list')
        elif pt == float:
               raise Exception('Given Value must be either tuple or list')
        elif pt == list:
               Pl = Pl
                # return Pl
        array_form = Pl
        unchecked = [True] * len(Pl)
        cyclic_form = []
        for i in range(len(Pl)):
            if unchecked[i]:
                cycle = []
                cycle.append(i)
                unchecked[i] = False
                j = i
                while unchecked[array_form[j]]:
                    j = array_form[j]
                    cycle.append(j)
                    unchecked[j] = False
                if len(cycle) > 1:
                    cyclic_form.append(cycle)
                    assert cycle == list(minlex(cycle, is_set=True))
        cyclic_form.sort()
        cyclic_form = cyclic_form[:]
        return cyclic_form
 
def transpositions(Pl):
        """
        Return the permutation decomposed into a list of transpositions.
        It is always possible to express a permutation as the product of
        transpositions, see [1]
        Examples
        ========
        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([[1, 2, 3], [0, 4, 5, 6, 7]])
        >>> t = p.transpositions()
        >>> t
        [(0, 7), (0, 6), (0, 5), (0, 4), (1, 3), (1, 2)]
        >>> print(''.join(str(c) for c in t))
        (0, 7)(0, 6)(0, 5)(0, 4)(1, 3)(1, 2)
        >>> Permutation.rmul(*[Permutation([ti], size=p.size) for ti in t]) == p
        True
        References
        ==========
        .. [1] https://en.wikipedia.org/wiki/Transposition_%28mathematics%29#Properties
        """
        al = [i for i in range(len(Pl))]
        if al == Pl is True:
               return [0, 0]    
        else:
               a = cyclic_form(Pl)
               res = []
               for x in a:
                   nx = len(x)
                   if nx == 2:
                       res.append(tuple(x))
                   elif nx > 2:
                       first = x[0]
                       for y in x[nx - 1:0:-1]:
                           res.append((first, y))
               return res      


def transpose(M):
       columns = len(M[0][:])
       rows = len(M)
       tmat = [[i*0.0 for i in range(rows)] for j in range(columns)]
       for i in range(rows):
              for j in range(columns):
                     tmat[j][i] = M[i][j]
       return tmat

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

def perm(n, r):
       pp = factorial(n)/factorial(n-r)
       return int(pp)

def comb(n, r):
       cc = factorial(n)/(factorial(n-r)*factorial(r))
       return cc

def Gamma(n):
       if n <= 0:
              return None
       else:
              return factorial(n-1)

def dot(v1, v2):
     return sum([x*y for x,y in zip(v1,v2)])

def matmul(A, B):
       acolumns = len(A[0][:])
       arows = len(A)
       bcolumns = len(B[0][:])
       brows = len(B)
       if acolumns == brows:
              nmat = [[i*0.0 for i in range(bcolumns)] for j in range(arows)]
              for i in range(arows):
                     Ar = A[i][:]
                     for j in range(bcolumns):
                            Bc = [B[i][j] for i in range(brows)]
                            Cij = dot(Ar, Bc) #sum([i*j for i, j in zip(Ar,Bc)])
                            nmat[i][j] = Cij
              return nmat
       
       elif acolumns != brows:
              raise Exception('Columns of matrix A ({}) needs to equal the Rows of Matrix B ({}) {} != {}'.format(acolumns, brows, acolumns, brows))
       
def kcycle(S, k):
       n = len(S)
       kn = factorial(n + 1)/(factorial(n - k + 1)*k)
       return kn
       

def multiply(n):
    total = 1
    for i in n:
        total *= i
    return total       

def sgnf(m):
       sgn = (-1)**m
       return sgn

def det(A):
       det = []
       colsl = len(A[0][:])
       p = [list(i) for i in list(permutations([pi for pi in range(colsl)]))]
       ts = []
       tns = []
       ss = []
       ais = []
       for pi in range(len(p)):
              tl = transpositions(p[pi])
              ts.append(tl)
              tn = len(transpositions(p[pi]))
              tns.append(tn)
       for i in tns:
              ss.append(sgnf(i))
       for fi in range(len(p)):
                            σ = [i + 1 for i in p[fi]]
                            sig = ss[fi]
                            for i, j in enumerate(σ):
                                   ai = A[i-1][j-1]
                                   ais.append(ai)
                            fin = sig*multiply(ais)
                            det.append(fin)
                            ais.clear()
       return sum(det)


def iden(n):
       mm = []
       for ni in range(n):
              mm.append([mi*0.0 for mi in range(n)])
       
       for nn in range(n):
              mm[nn][nn] = 1
       return mm

def zero(n):
       mm = []
       for ni in range(n):
              mm.append([mi*0.0 for mi in range(n)])
       
       for nn in range(n):
              mm[nn][nn] = 0.0
       return mm              

def rowred(A):
       rows = len(A)
       cols = len(A[0])
       itern = 0
       II = iden(rows)
       for i in range(cols):
              start = A[i][i]
              vals = [A[j][i]/start for j in range(1 + itern, rows)]
              nr = [[vals[v]*rv for rv in A[i][:]] for v in range(len(vals))]
              ni = [[vals[v]*rv for rv in II[i][:]] for v in range(len(vals))]
              rrows = [A[iv] for iv in range(1 + itern, cols)]
              rrowsi = [II[iv] for iv in range(1 + itern, cols)]
              # print(np.array(A))
              for g in range(len(rrows)):
                     vv = nr[g]
                     nn = rrows[g]
                     nnr = [i - j for i, j in zip(nn, vv)]
                     vvi = ni[g]
                     nni = rrowsi[g]
                     nnii = [i - j for i, j in zip(nni, vvi)]
                     A[g+1+itern] = nnr
                     II[g+1+itern] = nnii
              itern += 1
       return A, II

def solsys(A, II):
       rows = len(A)
       cols = len(A[0])
       itern = 0
       for i in range(cols):
              start = A[-1-i][-1-i]
              vals = [A[j][-1 - itern]/start for j in range(0, rows-itern)]
              nr = [[vals[v]*rv for rv in A[-1 - itern][:]] for v in range(len(vals))]
              rrows = [A[iv] for iv in range(0, cols-itern-1)]
              ni = [[vals[v]*rv for rv in II[-1 - itern][:]] for v in range(len(vals))]
              rrowsi = [II[iv] for iv in range(0, cols-itern-1)]
              for g in range(len(rrows)):
                     vv = nr[g]
                     nn = rrows[g]
                     nnr = [round(i - j, 5) for i, j in zip(nn, vv)]
                     vvi = ni[g]
                     nni = rrowsi[g]
                     nnii = [round(i - j, 5) for i, j in zip(nni, vvi)]
                     A[g] = nnr
                     II[g] = nnii
              itern += 1
       for i in range(rows):
              start = A[i][i]
              IIv = [iv/start for iv in II[i]]
              AAv = [av/start for av in A[i]]
              A[i] = AAv
              II[i] = IIv
       return A, II

# def LUf(A):
#        rows = len(A)
#        cols = len(A[0])
       
#        def U_mat(row, col):
              
#               BMat = [[0.0*j for j in range(col)] for i in range(row)]
              
#               for i in range(row):
                     
#                      for j in range(col):
                            
#                             if i == j:
#                                    strv = "U_{}{}".format(i, j)
#                                    BMat[i][j] = strv
#                             elif i > j:
#                                    BMat[i][j] = 0.0
#                             elif i < j:
#                                    strv = "U_{}{}".format(i, j)
#                                    BMat[i][j] = strv
#               return BMat
                            
                            
#        def L_mat(row, col):
              
#               BMatl = [[0.0*j for j in range(col)] for i in range(row)]
              
#               for i in range(row):
                     
#                      for j in range(col):
                            
#                             if i == j:
#                                    strv = "L_{}{}".format(i, j)
#                                    BMatl[i][j] = 1.0
#                             elif i < j:
#                                    BMatl[i][j] = 0.0
#                             elif i > j:
#                                    strv = "L_{}{}".format(i, j)
#                                    BMatl[i][j] = strv
#               return BMatl

def matmullu(A, B):
              acolumns = len(A[0][:])
              arows = len(A)
              bcolumns = len(B[0][:])
              brows = len(B)
              if acolumns == brows:
                     nmat = [[i*0.0 for i in range(bcolumns)] for j in range(arows)]
                     for i in range(arows):
                            Ar = A[i][:]
                            for j in range(bcolumns):
                                   Bc = [B[i][j] for i in range(brows)]
                                   print(Bc)
                                   Cc = [A[ia][j] for ia in range(brows)]
                                   print(Cc)
                                   Cij = sum(["{}*{}".format(Bc[i], Cc[i])  for i in range(acolumns)])
                                   print(Cij)
                                   nmat[i][j] = Cij
                     return nmat
              
              elif acolumns != brows:
                     raise Exception('Columns of matrix A ({}) needs to equal the Rows of Matrix B ({}) {} != {}'.format(acolumns, brows, acolumns, brows))
              

def LU_decomposition(A):
    """Perform LU decomposition using the Doolittle factorisation."""

    L = zero(len(A))
    U = zero(len(A))
    N = len(A)
    
    def uvals(Um, k, n):
            ulist = []
            for i in range(k):
                  uu = Um[i][n]
                  ulist.append(uu)
            return ulist
     
    def lvals(Lm, k, n):
            llist = []
            lu = Lm[k]
            lul = lu[0:k]
            return lul
       
    for k in range(N):
        L[k][k] = 1
        U[k][k] = (A[k][k] - dot(lvals(L, k, k), uvals(U, k, k))) / L[k][k]
        for j in range(k+1, N):
            U[k][j] = (A[k][j] - dot(lvals(L, k, k), uvals(U, j, j))) / L[k][k]
        for i in range(k+1, N):
            L[i][k] = (A[i][k] - dot(lvals(L, i, i), uvals(U, k, k))) / U[k][k]

    return L, U


def forward_sub(L, b):
    """Given a lower triangular matrix L and right-side vector b,
    compute the solution vector y solving Ly = b."""

    y = []
    for i in range(len(b)):
        y.append(b[i])
        for j in range(i):
            y[i]=y[i]-(L[i][j]*y[j])
        y[i] = y[i]/L[i][i]

    return y

def backward_sub(U, y):
    """Given a lower triangular matrix U and right-side vector y,
    compute the solution vector x solving Ux = y."""

    x = zero(len(y))

    for i in range(len(x), 0, -1):
      x[i-1] = (y[i-1] - dot(U[i-1][i:], x[i:])) / U[i-1][i-1]

    return x

def lu_solve(L, U, b):
    # Step 1: Solve Uy = b using forward substitution
    # Step 2: Solve Lx = y using backward substitution
    y = forward_sub(L,b)
    x = backward_sub(U,y)
    return x

def linear_solve(A, b):
    L, U = LU_decomposition(A)
    x = lu_solve(L, U,b)
    return x


def LUdecomp(a):
    n = len(a)
    for k in range(0,n-1):
        for i in range(k+1,n):
            if a[i,k] != 0.0:
                lam = a[i,k]/a[k,k]
                a[i,k+1:n] = a[i,k+1:n] - lam*a[k,k+1:n]
                a[i,k] = lam
    return a

def sparsity(A):
       rows = len(A)
       cols = len(A[0])
       Zs = 0
       Nzs = 0
       tot = 0
       for i in range(rows):
              for j in range(cols):
                     val_ij = A[i][j]
                     if val_ij == 0 or val_ij == 0.0:
                            Zs += 1
                            tot += 1
                     elif val_ij != 0 or val_ij != 0.0:
                            Nzs += 1
                            tot += 1
       return Zs, Nzs, tot, Zs/tot

def csc_groups(A):
       rows = len(A)
       cols = len(A[0])
       rows_l = []
       cols_l = []
       data_l = []
       Zs = []
       Nzs = []
       tot = 0
       for i in range(rows):
              for j in range(cols):
                     val_ij = A[i][j]
                     if val_ij == 0 or val_ij == 0.0:
                            pass
                     elif val_ij != 0 or val_ij != 0.0:
                            rows_l.append(int(i))
                            cols_l.append(int(j))
                            data_l.append(A[i][j])
       return rows_l, cols_l, data_l

Av, IIa = rowred(AM)
B, IF = solsys(Av, IIa)
# pp((np.array(IF)).tolist())
# print(det(A))
# pp(np.array(Av))
# pp(np.array(IIa))
# pp(np.array(B))
# pp(np.array(IF))

AA = np.array([[3, 2, 3],     [4, 6, 6],     [7, 4, 9]])
aa = transpose([[3, 4, 7], [2, 6, 4], [3, 6, 9]])
L, U = LU_decomposition(aa)
# LUf(AM)
# print(P)
# print(np.array(L))
# print(np.array(U))

P, LL, UU = la.lu(AA.T)
# print(np.array(P))
# print(np.array(LL))
# print(np.array(UU))

b = [6,-4,27]
As = [[3, 2, 3],     [4, 6, 6],     [7, 4, 9]]
LU = la.lu_factor(As)
x = la.lu_solve(LU, b)
print(linear_solve(As,b))
print(x)

ZS, NZS, TOT, sps = sparsity([[3, 0, 3],     [4, 0.0, 0],     [0, 4, 9]])
print(ZS, NZS, TOT, sps)

S6 = 6**0.5

P = [[13/3 + 7*S6/3, -23/3 - 22*S6/3, 10/3 + 5 * S6],
    [13/3 - 7*S6/3, -23/3 + 22*S6/3, 10/3 - 5 * S6],
    [1/3, -8/3, 10/3]]
print(P)
row = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
col = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])
data = np.array([3, 2, 3, 4, 6, 6, 7, 4, 9])
print(sc.sparse.csc_matrix((data,(row,col)), shape=(3,3)).todense() )
rowsl, colsl, datal = csc_groups(As)
# print(rowsl, colsl, datal)
print(sc.sparse.csc_matrix((datal,(rowsl,colsl)), shape=(3,3)).todense() )