# include <stdlib.h>
# include <stdio.h>
# define Long long long int
from copy import deepcopy
from typing import List
from sympy.functions.combinatorial.numbers import totient
from sympy.ntheory import factorint
from utils.tests import dump_args


# /* factors squarefree n, writes prime divisors of n to P in increasing order
#    writes phi(n) to phi */
# int factor64(Long n, Long P[10], Long * phi){ /* factors integers */
#   if(n == 0) return -1; else if(n == 1) return 0;
#   if(n % 2 == 0) return -1; Long i,q,q2,m;
#   for(i=0, q=3, q2=9, m=n, *phi=1; q2 <= m && q2 >= q; q2+=4*q+4, q+=2){
#     if(m % q == 0){ m /= q; if(m % q == 0) return -1; P[i++] = q; *phi *= q-1; }}
#   P[i++] = m; *phi *= m-1;
#   for(q=i-1,q2=0; q>q2; q--,q2++){ m=P[q]; P[q]=P[q2]; P[q2]=m; }
#   return i; }

class SqFreeFactors:
    """
    A class used to represent square free factors of an integer

    value: integer value
    totient: Euler's totient of the integer
    factors: the prime factors of the integer
    """
    value = None
    totient = None
    factors = None

    def __init__(self, value: int) -> None:
        self.value = value
        self.totient = totient(value)
        dict_factors = factorint(value)
        self.factors = [*dict_factors]
        self.factors.sort()

    def __repr__(self):
        return f"{self.factors}"

    def has_prime(self) -> bool:
        return self.value > 1

    def remove_last(self):
        prime = self.factors[-1]
        self.value //= prime
        self.totient //= prime - 1
        del self.factors[-1]
        return prime


@dump_args
def sps4(m: SqFreeFactors, e: int, alpha: bool):
    ms, es = deepcopy(m), e
    while ms.has_prime():
        p = ms.remove_last()
        sps4(ms, es, not alpha)
        es *= p


def cyclotomic(m: SqFreeFactors, alpha: bool):
    """Return a cyclotomic or inverse cyclotomic polynomial.

    Parameters
    ==========
    m : SqFreeFactors
        Square free odd integer
    alpha: bool
        True to compute the cyclotomic polynomial, otherwise compute the inverse cyclotomic polynomial
    """
    pass


'''/*SPS4
Input:
   -m,e : squarefree integers
   -phi_m = phi(m)
   -alpha : a boolean
   -Df : degree of input polynomial f
   -D : the least degree of a polynomial we know will occur at a later
      stage of computation
   -A : an array containing the terms of f up to degree
      floor( min(Df,D)/2 )
   -P : a list of distinct primes, in increasing order, that contains
      all the prime divisors of m
   -mark : an integer whose ith bit is 1 if P[i] divides m
Output:
    -If alpha is nonzero, SPS4 computes the terms of g=f*Phi_m(z^e),
       where Phi_m is the mth cyclotomic polynomial, up to degree
       floor( min(Dg,D)/2 ), where Dg is the degree of g.  The
       computed coefficients of g are written to array A.
    -If alpha=0, SPS4 instead computes the terms of g=Psi_m(z^e),
       where Psi_m is the mth inverse cyclotomic polynomial.
    -SPS4 returns Dg, the degree of g  */
'''


# Long SPS4(Long m, Long phi_m, Long e, char alpha, Long Df, Long D,
#  Long * A, Long P[10], int mark){
def sps4_old(m: int, phi_m: int, e: int, alpha: bool, Df: int, D: int, A: List[int], P: List[int], mark: int):
    #
    #    Long i,k,d;
    #    Long Dg = Df+e*(alpha ? phi_m : m - phi_m);
    Dg: int = Df + e * (phi_m if alpha else m - phi_m)

    #    Long D1 = Dg < D ? Dg : D;
    D1: int = Dg if Dg < D else Dg
    #    Long e1 = e, m1 = m, mark1=mark;
    e1: int = e
    m1: int = m
    mark1: int = mark


#
#    for(k=mark,i=0; k>0; k=k>>1,i++){
#       if(k&1){
#          m1/=P[i]; phi_m/=P[i]-1; mark1 -= (1<<i);
#          if(m1>1 || !alpha){
#             Df = SPS4(m1, phi_m, e1, !alpha, Df, D1, A, P, mark1); }
#          e1 *= P[i];
#    }}
#    if (!alpha) return Dg;
#
#    /* Generate higher-degree terms */
#    Long halfDf = Df>>1;  D1=D1>>1;
#    k=halfDf+1; i=Df-halfDf-1;
#    if(Df & 1){ while(i>=0 && k<=D1){ A[k] = -A[i]; i--; k++; }}
#    else { while(i>=0 && k<=D1){ A[k] = A[i]; i--; k++;}}
#
#    /* multiply by 1/(1-z^(me/p)) for p|m */
#    for(k=0, m=m*e; mark>0; k++,mark=mark>>1){
#       if(mark & 1){ d=m/P[k];
#          for(i=d; i<=D1; i++){ A[i] = A[i]+A[i-d]; }
#    }}
#
#    /* multiply by 1-z^(me) */
#    for(i=D1; i>=m; i--){ A[i] = A[i] - A[i-m]; }
#    return Dg;
# }


'''
Long cyclotomic(Long n, char alpha){
    Long deg,phi,i,P[10];
    printf("%lld & ",n);
    int omega = factor64(n,P,&phi);
    if(omega < 0){ printf("Input is not squarefree\n"); exit(0); }
    for(i=0; i<omega; i++) printf("(%lld)", P[i]);
    deg = alpha ? phi : n-phi;
    Long * A = malloc(sizeof(Long)*(deg/2+1));
    A[0] = alpha ? 1 : -1;
    SPS4(n, phi, 1, alpha, 0, deg, A, P, (1<<omega)-1);
    Long ht=0;
    for(i=0; i<=deg/2; i++){
        ht = ht>A[i] ? ht : A[i];
        ht = ht>-A[i] ? ht : -A[i];
    }
    printf(" & %lld\n",ht);
    free(A);
}
'''

'''
int main(nargs, args)
int nargs;
char *args[];
{
    Long n, alpha;
    if(nargs != 3){
        printf("Usage : ./a.out n l\n");
        printf("n : squarefree, positive odd number\n");
        printf("l : an integer\n");
        printf("Computes $\\Phi_n(z)$ if l is nonzero, and $\\Psi_n(z)$ otherwise.\n");
        return 0;
    }
    sscanf(args[1],"%lld",&n);
    sscanf(args[nargs-1], "%d", &alpha);
    cyclotomic(n,alpha);
    return 0;
}
'''
