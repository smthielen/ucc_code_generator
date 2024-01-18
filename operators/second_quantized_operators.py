#
#NOTE:From an older version, not used anymore. Kept only for completeness.
#
#All second-quantized operators
#
#SMT 5/2023
#


from sympy.physics.secondquant import (AntiSymmetricTensor, F, Fd, NO)
from sympy import (symbols, Rational, Dummy)


def sigma1(is_ucc, only_real):
    """
    Singles amplitude operator
    - is_ucc: UCC t1 amplitude operator if True, normal CC one if False
    - only_real: relevant for UCC, if True only real cluster amplitudes
    """
    
    i = symbols('i', below_fermi=True, cls=Dummy)
    a = symbols('a', above_fermi=True, cls=Dummy)
    sigma_ai = AntiSymmetricTensor('t1', (a,), (i,))
    sigma_ia = AntiSymmetricTensor('tc1', (i,), (a,))
    ai = NO(Fd(a)*F(i))
    ia = NO(Fd(i)*F(a))

    if is_ucc:
        if only_real:
            Sigma1 = sigma_ai*ai - sigma_ai*ia
        else:
            Sigma1 = sigma_ai*ai + sigma_ia*ia
    else:
        Sigma1 = sigma_ai*ai

    return Sigma1
   

def sigma2(is_ucc, only_real):
    """
    Doubles amplitude operator
    - is_ucc: UCC t2 amplitude operator if True, normal CC one if False
    - only_real: relevant for UCC, if True only real cluster amplitudes
    """

    i, j = symbols('i,j', below_fermi=True, cls=Dummy)
    a, b = symbols('a,b', above_fermi=True, cls=Dummy)
    sigma_abij = AntiSymmetricTensor('t2', (a, b), (i, j))
    sigma_ijab = AntiSymmetricTensor('tc2', (i, j), (a, b))
    abji = NO(Fd(a)*Fd(b)*F(j)*F(i))
    ijba = NO(Fd(i)*Fd(j)*F(b)*F(a))

    if is_ucc:
        if only_real:
            Sigma2 = Rational(1, 4)*sigma_abij*abji - Rational(1, 4)*sigma_abij*ijba
        else:
            Sigma2 = Rational(1, 4)*sigma_abij*abji + Rational(1, 4)*sigma_ijab*ijba
    else:
        Sigma2 = Rational(1, 4)*sigma_abij*abji
    
    return Sigma2


def sigma3(is_ucc, only_real):
    """
    Triples amplitude operator
    - is_ucc: UCC t3 amplitude operator if True, normal CC one if False
    - only_real: relevant for UCC, if True only real cluster amplitudes
    """

    i, j, k = symbols('i,j,k', below_fermi=True, cls=Dummy)
    a, b, c = symbols('a,b,c', above_fermi=True, cls=Dummy)
    sigma_abcijk = AntiSymmetricTensor('t3', (a, b, c), (i, j, k))
    sigma_ijkabc = AntiSymmetricTensor('tc3', (i, j, k), (a, b, c))
    abckji = NO(Fd(a)*Fd(b)*Fd(c)*F(k)*F(j)*F(i))
    ijkcba = NO(Fd(i)*Fd(j)*Fd(k)*F(c)*F(b)*F(a))

    if is_ucc:
        if only_real:
            Sigma3 = Rational(1, 36)*sigma_abcijk*abckji - Rational(1, 36)*sigma_abcijk*ijkcba
        else:
            Sigma3 = Rational(1, 36)*sigma_abcijk*abckji + Rational(1, 36)*sigma_ijkabc*ijkcba
    else:
        Sigma3 = Rational(1, 36)*sigma_abcijk*abckji
    
    return Sigma3


def sigma12(is_ucc, only_real):
    """
    Combined singles and doubles amplitude operator
    - is_ucc: UCC t1 and t2 amplitude operators if True, normal CC ones if False
    - only_real: relevant for UCC, if True only real cluster amplitudes
    """
    
    Sigma1 = sigma1(is_ucc, only_real)
    Sigma2 = sigma2(is_ucc, only_real)
    Sigma = Sigma1 + Sigma2
    return Sigma


def sigma123(is_ucc, only_real):
    """
    Combined singles, doubles and triples amplitude operator
    - is_ucc: UCC t1, t2 and t3 amplitude operators if True, normal CC ones if False
    - only_real: relevant for UCC, if True only real cluster amplitudes
    """
    
    Sigma1 = sigma1(is_ucc, only_real)
    Sigma2 = sigma2(is_ucc, only_real)
    Sigma3 = sigma3(is_ucc, only_real)
    Sigma = Sigma1 + Sigma2 + Sigma3
    return Sigma


def f_int():
    """
    HF orbital enegry contributions
    """

    p, q = symbols('p,q', cls=Dummy)
    f = AntiSymmetricTensor('f', (p,), (q,))
    pq = NO(Fd(p)*F(q))
    F_int = f*pq
    return F_int


def v_int(op_type):
    """
    exchange enegry contributions
    """

    #Rewrite this check properly....
    if(op_type != 'ord' and op_type != 'nd' and op_type != 'rest'):
        print('Incompatible setup of v_int!!')

    if (op_type == 'ord'):
        p, q, r, s = symbols('p,q,r,s', cls=Dummy)
        v = AntiSymmetricTensor('i', (p, q), (r, s))
        pqsr = NO(Fd(p)*Fd(q)*F(s)*F(r))
        V_int = Rational(1, 4)*v*pqsr

    if (op_type == 'nd'):
        i, j = symbols('i,j', below_fermi=True, cls=Dummy)
        a, b = symbols('a,b', above_fermi=True, cls=Dummy)
        v_abij = AntiSymmetricTensor('i', (a, b), (i, j)) 
        abji = NO(Fd(a)*Fd(b)*F(j)*F(i))
        v_ijab = AntiSymmetricTensor('i', (i, j), (a, b)) 
        ijba = NO(Fd(i)*Fd(j)*F(b)*F(a))
        V_int = Rational(1, 4)*(v_abij*abji + v_ijab*ijba)

    if (op_type == 'rest'):
        i, j, k, l = symbols('i,j,k,l', below_fermi=True, cls=Dummy)
        a, b, c, d = symbols('a,b,c,d', above_fermi=True, cls=Dummy)
        v_abcd = AntiSymmetricTensor('i', (a, b), (c, d)) 
        abdc = NO(Fd(a)*Fd(b)*F(d)*F(c))
        v_ijkl = AntiSymmetricTensor('i', (i, j), (k, l)) 
        ijlk = NO(Fd(i)*Fd(j)*F(l)*F(k))
        v_iabj = AntiSymmetricTensor('i', (i, a), (b, j)) 
        iajb = NO(Fd(i)*Fd(a)*F(j)*F(b))
        v_aibc = AntiSymmetricTensor('i', (a, i), (b, c)) 
        aicb = NO(Fd(a)*Fd(i)*F(c)*F(b))
        v_ijka = AntiSymmetricTensor('i', (i, j), (k, a)) 
        ijak = NO(Fd(i)*Fd(j)*F(a)*F(k))
        v_abci = AntiSymmetricTensor('i', (a, b), (c, i)) 
        abic = NO(Fd(a)*Fd(b)*F(i)*F(c))
        v_iajk = AntiSymmetricTensor('i', (i, a), (j, k)) 
        iakj = NO(Fd(i)*Fd(a)*F(k)*F(j))
        V_int = Rational(1, 4)*(v_abcd*abdc + v_ijkl*ijlk) + v_iabj*iajb + Rational(1, 2)*(v_aibc*aicb + v_ijka*ijak + v_abci*abic + v_iajk*iakj)
        
    return V_int


def f_v_int():
    """
    Combination of F and V for H_0
    """

    F_int = f_int()
    V_int = v_int('ord')
    F_V_int = F_int + V_int
    return F_V_int


def ti_ia():
    """
    Singles target indices i,a for amplitude equations
    """

    i = symbols('i', below_fermi=True)
    a = symbols('a', above_fermi=True)
    ia = NO(Fd(i)*F(a))
    return ia


def ti_bj():
    """
    Singles target indices j,b for ISR
    """

    l = symbols('l', below_fermi=True)
    d = symbols('d', above_fermi=True)
    bj = NO(Fd(d)*F(l))
    return bj


def ti_ijba():
    """
    Doubles target indices i,j,a,b for amplitude equations
    """

    i, j = symbols('i,j', below_fermi=True)
    a, b = symbols('a,b', above_fermi=True)
    ijba = NO(Fd(i)*Fd(j)*F(b)*F(a))
    return ijba


def ti_cdlk():
    """
    Doubles target indices k,l,c,d for ISR
    """

    k ,l = symbols('k,l', below_fermi=True)
    c, d = symbols('c,d', above_fermi=True)
    cdlk = NO(Fd(c)*Fd(d)*F(l)*F(k))
    return cdlk


def ti_ijkcba():
    """
    Triples target indices i,j,k,a,b,c for amplitude equations
    """

    i, j, k = symbols('i,j,k', below_fermi=True)
    a, b, c = symbols('a,b,c', above_fermi=True)
    ijkcba = NO(Fd(i)*Fd(j)*Fd(k)*F(c)*F(b)*F(a))
    return ijkcba
