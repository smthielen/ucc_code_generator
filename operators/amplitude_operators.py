#
#Amplitude operators
#
#->CC amplitudes contain "cc" in the name, i. e. "s2cc" would be
#normal Coupled Cluster while "s2" would be UCC
#The corresponding bool is set in the return_ampl_operator function
#
#SMT 2023-2024
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
    sigma_ia = AntiSymmetricTensor('tc1', (a,), (i,))
    ai = NO(Fd(a)*F(i))
    ia = NO(Fd(i)*F(a))

    if is_ucc:
        if only_real:
            Sigma1 = sigma_ai*ai - sigma_ai*ia
        else:
            Sigma1 = sigma_ai*ai - sigma_ia*ia
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
    sigma_ijab = AntiSymmetricTensor('tc2', (a, b), (i, j))
    abji = NO(Fd(a)*Fd(b)*F(j)*F(i))
    ijba = NO(Fd(i)*Fd(j)*F(b)*F(a))

    if is_ucc:
        if only_real:
            Sigma2 = Rational(1, 4)*sigma_abij*abji - Rational(1, 4)*sigma_abij*ijba
        else:
            Sigma2 = Rational(1, 4)*sigma_abij*abji - Rational(1, 4)*sigma_ijab*ijba
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
    sigma_ijkabc = AntiSymmetricTensor('tc3', (a, b, c), (i, j, k))
    abckji = NO(Fd(a)*Fd(b)*Fd(c)*F(k)*F(j)*F(i))
    ijkcba = NO(Fd(i)*Fd(j)*Fd(k)*F(c)*F(b)*F(a))

    if is_ucc:
        if only_real:
            Sigma3 = Rational(1, 36)*sigma_abcijk*abckji - Rational(1, 36)*sigma_abcijk*ijkcba
        else:
            Sigma3 = Rational(1, 36)*sigma_abcijk*abckji - Rational(1, 36)*sigma_ijkabc*ijkcba
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


def sigma13(is_ucc, only_real):
    """
    Combined singles and triples amplitude operator
    - is_ucc: UCC t1 and t3 amplitude operators if True, normal CC ones if False
    - only_real: relevant for UCC, if True only real cluster amplitudes
    """
    
    Sigma1 = sigma1(is_ucc, only_real)
    Sigma3 = sigma3(is_ucc, only_real)
    Sigma = Sigma1 + Sigma3
    return Sigma


amplitude_dict = {'s1': sigma1, 's2': sigma2, 's3': sigma3,
        's1cc': sigma1, 's2cc': sigma2, 's3cc': sigma3,
        's12': sigma12, 's13': sigma13, 's123': sigma123,
        's12cc': sigma12, 's13cc': sigma13, 's123cc': sigma123}


def return_ampl_operator(amp_type, only_real):
    if amp_type not in amplitude_dict:
        raise ValueError('Amplitude operator key does not exist!')
    if not isinstance(only_real, bool): 
        raise TypeError(f'Specification needs to be of boolean type!')
    
    is_ucc = True
    if (amp_type.find('cc') != -1):
        is_ucc = False

    return amplitude_dict[amp_type](is_ucc, only_real)
