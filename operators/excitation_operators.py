#
#Excitation operators
#
#Normal-ordered strings containig second-quantized operators
#to act on the bra and ket reference wavefunction, respectively
#
#SMT 2023-2024
#


from sympy.physics.secondquant import (AntiSymmetricTensor, F, Fd, NO)
from sympy import (symbols)


def ti_ia():
    """
    Singles target indices i,a for amplitude equations
    """

    i = symbols('i', below_fermi=True)
    a = symbols('a', above_fermi=True)
    ia = NO(Fd(i)*F(a))
    return ia, [[i], [a]]


def ti_bj():
    """
    Singles target indices j,b for ISR
    """

    j = symbols('j', below_fermi=True)
    b = symbols('b', above_fermi=True)
    bj = NO(Fd(b)*F(j))
    return bj, [[j], [b]]


def ti_bckj():
    """
    Doubles target indices for the sd block
    """

    j,k = symbols('j,k', below_fermi=True)
    b,c = symbols('b,c', above_fermi=True)
    bckj = NO(Fd(b)*Fd(c)*F(k)*F(j))
    return bckj, [[j, k], [b, c]]


def ti_bcdlkj():
    """
    Triples target indices for the st block
    """

    j,k,l = symbols('j,k,l', below_fermi=True)
    b,c,d = symbols('b,c,d', above_fermi=True)
    bcdlkj = NO(Fd(b)*Fd(c)*Fd(d)*F(l)*F(k)*F(j))
    return bcdlkj, [[j, k, l], [b, c, d]]


def ti_ijba():
    """
    Doubles target indices i,j,a,b for amplitude equations
    """

    i, j = symbols('i,j', below_fermi=True)
    a, b = symbols('a,b', above_fermi=True)
    ijba = NO(Fd(i)*Fd(j)*F(b)*F(a))
    return ijba, [[i,j], [a,b]]


def ti_ck():
    """
    Singles target indices for the ds block
    """

    k = symbols('k', below_fermi=True)
    c = symbols('c', above_fermi=True)
    ck = NO(Fd(c)*F(k))
    return ck, [[k], [c]]


def ti_cdlk():
    """
    Doubles target indices k,l,c,d for ISR
    """

    k ,l = symbols('k,l', below_fermi=True)
    c, d = symbols('c,d', above_fermi=True)
    cdlk = NO(Fd(c)*Fd(d)*F(l)*F(k))
    return cdlk, [[k,l], [c,d]]


def ti_ijkcba():
    """
    Triples target indices i,j,k,a,b,c for amplitude equations
    """

    i, j, k = symbols('i,j,k', below_fermi=True)
    a, b, c = symbols('a,b,c', above_fermi=True)
    ijkcba = NO(Fd(i)*Fd(j)*Fd(k)*F(c)*F(b)*F(a))
    return ijkcba, [[i,j,k], [a,b,c]]


def ti_dl():
    """
    Singles target indices for the ts block
    """

    l = symbols('l', below_fermi=True)
    d = symbols('d', above_fermi=True)
    dl = NO(Fd(d)*F(l))
    return dl, [[l], [d]]


def ti_cdemlk():
    """
    Triples target indices for the dt block
    """

    k,l,m = symbols('k,l,m', below_fermi=True)
    c,d,e = symbols('c,d,e', above_fermi=True)
    cdemlk = NO(Fd(c)*Fd(d)*Fd(e)*F(m)*F(l)*F(k))
    return cdemlk, [[k,l,m], [c,d,e]]


def ti_deml():
    """
    Doubles target indices for the td block
    """

    l,m = symbols('l,m', below_fermi=True)
    d,e = symbols('d,e', above_fermi=True)
    deml = NO(Fd(d)*Fd(e)*F(m)*F(l))
    return deml, [[l,m], [d,e]]


def ti_defnml():
    """
    Triples target indices for the tt block
    """

    l,m,n = symbols('l,m,n', below_fermi=True)
    d,e,f = symbols('d,e,f', above_fermi=True)
    defnml = NO(Fd(d)*Fd(e)*Fd(f)*F(n)*F(m)*F(l))
    return defnml, [[l,m,n], [d,e,f]]


lhs_excitation_dict = {'s': ti_ia, 'd': ti_ijba, 't': ti_ijkcba,
        'ss': ti_ia, 'sd': ti_ia, 'st': ti_ia,
        'ds': ti_ijba, 'dd': ti_ijba, 'dt': ti_ijba,
        'ts': ti_ijkcba, 'td': ti_ijkcba, 'tt': ti_ijkcba
        }


rhs_excitation_dict = {'ss': ti_bj, 'sd': ti_bckj, 'st': ti_bcdlkj,
        'ds': ti_ck, 'dd': ti_cdlk, 'dt': ti_cdemlk,
        'ts': ti_dl, 'td': ti_deml, 'tt': ti_defnml
        }


def return_excitation(exc_type):
    if exc_type not in excitation_dict:
        raise ValueError('No valid input for calculation: {exc_type}')
    
    #Initialize the lhs and rhs expressions as well the corresponding
    #lists of target indices: first list for occupied indices,
    #second one for virtuals
    Lhs = Rhs = 1
    Lhs_ti = Rhs_ti = [[], []]
    
    if exc_type == '00':
        pass
    elif exc_type in ['s','d','t']:
        Lhs = lhs_excitation_dict[exc_type]()[0]
        Lhs_ti = excitation_dict[exc_type]()[1]
    else:
        Lhs = lhs_excitation_dict[exc_type]()[0]
        Lhs_ti = lhs_excitation_dict[exc_type]()[1]
        Rhs = rhs_excitation_dict[exc_type]()[0]
        Rhs_ti = rhs_excitation_dict[exc_type]()[1]

    return Lhs, Rhs, Lhs_ti, Rhs_ti
