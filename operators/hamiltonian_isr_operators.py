#Hamiltonian operators
#SMT 5/2023


#TODO: write proper comments!!


from sympy.physics.secondquant import (AntiSymmetricTensor, F, Fd, NO)
from sympy import (symbols, Rational, Dummy)


def f_int():
    """
    HF orbital enegry contributions
    """

    p, q = symbols('p,q', cls=Dummy)
    f = AntiSymmetricTensor('f', (p,), (q,))
    pq = NO(Fd(p)*F(q))
    F_int = f*pq
    return F_int


def vo_int():
    """
    exchange enegry contributions
    """

    p, q, r, s = symbols('p,q,r,s', cls=Dummy)
    v = AntiSymmetricTensor('i', (p, q), (r, s))
    pqsr = NO(Fd(p)*Fd(q)*F(s)*F(r))
    V_int = Rational(1, 4)*v*pqsr
        
    return V_int

def vn_int():
    """
    exchange enegry contributions
    """

    i, j = symbols('i,j', below_fermi=True, cls=Dummy)
    a, b = symbols('a,b', above_fermi=True, cls=Dummy)
    v_abij = AntiSymmetricTensor('i', (a, b), (i, j)) 
    abji = NO(Fd(a)*Fd(b)*F(j)*F(i))
    v_ijab = AntiSymmetricTensor('i', (i, j), (a, b)) 
    ijba = NO(Fd(i)*Fd(j)*F(b)*F(a))
    V_int = Rational(1, 4)*(v_abij*abji + v_ijab*ijba)

    return V_int


def vr_int():
    """
    exchange enegry contributions
    """

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
    V_int = vo_int()
    F_V_int = F_int + V_int
    return F_V_int


def d_op():
    """
    One-particle property d
    """

    p, q = symbols('p,q', cls=Dummy)
    d = AntiSymmetricTensor('d', (p,), (q,))
    pq = NO(Fd(p)*F(q))
    D_op = d*pq
    return D_op


hamiltonian_isr_dict = {'f': f_int, 'vo': vo_int, 'vn': vn_int,
        'vr': vr_int, 'fv': f_v_int, 'd': d_op}


def return_ham_isr_fragment(ham_isr_type):
    if ham_isr_type not in hamiltonian_isr_dict:
        raise ValueError(f'The hamiltonian/ISR has not been set up correctly!\n'
                'Key does not exist!')
    return hamiltonian_isr_dict[ham_isr_type]()
