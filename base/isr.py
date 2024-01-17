#
#Hamiltonian setup
#
#->Specification for ordinary coupled cluster amplitudes (no
#  de-excitations) ca be done by adding 'cc' at the end of the 
#  amplitude operator string
#->Additionally, the only_real keyword specifies whether
#  UCC cluster amplitudes should only be real. When initializing
#  a CC Hamiltonian, this key must be also set, but is irrelevant
#  in the computation
#
#SMT 5/2023
#


#TODO: write proper comments!!


from base.commutator import *


class isr_ucc2():
    """
    UCC2 ISR
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsd'
        self.pertubation_order = 2
        self.comm_list = []
   
        #ISR zeroth-order D
        self.comm_list.append(no_commutator('d', 1.))
        #ISR first-order [D, Sigma]
        self.comm_list.append(commutator('d', 's12', 1.))
        #ISR second-order 0.5*[[D, Sigma], Sigma]
        self.comm_list.append(nested_commutator(['d', 's12', 's12'], ['ord'], 0.5))


class isr_ucc3():
    """
    UCC3 ISR
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsd'
        self.pertubation_order = 3
        self.comm_list = []
   
        #ISR zeroth-order D
        self.comm_list.append(no_commutator('d', 1.))
        #ISR first-order [D, Sigma]
        self.comm_list.append(commutator('d', 's12', 1.))
        #ISR second-order 0.5*[[D, Sigma], Sigma]
        self.comm_list.append(nested_commutator(['d', 's12', 's12'], ['ord'], 0.5))
        #ISR third-order 1/6*[[[D, Sigma], Sigma], Sigma]
        self.comm_list.append(nested_commutator(['d', 's2', 's2', 's2'], ['ord', 'ord'], 1./6.))


class isr_d():
    """
    only d....
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsdt'
        self.pertubation_order = 0
        self.comm_list = []
   
        #ISR zeroth-order D
        self.comm_list.append(no_commutator('d', 1.))


class isr_ucc4_1st():
    """
    UCC4 ISR up to 1st order
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsdt'
        self.pertubation_order = 1
        self.comm_list = []
   
        #ISR zeroth-order D
        self.comm_list.append(no_commutator('d', 1.))
        #ISR singly-nested [D, Sigma] up to 2nd order
        self.comm_list.append(commutator('d', 's2', 1.))


class isr_ucc4_2nd():
    """
    UCC4 ISR up to 2nd order
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsdt'
        self.pertubation_order = 2
        self.comm_list = []
   
        #ISR zeroth-order D
        self.comm_list.append(no_commutator('d', 1.))
        #ISR singly-nested [D, Sigma] up to 2nd order
        self.comm_list.append(commutator('d', 's123', 1.))
        #ISR doubly-nested 0.5*[[D, Sigma], Sigma] up to 4 order
        self.comm_list.append(nested_commutator(['d', 's2', 's2'], ['ord'], 0.5))


class isr_ucc4_3rd():
    """
    UCC4 ISR up to 3rd order
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsdt'
        self.pertubation_order = 3
        self.comm_list = []
   
        #ISR zeroth-order D
        self.comm_list.append(no_commutator('d', 1.))
        #ISR singly-nested [D, Sigma] up to 2nd order
        self.comm_list.append(commutator('d', 's123', 1.))
        #ISR doubly-nested 0.5*[[D, Sigma], Sigma] up to 3rd order
        self.comm_list.append(nested_commutator(['d', 's13', 's2'], ['ord'], 0.5))
        self.comm_list.append(nested_commutator(['d', 's2', 's13'], ['ord'], 0.5))
        #ISR triply-nested 1/6*[[[D, Sigma], Sigma], Sigma] up to 3rd order
        self.comm_list.append(nested_commutator(['d', 's2', 's2', 's2'], ['ord', 'ord'], 1./6.))


class isr_ucc4_double():
    """
    UCC4 ISR, the doubly nested commutators
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsdt'
        self.pertubation_order = 4
        self.comm_list = []
   
        #ISR zeroth-order D
        self.comm_list.append(no_commutator('d', 1.))
        #ISR singly-nested [D, Sigma] up to 2nd order
        self.comm_list.append(commutator('d', 's123', 1.))
        #ISR doubly-nested 0.5*[[D, Sigma], Sigma] up to 4th order
        self.comm_list.append(nested_commutator(['d', 's123', 's123'], ['ord'], 0.5))


class isr_ucc4_triple():
    """
    UCC4 ISR, the triply nested commutators are of the form s2*s2*s2 (3rd order),
    s1*s2*S2 and s3*s2*s2 (fourth order)
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsdt'
        self.pertubation_order = 4
        self.comm_list = []
   
        #ISR zeroth-order D
        self.comm_list.append(no_commutator('d', 1.))
        #ISR singly-nested [D, Sigma] up to 2nd order
        self.comm_list.append(commutator('d', 's123', 1.))
        #ISR doubly-nested 0.5*[[D, Sigma], Sigma] up to 4 order
        self.comm_list.append(nested_commutator(['d', 's123', 's123'], ['ord'], 0.5))
        #ISR triply-nested 1/6*[[[D, Sigma], Sigma], Sigma] up to 4th order
        self.comm_list.append(nested_commutator(['d', 's2', 's2', 's2'], ['ord', 'ord'], 1./6.))
        self.comm_list.append(nested_commutator(['d', 's13', 's2', 's2'], ['ord', 'ord'], 1./6.))
        self.comm_list.append(nested_commutator(['d', 's2', 's13', 's2'], ['ord', 'ord'], 1./6.))
        self.comm_list.append(nested_commutator(['d', 's2', 's2', 's13'], ['ord', 'ord'], 1./6.))


class isr_ucc4_quadruple():
    """
    UCC4 ISR, the quadruply nested commutator
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsdt'
        self.pertubation_order = 4
        self.comm_list = []
   
        #ISR quadruply-nested 1/24*[[[[D, Sigma], Sigma], Sigma], Sigma]
        self.comm_list.append(nested_commutator(['d', 's2', 's2', 's2', 's2'], ['ord', 'ord', 'ord'], 1./24.))


'''
class hamiltonian_ucc3():
    """
    UCC3 Hamiltonian
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsd'
        self.pertubation_order = 3
        self.comm_list = []
    
        #UCC Hbar_0 Terms F + V
        self.comm_list.append(no_commutator('fv'))
        #UCC Hbar_1 first term [F, Sigma]
        self.comm_list.append(commutator('f', 's12', 1.))
        #UCC Hbar_1 second term 0.5*[V, Sigma]
        self.comm_list.append(commutator('vo', 's12', 0.5))
        #UCC Hbar_1 third term 0.5*[V_R, Sigma]
        self.comm_list.append(commutator('vr', 's12', 0.5))
        
        #UCC Hbar_2 first term 1./12.*[[V_N, Sigma], Sigma]
        self.comm_list.append(nested_commutator(['vn', 's12', 's12'], ['ord'], 1./12.))
        #UCC Hbar_2 second term 1./4.*[[V, Sigma]_R, Sigma]
        self.comm_list.append(nested_commutator(['vo', 's12', 's12'], ['rest'], 1./4.))
        #UCC Hbar_2 third term 1./4.*[[V_R, Sigma]_R, Sigma]
        self.comm_list.append(nested_commutator(['vr', 's12', 's12'], ['rest'], 1./4.))

'''
