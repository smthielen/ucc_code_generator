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
#SMT 2023-2024
#
#TODO: write proper comments!!


from base.commutator import *


class hamiltonian_f():
    """
    just f...
    """
    
    def __init__(self, only_real):
        #only_real doesnt matter here....
        self.only_real = only_real
        #domain doesnt matter here....
        self.domain = 'uccsdt'
        self.pertubation_order = 0
        self.comm_list = []
    
        #UCC Hbar_0 Term F
        self.comm_list.append(no_commutator('f', 1.))


class hamiltonian_ucc2():
    """
    UCC2 Hamiltonian
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsd'
        self.pertubation_order = 2
        self.comm_list = []
   
        #UCC Hbar_0 Terms F + V
        self.comm_list.append(no_commutator('fv', 1.))
        #UCC Hbar_1 first term [F, Sigma]
        self.comm_list.append(commutator('f', 's12', 1.))
        #UCC Hbar_1 second term 0.5*[V, Sigma]
        self.comm_list.append(commutator('vo', 's12', 0.5))
        #UCC Hbar_1 third term 0.5*[V_R, Sigma]
        self.comm_list.append(commutator('vr', 's12', 0.5))


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
        self.comm_list.append(no_commutator('fv', 1.0))
        #UCC Hbar_1 first term [F, Sigma]
        self.comm_list.append(commutator('f', 's12', 1.0))
        #UCC Hbar_1 second term 0.5*[V, Sigma]
        self.comm_list.append(commutator('vo', 's12', 0.50))
        #UCC Hbar_1 third term 0.5*[V_R, Sigma]
        self.comm_list.append(commutator('vr', 's12', 0.50))
        
        #UCC Hbar_2 first term 1./12.*[[V_N, Sigma], Sigma]
        self.comm_list.append(nested_commutator(['vn', 's12', 's12'], ['ord'], 1.0/12.0))
        #UCC Hbar_2 second term 1./4.*[[V, Sigma]_R, Sigma]
        self.comm_list.append(nested_commutator(['vo', 's12', 's12'], ['rest'], 1.0/4.0))
        #UCC Hbar_2 third term 1./4.*[[V_R, Sigma]_R, Sigma]
        self.comm_list.append(nested_commutator(['vr', 's12', 's12'], ['rest'], 1.0/4.0))


class hamiltonian_ucc4_Hbar01_1st():
    """
    UCC4 version of the Hamiltonian for Hbar_0 + 1st order Hbar_1
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsdt'
        self.pertubation_order = 1
        self.comm_list = []
    
        #UCC Hbar_0 Terms F + V
        self.comm_list.append(no_commutator('fv', 1.))
        #UCC Hbar_1 first term [F, Sigma]
        self.comm_list.append(commutator('f', 's2', 1.))


class hamiltonian_ucc4_Hbar01_2nd():
    """
    UCC4 version of the Hamiltonian for Hbar_0 + 2nd order Hbar_1
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsdt'
        self.pertubation_order = 2
        self.comm_list = []
    
        #UCC Hbar_0 Terms F + V
        self.comm_list.append(no_commutator('fv', 1.))
        #UCC Hbar_1 first term [F, Sigma]
        self.comm_list.append(commutator('f', 's123', 1.))
        #UCC Hbar_1 second term 0.5*[V, Sigma]
        self.comm_list.append(commutator('vo', 's2', 0.5))
        #UCC Hbar_1 third term 0.5*[V_R, Sigma]
        self.comm_list.append(commutator('vr', 's2', 0.5))


class hamiltonian_ucc4_Hbar01():
    """
    UCC4 version of the Hamiltonian for Hbar_1 (terms up to 3rd order)
    Including the trivial Hbar_0
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsdt'
        self.pertubation_order = 3
        self.comm_list = []
    
        #UCC Hbar_0 Terms F + V
        self.comm_list.append(no_commutator('fv', 1.))
        #UCC Hbar_1 first term [F, Sigma]
        self.comm_list.append(commutator('f', 's123', 1.))
        #UCC Hbar_1 second term 0.5*[V, Sigma]
        self.comm_list.append(commutator('vo', 's123', 0.5))
        #UCC Hbar_1 third term 0.5*[V_R, Sigma]
        self.comm_list.append(commutator('vr', 's123', 0.5))
       

class hamiltonian_ucc4_Hbar2():

    """
    Hbar_2 part of the UCC4 Hamiltonian
    3rd order contributions from eri*s2*s2, 4th order contributions from eri*s1*s2 and eri*s2*s3
    Linearity of the commutator has been used to split up terms and reduce computational cost
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsdt'
        self.pertubation_order = 4
        self.comm_list = []

        #UCC Hbar_2 first term 1./12.*[[V_N, Sigma], Sigma]
        self.comm_list.append(nested_commutator(['vn', 's123', 's2'], ['ord'], 1./12.))
        self.comm_list.append(nested_commutator(['vn', 's2', 's13'], ['ord'], 1./12.))
        #UCC Hbar_2 second term 1./4.*[[V, Sigma]_R, Sigma]
        self.comm_list.append(nested_commutator(['vo', 's123', 's2'], ['rest'], 1./4.))
        self.comm_list.append(nested_commutator(['vo', 's2', 's13'], ['rest'], 1./4.))
        #UCC Hbar_2 third term 1./4.*[[V_R, Sigma]_R, Sigma]
        self.comm_list.append(nested_commutator(['vr', 's123', 's2'], ['rest'], 1./4.))
        self.comm_list.append(nested_commutator(['vr', 's2', 's13'], ['rest'], 1./4.))
        

class hamiltonian_ucc4_Hbar2_onlyT2():

    """
    Hbar_2 part of the UCC4 Hamiltonian with only the 3rd order contributions from eri*s2*s2
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsdt'
        self.pertubation_order = 3
        self.comm_list = []

        #UCC Hbar_2 first term 1./12.*[[V_N, Sigma], Sigma]
        self.comm_list.append(nested_commutator(['vn', 's2', 's2'], ['ord'], 1./12.))
        #UCC Hbar_2 second term 1./4.*[[V, Sigma]_R, Sigma]
        self.comm_list.append(nested_commutator(['vo', 's2', 's2'], ['rest'], 1./4.))
        #UCC Hbar_2 third term 1./4.*[[V_R, Sigma]_R, Sigma]
        self.comm_list.append(nested_commutator(['vr', 's2', 's2'], ['rest'], 1./4.))


class hamiltonian_ucc4_Hbar3():

    """
    Hbar_3 part of the UCC4 Hamiltonian
    Only 4th order contributions from eri*s2*s2*s2
    Linearity of the commutator has been used to split up terms and reduce computational cost
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsdt'
        self.pertubation_order = 4
        self.comm_list = []

        #UCC Hbar_3 first term 1./24.*[[[V_N, Sigma], Sigma]_R, Sigma]
        self.comm_list.append(nested_commutator(['vn', 's2', 's2', 's2'], ['ord','rest'], 1./24.))
        #UCC Hbar_3 second term 1./8.*[[[V_R, Sigma]_R, Sigma]_R, Sigma]
        self.comm_list.append(nested_commutator(['vr', 's2', 's2', 's2'], ['rest','rest'], 1./8.))
        #UCC Hbar_3 third term 1./8.*[[[V, Sigma]_R, Sigma]_R, Sigma]
        self.comm_list.append(nested_commutator(['vo', 's2', 's2', 's2'], ['rest','rest'], 1./8.))
        #UCC Hbar_3 fourth term -1./24.*[[[V, Sigma]_R, Sigma], Sigma]
        self.comm_list.append(nested_commutator(['vo', 's2', 's2', 's2'], ['rest','ord'], -1./24.))
        #UCC Hbar_3 fifth term -1./24.*[[[V_R, Sigma]_R, Sigma], Sigma]
        self.comm_list.append(nested_commutator(['vr', 's2', 's2', 's2'], ['rest','ord'], -1./24.))


class hamiltonian_ucc5_Hbar2():

    """
    Hbar_2 part of the UCC5 Hamiltonian
    5th order contributions from eri*s1*s1, eri*s1*s3, eri*s3*s1 and eri*s3*s3
    Linearity of the commutator has been used to split up terms and reduce computational cost
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsdt'
        self.pertubation_order = 5
        self.comm_list = []

        #UCC Hbar_2 first term 1./12.*[[V_N, Sigma], Sigma]
        self.comm_list.append(nested_commutator(['vn', 's13', 's13'], ['ord'], 1./12.))
        #UCC Hbar_2 second term 1./4.*[[V, Sigma]_R, Sigma]
        self.comm_list.append(nested_commutator(['vo', 's13', 's13'], ['rest'], 1./4.))
        #UCC Hbar_2 third term 1./4.*[[V_R, Sigma]_R, Sigma]
        self.comm_list.append(nested_commutator(['vr', 's13', 's13'], ['rest'], 1./4.))
