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

'''
class no_commutator:
    """
    Auxiliary class to handle expressions without a commutator 
    """

    def __init__(self, arg):
        self.arg = arg


class commutator:
    """
    Class for commutator expressions
    """

    def __init__(self, arg1 , arg2, prefactor):
        self.arg1 = arg1
        self.arg2 = arg2
        self.prefactor = prefactor


class nested_commutator:
    """
    Class for nested commutator expressions
    """

    def __init__(self, arg_list, comm_type_list, prefactor):
        self.arg_list = arg_list
        self.comm_type_list = comm_type_list   
        self.prefactor = prefactor

        #TODO: do this properly
        if ((len(self.arg_list)-2) != len(self.comm_type_list)):
            print('Incompatible setup for nested commutator!!')
'''


class hamiltonian_fv_uccsdt():
    """
    UCCSD F+V Hamiltonian up to 1st order
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsdt'
        self.pertubation_order = 1
        self.comm_list = []
   
        #UCC Hbar_0 Terms F + V
        self.comm_list.append(no_commutator('fv', 1.))
        #UCC Hbar_1 first term [F, Sigma]
        self.comm_list.append(commutator('f', 's123', 1.))


class hamiltonian_fv_uccsd():
    """
    UCCSD F+V Hamiltonian up to 1st order
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsd'
        self.pertubation_order = 1
        self.comm_list = []
   
        #UCC Hbar_0 Terms F + V
        self.comm_list.append(no_commutator('fv', 1.))
        #UCC Hbar_1 first term [F, Sigma]
        self.comm_list.append(commutator('f', 's12', 1.))


class hamiltonian_2nd_uccsdt():
    """
    UCC2 Hamiltonian
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsdt'
        self.pertubation_order = 2
        self.comm_list = []
   
        #UCC Hbar_0 Terms F + V
        self.comm_list.append(no_commutator('f', 1.))
        #UCC Hbar_1 first term [F, Sigma]
        #self.comm_list.append(commutator('f', 's123', 1.))
        #UCC Hbar_1 second term 0.5*[V, Sigma]
        #self.comm_list.append(commutator('vo', 's123', 0.5))
        #UCC Hbar_1 third term 0.5*[V_R, Sigma]
        #self.comm_list.append(commutator('vr', 's123', 0.5))


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
        
        #UCC Hbar_2 first term 1./12.*[[V_N, Sigma], Sigma]
        #self.comm_list.append(nested_commutator(['vn', 's1', 's2'], ['ord'], 1./12.))
        #self.comm_list.append(nested_commutator(['vn', 's2', 's1'], ['ord'], 1./12.))
        #self.comm_list.append(nested_commutator(['vn', 's2', 's2'], ['ord'], 1./12.))
        #self.comm_list.append(nested_commutator(['vn', 's2', 's3'], ['ord'], 1./12.))
        #self.comm_list.append(nested_commutator(['vn', 's3', 's2'], ['ord'], 1./12.))
        #UCC Hbar_2 second term 1./4.*[[V, Sigma]_R, Sigma]
        #self.comm_list.append(nested_commutator(['vo', 's1', 's2'], ['rest'], 1./4.))
        #self.comm_list.append(nested_commutator(['vo', 's2', 's1'], ['rest'], 1./4.))
        #self.comm_list.append(nested_commutator(['vo', 's2', 's2'], ['rest'], 1./4.))
        #self.comm_list.append(nested_commutator(['vo', 's2', 's3'], ['rest'], 1./4.))
        #self.comm_list.append(nested_commutator(['vo', 's3', 's2'], ['rest'], 1./4.))
        #UCC Hbar_2 third term 1./4.*[[V_R, Sigma]_R, Sigma]
        #self.comm_list.append(nested_commutator(['vr', 's1', 's2'], ['rest'], 1./4.))
        #self.comm_list.append(nested_commutator(['vr', 's2', 's1'], ['rest'], 1./4.))
        #self.comm_list.append(nested_commutator(['vr', 's2', 's2'], ['rest'], 1./4.))
        #self.comm_list.append(nested_commutator(['vr', 's2', 's3'], ['rest'], 1./4.))
        #self.comm_list.append(nested_commutator(['vr', 's3', 's2'], ['rest'], 1./4.))


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
        
        #UCC Hbar_2 first term 1./12.*[[V_N, Sigma], Sigma]
        #self.comm_list.append(nested_commutator(['vn', 's1', 's1'], ['ord'], 1./12.))
        #self.comm_list.append(nested_commutator(['vn', 's1', 's3'], ['ord'], 1./12.))
        #self.comm_list.append(nested_commutator(['vn', 's3', 's1'], ['ord'], 1./12.))
        #self.comm_list.append(nested_commutator(['vn', 's3', 's3'], ['ord'], 1./12.))
        #UCC Hbar_2 second term 1./4.*[[V, Sigma]_R, Sigma]
        #self.comm_list.append(nested_commutator(['vo', 's1', 's1'], ['rest'], 1./4.))
        #self.comm_list.append(nested_commutator(['vo', 's1', 's3'], ['rest'], 1./4.))
        #self.comm_list.append(nested_commutator(['vo', 's3', 's1'], ['rest'], 1./4.))
        #self.comm_list.append(nested_commutator(['vo', 's3', 's3'], ['rest'], 1./4.))
        #UCC Hbar_2 third term 1./4.*[[V_R, Sigma]_R, Sigma]
        #self.comm_list.append(nested_commutator(['vr', 's1', 's1'], ['rest'], 1./4.))
        #self.comm_list.append(nested_commutator(['vr', 's1', 's3'], ['rest'], 1./4.))
        #self.comm_list.append(nested_commutator(['vr', 's3', 's1'], ['rest'], 1./4.))
        #self.comm_list.append(nested_commutator(['vr', 's3', 's3'], ['rest'], 1./4.))


class hamiltonian_uccsdt_sd():
    """
    Part of the UCC4 Hamiltonian containing one singles and one doubles amplitude
    """

    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsd'
        self.pertubation_order = 3
        self.comm_list = []


class hamiltonian_1st_order_uccsdt():
    """
    UCC2 Hamiltonian
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsdt'
        self.pertubation_order = 1
        self.comm_list = []
   
        #UCC H_0 Terms F + V
        self.comm_list.append(no_commutator('fv'))
        #UCC H_1 first term [F, Sigma]
        self.comm_list.append(commutator('f', 's123', 1.))


class hamiltonian_2nd_order_uccsdt():
    """
    UCC2 Hamiltonian
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsdt'
        self.pertubation_order = 2
        self.comm_list = []
   
        #UCC H_0 Terms F + V
        self.comm_list.append(no_commutator('fv'))
        #UCC H_1 first term [F, Sigma]
        self.comm_list.append(commutator('f', 's123', 1.))
        #UCC H_1 second term 0.5*[V, Sigma]
        self.comm_list.append(commutator('vo', 's123', 0.5))
        #UCC H_1 third term 0.5*[V_R, Sigma]
        self.comm_list.append(commutator('vr', 's123', 0.5))


class ham_only_4th_e():
    """
    UCC2 Hamiltonian
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsdt'
        self.pertubation_order = 4
        self.comm_list = []

        self.comm_list.append(nested_commutator(['vn', 's2', 's2', 's2'], ['ord','rest'], 1./24.))
        self.comm_list.append(nested_commutator(['vr', 's2', 's2', 's2'], ['rest','rest'], 1./8.))
        self.comm_list.append(nested_commutator(['vo', 's2', 's2', 's2'], ['rest','rest'], 1./8.))
        self.comm_list.append(nested_commutator(['vo', 's2', 's2', 's2'], ['rest','ord'], -1./24.))
        self.comm_list.append(nested_commutator(['vr', 's2', 's2', 's2'], ['rest','ord'], -1./24.))


class hamiltonian_3rd_order_uccsdt():
    """
    UCC3 Hamiltonian
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        self.domain = 'uccsdt'
        self.pertubation_order = 3
        self.comm_list = []
    
        #UCC H_0 Terms F + V
        self.comm_list.append(no_commutator('fv'))
        #UCC H_1 first term [F, Sigma]
        self.comm_list.append(commutator('f', 's123', 1.))
        #UCC H_1 second term 0.5*[V, Sigma]
        self.comm_list.append(commutator('vo', 's123', 0.5))
        #UCC H_1 third term 0.5*[V_R, Sigma]
        self.comm_list.append(commutator('vr', 's123', 0.5))
        
        #UCC H_2 first term 1./12.*[[V_N, Sigma], Sigma]
        self.comm_list.append(nested_commutator(['vn', 's123', 's123'], ['ord'], 1./12.))
        self.comm_list.append(nested_commutator(['vo', 's123', 's123'], ['rest'], 1./4.))
        self.comm_list.append(nested_commutator(['vr', 's123', 's123'], ['rest'], 1./4.))


class hamiltonian_ucc4():
    """
    UCC4 Hamiltonian
    """
    
    def __init__(self, only_real):
        self.only_real = only_real
        
        self.domain = 'uccsdt'
        self.pertubation_order = 4
        self.comm_list = []
    

               
        #UCC H_0 Terms F + V
        self.comm_list.append(no_commutator(f_v))
        #UCC H_1 first term [F, Sigma]
        self.comm_list.append(commutator(f, s123_1, 1.))
        #UCC H_1 second term 0.5*[V, Sigma]
        self.comm_list.append(commutator(v_o_1, s123_2, 0.5))
        #UCC H_1 third term 0.5*[V_R, Sigma]
        self.comm_list.append(commutator(v_r_1, s123_3, 0.5))
    
        #UCC H_2
        #self.comm_list.append(nested_commutator([v_n_1, s123_4, s123_5], ['ord'], 1./12.))
        #self.comm_list.append(nested_commutator([v_o_2, s123_6, s123_7], ['rest'], 1./4.))
        #self.comm_list.append(nested_commutator([v_r_2, s123_8, s123_9], ['rest'], 1./4.))

        #UCC H_3
        #self.comm_list.append(nested_commutator([v_n_, s123_10, s123_11, s123_12], ['ord','rest'], 1./24.))
        self.comm_list.append(nested_commutator([v_n_2, s1_1, s2_1], ['ord'], 1./24.))
        self.comm_list.append(nested_commutator([v_n_3, s2_2, s1_2], ['ord'], 1./24.))
        self.comm_list.append(nested_commutator([v_n_4, s2_3, s2_4, s2_5], ['ord','rest'], 1./24.))
        self.comm_list.append(nested_commutator([v_n_5, s2_6, s3_1], ['ord'], 1./24.))
        self.comm_list.append(nested_commutator([v_n_6, s3_2, s2_7], ['ord'], 1./24.))
        #self.comm_list.append(nested_commutator([v_r__, s123_13, s123_14, s123_15], ['rest','rest'], 1./8.))
        self.comm_list.append(nested_commutator([v_r_3, s1_3, s2_8], ['rest'], 1./8.))
        self.comm_list.append(nested_commutator([v_r_4, s2_9, s1_4], ['rest'], 1./8.))
        self.comm_list.append(nested_commutator([v_r_5, s2_10, s2_11, s2_12], ['rest','rest'], 1./8.))
        self.comm_list.append(nested_commutator([v_r_6, s2_13, s3_3], ['rest'], 1./8.))
        self.comm_list.append(nested_commutator([v_r_7, s3_4, s2_14], ['rest'], 1./8.))
        #self.comm_list.append(nested_commutator([v_o__, s123_16, s123_17, s123_18], ['rest','rest'], 1./8.))
        self.comm_list.append(nested_commutator([v_o_3, s1_5, s2_15], ['rest'], 1./8.))
        self.comm_list.append(nested_commutator([v_o_4, s2_16, s1_6], ['rest'], 1./8.))
        
        self.comm_list.append(nested_commutator([v_o_5, s2_17, s2_18, s2_19], ['rest','rest'], 1./8.))
        
        self.comm_list.append(nested_commutator([v_o_6, s2_20, s3_5], ['rest'], 1./8.))
        self.comm_list.append(nested_commutator([v_o_7, s3_6, s2_21], ['rest'], 1./8.))
        #self.comm_list.append(nested_commutator([v_o___, s123_19, s123_20, s123_21], ['rest','ord'], -1./24.))
        self.comm_list.append(nested_commutator([v_o_8, s1_7, s2_22], ['rest'], -1./24.))
        self.comm_list.append(nested_commutator([v_o_9, s2_23, s1_8], ['rest'], -1./24.))
        self.comm_list.append(nested_commutator([v_o_10, s2_24, s2_25, s2_26], ['rest','ord'], -1./24.))
        self.comm_list.append(nested_commutator([v_o_11, s2_27, s3_7], ['rest'], -1./24.))
        self.comm_list.append(nested_commutator([v_o_12, s3_8, s2_28], ['rest'], -1./24.))
        #self.comm_list.append(nested_commutator([v_r___, s123_22, s123_23, s123_24], ['rest','ord'], -1./24.))
        self.comm_list.append(nested_commutator([v_r_8, s1_9, s2_29], ['rest'], -1./24.))
        self.comm_list.append(nested_commutator([v_r_9, s2_30, s1_10], ['rest'], -1./24.))
        self.comm_list.append(nested_commutator([v_r_10, s2_31, s2_32, s2_33], ['rest','ord'], -1./24.))
        self.comm_list.append(nested_commutator([v_r_11, s2_34, s3_9], ['rest'], -1./24.))
        self.comm_list.append(nested_commutator([v_r_12, s3_10, s2_35], ['rest'], -1./24.))
'''
#class hamiltonian_ccsd():
#    """
#    CCSD Hamiltonian
#    """
#    
#    comm_list = []
#
#    #UCC H_0 Terms F + V
#    comm_list.append(no_commutator(f_v_int('n')))
#    #UCC H_1 first term [F, Sigma]
#    comm_list.append(commutator(sigma12('n', True), sigma1('n', True), True, True, 2))
#    #UCC H_1 second term 0.5[V, Sigma]
#    comm_list.append(commutator(sigma12('n', True), sigma1('n', True), True, True, 2))
#    #UCC H_1 second term 0.5[V, Sigma]
#    comm_list.append(commutator(sigma12('n', True), sigma1('n', True), True, True, 2))
#
#    def __init__(self):
#        pass
'''
