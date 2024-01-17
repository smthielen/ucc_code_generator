#
#main file to start calculations
#
#SMT 5/2023
#

import numpy as np
import time

#For gs energy contribution
from sympy.physics.secondquant import (AntiSymmetricTensor, KroneckerDelta, wicks,
        NO, evaluate_deltas, substitute_dummies, Commutator,
        simplify_index_permutations, PermutationOperator, Symbol, Add, Mul) #simplify_index_permutations and PermutationOperator ??

#TODO: maybe substitute * here by only explicitly needed stuff?!
from base.hamiltonian import *
from base.isr import *
from base.evaluate import *

from utils.pretty_dummies import * 
from utils.term_list import *
from utils.assume_diag_fock import *
from utils.permutations_finder import *

from operators.excitation_operators import return_excitation
from operators.hamiltonian_isr_operators import *
from operators.amplitude_operators import return_ampl_operator

#True: only real cluster amplitudes
#UCC4 full hamiltonian + UCC4+s/t additional terms (from UCC5)
#hamiltonian_list = [hamiltonian_ucc4_Hbar01(False), hamiltonian_ucc4_Hbar2(False), hamiltonian_ucc4_Hbar3(False), hamiltonian_ucc5_Hbar2(False)]
hamiltonian_list = [hamiltonian_ucc3(False)]


#TODO:Hacked!! Do this properly!!
#diag_block = ['ss','dd','tt']

#Calculation type!!
calc_type = 'ss'

#valid inputs: see in the excitaion_dict in operators/excitation_operators.py
def main():

    for h in hamiltonian_list:

        #Main calculation
        ts = time.time()
        eq, lti, rti, pdd = evaluate_hamiltonian_isr(h, calc_type)
        te = time.time()
        print("\nMAIN CALCULATION TOOK TIME:")
        print(float(te-ts))


        occs = rti[0]
        virs = rti[1]

        occs = []
        virs = []

        print('permutations o/v:')
        print(occs,virs)


        #Permutations finder
        ts = time.time()
        perms = find_permutations(eq, occs, virs, pdd)
        te = time.time()
        print("\nPERMUTATIONS FINDER TOOK TIME:")
        print(float(te-ts))

        #term list
        ls = build_term_list(eq, perms)
        save_term_list(h, calc_type, ls, eq, float(te-ts))
        write_latex(h, calc_type, ls)
    
        '''
        vals = []
        for v1 in perms.values():
            for v2 in v1:
                for v3 in v2:
                    vals.append(v3)

        print(vals)

        for ii, i in enumerate(vals):
            for j in vals[ii+1:]:
                if (i == j):
                    print("SAME SAME!!")

        '''

        '''
        filter_out = []
        for key in perms:
            filter_out.append(perms[key])


        proper=[]
        for i in filter_out:
            for j in i:
                for k in j:
                    proper.append(k)


        eq_ = 0
        for i, ii in enumerate(eq.args):
            if i not in proper:
                eq_ += ii


        eq_ = substitute_dummies(-1.0 *eq_, new_indices=True, pretty_indices=pdd)
        eq_ = substitute_dummies(-1.0 *eq_, new_indices=True, pretty_indices=pdd)
        
        perms2 = find_permutations(eq_, occs, virs, pdd)
        #perms3 = find_permutations(eq, occs, virs, pdd)
        #perms4 = find_permutations(eq, occs, virs, pdd)
        #perms5 = find_permutations(eq, occs, virs, pdd)
        '''
        print('############################\n')
        for i, ii in enumerate(eq.args):
            print(i)
            print(ii)
    
        
        print(perms)
        '''
        for i, ii in enumerate(eq_.args):
            print(i)
            print(ii)
        print(perms2)
    
        #print(perms3)
        #print(perms4)
        #print(perms5)
        #for i in perms:
        #    print(i)

        #for i in p_eq.args:
        #    print(i)
        
        #Build term list
        #ts = time.time()
        #ls,ps = build_term_list(eq, [])
        #te = time.time()
        #print("\nTERM LIST BUILDER TOOK TIME:")
        #print(float(te-ts))
        #save_term_list(h, calc_type, ls, eq, float(te-ts))
        

        #i = lti[0][0] 
        #j = lti[0][1] 
        #a = lti[1][0] 
        #b = lti[1][1]

        #print(i)
        #print(j)
        #print(a)
        #print(b)

        #po = PermutationOperator(i,j,a)
        #pv = PermutationOperator(a,b)
        '''
        
'''
        lhs, rhs = return_excitation(calc_type)
        gs_e = wicks(lhs*gs_e*rhs,
                simplify_kronecker_deltas=True, keep_only_fully_contracted=True)
        if (gs_e != 0):
            gs_e = substitute_dummies(gs_e, new_indices=True,
                    pretty_indices=dummies)
        
        te=time.time()
        print("\nGS ENERGY CONTRIBUTION TOOK TIME:")
        print(float(te-ts))
        mvp_ = build_term_list(gs_e) 
        print("\nmvp:")
        print_term_list(mvp_)
    else:
        print("\nNO GS ENERGY CONTRIBUTION IN OFF-DIAGONAL BLOCK!!")

    #now, set explicit pertubation order for scheme; for >2, all gs energy contributions
    #will be discarded
    #h.pertubation_order = pert_order
 

    ts = time.time()
    eq = evaluate_hamiltonian(h, calc_type)
    eq = substitute_dummies(eq, new_indices=True,
            pretty_indices=dummies)
    res = eq - gs_e
    res = substitute_dummies(res, new_indices=True,
            pretty_indices=dummies)
    final_eq = build_term_list(res) 
    
    final_eq_df = return_diag_fock(final_eq)
    print("\nSECULAR MATRIX FOR THE "+str(calc_type)+" BLOCK IN ORDER "+str(pert_order))
    print("\nUSING HAMILTONIAN:"+str(type(h)))

    print_term_list(final_eq_df)
    te=time.time()
    print("\nSECULAR MATRIX TOOK TIME:")
    print(float(te-ts))
'''

if __name__ == "__main__":
    print('Starting (U)CC code generator ;)')
    main()
