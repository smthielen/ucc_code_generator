#
#main routine to evaluate the Hamiltonian
#
#SMT 5/2023
#


#TODO: Check dependancies!!


from sympy.physics.secondquant import (AntiSymmetricTensor, KroneckerDelta, wicks,
        NO, evaluate_deltas, substitute_dummies, Commutator,
        simplify_index_permutations, PermutationOperator, Symbol, Add, Mul) #simplify_index_permutations and PermutationOperator ??

from sympy import (symbols, Rational, latex, Dummy, Float)

from sympy.core import numbers

from base.commutator import *

from operators.excitation_operators import return_excitation
from operators.hamiltonian_isr_operators import return_ham_isr_fragment
from operators.amplitude_operators import return_ampl_operator

from utils.latex import write_latex
from utils.pretty_dummies import return_pretty_dummies
from utils.term_list import build_term_list


#canonical_order_dict = {
#    'occ': 'ijklmno',
#    'vir': 'abcdefg'
#}


#final_substitute_dict = {
#    'old_e': 'klmno',
#    'new_e': 'jklmn',
#    'old_t2': 'klmno',
#    'new_t2': 'jklmn'
#}


def check_order(obj, order): 
    """
    Checks order of a term
    """

    proper = 0

    if (obj == 0):
        print('No contribution')
    else:
        for x in obj.args:
            term_order = 0
            for y in x.args:
                if isinstance(y, AntiSymmetricTensor):
                    if (str(y.symbol) == 'i' or str(y.symbol) == 't2' or str(y.symbol) == 'tc2'):
                        term_order += 1
                    if (str(y.symbol) == 't1' or str(y.symbol) == 'tc1' or str(y.symbol) == 't3' or str(y.symbol) == 'tc3'):
                        term_order += 2
                if isinstance(y, Symbol):
                    if (str(y) == 'i'):
                        term_order += 1
            if (term_order <= order):
                proper += x

    return proper


def check_commutator_type(nested_comm, comm_type, domain):
    """
    Checks if a term is not allowed due to n or r operation
    """

    proper_comm = 0

    exc_lvl = 0
    if (domain == 'uccsd'):
        #print('Domain is UCCSD')
        exc_lvl = [1.0, 2.0]
    elif (domain == 'uccsdt'):
        exc_lvl = [1.0, 2.0, 3.0]

    #print('exc_lvl=')
    #print(exc_lvl)

    for x in nested_comm.args:
        num_creators = 0
        num_annihilators = 0
        for y in x.args:
            if (isinstance(y, NO)):
                num_creators = len(list(y.iter_q_creators()))
                num_annihilators = len(list(y.iter_q_annihilators()))
        if (num_creators == 0 and num_annihilators == 0):
            #This case could give a hint to disconnected diagrams, Must be taken care of this??
            print('Is this case actually to be taken care of??')
            print('a')
        elif ((num_creators == 0 and num_annihilators != 0) or (num_creators != 0 and num_annihilators == 0)):
            if (comm_type == 'nd' and ((0.5*num_creators in exc_lvl) or (0.5*num_annihilators in exc_lvl))):
                proper_comm += x
                print('aa')
            elif (comm_type == 'rest' and ((0.5*num_creators) not in exc_lvl) and (num_creators != 0)):
                proper_comm += x
                print('aaa')
            elif (comm_type == 'rest' and ((0.5*num_annihilators) not in exc_lvl) and (num_annihilators != 0)):
                proper_comm += x
                print('aaaa')
        
        elif (comm_type == 'rest'):
            proper_comm += x

    return proper_comm


def evaluate_hamiltonian_isr(hamiltonian_isr, calc_type):
    """
    Evaluation of the Hamiltonian
    """

    equation = 0
    lhs = return_excitation(calc_type)[0]
    rhs = return_excitation(calc_type)[1]
    lhs_ti = return_excitation(calc_type)[2]
    rhs_ti = return_excitation(calc_type)[3]
    print("\nlhs,rhs:")
    print(lhs, rhs)
    print("\nlhs_ti,rhs_ti:")
    print(lhs_ti, rhs_ti)
    
    pretty_dummies_dict = return_pretty_dummies(calc_type)

    for obj in hamiltonian_isr.comm_list:

        #lhs, rhs = return_excitation(calc_type)
        #print("\nlhs,rhs:")
        #print(lhs, rhs)
    
        #pretty_dummies_dict = return_pretty_dummies(calc_type)

        #First evaluate terms with no commuataor involved

        print(type(obj))
        if (isinstance(obj, no_commutator)):
            print('22')
            print("\nEvaluating term with no commutator...")
            #Include the factor of 1.0 to ensure a Mul object is formed
            expr = 1.0 * wicks(lhs*return_ham_isr_fragment(obj.arg)*rhs,
                    simplify_kronecker_deltas=True, keep_only_fully_contracted=True)
            equation += substitute_dummies(expr, new_indices=True,
                    pretty_indices=pretty_dummies_dict)
            print("\n...done!")

        #Now evaluate all terms that include commutators
        if (isinstance(obj, commutator)):
            print("\nEvaluating term with a single commutator...")
            comm = wicks(Commutator(return_ham_isr_fragment(obj.arg1),
                return_ampl_operator(obj.arg2, hamiltonian_isr.only_real)))
            expr = obj.prefactor*wicks(lhs*comm*rhs, simplify_kronecker_deltas=True,
                    keep_only_fully_contracted=True)
            equation += substitute_dummies(expr, new_indices=True,
                    pretty_indices=pretty_dummies_dict)
            print("\n...done!")
        
        #Finally, evaluate nested commutators
        if (isinstance(obj, nested_commutator)):
            print("\nEvaluating a nested commutator...")
            
            #Nested commutator intermediate storage
            nested_comm = 0

            i = 0
            #print("NESTED COMMUTATOR WITH OBJECTS:")
            #print(len(x.arg_list))
            
            while (i < (len(obj.arg_list))):

                #print('i='+str(i))
                #print('NESTED COMM:')
                #print(nested_comm)
                
                if (i == 0):
                    print("\nNesting lvl 0")
                    nested_comm = wicks(Commutator(return_ham_isr_fragment(obj.arg_list[i]),
                            return_ampl_operator(obj.arg_list[i+1], hamiltonian_isr.only_real)))
                    i+=2
                    
                else:
                    print("\nNesting lvl "+str(i-1))
                    #print('enter')
                    if (obj.comm_type_list[i-2] == 'ord'): 
                        nested_comm = wicks(Commutator(nested_comm,
                                return_ampl_operator(obj.arg_list[i], hamiltonian_isr.only_real)))
                        
                    else: 
                        nested_comm = check_commutator_type(nested_comm,
                                obj.comm_type_list[i-2], hamiltonian_isr.domain)
                        nested_comm = wicks(Commutator(nested_comm,
                                return_ampl_operator(obj.arg_list[i], hamiltonian_isr.only_real)))

                    i+=1
            expr = obj.prefactor*wicks(lhs*nested_comm*rhs, simplify_kronecker_deltas=True,
                           keep_only_fully_contracted=True)
            equation += substitute_dummies(expr, new_indices=True,
                    pretty_indices=pretty_dummies_dict)
            print("\n...done!")

    #print('pert_order:')
    #print(hamiltonian_isr.pertubation_order)
    
    #Check order of all terms, i. e. that no term with a higher order
    #than specified is taken into the final expression
    equation = check_order(equation, hamiltonian_isr.pertubation_order)
    

    #TODO: Does that work this way?!
    if (equation != 0):
        equation = substitute_dummies(equation, new_indices=True, 
                pretty_indices=pretty_dummies_dict)
    
    #print('\n',equation)

    #eq = build_term_list(equation)
    
    #print('\neq:')
    #for i in eq:
    #    print(i.tensors)
    #    print(i.indices)
    #    print(i.prefactor)
    #print(perm_eq)

    #write_latex(calc_type, dict_eq, perm_eq)

    return equation, lhs_ti, rhs_ti, pretty_dummies_dict
