#
#main routine to evaluate the Hamiltonian
#
#SMT 2023-2024
#


from sympy.physics.secondquant import (NO, Commutator, wicks,
        AntiSymmetricTensor, substitute_dummies)
from base.commutator import *
from operators.excitation_operators import return_excitation
from operators.hamiltonian_isr_operators import return_ham_isr_fragment
from operators.amplitude_operators import return_ampl_operator
from utils.pretty_dummies import return_pretty_dummies

#TODO:needed??
#from sympy.core import numbers


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
                    
                    #ERI and t2 amplitudes are first-order properties
                    if (str(y.symbol) == 'i' or str(y.symbol) == 't2' or 
                            str(y.symbol) == 'tc2'):
                        term_order += 1
                    
                    #t1 and t3 amplitudes are second-order properties
                    if (str(y.symbol) == 't1' or str(y.symbol) == 'tc1' or
                            str(y.symbol) == 't3' or str(y.symbol) == 'tc3'):
                        term_order += 2
                
                #TODO: Is this condition needed after all?
                #if isinstance(y, Symbol):
                #    if (str(y) == 'i'):
                #        term_order += 1
                #TODO

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
        exc_lvl = [1.0, 2.0]
    elif (domain == 'uccsdt'):
        exc_lvl = [1.0, 2.0, 3.0]

    for x in nested_comm.args:
        num_creators = 0
        num_annihilators = 0
        
        for y in x.args:
            if (isinstance(y, NO)):
                num_creators = len(list(y.iter_q_creators()))
                num_annihilators = len(list(y.iter_q_annihilators()))
       
        #Fully contracted expressions always vanish except
        #for the last commutator which is always of "ord" type.
        #Therefore, this type of expressions may be discarded here.
        if (num_creators == 0 and num_annihilators == 0):
            pass

        #At first, check for diagrans that have only open lines above or beneath
        #all vertices
        elif ((num_creators == 0 and num_annihilators != 0) or 
                (num_creators != 0 and num_annihilators == 0)):
            
            #Now, in the non-diagonal case, check if the excitation level of the
            #diagram is matching the underlying UCC domain
            if (comm_type == 'nd' and ((0.5*num_creators in exc_lvl) or
                    (0.5*num_annihilators in exc_lvl))):
                proper_comm += x

            #If the excitation level is not in the corresponding UCC domain,
            #this is contributing to the rest part
            elif (comm_type == 'rest' and ((0.5*num_creators) not in exc_lvl) and 
                    (num_creators != 0)):
                proper_comm += x
            elif (comm_type == 'rest' and ((0.5*num_annihilators) not in exc_lvl) and
                    (num_annihilators != 0)):
                proper_comm += x
        
        #Diagrams that have open lines above and beneath all vertices are
        #by definition always belonging to the rest part
        elif (comm_type == 'rest'):
            proper_comm += x

    return proper_comm


def evaluate(braket, calc_type):
    """
    Evaluation of the Hamiltonian
    """

    equation = 0
    pretty_dummies_dict = return_pretty_dummies(calc_type)

    #Get lhs and rhs excitation operators
    lhs = return_excitation(calc_type)[0]
    rhs = return_excitation(calc_type)[1]
    
    #Get also lhs and rhs target indices
    lhs_ti = return_excitation(calc_type)[2]
    rhs_ti = return_excitation(calc_type)[3]
    
    for obj in braket.comm_list:

        #First evaluate terms with no commuataor involved
        if (isinstance(obj, no_commutator)):
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
                    return_ampl_operator(obj.arg2, braket.only_real)))
            expr = obj.pref*wicks(lhs*comm*rhs, simplify_kronecker_deltas=True,
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
            
            while (i < (len(obj.arg_list))):
                if (i == 0):
                    print("\nNesting lvl 0")
                    nested_comm = wicks(Commutator(
                            return_ham_isr_fragment(obj.arg_list[i]),
                            return_ampl_operator(obj.arg_list[i+1],
                            braket.only_real)))
                    i+=2
                else:
                    print("\nNesting lvl "+str(i-1))
                    
                    #ordinary operator, no need for term checkup
                    if (obj.comm_type_list[i-2] == 'ord'): 
                        nested_comm = wicks(Commutator(nested_comm,
                                return_ampl_operator(obj.arg_list[i],
                                braket.only_real)))
                    
                    #ND or R part, check for corresponding contributions!
                    else: 
                        nested_comm = check_commutator_type(nested_comm,
                                obj.comm_type_list[i-2], braket.domain)
                        nested_comm = wicks(Commutator(nested_comm,
                                return_ampl_operator(obj.arg_list[i],
                                braket.only_real)))
                    i+=1

            expr = obj.pref*wicks(lhs*nested_comm*rhs,
                    simplify_kronecker_deltas=True,
                    keep_only_fully_contracted=True)
            equation += substitute_dummies(expr, new_indices=True,
                    pretty_indices=pretty_dummies_dict)
            print("\n...done!")

    #Check order of all terms, i. e. that no term with a higher order
    #than specified is taken into the final expression
    equation = check_order(equation, braket.pertubation_order)
    
    #Now do a final substitution of dummy indices
    if (equation != 0):
        equation = substitute_dummies(equation, new_indices=True, 
                pretty_indices=pretty_dummies_dict)
    
    return equation, lhs_ti, rhs_ti, pretty_dummies_dict
