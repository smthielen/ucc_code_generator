#
#routine to search a/b and i/j permutations in an equation
#TODO: apply also to a/b/c and i/j/k !!
#
#SMT 8/2023
#


import numpy
from sympy.physics.secondquant import (AntiSymmetricTensor, KroneckerDelta,
        substitute_dummies,
        simplify_index_permutations, PermutationOperator, Symbol, Add, Mul) #simplify_index_permutations and PermutationOperator ??


def return_perms(o, v):
    """
    Return perms according to target indices (in a very dirty manner....).
    """
    P = PermutationOperator
    
    perms = []
    perm_strings = []

    if (len(o) == 2 and len(v) == 2):
        #ijab for doubles
        perms.append([[P(o[0], o[1])], [P(v[0], v[1])], [P(o[0], o[1]), P(v[0], v[1])]]) 
        perm_strings.append(str(o[0])+str(o[1])+str(v[0])+str(v[1]))

    if (len(o) == 2):
        #ij for doubles
        perms.append([[P(o[0], o[1])]])
        perm_strings.append(str(o[0])+str(o[1]))

    if (len(v) == 2):
        #ab for doubles
        perms.append([[P(v[0], v[1])]]) 
        perm_strings.append(str(v[0])+str(v[1]))

    #TODO: This is such a dirty way to do it....maybe do mor generally?!
    if (len(o) == 3 and len(v) == 3):
        perms.append([[P(o[0], o[1])], [P(o[0], o[2])], [P(o[1], o[2])], 
                      [P(o[0], o[1]), P(o[0], o[2])], [P(o[0], o[2]), P(o[0], o[1])], 
                      [P(v[0], v[1])], [P(v[0], v[2])], [P(v[1], v[2])], 
                      [P(v[0], v[1]), P(v[0], v[2])], [P(v[0], v[2]), P(v[0], v[1])],
                      [P(o[0], o[1]), P(v[0], v[1])], [P(o[0], o[2]), P(v[0], v[1])], 
                      [P(o[1], o[2]), P(v[0], v[1])],
                      [P(o[0], o[1]), P(o[0], o[2]), P(v[0], v[1])], 
                      [P(o[0], o[2]), P(o[0], o[1]), P(v[0], v[1])],
                      [P(o[0], o[1]), P(v[0], v[2])], [P(o[0], o[2]), P(v[0], v[2])],
                      [P(o[1], o[2]), P(v[0], v[2])], 
                      [P(o[0], o[1]), P(o[0], o[2]), P(v[0], v[2])],
                      [P(o[0], o[2]), P(o[0], o[1]), P(v[0], v[2])], 
                      [P(o[0], o[1]), P(v[1], v[2])], [P(o[0], o[2]), P(v[1], v[2])],
                      [P(o[1], o[2]), P(v[1], v[2])], 
                      [P(o[0], o[1]), P(o[0], o[2]), P(v[1], v[2])],
                      [P(o[0], o[2]), P(o[0], o[1]), P(v[1], v[2])], 
                      [P(o[0], o[1]), P(v[0], v[1]), P(v[0], v[2])],
                      [P(o[0], o[2]), P(v[0], v[1]), P(v[0], v[2])], 
                      [P(o[1], o[2]), P(v[0], v[1]), P(v[0], v[2])], 
                      [P(o[0], o[1]), P(o[0], o[2]), P(v[0], v[1]), P(v[0], v[2])], 
                      [P(o[0], o[2]), P(o[0], o[1]), P(v[0], v[1]), P(v[0], v[2])],
                      [P(o[0], o[1]), P(v[0], v[2]), P(v[0], v[1])],
                      [P(o[0], o[2]), P(v[0], v[2]), P(v[0], v[1])], 
                      [P(o[1], o[2]), P(v[0], v[2]), P(v[0], v[1])], 
                      [P(o[0], o[1]), P(o[0], o[2]), P(v[0], v[2]), P(v[0], v[1])], 
                      [P(o[0], o[2]), P(o[0], o[1]), P(v[0], v[2]), P(v[0], v[1])]])
        perm_strings.append(str(o[0])+str(o[1])+str(o[2])+str(v[0])+str(v[1])+str(v[2]))

    if (len(o) == 3):
        #ijk for triples
        perms.append([[P(o[0], o[1])], [P(o[0], o[2])], [P(o[1], o[2])], 
            [P(o[0], o[1]), P(o[0], o[2])], [P(o[0], o[2]), P(o[0], o[1])]])
        perm_strings.append(str(o[0])+str(o[1])+str(o[2])) 
    
    if (len(v) == 3):
        #abc for triples
        perms.append([[P(v[0], v[1])], [P(v[0], v[2])], [P(v[1], v[2])], 
            [P(v[0], v[1]), P(v[0], v[2])], [P(v[0], v[2]), P(v[0], v[1])]])
        perm_strings.append(str(v[0])+str(v[1])+str(v[2])) 

    if (len(o) == 3 and len(v) == 3):
        #ijab for triplets
        perms.append([[P(o[0], o[1])], [P(v[0], v[1])], [P(o[0], o[1]), P(v[0], v[1])]]) 
        perm_strings.append(str(o[0])+str(o[1])+str(v[0])+str(v[1]))
        #ikab for triplets
        perms.append([[P(o[0], o[2])], [P(v[0], v[1])], [P(o[0], o[2]), P(v[0], v[1])]]) 
        perm_strings.append(str(o[0])+str(o[2])+str(v[0])+str(v[1]))
        #jkab for triplets
        perms.append([[P(o[1], o[2])], [P(v[0], v[1])], [P(o[1], o[2]), P(v[0], v[1])]]) 
        perm_strings.append(str(o[1])+str(o[2])+str(v[0])+str(v[1]))
        #ijac for triplets
        perms.append([[P(o[0], o[1])], [P(v[0], v[2])], [P(o[0], o[1]), P(v[0], v[2])]]) 
        perm_strings.append(str(o[0])+str(o[1])+str(v[0])+str(v[2]))
        #ikac for triplets
        perms.append([[P(o[0], o[2])], [P(v[0], v[2])], [P(o[0], o[2]), P(v[0], v[2])]]) 
        perm_strings.append(str(o[0])+str(o[2])+str(v[0])+str(v[2]))
        #jkac for triplets
        perms.append([[P(o[1], o[2])], [P(v[0], v[2])], [P(o[1], o[2]), P(v[0], v[2])]]) 
        perm_strings.append(str(o[1])+str(o[2])+str(v[0])+str(v[2]))
        #ijbc for triplets
        perms.append([[P(o[0], o[1])], [P(v[1], v[2])], [P(o[0], o[1]), P(v[1], v[2])]]) 
        perm_strings.append(str(o[0])+str(o[1])+str(v[1])+str(v[2]))
        #ikbc for triplets
        perms.append([[P(o[0], o[2])], [P(v[1], v[2])], [P(o[0], o[2]), P(v[1], v[2])]]) 
        perm_strings.append(str(o[0])+str(o[2])+str(v[1])+str(v[2]))
        #jkbc for triplets
        perms.append([[P(o[1], o[2])], [P(v[1], v[2])], [P(o[1], o[2]), P(v[1], v[2])]]) 
        perm_strings.append(str(o[1])+str(o[2])+str(v[1])+str(v[2]))

    if (len(o) == 3):
        #ij for triples
        perms.append([[P(o[0], o[1])]])
        perm_strings.append(str(o[0])+str(o[1]))
        #ik for triples
        perms.append([[P(o[0], o[2])]])
        perm_strings.append(str(o[0])+str(o[2]))
        #jk for triples
        perms.append([[P(o[1], o[2])]])
        perm_strings.append(str(o[1])+str(o[2]))

    if (len(v) == 3):
        #ab for triples
        perms.append([[P(v[0], v[1])]]) 
        perm_strings.append(str(v[0])+str(v[1]))
        #ac for triples
        perms.append([[P(v[0], v[2])]]) 
        perm_strings.append(str(v[0])+str(v[2]))
        #bc for triples
        perms.append([[P(v[1], v[2])]]) 
        perm_strings.append(str(v[1])+str(v[2]))

    return perms, perm_strings


def check_same_type(expr1, expr2):
    """
    Checks if two expressions in the sympy Mul object contain the same
    tensors. This is a requirement for the permutation of target indices!
    """
    t1 = []
    t2 = []

    #TODO: Maybe also check the prefactor here?! Or not deeded I guess....
    
    for i in expr1.args:
        if isinstance(i, KroneckerDelta):
            t1.append('d')
        elif isinstance(i, AntiSymmetricTensor):
            t1.append(str(i.symbol))
            print(t1)

    for i in expr2.args:
        if isinstance(i, KroneckerDelta):
            t2.append('d')
        elif isinstance(i, AntiSymmetricTensor):
            t2.append(str(i.symbol))

    if (t1 == t2):
        return True
    else:
        return False


def find_permutations(eq, occs, virs, pdd):
    """
    Find permutations of target indices pairs in lti and rti.
    Use the pretty_dummy_dict pdd to substitute dummies IN EACH STEP!
    Save found permutations in the permutation_map.
    """

    #All possible permutations
    perms, perm_strings = return_perms(occs, virs)
    permutations_dict = {}

    for s in perm_strings:
        permutations_dict[s] = []

    #Array to check already found terms
    print(eq)

    if (eq == 0):
        print("0: NOTHING TO LOOK FOR PERMUTATIONS!!!!")
        return permutations_dict

    already_found_term = numpy.zeros(len(eq.args), dtype=bool)

    #go over all permutations as they are defined in the perms array
    #the ordering is crucial since the individual permutations have
    #common operations on the same terms. The permutations are checked
    #in their order of complexity (i. e. number of terms summarized).
    for np, permutation in enumerate(perms):
   
        print("PERMUTATION:")
        print(permutation)

        #loop over all pairs of terms
        for a1, arg1 in enumerate(eq.args):
            print("ARG1:")
            print(arg1)
            if (arg1 == 0):
                continue
            if (already_found_term[a1] == True):
                continue

            #check if all 3 (35) terms for a full o/v permutation are present
            check_perms = []
            
            #Array to check already found permutations
            already_found_perm = numpy.zeros(len(permutation), dtype=bool)

            for a2, arg2 in enumerate(eq.args[a1+1:]):
                if (arg2 == 0):
                    continue
                if (already_found_term[a1+1+a2] == True):
                    continue
                if not check_same_type(arg1, arg2):
                    continue

                #now throw all wanted permutations on the tensor in loop 2
                for pp, p in enumerate(permutation):

                    if (already_found_perm[pp] == True):
                        continue

                    trial = arg2
                    for q in p:
                        trial = q.get_permuted(trial)
                        trial = substitute_dummies(trial, new_indices=True, pretty_indices=pdd)

                    trial = trial - arg1
                    trial = substitute_dummies(trial, new_indices=True, pretty_indices=pdd)

                    #We found a match!
                    if (trial == 0):
                        check_perms.append(a1+1+a2)
                        already_found_perm[pp] = True
                        break
                    elif (isinstance(trial, Mul)):
                        if (trial.args[0] < 0.000000001):
                            check_perms.append(a1+1+a2)
                            already_found_perm[pp] = True
                            break
                
                #We found all 3/35 necessary matches!!! full o/v permutation found
                if (len(check_perms) == len(permutation)):
                    print('Found a full symmetry!')
                    check_perms.insert(0, a1)
                    permutations_dict[perm_strings[np]] += [check_perms]
                    for i in check_perms:
                        already_found_term[i] = True
                    break

    
    return permutations_dict
