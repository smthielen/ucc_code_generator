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


def return_perms(occs, virs):
    """
    Return perms according to target indices (in a very dirty manner....).
    """
    PO = PermutationOperator
    
    perms = []
    perm_strings = []

    #TODO: This is such a dirty way to do it....maybe do mor generally?!
    if (len(occs) == 2 and len(virs) == 2):
        #ijab for doubles
        perms.append([[PO(occs[0], occs[1])], [PO(virs[0], virs[1])], [PO(occs[0], occs[1]), PO(virs[0], virs[1])]]) 
        perm_strings.append(str(occs[0])+str(occs[1])+str(virs[0])+str(virs[1]))

    if (len(occs) == 2):
        #ij for doubles
        perms.append([[PO(occs[0], occs[1])]])
        perm_strings.append(str(occs[0])+str(occs[1]))

    if (len(virs) == 2):
        #ab for doubles
        perms.append([[PO(virs[0], virs[1])]]) 
        perm_strings.append(str(virs[0])+str(virs[1]))

    if (len(occs) == 3 and len(virs) == 3):
        perms.append([[PO(occs[0], occs[1])], [PO(occs[0], occs[2])], [PO(occs[1], occs[2])], [PO(occs[0], occs[1]), PO(occs[0], occs[2])], 
                      [PO(occs[0], occs[2]), PO(occs[0], occs[1])], [PO(virs[0], virs[1])], [PO(virs[0], virs[2])], [PO(virs[1], virs[2])], 
                      [PO(virs[0], virs[1]), PO(virs[0], virs[2])], [PO(virs[0], virs[2]), PO(virs[0], virs[1])], 
                      [PO(occs[0], occs[1]), PO(virs[0], virs[1])], [PO(occs[0], occs[2]), PO(virs[0], virs[1])], 
                      [PO(occs[1], occs[2]), PO(virs[0], virs[1])], [PO(occs[0], occs[1]), PO(occs[0], occs[2]), PO(virs[0], virs[1])], 
                      [PO(occs[0], occs[2]), PO(occs[0], occs[1]), PO(virs[0], virs[1])], [PO(occs[0], occs[1]), PO(virs[0], virs[2])], 
                      [PO(occs[0], occs[2]), PO(virs[0], virs[2])], [PO(occs[1], occs[2]), PO(virs[0], virs[2])], 
                      [PO(occs[0], occs[1]), PO(occs[0], occs[2]), PO(virs[0], virs[2])], [PO(occs[0], occs[2]), PO(occs[0], occs[1]), PO(virs[0], virs[2])], 
                      [PO(occs[0], occs[1]), PO(virs[1], virs[2])], [PO(occs[0], occs[2]), PO(virs[1], virs[2])], [PO(occs[1], occs[2]), PO(virs[1], virs[2])], 
                      [PO(occs[0], occs[1]), PO(occs[0], occs[2]), PO(virs[1], virs[2])], [PO(occs[0], occs[2]), PO(occs[0], occs[1]), PO(virs[1], virs[2])], 
                      [PO(occs[0], occs[1]), PO(virs[0], virs[1]), PO(virs[0], virs[2])], [PO(occs[0], occs[2]), PO(virs[0], virs[1]), PO(virs[0], virs[2])], 
                      [PO(occs[1], occs[2]), PO(virs[0], virs[1]), PO(virs[0], virs[2])], 
                      [PO(occs[0], occs[1]), PO(occs[0], occs[2]), PO(virs[0], virs[1]), PO(virs[0], virs[2])], 
                      [PO(occs[0], occs[2]), PO(occs[0], occs[1]), PO(virs[0], virs[1]), PO(virs[0], virs[2])],
                      [PO(occs[0], occs[1]), PO(virs[0], virs[2]), PO(virs[0], virs[1])], [PO(occs[0], occs[2]), PO(virs[0], virs[2]), PO(virs[0], virs[1])], 
                      [PO(occs[1], occs[2]), PO(virs[0], virs[2]), PO(virs[0], virs[1])], 
                      [PO(occs[0], occs[1]), PO(occs[0], occs[2]), PO(virs[0], virs[2]), PO(virs[0], virs[1])], 
                      [PO(occs[0], occs[2]), PO(occs[0], occs[1]), PO(virs[0], virs[2]), PO(virs[0], virs[1])]])
        perm_strings.append(str(occs[0])+str(occs[1])+str(occs[2])+str(virs[0])+str(virs[1])+str(virs[2]))

    if (len(occs) == 3):
        #ijk for triples
        perms.append([[PO(occs[0], occs[1])], [PO(occs[0], occs[2])], [PO(occs[1], occs[2])], 
            [PO(occs[0], occs[1]), PO(occs[0], occs[2])], [PO(occs[0], occs[2]), PO(occs[0], occs[1])]])
        perm_strings.append(str(occs[0])+str(occs[1])+str(occs[2])) 
    
    if (len(virs) == 3):
        #abc for triples
        perms.append([[PO(virs[0], virs[1])], [PO(virs[0], virs[2])], [PO(virs[1], virs[2])], 
            [PO(virs[0], virs[1]), PO(virs[0], virs[2])], [PO(virs[0], virs[2]), PO(virs[0], virs[1])]])
        perm_strings.append(str(virs[0])+str(virs[1])+str(virs[2])) 

    if (len(occs) == 3 and len(virs) == 3):
        #ijab for triplets
        perms.append([[PO(occs[0], occs[1])], [PO(virs[0], virs[1])], [PO(occs[0], occs[1]), PO(virs[0], virs[1])]]) 
        perm_strings.append(str(occs[0])+str(occs[1])+str(virs[0])+str(virs[1]))
        #ikab for triplets
        perms.append([[PO(occs[0], occs[2])], [PO(virs[0], virs[1])], [PO(occs[0], occs[2]), PO(virs[0], virs[1])]]) 
        perm_strings.append(str(occs[0])+str(occs[2])+str(virs[0])+str(virs[1]))
        #jkab for triplets
        perms.append([[PO(occs[1], occs[2])], [PO(virs[0], virs[1])], [PO(occs[1], occs[2]), PO(virs[0], virs[1])]]) 
        perm_strings.append(str(occs[1])+str(occs[2])+str(virs[0])+str(virs[1]))
        #ijac for triplets
        perms.append([[PO(occs[0], occs[1])], [PO(virs[0], virs[2])], [PO(occs[0], occs[1]), PO(virs[0], virs[2])]]) 
        perm_strings.append(str(occs[0])+str(occs[1])+str(virs[0])+str(virs[2]))
        #ikac for triplets
        perms.append([[PO(occs[0], occs[2])], [PO(virs[0], virs[2])], [PO(occs[0], occs[2]), PO(virs[0], virs[2])]]) 
        perm_strings.append(str(occs[0])+str(occs[2])+str(virs[0])+str(virs[2]))
        #jkac for triplets
        perms.append([[PO(occs[1], occs[2])], [PO(virs[0], virs[2])], [PO(occs[1], occs[2]), PO(virs[0], virs[2])]]) 
        perm_strings.append(str(occs[1])+str(occs[2])+str(virs[0])+str(virs[2]))
        #ijbc for triplets
        perms.append([[PO(occs[0], occs[1])], [PO(virs[1], virs[2])], [PO(occs[0], occs[1]), PO(virs[1], virs[2])]]) 
        perm_strings.append(str(occs[0])+str(occs[1])+str(virs[1])+str(virs[2]))
        #ikbc for triplets
        perms.append([[PO(occs[0], occs[2])], [PO(virs[1], virs[2])], [PO(occs[0], occs[2]), PO(virs[1], virs[2])]]) 
        perm_strings.append(str(occs[0])+str(occs[2])+str(virs[1])+str(virs[2]))
        #jkbc for triplets
        perms.append([[PO(occs[1], occs[2])], [PO(virs[1], virs[2])], [PO(occs[1], occs[2]), PO(virs[1], virs[2])]]) 
        perm_strings.append(str(occs[1])+str(occs[2])+str(virs[1])+str(virs[2]))

    if (len(occs) == 3):
        #ij for triples
        perms.append([[PO(occs[0], occs[1])]])
        perm_strings.append(str(occs[0])+str(occs[1]))
        #ik for triples
        perms.append([[PO(occs[0], occs[2])]])
        perm_strings.append(str(occs[0])+str(occs[2]))
        #jk for triples
        perms.append([[PO(occs[1], occs[2])]])
        perm_strings.append(str(occs[1])+str(occs[2]))

    if (len(virs) == 3):
        #ab for triples
        perms.append([[PO(virs[0], virs[1])]]) 
        perm_strings.append(str(virs[0])+str(virs[1]))
        #ac for triples
        perms.append([[PO(virs[0], virs[2])]]) 
        perm_strings.append(str(virs[0])+str(virs[2]))
        #bc for triples
        perms.append([[PO(virs[1], virs[2])]]) 
        perm_strings.append(str(virs[1])+str(virs[2]))


    print(perms)
    for i in perms:
        print(i)

    print(perm_strings)


    return perms, perm_strings

'''
def own_return_perms(occs, virs):
    """
    Return perms according to target indices (in a very dirty manner....).
    """
    PO = PermutationOperator

    perms_o = []
    perms_v = []
    perms_ov = []

    #TODO: This is such a dirty way to do it....maybe do mor generally?!
    if (len(occs) == 2):
        perms_o.append([[occs[0], occs[1]]]) 
        perms_ov.append([[occs[0], occs[1]]]) 
    elif (len(occs) == 3):
        perms_o.append([[occs[0], occs[1]]])
        perms_o.append([[occs[0], occs[2]]])
        perms_o.append([[occs[1], occs[2]]])
    if (len(virs) == 2):
        perms_v.append([[virs[0], virs[1]]]) 
        perms_ov.append([[virs[0], virs[1]]])
        #Mixed term!!
        perms_ov.append([[occs[0], occs[1]], [virs[0], virs[1]]]) 
    elif (len(virs) == 3):
        perms_v.append([[virs[0], virs[1]]]) 
        perms_v.append([[virs[0], virs[2]]]) 
        perms_v.append([[virs[1], virs[2]]]) 

    return perms_o, perms_v, perms_ov 
'''

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

    print(t1)
    print(t2)
    #t1 = t1.sort()
    #t2 = t2.sort()

    if (t1 == t2):
        print(t1)
        print(t2)
        print("SAME TYPE CONFIRMED")
        return True
    else:
        return False

'''
def own_permutation_operator(perm, eq, pdd):
    """
    self-written permutation operator PO
    """
   
    permuted_tensor = 1.0

    #leave the prefactor out
    for arg in eq.args[1:]:

        if isinstance(arg, KroneckerDelta):
            nd = arg.args
            for d, dd in enumerate(arg.args):
                if (dd == perm[0]):
                    nd[d] = perm[1] 
                elif (dd == perm[1]):
                    nd[d] = perm[0]
           
            new_delta = KroneckerDelta((nd[0], nd(1)))
            permuted_tensor *= new_delta
            permuted_tensor = substitute_dummies(permuted_tensor, new_indices=True, pretty_indices=pdd)

        elif isinstance(arg, AntiSymmetricTensor):
            nu = []
            nl = []
            for i in enumerate(arg.upper.args):
                print(i)
                nu.append(i[1])
            print(nu)
            for i in enumerate(arg.lower.args):
                print(i)
                nl.append(i[1])
            print(nl)
            print(arg)
            print(nu)
            print(nl)
            for u, uu in enumerate(arg.upper.args):
                if (uu == perm[0]):
                    nu[u] = perm[1] 
                elif (uu == perm[1]):
                    nu[u] = perm[0] 
            for l, ll in enumerate(arg.lower.args):
                if (ll == perm[0]):
                    nl[l] = perm[1] 
                elif (ll == perm[1]):
                    nl[l] = perm[0] 

            new_tensor = 1.0
            if len(nu) == 1:
                new_tensor = AntiSymmetricTensor(str(arg.symbol), (nu[0], ), (nl[0], ))
            elif len(nu) == 2:
                new_tensor = AntiSymmetricTensor(str(arg.symbol), (nu[0], nu[1]), (nl[0], nl[1]))
            elif len(nu) == 3:
                new_tensor = AntiSymmetricTensor(str(arg.symbol), (nu[0], nu[1], nu[2]), (nl[0], nl[1], nl[2]))

            permuted_tensor *= new_tensor
            #permuted_tensor = substitute_dummies(permuted_tensor, new_indices=True, pretty_indices=pdd)

    #TODO: Proper? Prefactor always at front in Mul right?!
    permuted_tensor *= -1.0*eq.args[0]

    #y not one more time....the more the merrier
    permuted_tensor = substitute_dummies(permuted_tensor, new_indices=True, pretty_indices=pdd)

    print('in eq:')
    print(eq)
    print('out eq:')
    print(permuted_tensor)
    
    return permuted_tensor
'''
        


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
                    print('arg2:')
                    print(arg2)
                    for q in p:
                        trial = q.get_permuted(trial)
                        trial = substitute_dummies(trial, new_indices=True, pretty_indices=pdd)

                    trial = trial - arg1
                    trial = substitute_dummies(trial, new_indices=True, pretty_indices=pdd)
                    print("\n")
                    print('combination '+str(a1)+','+str(a2))
                    print(trial)
                    print("\n") 

                    #We found a match!
                    if (trial == 0):
                        print('match!!')
                        check_perms.append(a1+1+a2)
                        already_found_perm[pp] = True
                        break
                    elif (isinstance(trial, Mul)):
                        if (trial.args[0] < 0.000000001):
                            print('match!!')
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


