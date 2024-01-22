#
#NOTE: just a bunch of unfnished stuff....
#
#
#
#
#SMT 2023-2024
#


from term_list import term_expression, common_tensor_expression, term_list


#TODO: Pretty badly hacked, but works fine, i guess....
occ_i = 'ijklmno'
vir_i = 'abcdefgh'
all_i = 'ijklmnoabcdefgh'



def swap_indices(pref, indices, perms):
    """
    Swap indices according to indexes
    IMPORTANT: This function is desigend such that one permutations
    always yields a sign chage of -1 and two permutations +1
    """
    proper_indices = indices
    for p in perms:
        proper_indeces[p[0]] = indices[p[1]]
        proper_indeces[p[1]] = indices[p[0]]
        pref *= -1

    return pref, proper_indices


def return_summ_indes_swap(indices):
    """

    """
    if len(indices) == 2:
        p2 = [[0,1]]
        return p2
    elif len(indices) == 3:
        p3 = [[0,1,2], ]


def smallest_first(pref, indices, obj):
    """
    Smallest index at the front
    """
    if (obj == 'i' or obj == 't2' or obj == 't2c'):
        #first two indices
        if (all_i.find(indices[0]) > all_i.find(indices[1])):
            pref, indices = swap_indices(pref, indices, [0, 1])
        #last two indices
        if (all_i.find(indices[2]) > all_i.find(indices[3])):
            pref, indices = swap_indices(pref, indices, [2, 3])
    elif (obj == 't3' or obj =='t3c'):
        #first three indices
        if (all_i.find(indices[0]) > all_i.find(indices[1])):
            pref, indices = swap_indices(pref, indices, [0, 1])
        if (all_i.find(indices[1]) > all_i.find(indices[2])):
            pref, indices = swap_indices(pref, indices, [1, 2])
        if (all_i.find(indices[0]) > all_i.find(indices[1])):
            pref, indices = swap_indices(pref, indices, [0, 1])
        #lastst three indices
        if (all_i.find(indices[3]) > all_i.find(indices[4])):
            pref, indices = swap_indices(pref, indices, [3, 4])
        if (all_i.find(indices[4]) > all_i.find(indices[5])):
            pref, indices = swap_indices(pref, indices, [4, 5])
        if (all_i.find(indices[3]) > all_i.find(indices[4])):
            pref, indices = swap_indices(pref, indices, [3, 4])
    #TODO: needed or not!??
    elif (obj == 'delta'):
        if (all_i.find(indices[0]) > all_i.find(indices[1])):
            pref, indices = swap_indices(pref, indices, [0, 1])
            #for a delta, the prefactor does not change
            pref *= -1

    return pref, indices


def canonical_eris(tensor_list):
    """
    Bring ERIs in canonical form:
    OOOV, OOVV, OVOV, OVVV
    """
    for i in tensor_list.expressions:
        for j in i.common_tensors:
            for k, kk in enumerate(j.tensors):
                if kk == 'i':
                    raw_ind = ''
                    for ind in j.indices[k]:
                        if ind in occ_i:
                            raw_ind.join('o')
                        elif ind in vir_i:
                            raw_ind.join('v')
                    #Now, go through all cases that need a major index chage
                    if (raw_ind == 'oovo'): # to ooov
                        j.prefactor, j.indices[k] = swap_indices(j.prefactor, j.indices[k], [2, 3])
                    elif (raw_ind == 'ovoo'): # to ooov
                        j.prefactor, j.indices[k] = swap_indices(j.prefactor, j.indices[k], [[0, 2], [1, 3]])
                    elif (raw_ind == 'vooo'): # to ooov
                        j.prefactor, j.indices[k] = swap_indices(j.prefactor, j.indices[k], [[0, 3], [1, 2]])
                    elif (raw_ind == 'ovvo'): # to ovov
                        j.prefactor, j.indices[k] = swap_indices(j.prefactor, j.indices[k], [2, 3])
                    elif (raw_ind == 'voov'): # to ovov
                        j.prefactor, j.indices[k] = swap_indices(j.prefactor, j.indices[k], [0, 1])
                    elif (raw_ind == 'vovo'): # to ovov
                        j.prefactor, j.indices[k] = swap_indices(j.prefactor, j.indices[k], [0, 1], [2, 3])
                    elif (raw_ind == 'vvoo'): # to oovv
                        j.prefactor, j.indices[k] = swap_indices(j.prefactor, j.indices[k], [[0, 2], [1, 3]])
                    elif (raw_ind == 'vovv'): # to ovvv
                        j.prefactor, j.indices[k] = swap_indices(j.prefactor, j.indices[k], [0, 1])
                    elif (raw_ind == 'vvov'): # to ovvv
                        j.prefactor, j.indices[k] = swap_indices(j.prefactor, j.indices[k], [[0, 2], [1, 3]])
                    elif (raw_ind == 'vvvo'): # to ovvv
                        j.prefactor, j.indices[k] = swap_indices(j.prefactor, j.indices[k], [[0, 3], [1, 2]])
                    #Smallest at the front!
                    j.prefactor, j.indices[k] = smallest_first(j.prefactor, j.indices[k], 'i')
    
    return tensor_list




def check_eri_type(t1, t2):
    """
    Checks if two terms possess the same type of ERI
    """
    t1_o1 = t1_o2 = t1_v1 = t1_v2 = 0
    t2_o1 = t2_o2 = t2_v1 = t2_v2 = 0
    #TODO: This is done pretty dirty....maybe do it better!?
    #TODO: Maybe combine it with the canonical ERI function....
    for i, ii in enumerate(t1.tensors):
        #ERI detected
        if ii == 'i':
            #Count occs and virs and check where in ERI
            for c, cc in enumerate(t1.indices[i]):
                if (occ_i.find(cc) != -1):
                    if (c > 2):
                        t1_o1 += 1
                    else:
                        t1_o2 += 1
                elif (vir_i.find(c) != -1):
                    if (c > 2):
                        t1_v1 += 1
                    else:
                        t1_v2 += 1
    for i, ii in enumerate(t2.tensors):
        #ERI detected
        if ii == 'i':
            #Count occs and virs and check where in ERI
            for c, cc in enumerate(t2.indices[i]):
                if (occ_i.find(cc) != -1):
                    if (c > 2):
                        t2_o1 += 1
                    else:
                        t2_o2 += 1
                elif (vir_i.find(c) != -1):
                    if (c > 2):
                        t2_v1 += 1
                    else:
                        t2_v2 += 1
    if (t1_o1 == t2_o1 and t1_o2 == t2_o2 and t1_v1 == t2_v1 and t1_v2 == t2_v2):
        return True
    else:
        return False


def check_ti_positions(t1, t2, pi):
    """
    Check if target indeces allow for a permutation, i. e. if they are placed
    at the correct tensors to be interchanged
    """
    i1_t1 = i1_t2 = i2_t1 = i2_t2 = ''

    #Look for tensors that have the target indices
    for i, ii in enumerate(t1.indices):
        if (ii.find(pi[0]) != -1):
            i1_t1 = t1.tensors[i]
        if (ii.find(pi[1]) != -1):
            i2_t1 = t1.tensors[i]
    for i, ii in enumerate(t2.indices):
        if (ii.find(pi[0]) != -1):
            i1_t2 = t2.tensors[i]
        if (ii.find(pi[1]) != -1):
            i2_t2 = t2.tensors[i]

    #This is the requirement for a i/j, a/b etc. permutational symmetry
    if (i1_t1 == i2_t2 and i1_t2 == i2_t1):
        return True
    else:
        return False


def get_summ_ind(expr):
    """
    Find all occ/vir summation indices in an expression
    """

    all_indices = ''
    for i in expr.indices:
        all_indices.append(i)

    occs = ''
    virs = ''

    for i, ii in enumerate(all_indices):
        for jj in all_indices[i+1:]:

            o1 = all_occ.find(ii)
            o2 = all_occ.find(jj)
            
            v1 = all_vir.find(ij) 
            v2 = all_vir.find(jj) 

            if (o1 != -1 and o2 != -1):
                occs.append(i)
            elif (v1 != -1 and v2 != -1):
                virs.append(i)

    return occs, virs


def check_for_permutations(common_t, perm):
    """
    Check for all permutations occuring in a common_tensor object
    """

    list_of_perms = []
    #TODO: Does this lead to a loss of information if one does not all checks everytime?!
    #auxiliary array to bookkeep all checked terms in the common_tensor object
    #already_checked = np.zeros(len(common_t))

    for i, ii in enumerate(common_t):
        occs_ii, virs_ii = get_summ_ind(ii) 

        for j, jj in enumerate(common_t[i:]):
            occs_jj, virs_jj = get_summ_ind(jj) 
            
            #BASIC CHECK
            #Are the terms already checked?
            #if (already_checked[i] != 0 or already_checked[j] != 0):
            #    continue
           
            #SPECIAL CHECKS
            #Check for same numbers of occ/vir indices
            if not ((len(occs_ii) == len(occs_jj)) and (len(virs_ii) == len(virs_jj))):
                continue
            #Check if ERIs have the same type
            if not (check_eri_type(ii, jj)):
                continue
            #Check if the positions of the target indeces allow for a permutation
            if not (check_ti_positions(ii, jj, perm)):
                continue
            
            #TODO: Fill in the actual stoff here....
            #Now do the dirty work, namely build up all perms....:(
            #Start with ERI ckeck


            #Now do all index permutations




    return list_of_perms



def find_permutations(tensor_list, perm_list):
    """
    Find permutations....
    """
    list_of_perms = []
    
    #Loop over all specified permutations
    for p in perm_list:
        for ii, i in enumerate(tensor_list.expressions):
            for ij, j in enumerate(i.common_tensors):
                check_for_permutations(j, p)

                #TODO: Fill in the actual stuff here....

    return list_of_perms


occ = 'abcdef'
vir = 'ijklmn'

def return_diag_fock(expr):

    diag_expr = []
    for el in expr:
        if 'f' in el.tensors:
            for i, j in zip(el.tensors, el.indices):
                if i == 'f':
                    if ((occ.find(j[0]) !=-1) and (occ.find(j[1]) !=-1)):
                        diag_expr.append(el)
                    elif ((vir.find(j[0]) !=-1) and (vir.find(j[1]) !=-1)):
                        diag_expr.append(el)
        else:
            diag_expr.append(el)
    
    return diag_expr
