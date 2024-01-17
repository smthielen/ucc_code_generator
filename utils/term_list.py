#
#Handling of sympy expressions after the evaluate routine
#
#-> term_expression: base class to contain all tensors, indices and the prefactor of a single term
#-> common_tensor_expression: class to summarize terms with common tensors, i. e. all eri*t2 terms,
#   eri*t1 terms, for instance, woild get a new class
#-> term_list: list of common_tensor_expression objects
#-> build_term_list: build a term_list object using the original sympy expression from evaluate
#
#
#SMT 5/2023 - 8/2023
#


from sympy.physics.secondquant import (AntiSymmetricTensor, KroneckerDelta, wicks,
        NO, evaluate_deltas, substitute_dummies, Commutator,
        simplify_index_permutations, PermutationOperator, Symbol, Add, Mul) #simplify_index_permutations and PermutationOperator ??

import numpy as np


class term_expression:
    """
    Class for expressions in the term list
    """

    def __init__(self):
        self.tensors = []
        self.indices = []
        self.prefactor = 0
        self.permutations = ''

    def add_tensor(self, tensor, indices):
        self.tensors.append(tensor)
        self.indices.append(indices)

    def add_prefactor(self, prefactor):
        self.prefactor = prefactor

    def add_permutations(self, permutations):
        self.permutations = permutations


class common_tensor_expression:
    """
    Class for expressions with common tensors
    """

    def __init__(self):
        self.common_tensors = []

    def add_common_tensor(self, expr):
        if not (isinstance(expr, term_expression) or isinstance(expr, Mul)):
            raise TypeError(f'only term_expressions can be added to a common_tensor_expression!')
        else:
            self.common_tensors.append(expr)


class term_list:
    """
    Class for the term list
    """

    def __init__(self):
        self.expressions = []

    def add_expression(self, expr):
        if not isinstance(expr, common_tensor_expression):
            raise TypeError(f'only common_tensor_expressions can be added to a term_list!')
        else:
            self.expressions.append(expr)


def get_key(val, dic):
   
    for key, value in dic.items():
        if val == value:
            return key
 
    return "key doesn't exist"


def build_term_list(expr, perm_dict):
    """
    Function to build a term_list. Takes the sympy expression and list of desired permutations as an input
    """

    #final list to be returned
    tensor_list = term_list()

    if (expr == 0):
        return tensor_list    

    #note all terms that are already taken care of
    already_checked = np.zeros(len(expr.args), dtype=bool) 

    #Check if its not a zero or something....
    if (isinstance(expr, Add)):
        for xx, x in enumerate(expr.args):

            #Check if already in the term
            if (already_checked[xx] == True):
                continue

            te = term_expression()
            if (isinstance(x, Mul)):
                for y in x.args:
                    indices = ''
                    if (isinstance(y, KroneckerDelta)):
                        indices += str(y.args[0])
                        indices += str(y.args[1])
                        te.add_tensor('delta', indices)
                    elif (isinstance(y, AntiSymmetricTensor)):
                        for u in y.upper:
                            indices += u.name
                        for l in y.lower:
                            indices += l.name
                        #print(str(y.symbol))
                        te.add_tensor(str(y.symbol), indices)
                    #else:
                    #    print(',.')
                    #    print(y)
                te.add_prefactor(x.args[0])
               
                perm_operator = ''

                #check if involved in a permutation
                for perm_types in perm_dict.values():
                    for perms in perm_types:
                        if xx in perms:
                            print(xx)
                            print(perms)
                            print(get_key(perm_types, perm_dict))
                            perm_key = get_key(perm_types, perm_dict)
                            perm_operator += str(perm_key)
                            te.add_permutations(perm_operator)
                            for i in perms:
                                already_checked[i] = True

                #not involved in a permutation
                if (te.permutations == ''):
                    te.add_permutations('none')

                if not tensor_list.expressions:
                    cnew = common_tensor_expression()
                    cnew.add_common_tensor(te)
                    tensor_list.add_expression(cnew)

                else:
                    added = False
                    for c in tensor_list.expressions:
                        if c.common_tensors[0].tensors == te.tensors :
                            c.add_common_tensor(te)
                            added = True
                            break
                    if added == False:
                        cnew = common_tensor_expression()
                        cnew.add_common_tensor(te)
                        tensor_list.add_expression(cnew)

            
    
    else:
        print("Nothing to permute in a single term....")


    #Now, if wished, check for permutations
    #permutation_list = []
    
    #if not perms:
    #    print("Did not ask for search of permutations!")
    '''
    else:

        for p in perms:
            ti1, ti2 = p[0], p[1]
            print(ti1)
            print(ti2)

            permutation_list_c = []
            for nc, c in enumerate(tensor_list.expressions):

                for ni, i in enumerate(c.common_tensors):
                    ti1_pos1 = ti2_pos1 = -1
                    
                    for idx_pos, idx in enumerate(i.indices):
                        if idx.find(ti1) != -1:
                            ti1_pos1 = idx_pos
                        if idx.find(ti2) != -1:
                            ti2_pos1 = idx_pos

                    for nj, j in enumerate(c.common_tensors):
                        ti1_pos2 = ti2_pos2 = -1
                   
                        if nj > ni:
                            for idx_pos, idx in enumerate(j.indices):
                                print(idx_pos)
                                print(idx)
                                if idx.find(ti1) != -1:
                                    ti1_pos2 = idx_pos
                                if idx.find(ti2) != -1:
                                    ti2_pos2 = idx_pos

                            if ((i.tensors[ti1_pos1] == j.tensors[ti2_pos2]) and (i.tensors[ti2_pos1] == j.tensors[ti1_pos2])):
                                print('heyja')
                                print(ti1_pos1)
                                print(ti2_pos2)
                                print(ti2_pos1)
                                print(ti1_pos2)
                                print(i.tensors)
                                print(i.indices)
                                print(j.tensors)
                                print(j.indices)
                                rest_ti1_pos1= ''.join(sorted(i.indices[ti1_pos1].replace(ti1,'')))
                                print(rest_ti1_pos1)
                                rest_ti2_pos1= ''.join(sorted(i.indices[ti2_pos1].replace(ti2,'')))
                                print(rest_ti2_pos1)
                                
                                rest_ti1_pos2= ''.join(sorted(j.indices[ti1_pos2].replace(ti1,'')))
                                print(rest_ti1_pos2)
                                rest_ti2_pos2= ''.join(sorted(j.indices[ti2_pos2].replace(ti2,'')))
                                print(rest_ti2_pos2)

                                if (rest_ti1_pos1 == rest_ti2_pos2) and (rest_ti1_pos2 == rest_ti2_pos1):
                                    print('BINGO!!')
                                    permutation_list_c.append([ti1+ti2, nc+1, ni+1, nj+1])

            permutation_list.append(permutation_list_c)

            #TODO: Is an else case needed after all?!
            #else:
            #elif (isinstance(x, AntiSymmetricTensor)):
            #    print('.')
            #    print(x.symbol)
    '''
    return tensor_list #, permutation_list


def print_term_list(ls):
    """
    Simple print function for a term_list to be printed in the terminal
    """

    for ni, i in enumerate(ls.expressions):
        for nj, j in enumerate(i.common_tensors):
            print("OBJECT "+str(ni)+","+str(nj)+":\n")
            print(j.permutations)
            print(j.tensors)
            print(j.indices)
            print(j.prefactor)


def save_term_list(obj,cal,ls,expr,time):
    """
    Function to save term list to a file
    """

    obj_str = str(type(obj))
    p1=obj_str.rfind('.')+1
    p2=obj_str.rfind('\'')
    obj_str = obj_str[p1:p2]

    f = open(obj_str+'__'+cal,'w')
    f.write('Object: '+obj_str+'\n')
    f.write('Calculation: '+cal+'\n')
    f.write('EQUATION:\n')
    f.write('\n')

    sets = 0
    terms = 0

    for ni, i in enumerate(ls.expressions):
        sets+=1
        for nj, j in enumerate(i.common_tensors):
            terms+=1
            f.write(str("OBJECT "+str(ni+1)+","+str(nj+1)+":\n"))
            if not (j.permutations == ''):
                f.write(str(j.permutations)+'\n')
            f.write(str(j.tensors)+'\n')
            f.write(str(j.indices)+'\n')
            f.write(str(j.prefactor)+'\n')
        f.write('\n')

    f.write('NUMBER OF COMMON TERM SETS:'+str(sets)+'\n')
    f.write('NUMBER OF OVERALL TERMS:'+str(terms)+'\n')
    if (expr != 0):
        f.write('NUMBER OF OVERALL TERMS WITHOUT SYMMETRIZATION:'+str(len(expr.args))+'\n')

    #if (expr == 0 and terms == 0):
    #    f.write('SAME NUMBER AS IN SYMPY EXPRESSION! VERY GOODY!\n')
    #elif (len(expr.args) == terms):
    #    f.write('SAME NUMBER AS IN SYMPY EXPRESSION! VERY GOODY!\n')
    #else:
    #    ('NOT SAME TERM NUMBER AS IN SYMPY EXPRESSION. EVERYTHING IS WRONG!!\n')
    
    f.write('TOOK TIME:'+str(time))
    f.close() 
