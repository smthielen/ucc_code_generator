#
#routine to search a/b and i/j permutations in an equation
#TODO: apply also to a/b/c and i/j/k !!
#
#SMT 5/2023
#


def search_permutations(expr):

    tensor_dict = []

    if (isinstance(expr, Add)):
        for x in expr.args:
            dict_element = {}
            if (isinstance(x, Mul)):
                #print(x)
                #dict_element.update({'pref' : str(x.args[0])})
                for y in x.args:
                    if (isinstance(y, AntiSymmetricTensor)):
                        key = str(y.symbol)
                        value = ''
                        for u in y.upper:
                            value += u.name
                        for l in y.lower:
                            value += l.name
                        if (key in dict_element):
                            dict_element[key] += [value]
                        else:
                            dict_element.update({key : [value]})
                    else:
                        print('.')
                    #    print(type(y))
                tensor_dict.append([dict_element, x.args[0]])
            elif (isinstance(x, AntiSymmetricTensor)):
                key = str(x.symbol)
                value = ''
                for u in x.upper:
                    value += u.name
                for l in x.lower:
                    value += l.name
                if (key in dict_element):
                    dict_element[key] += [value]
                else:
                    dict_element.update({key : [value]})
                #Is this 1.0 always true?!
                tensor_dict.append([dict_element, 1.0])
    #print(tensor_dict)
    permutations = []

    length = len(tensor_dict)
    print('length of tensor_dict:')
    print(length)
    print('!tensor_dict!:')
    print(tensor_dict)

    iter1 = 0
    while (iter1 < length):
        #ind1 = ''
        #for i in tensor_dict[iter1][0].values():
        #    ind1 += i
        #print(ind1)
        iter2 = iter1 + 1
        while (iter2 < length):
            #ind2_ij = ''
            #ind2_ab = ''
            #First check if terms are compatible, i. e. the keys are the same
            if (tensor_dict[iter1][0].keys() == tensor_dict[iter2][0].keys()):
                #print('same')
                #ind2 = ''
                t_i1 = t_j1 = t_a1 = t_b1 = ''
                t_i2 = t_j2 = t_a2 = t_b2 = ''

                for k,v in tensor_dict[iter1][0].items():
                    for w in v:
                        if (w.find('i') != -1):
                            t_i1 = k
                        if (w.find('j') != -1):
                            t_j1 = k
                        if (w.find('a') != -1):
                            t_a1 = k
                        if (w.find('b') != -1):
                            t_b1 = k

                for k,v in tensor_dict[iter2][0].items():
                    for w in v:
                        if (w.find('i') != -1):
                            t_i2 = k
                        if (w.find('j') != -1):
                            t_j2 = k
                        if (w.find('a') != -1):
                            t_a2 = k
                        if (w.find('b') != -1):
                            t_b2 = k
                '''
                t_i1 = [k for k,v in tensor_dict[iter1][0].items() if v.find('i') != -1][0]
                t_j1 = [k for k,v in tensor_dict[iter1][0].items() if v.find('j') != -1][0]
                t_i2 = [k for k,v in tensor_dict[iter2][0].items() if v.find('i') != -1][0]
                t_j2 = [k for k,v in tensor_dict[iter2][0].items() if v.find('j') != -1][0]
                
                t_a1 = [k for k,v in tensor_dict[iter1][0].items() if v.find('a') != -1][0]
                t_b1 = [k for k,v in tensor_dict[iter1][0].items() if v.find('b') != -1][0]
                t_a2 = [k for k,v in tensor_dict[iter2][0].items() if v.find('a') != -1][0]
                t_b2 = [k for k,v in tensor_dict[iter2][0].items() if v.find('b') != -1][0]
                '''
                #Check for ij permutation, i. e. a and b must have the same key
                if ((t_a1 == t_a2 and t_b1 == t_b2) and (t_i1 == t_j2 and t_i2 == t_j1)):
                    
                    #Check if 1 and 2 are different tensor objects, i. e. have different keys
                    #The second condition, however, should be negligible
                    if (t_i1 != t_j1 and t_i2 != t_j2):
                        #print("Found ij symmetry for terms "+str(iter1)+","+str(iter2)+" for tensors of different types: "+print(t_i1)+","+print(t_j1))
                        print("Found ij symmetry for terms "+str(iter1)+","+str(iter2)+" for tensors of different types")
                        permutations.append([iter1, iter2, 'ij'])
                    
                    #Both 1 and 2 are the same tensor object, i. e. have the same key
                    #Just to be sure, check once again, although not necessary
                    if (t_i1 == t_j1 and t_i2 == t_j2):
                        dict1 = tensor_dict[iter1][0]
                        dict1_v = dict1[t_i1]
                        dict2 = tensor_dict[iter2][0]
                        dict2_v = dict1[t_i2]
                        
                        #Check if a and b occur in the same value, remember they have the same key
                        a1_pos = 0
                        b1_pos = 0
                        for idx in dict1_v:
                            if (idx.find('a') != -1):
                                a1_pos = dict1_v.index(idx)
                            if (idx.find('b') != -1):
                                b1_pos = dict1_v.index(idx)
                        a2_pos = 0
                        b2_pos = 0
                        for idx in dict2_v:
                            if (idx.find('a') != -1):
                                a2_pos = dict2_v.index(idx)
                            if (idx.find('b') != -1):
                                b2_pos = dict2_v.index(idx)

                        if (a1_pos == b1_pos and a2_pos == b2_pos):
                            print("Found ij symmetry for terms "+str(iter1)+","+str(iter2)+" for tensors of SAME type")#: "+print(t_i1))
                            permutations.append([iter1, iter2, 'ij'])
                        if (a1_pos != b1_pos and a2_pos != b2_pos):
                            print("Found ij symmetry for terms "+str(iter1)+","+str(iter2)+" for tensors of SAME type")#+print(t_i1))
                            permutations.append([iter1, iter2, 'ij'])

            iter2 += 1
        iter1 += 1

    for i in tensor_dict:
        i[1] = float(i[1])

    #print(tensor_dict)

    return tensor_dict, permutations
