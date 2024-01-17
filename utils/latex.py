#
#Latex output for equations
#
#SMT 5/2023
#


import numpy
import subprocess
import fractions


def lhs(t):
    lhs_str = ''
    #if (t == 'd'):
    #    lhs_str += '& 0 \ '+'\\'+'overset{!}{=} \ '+'\\'+'braket{'+'\\'+'Phi^{ab}_{ij}|'+'\\'+'bar{H}|'+'\\'+'Phi_0} \ = \ '
    #elif (t == 't1'):
    #    lhs_str += 't^{a}_{i}'
    #else:
    lhs_str += str(t)
    lhs_str += '\\' + ' ='
    return lhs_str


def prefactor(o, f):
    pref_str = ''
    pref = fractions.Fraction(str(o.prefactor)).limit_denominator(1000)
    if (pref.numerator < 0):
        pref_str += '-'
    else:
        if not f:
            pref_str += '+'
    if (pref.denominator != 1):
        pref_str += '\\'+'frac{'+str(abs(pref.numerator))+'}{'+str(pref.denominator)+'}'
    return pref_str


def permutation(o):
    perm_str = ''
    perm = o.permutations
    if (perm != 'none'):
        perm_str += '\\'+'mathcal{P}^-_{('+str(perm)+')}'
    return perm_str


def indices(o):
    l = int(len(o)*0.5)
    upper = str(o[:l])
    lower = str(o[l:])
    return upper, lower


def tensor(t, ind):
    if (t == 'delta'):
        tensor_str = '\\'+'delta_{'+ind[0]+ind[1]+'}'
    elif (t.find('t') != -1):
        tensor_str = '\\'+'sigma^{'+ind[0]
        if (t.find('c') != -1):
            tensor_str += '*'
        tensor_str += '}_{'+ind[1]+'}'
    elif (t == 'i'):
        tensor_str = '\\'+'braket{'+ind[0]+'||'+ind[1]+'}'
    elif (t == 'f'):
        tensor_str = 'f_{'+ind[0]+ind[1]+'}'
    elif (t == 'd'):
        tensor_str = 'd_{'+ind[0]+ind[1]+'}'
    return tensor_str


def write_latex(obj, typ, terms):
    """
    Generate Latex output
    """
    
    obj_str = str(type(obj))
    p1=obj_str.rfind('.')+1
    p2=obj_str.rfind('\'')
    obj_str = obj_str[p1:p2]

    file = open(obj_str+'__'+typ+'.tex', 'w')
    file.write('\\'+'documentclass{article}'+'\n')
    file.write('\\'+'usepackage[a4paper]{geometry}'+'\n')
    file.write('\\'+'usepackage{amsmath, amssymb}'+'\n')
    file.write('\\'+'usepackage{xcolor}'+'\n')
    file.write('\\'+'usepackage{braket}'+'\n')
    file.write('\n')
    file.write('\\'+'begin{document}'+'\n')

    #equation
    file.write('\\'+'begin{align*}'+'\n')
    #LHS of equation
    file.write(lhs(typ)+'\\'+'\\'+'\n')

    #now loop over terms in RHS
    first_line = True
    for common in terms.expressions:
        #file.write('\\'+'\\'+'\n')
        #file.write('Tensor product!!'+'\\'+'\\'+'\n')
        #file.write('\\'+'\\'+'\n')

        line_str = '& '+'\\'+' '+'\\'+' '
        
        tensor_count = 0
        for c in common.common_tensors:
            print(c)

            term_str = ''
            term_str += prefactor(c, first_line)
            term_str += ' '
            term_str += permutation(c)
            term_str += ' '
            
            for tn, t in enumerate(c.tensors):
                inds = indices(c.indices[tn])
                term_str += tensor(t, inds)
                term_str += ' '

            line_str += term_str
            print(line_str)
                
            first_line = False
            tensor_count += 1
            
            #if (tensor_count == 4):
            if (tensor_count == 3):
                file.write(line_str+'\\'+'\\'+'\n')
                line_str = '& '+'\\'+' '+'\\'+' '+'\\'+' '+'\\'+' '
                tensor_count = 0

        file.write(line_str+'\\'+'\\'+'\n')
    
    file.write('\\'+'end{align*}'+'\n')
    file.write('\\'+'end{document}'+'\n')
    file.close()
