#
#Main file to start calculations
#
#SMT 2023-2024
#


import time
from base.hamiltonian import *
from base.isr import *
from base.evaluate import *
from utils.pretty_dummies import * 
from utils.term_list import *
from utils.permutations_finder import *
from utils.latex import *


#True: only real cluster amplitudes
hamiltonian_list = [hamiltonian_ucc2(False)]

#Calculation type
calc_type = 'd'

def main():

    for h in hamiltonian_list:

        #Main calculation
        ts = time.time()
        eq, lti, rti, pdd = evaluate(h, calc_type)
        te = time.time()
        print("\nMAIN CALCULATION TOOK TIME:")
        print(float(te-ts))

        #Specify if left or right hand side target indices
        #are chosen for the permutations finder
        occs = rti[0]
        virs = rti[1]

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
    

if __name__ == "__main__":
    print('Starting (U)CC code generator')
    main()
