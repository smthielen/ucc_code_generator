#
#Commutator class definition
#
#-> no_commutator: as it says, no commutator, i. e. F and V part
#-> commutator: single commutator, no nesting
#-> nested_commutator: singly or higer nested commutator
#
#SMT 5/2023 - 8/2023
#


class no_commutator:
    """
    Auxiliary class to handle expressions without a commutator 
    """

    def __init__(self, arg, prefactor):
        self.arg = arg
        self.prefactor = prefactor


class commutator:
    """
    Class for commutator expressions
    """

    def __init__(self, arg1 , arg2, prefactor):
        self.arg1 = arg1
        self.arg2 = arg2
        self.prefactor = prefactor


class nested_commutator:
    """
    Class for nested commutator expressions
    """

    def __init__(self, arg_list, comm_type_list, prefactor):
        self.arg_list = arg_list
        self.comm_type_list = comm_type_list   
        self.prefactor = prefactor

        #TODO: do this properly
        if ((len(self.arg_list)-2) != len(self.comm_type_list)):
            print('Incompatible setup for nested commutator!!')
