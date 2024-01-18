#
#Commutator class definition
#
#-> no_commutator: as it says, no commutator, i. e. F and V part
#-> commutator: single commutator, no nesting
#-> nested_commutator: singly or higer nested commutator
#
#SMT 2023-2024
#


class no_commutator:
    """
    Auxiliary class to handle expressions without a commutator
    """

    def __init__(self, arg, pref):
        self.arg = arg
        self.pref = pref


class commutator:
    """
    Class for commutator expressions
    """

    def __init__(self, arg1 , arg2, pref):
        self.arg1 = arg1
        self.arg2 = arg2
        self.pref = pref


class nested_commutator:
    """
    Class for nested commutator expressions
    """

    def __init__(self, arg_list, comm_type_list, pref):
        self.arg_list = arg_list
        self.comm_type_list = comm_type_list   
        self.pref = pref

        #check for proper setup of nested commutators 
        if ((len(self.arg_list)-2) != len(self.comm_type_list)):
            raise ValueError('No valid set up of nested commutator: '\
                    'Commutator type list ({len(self.comm_type_list)} entries) '\
                    'must have two less entries than argument list '\
                    '({len(self.arg_list)} entries).')
