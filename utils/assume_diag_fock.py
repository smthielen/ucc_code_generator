

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
