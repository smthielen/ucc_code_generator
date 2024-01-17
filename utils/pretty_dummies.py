#
#summary of pretty dummies dicts
#
#SMT 5/2023
#


pretty_dummies_00 = {
    'above': 'abcdefg',
    'below': 'ijklmno',
    'general': 'pqrstu'
}


pretty_dummies_s = {
    'above': 'bcdefg',
    'below': 'jklmno',
    'general': 'pqrstu'
}


pretty_dummies_d = {
    'above': 'cdefg',
    'below': 'klmno',
    'general': 'pqrstu'
}


pretty_dummies_t = {
    'above': 'defgh',
    'below': 'lmnov',
    'general': 'pqrstu'
}


pretty_dummies_sd = {
    'above': 'defg',
    'below': 'lmno',
    'general': 'pqrstu'
}


pretty_dummies_dd = {
    'above': 'efgx',
    'below': 'mnoy',
    'general': 'pqrstu'
}


pretty_dummies_dt = {
    'above': 'fgx',
    'below': 'noy',
    'general': 'pqrstu'
}


pretty_dummies_tt = {
    'above': 'gh',
    'below': 'oy',
    'general': 'pqrstu'
}


pretty_dummies_dict = {'00': pretty_dummies_00, 's': pretty_dummies_s,
                       'd': pretty_dummies_d, 't': pretty_dummies_t,
                       'ss': pretty_dummies_d,
                       'sd': pretty_dummies_sd, 'ds': pretty_dummies_sd,
                       'dd': pretty_dummies_dd, 'dt': pretty_dummies_dt,
                       'td': pretty_dummies_dt, 'ts': pretty_dummies_dd,
                       'st': pretty_dummies_dd, 'tt': pretty_dummies_tt}

def return_pretty_dummies(t):
    if t not in pretty_dummies_dict:
        raise ValueError('No valid input for calculation: {exc_type}')
    return pretty_dummies_dict[t]
