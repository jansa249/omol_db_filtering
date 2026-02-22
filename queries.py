"""
Dictionary for SMARTS queries for interesting molecule substructures.
Every entry here is used by other scripts.
Custom entries can be added after 'phosphate'.
"""
QUERY_DICT = {
    'adenine' : 'N1~C~N~C2~C1~N~C~N~C2~[NH2]',  # mandatory
    'guanine' : 'N1~C~N~C2~C1~N~C(~N)~N~C2~O',  # mandatory
    'uracthym' : 'N1~C~C~C(~O)~N~C1(~O)',       # mandatory
    'cytosine' : 'N1~C~C~C(~[NH2])~N~C1(~O)',   # mandatory
    'dribose_n' : 'O~C~C1~C(~O)~C~C(~N)~O1',    # mandatory
    'dribose' : 'C~C1~C(~O)~C~C~O1',            # mandatory
    'ribose' : 'C~C1~C(~O)~C(~O)~C~O1',         # mandatory
    'phosphate' : 'P(~O)(~O)(~O)(~O)',          # mandatory
}

