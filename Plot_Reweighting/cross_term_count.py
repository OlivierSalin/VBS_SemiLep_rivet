import numpy as np
import math
import json
import re
import itertools

all_operators=["FM0","FM1","FM2","FM3","FM4","FM5","FM7","FM8","FM9",
            "FS0","FS1","FS2",
            "FT0","FT1","FT2","FT3","FT4","FT5","FT6","FT7"]


cross_terms_ = [f"{op1}vs{op2}" for op1, op2 in itertools.combinations(all_operators, 2) if op1[:2] == op2[:2]]

# Count cross-terms per operator type
operator_types = {'FM': 0, 'FS': 0, 'FT': 0}
for ct in cross_terms_:
    if ct.startswith('FM'):
        operator_types['FM'] += 1
    elif ct.startswith('FS'):
        operator_types['FS'] += 1
    elif ct.startswith('FT'):
        operator_types['FT'] += 1
print('Number of cross-terms per operator type:')
for op_type, count in operator_types.items():
    print(f"{op_type}: {count}")
