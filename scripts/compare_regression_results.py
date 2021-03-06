#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 
    Compare expected actual regression test output.
    Author: Lorne Whiteway.
"""


# See LW's notes p. GL154

def regression_test_analysis(file1, file2):
    import healpy as hp
    import numpy as np
    import sys
    
    m1 = hp.read_map(file1, verbose=False)
    m2 = hp.read_map(file2, verbose=False)
    
    if np.max(np.abs(m1-m2)) > 1e-8:
        print("Files {} and {} DIFFER".format(file1, file2))
        return 1
    else:
        print("Files {} and {} are the same (up to tolerance)".format(file1, file2))
        return 0
    

if __name__ == '__main__':
    import sys
    import traceback
    
    res = 1
    
    try:
        if len(sys.argv) != 3:
            print("Usage: compare_regression_results.py results_filename_1 results_filename_2")
        else:
            res = regression_test_analysis(sys.argv[1], sys.argv[2])
            
        sys.exit(res)
        
    except Exception as err:

        print(traceback.format_exc())
        sys.exit(res)
        
