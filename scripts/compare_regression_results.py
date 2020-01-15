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
        print("Files differ")
        return 1
    else:
        print("Files are the same")
        return 0
    

if __name__ == '__main__':
    import sys
    
    res = 1

    if len(sys.argv) != 3:
        print("Usage: compare_regression_results.py results_filename_1 results_filename_2")
    else:
        res = regression_test_analysis(sys.argv[1], sys.argv[2])
        
    sys.exit(res)
    
    
    
