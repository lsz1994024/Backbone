# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 10:40:07 2021

@author: slaiad
"""

import pandas as pd
import re

import numpy as np

def foo(seq):
    for aa in seq:
        print(aa)
    return 1

if __name__ == '__main__':
    
    dtype = [('seq', np.unicode_, 6), ('pcMass', float)]
    pepsWithMass = np.array([], dtype)
    
    peps = ['DSA','DSA','DWQ']
    a = [tuple([pep, 2]) for pep in range(5)]
    # a = [[i, i+3.2] for i in range(5)]
    b = np.array([[i, i+3.2] for i in range(5)], dtype)
    pepsWithMass = np.array(a, dtype)
    
    
  