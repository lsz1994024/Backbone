# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 10:40:07 2021

@author: slaiad
"""

import pandas as pd
import re

import numpy as np
from time import time
def foo(seq):
    for aa in seq:
        print(aa)
    return 1

if __name__ == '__main__':  
    
    
    tag = 'SAW'
    
    seq = 'DWSADFSAJPOWQDJPINFADSKMSAWFDPOKFDAOIJOPSDJDSA'
    
    t1 = time()
    a = tag in seq
    
    t2 = time()
    
    b = re.match(tag, seq)
    t3 = time()
    
    print(t2-t1)
    print(t3-t2)
    
    
  