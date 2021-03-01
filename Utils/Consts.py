# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 17:07:55 2020

@author: slaiad
"""

import numpy as np

from pyteomics import mass
from math import sqrt, exp, pi, cos

MAX_TAG_LEN = 20
MIN_TAG_LEN = 3
NUM_PERMUTATION = 8

ATOM_MASS = \
    { #mono  https://www.unimod.org/masses.html
      'C' : 12,
      'N' : 14.003074,
      'O' : 15.99491463,
      'S' : 31.9720707,
      'P' : 30.973762,
      'H' : 1.007825035,
      'D' : 2.01410,
      'e' : 1.007825035/(1836+1),
      'Z' : 1.007825035*1836/(1836+1)
      }
    
H2O = ATOM_MASS['H']*2 + ATOM_MASS['O']
NH3 = ATOM_MASS['H']*3 + ATOM_MASS['N']
# def LOG()
AA = ['G', 'A', 'S', 'P','V','T','C','I','N','D','Q','K','E','M','H','F','R','Y','W']



ATOM_VAN_AREA = \
{
  'C' : 4*pi*1.7**2,
  'N' : 4*pi*1.55**2,
  'O' : 4*pi*1.52**2,
  'S' : 4*pi*1.85**2,
  'P' : 4*pi*1.9**2
}

    



    
FONT = {'family' : 'Times New Roman',
'weight' : 'normal',
'size'   : 18,
}


    
AA_S_NAME = \
    {
     'ALA' : 'A',
     'ARG' : 'R',
     'ASN' : 'N',
     'ASP' : 'D',
     'CYS' : 'C',
     'GLU' : 'E',
     'GLN' : 'Q',
     'GLY' : 'G',
     'HIS' : 'H',
     # 'HYP' : 'O',
     'ILE' : 'I',
     'LEU' : 'I', # NOTICE!
     'LYS' : 'K',
     'MET' : 'M',
     'PHE' : 'F',
     'PRO' : 'P',
     # 'GLP' : 'U',
     'SER' : 'S',
     'THR' : 'T',
     'TRP' : 'W',   
     'TYR' : 'Y',
     'VAL' : 'V'
     }
    
NOT_R_TOPO = ['CA', 'C', 'N', 'O']
    

# print(EXTRACT_AA_TOL)

AA_RES_MASS = {}
for i in AA:
    AA_RES_MASS[i] = mass.std_aa_mass[i]
    
#PXD013040  bsa TREATED 1     it works good for PXD022999
# AA_RES_MASS['m'] = AA_RES_MASS['M'] + 15.99491463
# AA_RES_MASS['C'] += 57.021464



#PXD018758 ACETYL
# AA_RES_MASS['m'] = AA_RES_MASS['M'] + 15.995
# AA_RES_MASS['c'] = AA_RES_MASS['C'] + 57.021464
# AA_RES_MASS['k'] = AA_RES_MASS['K'] + 42.0106

PRECISION = 3

dtype = [('aaName', np.unicode_, 1), ('resMass', float)]
RAW_AA_INFO = np.array([],dtype)

for aa in AA_RES_MASS:
    aaInfo = np.array([(aa, AA_RES_MASS[aa])],dtype)
    RAW_AA_INFO = np.concatenate((RAW_AA_INFO, aaInfo), axis = 0)
    
AA_INFO = np.sort(RAW_AA_INFO, order = ['resMass'])

if __name__ == '__main__':
    a = 'QFASQANVVGPWIQTK'[::-1]
    mass = []
    temp = AA_RES_MASS['Q'] + ATOM_MASS['H']*2 + ATOM_MASS['O']
    print(temp)
    for i in range(1,10):
        temp += AA_RES_MASS[a[i]]
        # print(a[i], i)
    
        mass.append(temp)
    
    # mass = 0
    # for i in a:
    #     mass += AA_RES_MASS[i]
    print(mass)

