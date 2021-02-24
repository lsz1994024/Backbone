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
    

# def LOG()
MU = 0
def gaussian(x, SIGMA2):
    
    f = 1/sqrt(2*pi*SIGMA2)*exp(-(x-MU)**2/2/SIGMA2)
    return f

# x = np.linspace(0, EXTRACT_AA_TOL, 50)
# y = []
# for i in x:
#     # y.append(gaussian(i, 0.05)/gaussian(0,0.05))
#     y.append(0.5*cos(pi/EXTRACT_AA_TOL*i) + 0.5)
# import matplotlib.pyplot as plt
# plt.plot(x,y)
ATOM_VAN_AREA = \
{
  'C' : 4*pi*1.7**2,
  'N' : 4*pi*1.55**2,
  'O' : 4*pi*1.52**2,
  'S' : 4*pi*1.85**2,
  'P' : 4*pi*1.9**2
}

    
R_PART_NORM = \
{
'G' : 1,
'A' : ATOM_VAN_AREA['C']*1,
'S' : ATOM_VAN_AREA['C']*1 + ATOM_VAN_AREA['N']*0 + ATOM_VAN_AREA['O']*1,
'P' : ATOM_VAN_AREA['C']*3 + ATOM_VAN_AREA['N']*0 + ATOM_VAN_AREA['O']*0,
'V' : ATOM_VAN_AREA['C']*3 + ATOM_VAN_AREA['N']*0 + ATOM_VAN_AREA['O']*0,
'T' : ATOM_VAN_AREA['C']*2 + ATOM_VAN_AREA['N']*0 + ATOM_VAN_AREA['O']*1,
'C' : ATOM_VAN_AREA['C']*1 + ATOM_VAN_AREA['N']*0 + ATOM_VAN_AREA['O']*0 + ATOM_VAN_AREA['S']*1,
'O' : ATOM_VAN_AREA['C']*3 + ATOM_VAN_AREA['N']*0 + ATOM_VAN_AREA['O']*1,
'I' : ATOM_VAN_AREA['C']*4 + ATOM_VAN_AREA['N']*0 + ATOM_VAN_AREA['O']*0,
'N' : ATOM_VAN_AREA['C']*2 + ATOM_VAN_AREA['N']*1 + ATOM_VAN_AREA['O']*1,
'D' : ATOM_VAN_AREA['C']*2 + ATOM_VAN_AREA['N']*0 + ATOM_VAN_AREA['O']*2,
'U' : ATOM_VAN_AREA['C']*3 + ATOM_VAN_AREA['N']*0 + ATOM_VAN_AREA['O']*1,    
'Q' : ATOM_VAN_AREA['C']*3 + ATOM_VAN_AREA['N']*1 + ATOM_VAN_AREA['O']*1,
'K' : ATOM_VAN_AREA['C']*4 + ATOM_VAN_AREA['N']*1 + ATOM_VAN_AREA['O']*0,
'E' : ATOM_VAN_AREA['C']*3 + ATOM_VAN_AREA['N']*0 + ATOM_VAN_AREA['O']*2,
'M' : ATOM_VAN_AREA['C']*3 + ATOM_VAN_AREA['N']*0 + ATOM_VAN_AREA['O']*0 + ATOM_VAN_AREA['S']*1,
'H' : ATOM_VAN_AREA['C']*4 + ATOM_VAN_AREA['N']*2 + ATOM_VAN_AREA['O']*0,
'F' : ATOM_VAN_AREA['C']*7 + ATOM_VAN_AREA['N']*0 + ATOM_VAN_AREA['O']*0,
'R' : ATOM_VAN_AREA['C']*4 + ATOM_VAN_AREA['N']*3 + ATOM_VAN_AREA['O']*0,
'Y' : ATOM_VAN_AREA['C']*7 + ATOM_VAN_AREA['N']*0 + ATOM_VAN_AREA['O']*1,
'W' : ATOM_VAN_AREA['C']*9 + ATOM_VAN_AREA['N']*1 + ATOM_VAN_AREA['O']*0 
}


#from wiki https://zh.wikipedia.org/wiki/%E6%B0%A8%E5%9F%BA%E9%85%B8
AA_MASS = \
    {
    'G' : 75.06714,
    'A' : 89.09404,
    'S' : 105.09344,
    'P' : 115.13194,
    'V' : 117.14784,
    'T' : 119.12034,
    'C' : 121.15404,
    # 'O' : 131.13, no data
    'I' : 131.17464, #same as L
    'N' : 132.11904,
    'D' : 133.10384,
    # 'U' : 139.11,    no data
    'Q' : 146.14594,
    'K' : 146.18934,
    'E' : 147.13074,
    'M' : 149.20784,
    'H' : 155.15634,
    'F' : 165.19184,
    'R' : 174.20274,
    'Y' : 181.19124,
    'W' : 204.22844
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

AA = [key for key in AA_MASS]
AA_RES_MASS = {}
for i in AA:
    AA_RES_MASS[i] = mass.std_aa_mass[i]
    
#PXD013040  bsa TREATED 1     it works good for PXD022999
AA_RES_MASS['m'] = AA_RES_MASS['M'] + 15.99491463
AA_RES_MASS['C'] += 57.021464



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
# aaResMassTolRangesUpper = np.array([(AA_RES_MASS[key] + EXTRACT_AA_TOL) for key in AA_RES_MASS])
# aaResMassTolRangesLower = np.array([(AA_RES_MASS[key] - EXTRACT_AA_TOL) for key in AA_RES_MASS])
# aaKeys = np.array(list(AA_RES_MASS.keys()))

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

