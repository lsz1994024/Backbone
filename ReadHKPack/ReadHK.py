# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 15:02:09 2021

@author: slaiad
"""
from ReadTagsPack.consts import MIN_HK_SPEC_LEN

def readFromHK(path):
    spectra = {}
    
    with open(path, 'r') as hk:
        line = hk.readline()
        
        lastSign = 'S'
        lastScanNo = 0
        while line is not None and line != '':
            
            curSign = line[0]
            if lastSign == 'S' and curSign == 'S':
                lastScanNo = int(line.split()[1])
                
                line = hk.readline()
                continue
            
            if lastSign == 'S' and curSign == 'P': 
                # newSpec = [lastScanNo]
                newSpec = []
                newSpec.append(float(line.split()[1]))
                lastSign = 'P'
                
                line = hk.readline()
                continue
            
            if lastSign == 'P' and curSign == 'P':
                newSpec.append(float(line.split()[1]))
                
                line = hk.readline()
                continue
            
            if lastSign == 'P' and curSign == 'S':
                if len(newSpec) >= MIN_HK_SPEC_LEN:
                    spectra[lastScanNo] = sorted(newSpec)
                lastScanNo = int(line.split()[1])
                lastSign = 'S'
                
                line = hk.readline()
                continue
    
    return spectra



if __name__ == '__main__':
    filePath = 'G:/Dataset/PXD022999/191122_MK_SIO13_P2-GM1.txt'
    spectra = readFromHK(filePath)
