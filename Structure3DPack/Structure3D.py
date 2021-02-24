# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 22:15:47 2020

@author: slaiad
"""

import numpy as np
from ReadTagsPack.consts import PRECISION, AA_S_NAME
import os
import pandas as pd
import concurrent.futures as cf

def mergeAtoms2AA(atomCoors):
    startI = 0
    endI = 0
    dtype = [('aaId', int), ('x', float), ('y', float), ('z', float)]
    aaCoors = np.array([], dtype)  
        
    while endI != len(atomCoors['aaId']):
        if atomCoors['aaId'][endI + 1] != atomCoors['aaId'][endI]:
            xAve = np.average(atomCoors['x'][startI : endI+1])
            yAve = np.average(atomCoors['y'][startI : endI+1])
            zAve = np.average(atomCoors['z'][startI : endI+1])
            aaCoor = np.array([(atomCoors['aaId'][startI], round(xAve, PRECISION), round(yAve, PRECISION), round(zAve, PRECISION))], dtype)
            aaCoors = np.concatenate((aaCoors, aaCoor), axis = 0)
            
            startI = endI + 1
            endI = endI + 1


        else:
            endI = endI + 1
            if endI == len(atomCoors['aaId']) - 1:
                xAve = np.average(atomCoors['x'][startI : endI+1])
                yAve = np.average(atomCoors['y'][startI : endI+1])
                zAve = np.average(atomCoors['z'][startI : endI+1])
                aaCoor = np.array([(atomCoors['aaId'][startI], round(xAve, PRECISION), round(yAve, PRECISION), round(zAve, PRECISION))], dtype)
                aaCoors = np.concatenate((aaCoors, aaCoor), axis = 0)
                break
            
    return aaCoors

def writeToCsv(aaCoors, file):
    dataFrame = pd.DataFrame({'aaId' : aaCoors['aaId'], 'x' : aaCoors['x'], 'y' : aaCoors['y'], 'z' : aaCoors['z']})
    
    toFile = '../TempData/ConvertedAACoors/%s.csv' % file.split('.')[0] 
    dataFrame.to_csv(toFile, index = False, header = False, sep = ',')
    
def convertOneFile(file):
    fullFilePath = r'D:\oneDrive\OneDrive - HKUST Connect\structure\%s' % file
    with open(fullFilePath, 'r') as fid:
        line = fid.readline()

        dtype = [('aaId', int), ('x', float), ('y', float), ('z', float)]
        atomCoors = np.array([], dtype)          
        
        while line is not None and line != '':
            
            line = line.split()
            
            if len(line) < 10 or line[0] != 'ATOM':
                line = fid.readline()
                continue
            
            if line[5] not in AA_S_NAME:
                line = fid.readline()
                continue
            
            #ATOM   4364 C  CB  . TYR D 2 145 ? 73.443  56.078 -0.943  1.00 17.36  ? 145 TYR D CB  1 
            newAA = np.array([(line[8], line[10], line[11], line[12])], dtype)
            atomCoors = np.concatenate((atomCoors, newAA), axis = 0)
            
            line = fid.readline()
        aaCoors = mergeAtoms2AA(atomCoors)
        
        writeToCsv(aaCoors, file)
    
def convertCifs2FakeImgs(files):
    #from cif file, convert ATOM rows into fake 3D images
    
    errorFiles = []
    for file in files:
        if os.path.exists('../TempData/ConvertedAACoors/%s.csv' % file) :
            continue
        
        try:
            convertOneFile(file)
        except:
            errorFiles.append(file)
        else:
            continue  
    return errorFiles
        
        

if __name__ == '__main__':
    filesFolder = 'D:\oneDrive\OneDrive - HKUST Connect\structure'
    files = os.listdir(filesFolder)
    files = [x for x in files if x[-3:] == 'cif' and len(x) == 8]
    # files = ['1abi.cif']
    # totalErrorList = convertCifs2FakeImgs(files)
    
    lenFiles = len(files)
    inputList = []
    # totalList = []
    for i in range(40):
        startI = i*1200
        endI = (i+1) * 1200
        if endI > lenFiles:
            endI = lenFiles
        inputList.append(files[startI : endI])
        # totalList.extend(files[startI : endI])
    
    totalErrorList = []
    with cf.ProcessPoolExecutor(max_workers = 8) as executor:
        results = executor.map(convertCifs2FakeImgs, inputList)
    # futures = [executor.submit(getBucket, hashMat)]
    # for future in cf.as_completed(futures):
    #     totalBucketSet = future.result()
    for result in results:
        totalErrorList.extend(result)
    
    # freesasa.ResidueArea()
    