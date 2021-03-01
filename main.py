# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 20:03:35 2020

@author: slaiad
"""

#%% import
from ProcessDbPack.ParsePep import getAllProts, getAllPepsForParallel
from ProcessSpecPack.ParseMsms import readPeaksFromXML, readTagsFromMS2
from ProcessSpecPack.PreProcess import deIsotope
import concurrent.futures as cf
from itertools import repeat
import pandas as pd
import numpy as np
from Utils.Funcs import binarySearch
from ScorePack.Match import getPepCand, cleanUpTags, findPtm
from Utils.Funcs import divideInputList
from time import time
from Parameters import *

#test
from ProcessDbPack.ParsePep import calcuSeqMass
#% functions

def getFeasiblePep(pcMass, allPeps, massShift):
    lb = binarySearch(allPeps['pcMass'], 0, len(allPeps)-1, pcMass - massShift)
    ub = binarySearch(allPeps['pcMass'], 0, len(allPeps)-1, pcMass + 1)
    return list(allPeps['seq'][lb : ub])
    
def searchSome(scanNoList):
    global specDict, allPeps
    psms = []
    # i = 0
    for scanNo in scanNoList:
            
        pcMass = specDict[scanNo][0]
        rawSpec = specDict[scanNo][1]
        mzs = deIsotope(rawSpec)
        
        if len(mzs) == 0:
            continue
        
        
        tags = readTagsFromMS2(mzs)
        print(mzs)
        # print('finish tag', tags)
            
        reliableTags = []
        for tag in tags:
            
            if tag[1] < TAG_SCORE_THRES:
                break
            
            reliableTags.append(tuple(tag[0:3]))
                
            if len(reliableTags) == MAX_TAGS_NUM:
                break
            
        if len(reliableTags) == 0:
            continue
        
        # print(reliableTags)
        cleanUpTags(reliableTags)
        
        feasiblePeps = getFeasiblePep(pcMass, allPeps, MASS_SHIFT)
        
        # print(reliableTags)
        ptm = findPtm(mzs, reliableTags, feasiblePeps, pcMass, MS1_TOL)
        psms.extend(ptm)  
        print('scan %d ' % scanNo, psms)
    return psms
    
    
def doSearch(specDict, allPeps):
    
    # tagDict = {}
    # pepCandDict = {}
    scanTagCand = {}
    
    scanNoList = [key for key in specDict]
    
    scanNoList = [2608]#1844,2462,2608,13540,18113,37531,40813,40972,44676,46721,47401,47547,47607] #test
    #1844,2462,2608, C   47401 no   other M
    dividedScanList = divideInputList(scanNoList, 14)
    
    print('len of scanNo list ', len(scanNoList))
    allPsms = []
    print("lsz start computation")
    with cf.ProcessPoolExecutor(max_workers = 14) as executor:
        # results = executor.map(searchOne, scanNoList)
        for result in executor.map(searchSome, dividedScanList):
        # for result in results:
            # if len(result[2]) != 0:
            allPsms.extend(result)
                # print(result)
                # break
                # scanTagCand[result[0]] = [result[1], result[2]]
    print('lsz finished computation')

    return allPsms

def parseDB(protList):
    dividedList = divideInputList(protList, 100)
    
    dtype = [('seq', np.unicode_, MAX_LENGTH_OF_PEP), ('pcMass', float)]
    allPeps = np.array([], dtype)
    
    with cf.ProcessPoolExecutor(max_workers = 88) as executor:
        for result in executor.map(getAllPepsForParallel, dividedList):
            allPeps = np.concatenate((allPeps, result), axis = 0)

    allPeps = np.unique(allPeps)        
    allPeps = (np.sort(allPeps, order = ['pcMass']))
    
    return allPeps

#%% Read DB and Spectra
if __name__ == '__main__':

    
    protList = getAllProts(DB_PATH)
    # allPeps = getAllPeps(protList)
    t2 = time()
    
    
    allPeps = parseDB(protList)
    
    t3 = time()        
    print("finishe parsing database ", t3-t2)
    print("Number of all peps = ", len(allPeps))
    
    specDict = readPeaksFromXML(SPEC_PATH)
    
    #%% doSearch
    # runcell('import', '/home/slaiad/Code/TagTree/main.py')
    from Parameters import *
    print("start Searching")
    
    t4 = time()
    allPsms = doSearch(specDict, allPeps)
    t5 = time()
    
    print("finish searching ", t5-t4)
    from Utils.Funcs import sendEmail
    # sendEmail()
    #% write
    # dataPsms = pd.DataFrame({'scanNo' : [psm[0] for psm in allPsms],
    #                                 'tags'   : [psm[1] for psm in allPsms],
    #                                 'seqs'   : [psm[2] for psm in allPsms]})#, orient='scanNo',columns=['Tags', 'Candidates'])
    # filePath = pd.ExcelWriter(r'TempData/dataPsms.xlsx')
    # dataPsms.to_excel(r'TempData/dataPsms.xlsx')
    
    #%% verification
    dataMascot = pd.read_excel('/home/slaiad/Code/TagTree/testData/MK_SIO13_P2_GM1.xlsx')
    
    masDict = {}
    # cometModDict = {}
    for i in range(len(dataMascot['scanNo'])):
        scanNo = dataMascot['scanNo'][i]
        pep = dataMascot['pep'][i]
        if scanNo in masDict:
            masDict[scanNo].append(pep)
        else:
            masDict[scanNo] = []     
            masDict[scanNo].append(pep)
        
    result = []
    numSameMascot = 0
    numWrong      = 0
    numNoInMasCot = 0
    numFirstRight = 0
    numTwoRight   = 0
    for psm in allPsms:
        scanNo = psm[0]
        tags = psm[1]
        pepCand = psm[2]
        lenFeasPeps = psm[3]
        
        if scanNo not in masDict:
            numNoInMasCot += 1
            result.append([scanNo, pepCand, tags, [], 'notMas', lenFeasPeps, 'notMas'])
            continue
        
        pepsTrue = [p.replace("L", "I") for p in masDict[scanNo]]
        
        if len(pepCand) == 0:
            numWrong += 1
            result.append([scanNo, pepCand, tags, pepsTrue, False, lenFeasPeps, False])
            continue
        
        
        pepsTrueUpper = [p.upper() for p in pepsTrue]
        pepCandUpper = [cand.upper() for cand in pepCand]
        
        
        if set(pepsTrueUpper) & set(pepCandUpper) != set():
            numSameMascot += len(set(pepsTrueUpper) & set(pepCandUpper))
            
            result.append([scanNo, pepCand, tags, pepsTrue, True, lenFeasPeps, pepCandUpper[0] in pepsTrueUpper])
            if pepCandUpper[0] in pepsTrueUpper:
                numFirstRight += 1
            
            if (len(pepCandUpper) >= 2) and (pepCandUpper[0] in pepsTrueUpper or pepCandUpper[1] in pepsTrueUpper):
                numTwoRight += 1
                
        else:
            numWrong += 1
            result.append([scanNo, pepCand, tags, pepsTrue, False, lenFeasPeps, False])
            
    print('First Right', numFirstRight)
    print('Two right', numTwoRight)
    print('Contains right', numSameMascot)
    print('Wrong', numWrong)
    print('Not in Mascot', numNoInMasCot)
    print('mascot', len(dataMascot['scanNo']))
    
    #%
    dataVeri = pd.DataFrame({'scanNo'     : [res[0] for res in result],
                             'pepCand'    : [res[1] for res in result],
                             'tags'       : [res[2] for res in result],
                             'truth'      : [res[3] for res in result],
                             'inOrnot'    : [res[4] for res in result],
                             'numFesiPep' : [res[5] for res in result],
                             'firstRight' : [res[6] for res in result]
                             })#, orient='scanNo',columns=['Tags', 'Candidates'])
    # filePath = pd.ExcelWriter(r'TempData/dataVeri.xlsx')
    dataVeri.to_excel(r'TempData/tag1.6_AA0.5_pcTol10ppm_simiXreliabilityTop10LCS_QeGA_nTerm42p011_varOxiM15p99_tagClustered.xlsx')
    
    
    
   

