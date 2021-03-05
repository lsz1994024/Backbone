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
from Utils.Funcs import binarySearch, binarySearchList
from ScorePack.Match import getPepCand, cleanUpTags, findPtm
from Utils.Funcs import divideInputList
from time import time
from Parameters import *
from suffix_trees import STree
import re

#test
from ProcessDbPack.ParsePep import calcuSeqMass
#% functions

def getFeasiblePepTest(pcMass, allPeps, massShift):
    lb = binarySearch(allPeps['pcMass'], 0, len(allPeps)-1, pcMass - massShift)
    ub = binarySearch(allPeps['pcMass'], 0, len(allPeps)-1, pcMass + 1)
    return list(allPeps['seq'][lb : ub])

def getFeasiblePepIndex(tag):
    global bigString, wordStarts, lenWordStarts
    # print(len(bigString))
    forwTagIndexs = [m.start() for m in re.finditer(tag, bigString)]
    backTagIndexs = [m.start() for m in re.finditer(tag[::-1], bigString)]
    forwPepIndexs = binarySearchList(wordStarts, 0, lenWordStarts, forwTagIndexs)
    backTagIndexs = binarySearchList(wordStarts, 0, lenWordStarts, backTagIndexs)
    return list(set(forwPepIndexs+backTagIndexs))
    
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
        # print(mzs)
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
        # cleanUpTags(reliableTags)
        
        pepCandDict = dict.fromkeys([tag[0] for tag in reliableTags], [])
        for tag in pepCandDict:
            pepIndexs = getFeasiblePepIndex(tag)
            pepCandDict[tag] = [allPeps[i] for i in pepIndexs]
            
        # pepIndexs = getFeasiblePepIndex(tag)
        # print(len(feasiblePeps))
        
        # print(reliableTags)
        ptm = findPtm(mzs, reliableTags, pepCandDict, pcMass, MS1_TOL)
        for p in ptm:
            psms.append([scanNo, p])  
    return psms
    
    
def doSearch(specDict, allPeps):
    
    # tagDict = {}
    # pepCandDict = {}
    scanTagCand = {}
    
    scanNoList = [key for key in specDict]
    
    # scanNoList = [47780, 1844, 2608]#1844,2462,2608,13540,18113,37531,40813,40972,44676,46721,47401,47547,47607] #test
    #1844,2462,2608, C   47401 no   other M
    dividedScanList = divideInputList(scanNoList, 2000)
    
    print('len of scanNo list ', len(scanNoList))
    allPsms = []
    print("lsz start computation")
    with cf.ProcessPoolExecutor(max_workers = 88) as executor:
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
    
    allPeps = []
    
    with cf.ProcessPoolExecutor(max_workers = 88) as executor:
        for result in executor.map(getAllPepsForParallel, dividedList):
            allPeps.extend(result)

    allPeps = list(set(allPeps))       
    
    return allPeps

def buildOneTree(pepList):
    global allPeps
    return STree.STree(list(allPeps[pepList[0]:pepList[-1]]))
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
    #%% build trees
    
    bigString = '.'.join(allPeps)
    wordStarts = [0]
    last = 0
    for pep in allPeps:
        # print(last)
        last += len(pep)+1
        wordStarts.append(last)
    lenWordStarts = len(wordStarts)
    #%% doSearch
    # runcell('import', '/home/slaiad/Code/TagTree/main.py')
    from Parameters import *
    print("start Searching")
    
    t4 = time()


    allPsms = doSearch(specDict, allPeps)
    t5 = time()
    
    print("finish searching ", t5-t4)
    from Utils.Funcs import sendEmail
    sendEmail()
    
    #% write
    # psms.append( tuple((tag, pep, pcMassDiff, modPosLb, modPosRb, Lbs, Rbs, numMod)) )
    dataPsms = pd.DataFrame({'scanNo'   : [psm[0] for psm in allPsms],
                             'tags'     : [psm[1][0] for psm in allPsms],
                             'pep'      : [psm[1][1] for psm in allPsms],
                             'massdiff' : [psm[1][2] for psm in allPsms], 
                             'modPosLb' : [psm[1][3] for psm in allPsms],
                             'modPosRb' : [psm[1][4] for psm in allPsms],
                             'lbs'      : [psm[1][5] for psm in allPsms],
                             'rbs'      : [psm[1][6] for psm in allPsms],
                             'numMod'   : [psm[1][7] for psm in allPsms]})
    filePath = pd.ExcelWriter(r'TempData/ptmRevise.xlsx')
    dataPsms.to_excel(r'TempData/ptmRevise.xlsx')

    #%% verification
    # dataMascot = pd.read_excel('/home/slaiad/Code/TagTree/testData/MK_SIO13_P2_GM1.xlsx')
    
    # masDict = {}
    # # cometModDict = {}
    # for i in range(len(dataMascot['scanNo'])):
    #     scanNo = dataMascot['scanNo'][i]
    #     pep = dataMascot['pep'][i]
    #     if scanNo in masDict:
    #         masDict[scanNo].append(pep)
    #     else:
    #         masDict[scanNo] = []     
    #         masDict[scanNo].append(pep)
        
    # result = []
    # numSameMascot = 0
    # numWrong      = 0
    # numNoInMasCot = 0
    # numFirstRight = 0
    # numTwoRight   = 0
    # for psm in allPsms:
    #     scanNo = psm[0]
    #     tags = psm[1]
    #     pepCand = psm[2]
    #     lenFeasPeps = psm[3]
        
    #     if scanNo not in masDict:
    #         numNoInMasCot += 1
    #         result.append([scanNo, pepCand, tags, [], 'notMas', lenFeasPeps, 'notMas'])
    #         continue
        
    #     pepsTrue = [p.replace("L", "I") for p in masDict[scanNo]]
        
    #     if len(pepCand) == 0:
    #         numWrong += 1
    #         result.append([scanNo, pepCand, tags, pepsTrue, False, lenFeasPeps, False])
    #         continue
        
        
    #     pepsTrueUpper = [p.upper() for p in pepsTrue]
    #     pepCandUpper = [cand.upper() for cand in pepCand]
        
        
    #     if set(pepsTrueUpper) & set(pepCandUpper) != set():
    #         numSameMascot += len(set(pepsTrueUpper) & set(pepCandUpper))
            
    #         result.append([scanNo, pepCand, tags, pepsTrue, True, lenFeasPeps, pepCandUpper[0] in pepsTrueUpper])
    #         if pepCandUpper[0] in pepsTrueUpper:
    #             numFirstRight += 1
            
    #         if (len(pepCandUpper) >= 2) and (pepCandUpper[0] in pepsTrueUpper or pepCandUpper[1] in pepsTrueUpper):
    #             numTwoRight += 1
                
    #     else:
    #         numWrong += 1
    #         result.append([scanNo, pepCand, tags, pepsTrue, False, lenFeasPeps, False])
            
    # print('First Right', numFirstRight)
    # print('Two right', numTwoRight)
    # print('Contains right', numSameMascot)
    # print('Wrong', numWrong)
    # print('Not in Mascot', numNoInMasCot)
    # print('mascot', len(dataMascot['scanNo']))
    
    # #%
    # dataVeri = pd.DataFrame({'scanNo'     : [res[0] for res in result],
    #                          'pepCand'    : [res[1] for res in result],
    #                          'tags'       : [res[2] for res in result],
    #                          'truth'      : [res[3] for res in result],
    #                          'inOrnot'    : [res[4] for res in result],
    #                          'numFesiPep' : [res[5] for res in result],
    #                          'firstRight' : [res[6] for res in result]
    #                          })#, orient='scanNo',columns=['Tags', 'Candidates'])
    # # filePath = pd.ExcelWriter(r'TempData/dataVeri.xlsx')
    # dataVeri.to_excel(r'TempData/tag1.6_AA0.5_pcTol10ppm_simiXreliabilityTop10LCS_QeGA_nTerm42p011_varOxiM15p99_tagClustered.xlsx')
    
    
    # #%%test for suffix tree
    # from suffix_trees import STree
    # peps = list(allPeps['seq'][0:10000])
    
    # t1 = time()
    # st = STree.STree(peps)
    # t2 = time()
    # print('build the tree', t2-t1)
    # #%%
    # tag = 'GTA'
    # t22 = time()
    # index = st.find_all(tag)
    # t3 = time()
    # print('find tags', t3-t22)
    # wordStarts = []
    # wordStarts = st.word_starts
    # wordStarts.append(len(peps[-1]) + wordStarts[-1])
    # cand = []
    # idds = []
    # for i in index:
    #     # print(i)
    #     # print(peps[binarySearch(wordStarts, 0, len(wordStarts), i)])
    #     try:
    #         cand.append(peps[binarySearch(wordStarts, 0, len(wordStarts), i)])
    #         idds.append(binarySearch(wordStarts, 0, len(wordStarts), i))
    #     except:
    #         print('error', i)
    # t4 = time()
    # print('normal binasea', t4-t3)
    
    # ids = binarySearchList(wordStarts, 0, len(wordStarts), sorted(list(index)))
    # tt4 = time()
    # print('list binasea', tt4-t4)
    # print(set(idds) == set(ids))
    # print('len cand', len(cand))
    # # for i in index:
    # #     print(i, binarySearch(wordStarts, 0, len(wordStarts), i), peps[binarySearch(wordStarts, 0, len(wordStarts), i)])    
    # # %%
    # peps = list(allPeps['seq'])
    # tt1 = time()
    # candV2 = []
    # tag = 'GTHR'
    # for pep in peps:
    #     if tag in pep:
    #         candV2.append(pep)
    # tt2 = time()
    # print('traverse ', tt2-tt1)  
    # print(len(candV2))
    
    # bigString = '.'.join(peps)
    # wordStarts = [0]
    # last = 0
    # for pep in allPeps['seq']:
    #     # print(last)
    #     last += len(pep)+1
    #     wordStarts.append(last)
    
    # tag = 'DEF'
    # tt3 = time()
    # a = [m.start() for m in re.finditer(tag, bigString)]
    # ids = binarySearchList(wordStarts, 0, len(wordStarts), a)
    # tt4 = time()
    # print('join', tt4-tt3)
    
    
    
    
   

