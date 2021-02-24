#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 10:28:26 2021

@author: slaiad@ust.hk
"""

# -*- coding: utf-8 -*-

from Utils.Consts import MIN_TAG_LEN, MAX_TAG_LEN, AA_INFO, AA_RES_MASS
from Parameters import EXTRACT_AA_TOL
from pyteomics import mgf
from math import pi, cos, sqrt
from pyteomics import mzxml

def readPeaksFromXML(mzXmlPath):
    with mzxml.read(mzXmlPath) as spectra:
        specDict = {}
        for spectrum in spectra:
            # print(spectrum)
            # break
            specDict[int(spectrum['num'])] = [spectrum['precursorMz'][0]['precursorMz']*spectrum['precursorMz'][0]['precursorCharge'],
                                              [peak for peak in zip(spectrum['m/z array'], spectrum['intensity array'])]]
                                              
    return specDict

def scoreTag(massDif):
    s = 0.5*cos(pi/EXTRACT_AA_TOL*massDif) + 0.5
    return s

def findAA(deltaMz, lenCurTags):
    realAA = ''
    score = 0
    # realLow = 0
    for i in range(len(AA_INFO)):
        if deltaMz >= AA_INFO['resMass'][i] - EXTRACT_AA_TOL:
            # realLow = AA_INFO['resMass'][i] - EXTRACT_AA_TOL
            if i == len(AA_INFO) - 1:#if the end, dont go on judge
                realAA = AA_INFO['aaName'][i]
                score = scoreTag(abs(deltaMz - AA_INFO['resMass'][i]))
                break
        else:
            if deltaMz <= AA_INFO['resMass'][i - 1] + EXTRACT_AA_TOL:
                realAA = AA_INFO['aaName'][i - 1]
                score = scoreTag(abs(deltaMz - AA_INFO['resMass'][i - 1]))
            break
    
    if score < 0.7 and lenCurTags == 0:# if it is the start aa, it has to be precise
        realAA = ''
        
    if score < 0.2 :
        realAA = ''
    
    return realAA, round(score, 7)

def extractTagByCheckAA(mzs, startId, currentTag, currentScore, aaCanBeHead, extractedTags, lastTag, lastScore, involvedIds, involvedScores, needRevise): #maybe no return
    #last tag is for back-recursing, make currentTag back to lastTag
    tagModified = False  #for this start Id, no new AA has been added to tag yet
    lenCurTags = len(currentTag)  #note down origin length
    
    candidatesId = getCandidatesForPeak(mzs, startId)#take followId as new start
    if len(candidatesId) == 0 and currentTag != '': #if no candi, then add currentTag, and return last
        if len(currentTag) >= MIN_TAG_LEN:
            extractedTags.append([currentTag, round(currentScore, 7), involvedIds, involvedScores])
            # print(currentTag)
        return lastTag, round(lastScore, 7), involvedIds.rsplit(' ', 1)[0], involvedScores.rsplit(' ', 1)[0]
    # print('sId ', startId, ' can ids ', candidatesId)
    for candId in candidatesId:
        # print(candId, mzs[candId], startId, mzs[startId]) 
        realAA, score = findAA(mzs[candId] - mzs[startId], lenCurTags)

            
        if realAA != '':#recusion continues
            aaCanBeHead[candId] = False #do not use this as start outside 
            oriCurrentTag = currentTag
            oriCurrentScore = currentScore
            currentTag += realAA  #new AA found, Tag extends
            currentScore += score
            tagModified = True
            involvedIds += (' ' + str(candId))
            involvedScores += (' ' + str(score))
            
            # if (score - 0.9311951) < 1e-6 :
                # print(realAA)
                # print(currentTag)
                # print(candId, mzs[candId], startId, mzs[startId])
                # print(candId,candidatesId[-1])
            # if currentTag == 'DYY':
            #     print(currentScore)
            if needRevise:
                mzs[candId] = mzs[startId] + AA_RES_MASS[realAA]
                
            if len(currentTag) == MAX_TAG_LEN:
                extractedTags.append([currentTag, round(currentScore, 7), involvedIds, involvedScores])
                # print(currentTag)
                return lastTag, round(lastScore, 7), involvedIds.rsplit(' ', 1)[0], involvedScores.rsplit(' ', 1)[0]
            
            currentTag, currentScore, involvedIds, involvedScores = extractTagByCheckAA(mzs, candId, currentTag, currentScore, aaCanBeHead, extractedTags, oriCurrentTag, oriCurrentScore, involvedIds, involvedScores, needRevise)
            if candId == candidatesId[-1]:
                return lastTag, round(lastScore, 7), involvedIds.rsplit(' ', 1)[0], involvedScores.rsplit(' ', 1)[0]
        else: # recursion stops
                
            if candId == candidatesId[-1]:#if its last candi then may have to add and return
                if (len(currentTag) == lenCurTags) and currentTag != '' \
                    and tagModified == False and len(currentTag) >= MIN_TAG_LEN: #all candi not working
                    # if (score - 0.5883309) < 1e-7 and startId == 92 :
                    #     print(realAA)
                    #     print(currentTag)
                    #     print(candId, mzs[candId], startId, mzs[startId])
                    #     print(candId,candidatesId[-1])
                    #     print(tagModified)
                    #     print(len(currentTag), lenCurTags)
                    extractedTags.append([currentTag, round(currentScore, 7), involvedIds, involvedScores])
                    # print(currentTag)
                return lastTag, round(lastScore, 7), involvedIds.rsplit(' ', 1)[0], involvedScores.rsplit(' ', 1)[0]
   
def getCandidatesForPeak(mzs, startId):
    mz = mzs[startId]
    indexes = []
    
    for k in range(startId + 1, len(mzs)):
        deltaMz = mzs[k] - mz
        
        if deltaMz < AA_INFO['resMass'][0] - EXTRACT_AA_TOL:
            continue
        if deltaMz > AA_INFO['resMass'][-1] + EXTRACT_AA_TOL:
            break
        indexes.append(k)
    # print(len(indexes))
    return indexes
    
def normalizeScore(tagsTotal):
    for tag in tagsTotal:
        tag[1] = round(tag[1]/sqrt(len(tag[0])), 7)

    return sorted(tagsTotal, key = (lambda x:x[1]), reverse = True)
    

def readTagsFromMS2(mzs):

    tagsTotal = []

    # filter the noise  todo
    aaCanBeHead = [True]*len(mzs)

    for startId in range(len(mzs) - 1):
        # if aaCanBeHead[startId] == True: #peaks that have been added wont be head later
        if 1: #peaks that have been added wont be head later    
            #get the candidates for start peak, it is very different from the following peaks
            needRevise = False
            if startId == 0 or startId == 1:
                needRevise = True
            extractedTags = []
            extractTagByCheckAA(mzs, startId, '', 0, aaCanBeHead, extractedTags, '', 0, str(startId), str(0), needRevise)
            
            tagsTotal.extend(extractedTags)
            
    normedTags = normalizeScore(tagsTotal)
    return normedTags               

#%%MAIN
if __name__ == '__main__':
    mgfFilePath = 'G:/Database/synthetic/TUM_first_pool_54_01_01_2xIT_2xHCD-1h-R2-tryptic/01650b_BD2-TUM_first_pool_54_01_01-2xIT_2xHCD-1h-R2.mgf'
    mgfFile = mgf.read(mgfFilePath)

    i = 0
    scanIndexes = [49475]
    totalTags = []
    for scanIndex in scanIndexes:
        scanIndex = 47430
        
        ms2 = mgfFile[scanIndex]
        
        tags = []
        tags = readTagsFromMS2(ms2['m/z array'])
        totalTags.extend(tags)
        break
        

    for tag in totalTags:
        tag[1] = round(tag[1]/len(tag[0]), 7)

    totalTags = sorted(totalTags, key = (lambda x:x[1]), reverse = True)
    
    
    
    
