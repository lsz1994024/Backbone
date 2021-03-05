#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 09:53:32 2021

@author: slaiad@ust.hk
"""

from fuzzywuzzy import fuzz
from ProcessDbPack.ParsePep import getTheoPeaks, calcuSeqMass
import numpy as np
from Parameters import NUM_EXTRA_MOD_PEAKS
from Utils.Consts import ATOM_MASS

def getIonCoverage(ions, mzs, massTol, lenTag):
    numMatch = 0
    for mz in mzs:
        minDiff = 100
        for ion in ions:
            if abs(mz - ion) < minDiff:
                minDiff = abs(mz - ion)
        
        if minDiff <= massTol:
            numMatch += 1
            # print(mz)
    # print(numMatch - 2, len(mzs) - lenTag - 2)
    # return float(numMatch - 2)/(len(mzs) - lenTag - 2)
    return (numMatch - 2)
    
def getModIons(originIons, mzs, massTol):
    
    numMod = 0
    for mz in mzs:
        for l in range(len(originIons)):
            if abs(mz - originIons[l]) < massTol:
                numMod += 1
                
    return numMod

def getModIonsV2(modInL, oriLIonsForMod, oriRIonsForMod, oriLIonsForCov, oriRIonsForCov, mzs, massTol, pcMassDiff, modPosL, modPosR):
    numModL = 0
    numModR = 0
    # lenMod = modPosR - modPosL
    # maxL = lenMod
    # minR = lenMod
    L = modPosL
    R = modPosR
    
    Lb = [L]
    Rb = [R]
    Lbsta = 'o'
    Rbsta = 'o'
    
    
    # print('minR, maxL',minR, maxL)
    # print(pcMassDiff)
    for mz in mzs:
        #mod
        for l in range(len(oriLIonsForMod)):
            # print(abs(eIon - oriLIons[l] - pcMassDiff) , massTol)
            if abs(mz - oriLIonsForMod[l][0] - pcMassDiff) < massTol or abs(mz - oriLIonsForMod[l][1] - pcMassDiff) < massTol or abs(mz - oriLIonsForMod[l][2] - pcMassDiff) < massTol:
                numModL += 1
                if modPosL + l < R:
                    R = modPosL + l
                    Rb.append(R)
                    Rbsta += 'm'
            if abs(mz - oriLIonsForMod[l][0]) < massTol or abs(mz - oriLIonsForMod[l][1]) < massTol or abs(mz - oriLIonsForMod[l][2]) < massTol:
                numModL += 1
                if modPosL + l + 1 > L:
                    L = modPosL + l + 1
                    Lb.append(L)
                    Lbsta += 'c'
                    
        for r in range(len(oriRIonsForMod)):
            if abs(mz - oriRIonsForMod[r][0] - pcMassDiff) < massTol or abs(mz - oriRIonsForMod[r][1] - pcMassDiff) < massTol or abs(mz - oriRIonsForMod[r][2] - pcMassDiff) < massTol:
                numModR += 1
                if modPosR - r > L:
                    L = modPosR - r
                    Lb.append(L)
                    Lbsta += 'm'
            if abs(mz - oriRIonsForMod[r][0]) < massTol or abs(mz - oriRIonsForMod[r][1]) < massTol or abs(mz - oriRIonsForMod[r][2]) < massTol:
                numModR += 1
                if modPosR - r - 1 < R:
                    R = modPosR - r - 1
                    Rb.append(R)
                    Rbsta += 'c'
                    
        #cov          
        for l in range(len(oriLIonsForCov)):
            # print(abs(eIon - oriLIons[l] - pcMassDiff) , massTol)
            if abs(mz - oriLIonsForCov[l][0] - (modInL == 1)*pcMassDiff) < massTol or abs(mz - oriLIonsForCov[l][1] - (modInL == 1)*pcMassDiff) < massTol or abs(mz - oriLIonsForCov[l][2] - (modInL == 1)*pcMassDiff) < massTol:
                numModL += 1
                    
        for r in range(len(oriRIonsForCov)):
            if abs(mz - oriRIonsForCov[r][0] - (modInL == 0)*pcMassDiff) < massTol or abs(mz - oriRIonsForCov[r][1] - (modInL == 0)*pcMassDiff) < massTol or abs(mz - oriRIonsForCov[r][2] - (modInL == 0)*pcMassDiff) < massTol:
                numModR += 1
                    
    # print('minR, maxL',minR, maxL)
    # print(modPosL, modPosR)
    # print((modPosL + (lenMod - min(maxL, lenMod))), (modPosR - (lenMod - minR)))
    return numModL + numModR, Lb, Rb, Lbsta, Rbsta
          
def findPtm(mzs, reliableTags, pepCandDict, pcMass, tol):
    massTol = pcMass*tol
    # additionalTags = []
    # for reliableTag in reliableTags:
    #     if reliableTag[0].count('Q') == 1:
    #         additionalTags.append(tuple((reliableTag[0].replace('Q', 'GA'), reliableTag[1], reliableTag[2])))
    #         additionalTags.append(tuple((reliableTag[0].replace('Q', 'AG'), reliableTag[1], reliableTag[2])))
    # reliableTags.extend(additionalTags)
    
    psms = []
    # print(reliableTags)
    for reliableTag in reliableTags:
        tag = reliableTag[0]
        lOfTag = len(tag)
        # if tag != 'DAISAETGT':
        #     continue
        # print(tag)
        # print('it RHF now')
        index = [int(i) for i in reliableTag[2].split()]
        
        pepCand = pepCandDict[tag]
        # print(pepCand)
        if len(pepCand) == 0:
            continue
        
        # print(pepCand)
        for pep in pepCand:
            lOfPep = len(pep)
            # if pep!= 'YHTVNGHNCEVRK':
            #     continue
            # print(pep)
            # print('GSCFHR')
            # pcMass = 3076.39902
            pcMassDiff = pcMass - calcuSeqMass(pep) 
            
            
            bPeaks, yPeaks, bH2OPeaks, yH2OPeaks, bNH3Peaks, yNH3Peaks, theoPeaks = getTheoPeaks(pep)
            
            isPepReversed = False
            
            if pep.find(tag) != -1:
                alignPos = pep.find(tag)
                lPeaks = bPeaks.copy()
                rPeaks = yPeaks.copy()
                lH2OPeaks = bH2OPeaks.copy()
                rH2OPeaks = yH2OPeaks.copy()
                lNH3Peaks = bNH3Peaks.copy()
                rNH3Peaks = yNH3Peaks.copy()
                
            else:
                isPepReversed = True
                pep = pep[::-1]
                alignPos = pep.find(tag)
                
                lPeaks = yPeaks.copy()
                rPeaks = bPeaks.copy()
                lH2OPeaks = yH2OPeaks.copy()
                rH2OPeaks = bH2OPeaks.copy()
                lNH3Peaks = yNH3Peaks.copy()
                rNH3Peaks = bNH3Peaks.copy()
                
            lIonMassTheo = np.array([lPeaks[i] for i in range(alignPos, alignPos + lOfTag + 1)]) #range(i, j) has i, does not have j
            lIonMassExpe = np.array([mzs[i] for i in index])
            # print(lIonMassTheo)
            if abs(pcMassDiff) < massTol:
                otherMzs = list(set(mzs) - set(lIonMassExpe))
                originIons = lPeaks + lH2OPeaks + lNH3Peaks + rPeaks + rH2OPeaks + rNH3Peaks
                numMod = getModIons(originIons, otherMzs, massTol)
                psms.append( tuple((tag, pep, pcMassDiff, '-', '-', '-', '-', numMod)) ) 
                continue
            # print(lIonMassExpe)
            alignMassDiff = np.average(lIonMassExpe - lIonMassTheo)
            # print(pcMassDiff)
            # print(abs(pcMassDiff - alignMassDiff) < massTol)
            modPosL = 1 #initial pos start from 1
            modPosR = lOfPep #initial pos
            modInL = 0
            
            if abs(pcMassDiff - alignMassDiff) < massTol: #very confident, almost do not need further validation
                #mod pos is in left
                
                modPosL = 1
                modPosR = alignPos
                modInL = 1  #Lions for cov should be compare to modified
                # print(modPosL, modPosR)
                #try to narrow down the mod pos range

                oriLIonsForMod = [i for i in zip(lPeaks[modPosL : modPosR], lH2OPeaks[modPosL : modPosR], lNH3Peaks[modPosL : modPosR] )]
                oriRIonsForMod = [i for i in zip(rPeaks[lOfPep - modPosR + 1 : ], rH2OPeaks[lOfPep - modPosR + 1 : ], rNH3Peaks[lOfPep - modPosR + 1 : ]  )]
                #for cov, they are used only to calcu ion cov, only oriMass should be used
                oriLIonsForCov = [i for i in zip(lPeaks[modPosR + lOfTag + 1 : ], lH2OPeaks[modPosR + lOfTag + 1 : ], lNH3Peaks[modPosR + lOfTag + 1 : ] )]
                oriRIonsForCov = [i for i in zip(rPeaks[1 : lOfPep - modPosR], rH2OPeaks[1 : lOfPep - modPosR], rNH3Peaks[1 : lOfPep - modPosR]  )]

                #
                otherMzs = list(set(mzs) - set(lIonMassExpe))
                numMod, modPosLb, modPosRb, Lbs, Rbs = getModIonsV2(modInL, oriLIonsForMod, oriRIonsForMod, oriLIonsForCov, oriRIonsForCov, otherMzs, 0.02, pcMassDiff, modPosL, modPosR)
                # if numMod <= 5:
                #     continue
                if numMod >= lOfPep or min(modPosRb) == max(modPosLb) <= 2:
                    if isPepReversed:
                        psms.append( tuple((tag, pep[::-1], pcMassDiff, lOfPep + 1 - modPosRb, lOfPep + 1 - modPosLb, Rbs, Lbs, numMod)) )
                    else:
                        psms.append( tuple((tag, pep, pcMassDiff, modPosLb, modPosRb, Lbs, Rbs, numMod)) )
                    
            elif abs(alignMassDiff) < massTol:
                #mod pos is in right
                # print('here')
                modPosL = alignPos + lOfTag + 1
                modPosR = lOfPep
                # print(modPosL, modPosR)
                
                oriLIonsForMod = [i for i in zip(lPeaks[modPosL : ], lH2OPeaks[modPosL : ], lNH3Peaks[modPosL : ] )]
                
                oriRIonsForMod = [i for i in zip(rPeaks[1 : lOfPep - modPosL + 1], rH2OPeaks[1 : lOfPep - modPosL + 1], rNH3Peaks[1 : lOfPep - modPosL + 1]  )]
                #for cov, they are used only to calcu ion cov, only oriMass should be used
                oriLIonsForCov = [i for i in zip(lPeaks[1 : modPosL - lOfTag - 1], lH2OPeaks[1 : modPosL - lOfTag - 1], lNH3Peaks[1 : modPosL - lOfTag - 1] )]
                if 1 > modPosL - lOfTag - 1:
                    oriLIonsForCov = []
                oriRIonsForCov = [i for i in zip(rPeaks[lOfPep - modPosL + 1 + 1 : ], rH2OPeaks[lOfPep - modPosL + 1 + 1 : ], rNH3Peaks[lOfPep - modPosL + 1 + 1 : ]  )]

                otherMzs = list(set(mzs) - set(lIonMassExpe))
                numMod, modPosLb, modPosRb, Lbs, Rbs = getModIonsV2(modInL, oriLIonsForMod, oriRIonsForMod, oriLIonsForCov, oriRIonsForCov, otherMzs, 0.02, pcMassDiff, modPosL, modPosR)
                # print('hasMod', hasMod)
                # print(numMod, modPosLb, modPosRb)
                # if numMod <= 5:
                #     continue
                if numMod >= lOfPep or (min(modPosRb) == max(modPosLb)) :
                    if isPepReversed:
                        psms.append( tuple((tag, pep[::-1], pcMassDiff, lOfPep + 1 - modPosRb, lOfPep + 1 - modPosLb, Rbs, Lbs, numMod)) )
                    else:
                        psms.append( tuple((tag, pep, pcMassDiff, modPosLb, modPosRb, Lbs, Rbs, numMod)) )
                else:
                    continue
            # else:
            #     print('todo more than one modification')
                
    return psms

def lcs(s1, s2): 
	m=[[0 for i in range(len(s2)+1)]  for j in range(len(s1)+1)]  #生成0矩阵，为方便后续计算，比字符串长度多了一列
	mmax=0   #最长匹配的长度
	p=0  #最长匹配对应在s1中的最后一位
	for i in range(len(s1)):
		for j in range(len(s2)):
			if s1[i]==s2[j]:
				m[i+1][j+1]=m[i][j]+1
				if m[i+1][j+1]>mmax:
					mmax=m[i+1][j+1]
					p=i+1
	return s1[p-mmax:p],mmax   #返回最长子串及其长度
 

def cleanUpTags(reliableTags):
    cleanedTags = []
    print('reliableTags', reliableTags)
    print(len(reliableTags))
    for i in range(len(reliableTags)):
        print('i', i)
        a = reliableTags[i]
        print('after')
        if len(reliableTags[i][0]) <= 6:
            continue
        for j in range(len(reliableTags) - 1, i, -1):
            # print(i,j)
            if fuzz.ratio(reliableTags[i][0], reliableTags[j][0]) >= 60:
                reliableTags.pop(j)
                
                
def simi(str1, str2):
    if str1 in str2 or str1[::-1] in str2:
        return 100 - min((len(str2) - len(str1)), 15)/2
    return (  max(fuzz.partial_ratio(str1, str2), fuzz.partial_ratio(str1[::-1], str2))
            + len(str1)/len(str2)*2
            - len(str1)/len(lcs(str1, str2))*3)


def getPepCand(tag, feasiblePeps):
    
    pepCands = []
    
    for pep in feasiblePeps:
        if tag in pep or tag[::-1] in pep:# exact containing
            pepCands.append(pep)
        
    return pepCands
    
if __name__ == '__main__':
    
    s1 = 'TKKAS'
    
    s2 = 'ISAKPPAK' #YIDIPKmIDAEDIVGTARPDEK
          # DDPVTNINNAFEVAEKYIDIPK
    tags = ['KAEMSAAIVG','KSIQNVI', 'KAEKNVI','KPTQM', 'SANM', 'SGQM', 'VAQM', 'mSAA', 'QNM', 'IVM', 'KSIGANVI', 'KSIAGNVI', 'KPTGAM', 'KPTAGM', 'SGGAM', 'SGAGM', 'VAGAM', 'VAAGM', 'GANM', 'AGNM']
    peps = ['VGVIAASMEAK', 'EKGVAASSAQK', 'IQAVASMVEK', 'CKWVNQIK', 'SEIRNISEK', 'VISGPmEKAK', 'QIGSMVEIAK', 'GIAFAEIQAR', 'EGIAFRPASK', 'TVGIPTAmAAK'] # print(fuzz.partial_ratio(s1,s2))
    # print(fuzz.ratio(s1,s2))
    
    # print(simi(s1,s2))
    from suffix_trees import STree
    st = STree.STree(tags)
    # st.build('KFWPMNASPEI')
    a = st.find_all('NVI')
    
    # print(st.word_starts)
    
    # for i in tags:
    #     print(len(i))
        
    index = [15, 23, 66, 75]
    pos = [0, 11, 19, 27, 33, 38, 43, 48, 53, 57, 61, 70, 79, 86, 93, 99, 105, 111, 117, 122] 
    from Utils.Funcs import binarySearch
    
    for i in index:
        
        print(i, binarySearch(pos, 0, len(pos), i), tags[binarySearch(pos, 0, len(pos), i)])
    # print(lcs(s1,s2))
    # for tag in tags:
    #     for pep in peps:
    #         if simi(tag, pep)>70:
    #             print(simi(tag, pep))
    #             print(tag, pep)
    
