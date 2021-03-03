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
    
def getModIons(oriLIons, oriRIons, mzs, massTol, pcMassDiff, modPosL, modPosR):
    
    # print(massTol)
    # print(mzs)
    numModL = 0
    numModR = 0
    lenMod = modPosR - modPosL
    maxL = lenMod
    minR = lenMod
    print('minR, maxL',minR, maxL)
    # print(pcMassDiff)
    for mz in mzs:
        for l in range(len(oriLIons)):
            # print(abs(eIon - oriLIons[l] - pcMassDiff) , massTol)
            if abs(mz - oriLIons[l][0] - pcMassDiff) < massTol or abs(mz - oriLIons[l][1] - pcMassDiff) < massTol or abs(mz - oriLIons[l][2] - pcMassDiff) < massTol:
                # print('mz', mz, oriLIons[l][0])
                if numModL == 0:
                    minR = l  #the first lIon corrects the max right mod pos
                    print('l',l)
                numModL += 1
                
        for r in range(len(oriRIons)):
            if abs(mz - oriRIons[r][0] - pcMassDiff) < massTol or abs(mz - oriRIons[r][1] - pcMassDiff) < massTol or abs(mz - oriRIons[r][2] - pcMassDiff) < massTol:
                # print('mz', mz, oriRIons[r][0])
                if numModR == 0:
                    maxL = r  #the first lIon corrects the max right mod pos
                    print('r',r)
                numModR += 1
    print('minR, maxL',minR, maxL)
    print(modPosL, modPosR)
    print((modPosL + (lenMod - min(maxL, lenMod))), (modPosR - (lenMod - minR)))
    return numModL + numModR, (modPosL + (lenMod - min(maxL, lenMod))), (modPosR - (lenMod - minR))

def getModIonsV2(oriLIons, oriRIons, mzs, massTol, pcMassDiff, modPosL, modPosR):
    numModL = 0
    numModR = 0
    lenMod = modPosR - modPosL
    maxL = lenMod
    minR = lenMod
    print('minR, maxL',minR, maxL)
    # print(pcMassDiff)
    for mz in mzs:
        for l in range(len(oriLIons)):
            # print(abs(eIon - oriLIons[l] - pcMassDiff) , massTol)
            if abs(mz - oriLIons[l][0] - pcMassDiff) < massTol or abs(mz - oriLIons[l][1] - pcMassDiff) < massTol or abs(mz - oriLIons[l][2] - pcMassDiff) < massTol:
                # print('mz', mz, oriLIons[l][0])
                if numModL == 0:
                    minR = l  #the first lIon corrects the max right mod pos
                    print('l',l)
                numModL += 1
                
        for r in range(len(oriRIons)):
            if abs(mz - oriRIons[r][0] - pcMassDiff) < massTol or abs(mz - oriRIons[r][1] - pcMassDiff) < massTol or abs(mz - oriRIons[r][2] - pcMassDiff) < massTol:
                # print('mz', mz, oriRIons[r][0])
                if numModR == 0:
                    maxL = r  #the first lIon corrects the max right mod pos
                    print('r',r)
                numModR += 1
    print('minR, maxL',minR, maxL)
    print(modPosL, modPosR)
    print((modPosL + (lenMod - min(maxL, lenMod))), (modPosR - (lenMod - minR)))
    return numModL + numModR, (modPosL + (lenMod - min(maxL, lenMod))), (modPosR - (lenMod - minR))
          
def findPtm(mzs, reliableTags, feasiblePeps, pcMass, tol):
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
        
        if tag != 'GVA':
            continue
        # print(tag)
        # print('it RHF now')
        index = [int(i) for i in reliableTag[2].split()]
        
        pepCand = getPepCand(tag, feasiblePeps)
        # print(pepCand)
        if len(pepCand) == 0:
            continue
        
        # print(pepCand)
        for pep in pepCand:
            if pep!= 'SKSYYICTSISTPAIGAGGSGSTGGAVGGK':
                continue
            # print(pep)
            # print('GSCFHR')
            # pcMass = 2755.258270771486
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
                
            lIonMassTheo = np.array([lPeaks[i] for i in range(alignPos, alignPos + len(tag) + 1)]) #range(i, j) has i, does not have j
            lIonMassExpe = np.array([mzs[i] for i in index])
            # print(lIonMassTheo)
            # print(lIonMassExpe)
            alignMassDiff = np.average(lIonMassExpe - lIonMassTheo)
            print(pcMassDiff)
            print(abs(pcMassDiff - alignMassDiff) < massTol)
            modPosL = 1 #initial pos start from 1
            modPosR = len(pep) #initial pos
            if abs(pcMassDiff - alignMassDiff) < massTol: #very confident, almost do not need further validation
                #mod pos is in left
                
                modPosL = 1
                modPosR = alignPos
                print(modPosL, modPosR)
                #try to narrow down the mod pos range
                oriLIons = [i for i in zip(lPeaks[modPosL + len(tag) + 1 : ], lH2OPeaks[modPosL + len(tag) + 1 : ], lNH3Peaks[modPosL + len(tag) + 1 : ] )]
                oriRIons = [i for i in zip(rPeaks[1 : len(pep) - modPosR], rH2OPeaks[1 : len(pep) - modPosR], rNH3Peaks[1 : len(pep) - modPosR]  )]
                otherMzs = list(set(mzs) - set(lIonMassExpe))
                getModIonsV2(oriLIons, oriRIons, otherMzs, 0.02, pcMassDiff, modPosL, modPosR)
                if isPepReversed:
                    psms.append( tuple((tag, pep[::-1], pcMassDiff, len(pep) + 1 - modPosR, len(pep) + 1 - modPosL)) )
                else:
                    psms.append( tuple((tag, pep, pcMassDiff, modPosL, modPosR)) )
                    
            elif abs(alignMassDiff) < massTol:
                #mod pos is in right
                # print('here')
                modPosL = alignPos + len(tag) + 1
                modPosR = len(pep)
                print(modPosL, modPosR)
                oriLIons = [i for i in zip(lPeaks[modPosL : ], lH2OPeaks[modPosL : ], lNH3Peaks[modPosL : ] )]
                oriRIons = [i for i in zip(rPeaks[1 : ], rH2OPeaks[1 : ], rNH3Peaks[1 : ]  )]
                otherMzs = list(set(mzs) - set(lIonMassExpe))
                numMod, modPosL, modPosR = getModIons(oriLIons, oriRIons, otherMzs, 0.02, pcMassDiff, modPosL, modPosR)
                # print('hasMod', hasMod)
                print(modPosL, modPosR)
                if numMod >= NUM_EXTRA_MOD_PEAKS:
                    if isPepReversed:
                        psms.append( tuple((tag, pep[::-1], pcMassDiff, len(pep) + 1 - modPosR, len(pep) + 1 - modPosL, numMod)) )
                    else:
                        psms.append( tuple((tag, pep, pcMassDiff, modPosL, modPosR, numMod)) )
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
    
    print(simi(s1,s2))
    
    # print(lcs(s1,s2))
    # for tag in tags:
    #     for pep in peps:
    #         if simi(tag, pep)>70:
    #             print(simi(tag, pep))
    #             print(tag, pep)
    
