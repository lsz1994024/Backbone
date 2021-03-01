#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 09:53:32 2021

@author: slaiad@ust.hk
"""

from fuzzywuzzy import fuzz
from ProcessDbPack.ParsePep import getTheoPeaks, calcuSeqMass

def getIonCoverage(ions, mzs, massTol, lenTag):
    numMatch = 0
    for mz in mzs:
        minDiff = 100
        for ion in ions:
            if abs(mz - ion) < minDiff:
                minDiff = abs(mz - ion)
        
        if minDiff <= massTol:
            numMatch += 1
            print(mz)
    print(numMatch - 2, len(mzs) - lenTag - 2)
    # return float(numMatch - 2)/(len(mzs) - lenTag - 2)
    return (numMatch - 2)
    
    
    
    
    
def findPtm(mzs, reliableTags, feasiblePeps, pcMass, tol):
    
    additionalTags = []
    for reliableTag in reliableTags:
        if reliableTag[0].count('Q') == 1:
            additionalTags.append(tuple((reliableTag[0].replace('Q', 'GA'), reliableTag[1], reliableTag[2])))
            additionalTags.append(tuple((reliableTag[0].replace('Q', 'AG'), reliableTag[1], reliableTag[2])))
    reliableTags.extend(additionalTags)
    
    psms = []
    print(reliableTags)
    for reliableTag in reliableTags:
        tag = reliableTag[0]
        print(tag)
        
        # if tag != 'KRVE':
        #     continue
        # print('it RHF now')
        index = [int(i) for i in reliableTag[2].split()]
        
        pepCand = getPepCand(tag, feasiblePeps)
        # print(pepCand)
        if len(pepCand) == 0:
            continue
        
        # print(pepCand)
        for pep in pepCand:
            # print(pep)
            # if pep!= 'GSCFHR':
            #     continue
            # print('GSCFHR')
            pcMassDiff = pcMass - calcuSeqMass(pep)  
            
            bPeaks, yPeaks, bH2OPeaks, yH2OPeaks, bNH3Peaks, yNH3Peaks, theoPeaks = getTheoPeaks(pep)
            # alignPos
            if pep.find(tag) != -1: # left to right, use b ions
                alignPos = pep.find(tag)
                print(tag, pep)
                modPos = alignPos # index starts from 1, for those life science people
                print(alignPos, 'N')
                bIonMassTheo = [bPeaks[i] for i in range(alignPos, alignPos + len(tag) + 1)] #range(i, j) has i, does not have j
                bIonMassExpe = [mzs[i] for i in index]
                
                # print([yIonMassExpe[i] - yIonMassTheo[i] for i in range(len(yIonMassTheo))])
                aveMassDiff = sum([bIonMassExpe[i] - bIonMassTheo[i] for i in range(len(bIonMassTheo))])/len(bIonMassTheo)
                # print(yIonMassTheo)
                # print(yIonMassExpe)
                # print(aveMassDiff)
                # print(abs(pcMassDiff - aveMassDiff) < pcMass*tol)
                ions = []
                
                if abs(aveMassDiff) < pcMass*tol:
                    modPos = -1
                    y1 = []
                    y2 = []
                    y3 = []
                    for i in range(len(pep) - alignPos + 1, len(pep) + 1):
                        y1 = yPeaks
                        y2 = yH2OPeaks
                        y3 = yNH3Peaks  #this is not yi ion, just copys of these 3 ionmass list
                        
                        y1[i] += aveMassDiff
                        y2[i] += aveMassDiff
                        y3[i] += aveMassDiff
                    ions.extend(y1)
                    ions.extend(y2)
                    ions.extend(y3)
                    ions.extend(bPeaks)
                    ions.extend(bH2OPeaks)
                    ions.extend(bNH3Peaks)
                    
                elif abs(pcMassDiff - aveMassDiff) < pcMass*tol:
                    b1 = []
                    b2 = []
                    b3 = []
                    for i in range(alignPos + len(tag) + 1, len(pep) + 1):
                        b1 = bPeaks
                        b2 = bH2OPeaks
                        b3 = bNH3Peaks  #this is not bi ion, just copys of these 3 ionmass list
                        
                        b1[i] += aveMassDiff
                        b2[i] += aveMassDiff
                        b3[i] += aveMassDiff
                    ions.extend(b1)
                    ions.extend(b2)
                    ions.extend(b3)
                    ions.extend(yPeaks)
                    ions.extend(yH2OPeaks)
                    ions.extend(yNH3Peaks)
                        
                else:
                    #check peak coverage
                    continue  

                cov = getIonCoverage(ions, mzs, pcMass*tol, len(tag))
                print('cov', cov)
                if cov > 0.1:
                    psms.append( tuple((tag, pep, pcMassDiff, modPos)) )
                
                
                
            else:                   # right to left, use y ions
                print(tag, pep)
                alignPos = pep[::-1].find(tag)
                print(alignPos, 'C')
                modPos = len(pep) - alignPos + 1 # index starts from 1, for those life science people
                yIonMassTheo = [yPeaks[i] for i in range(alignPos, alignPos + len(tag) + 1)] #range(i, j) has i, does not have j
                yIonMassExpe = [mzs[i] for i in index]
                
                # print([yIonMassExpe[i] - yIonMassTheo[i] for i in range(len(yIonMassTheo))])
                aveMassDiff = sum([yIonMassExpe[i] - yIonMassTheo[i] for i in range(len(yIonMassTheo))])/len(yIonMassTheo)
                # print(yIonMassTheo)
                # print(yIonMassExpe)
                # print(aveMassDiff)
                # print(abs(pcMassDiff - aveMassDiff) < pcMass*tol)
                ions = []
                b1 = list(bPeaks)
                b2 = list(bH2OPeaks)
                b3 = list(bNH3Peaks)
                y1 = list(yPeaks)
                y2 = list(yH2OPeaks)
                y3 = list(yNH3Peaks)
                for i in range(len(pep)):
                    b1[i] += pcMassDiff
                    b2[i] += pcMassDiff
                    b3[i] += pcMassDiff
                    y1[i] += pcMassDiff
                    y2[i] += pcMassDiff
                    y3[i] += pcMassDiff
                ions.extend(b1)
                ions.extend(b2)
                ions.extend(b3)
                ions.extend(yPeaks)
                ions.extend(yH2OPeaks)
                ions.extend(yNH3Peaks)
                
                ions.extend(y1)
                ions.extend(y2)
                ions.extend(y3)
                ions.extend(bPeaks)
                ions.extend(bH2OPeaks)
                ions.extend(bNH3Peaks)
                
                #start
                # if abs(aveMassDiff) < pcMass*tol:
                #     # print('here')
                #     modPos = -1
                #     b1 = []
                #     b2 = []
                #     b3 = []
                #     for i in range(len(pep) - alignPos + 1, len(pep) + 1):
                #         # print('here111')
                #         b1 = bPeaks
                #         b2 = bH2OPeaks
                #         b3 = bNH3Peaks  #this is not bi ion, just copys of these 3 ionmass list
                        
                #         b1[i] += aveMassDiff
                #         b2[i] += aveMassDiff
                #         b3[i] += aveMassDiff
                #     ions.extend(b1)
                #     ions.extend(b2)
                #     ions.extend(b3)
                #     ions.extend(yPeaks)
                #     ions.extend(yH2OPeaks)
                #     ions.extend(yNH3Peaks)
                #     # print(ions)
                    
                # elif abs(pcMassDiff - aveMassDiff) < pcMass*tol:
                #     y1 = []
                #     y2 = []
                #     y3 = []
                #     for i in range(alignPos + len(tag) + 1, len(pep) + 1):
                #         y1 = yPeaks
                #         y2 = yH2OPeaks
                #         y3 = yNH3Peaks
                        
                #         y1[i] += pcMassDiff
                #         y2[i] += pcMassDiff
                #         y3[i] += pcMassDiff
                #     ions.extend(y1)
                #     ions.extend(y2)
                #     ions.extend(y3)
                #     ions.extend(bPeaks)
                #     ions.extend(bH2OPeaks)
                #     ions.extend(bNH3Peaks)
                        
                # else:
                #     #check peak coverage
                #     continue  
                #end
                
                # print(ions)
                cov = getIonCoverage(list(set(ions)), mzs, pcMass*tol, len(tag))
                print('cov', cov)
                print('pcMassDiff',pcMassDiff)
                print('aveMassDiff',aveMassDiff)
                if cov > 10:
                    psms.append( tuple((tag, pep, pcMassDiff, modPos)) )
                
            
    # print(psms)
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
    # cleanTags = []
    for i in range(len(reliableTags)):
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
    
