#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 09:53:32 2021

@author: slaiad@ust.hk
"""

from fuzzywuzzy import fuzz
import numpy as np
from Parameters import MAX_LENGTH_OF_PEP
from collections import Counter

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


def getPepCand(reliableTags, feasiblePeps):
    
    dtype = [('seq', np.unicode_, MAX_LENGTH_OF_PEP), ('simi', float)]
    pepCands = np.array([], dtype)
    
    additionalTags = []
    
    for reliableTag in reliableTags:
        
        tag = reliableTag[0]
        reliability = reliableTag[1]
        
        if tag.count('Q') == 1:
            additionalTags.append( tuple((tag.replace('Q','GA'), reliability)) )
            additionalTags.append( tuple((tag.replace('Q','AG'), reliability)) )
    
    reliableTags.extend(additionalTags)
    for reliableTag in reliableTags:
            
        tag = reliableTag[0]
        reliability = reliableTag[1]
        # simiFuzz = np.array([], dtype)
        
        for pep in feasiblePeps:
            # if tag[0] in pep or tag[0][::-1] in pep:# exact containing
            similarity = simi(tag, pep)
            simiAndReliability = similarity*reliability
            # if pep == 'MYSYPARVPPPPPIAR':
            #     print(tag, similarity)
            # if pep == 'IAHYNKR':
            #     print(tag, similarity)
            if similarity >= 60:
                pepCands = np.concatenate((pepCands, np.array([(pep, simiAndReliability)], dtype)), axis = 0)
        #pepCands.extend(list(sortSimiFuzz['seq'][ : - min(20, len(sortSimiFuzz) + 1) : -1]))
        # print(simiFuzz)       
    sortPepCands = (np.sort(pepCands, order = 'simi'))[::-1]
    # print(sortPepCands)
    # print(list(sortPepCands['seq'][0 : min(10, len(sortPepCands))]))
        
    top10 = sortPepCands['seq'][0 : min(10, len(sortPepCands))]
    # if len(top10) == len(set(top10)):
    return list(top10)
    
    # freq = Counter(top10)
    # return [pep for pep in freq if freq[pep] > 1]
    # return list(sortPepCands['seq'][ - min(10, len(sortPepCands)) : ])


if __name__ == '__main__':
    
    # print(simi('QIWLEJSBS', 'IQWLEJSBS'))
    
    # [('KPSPAADI', 2.6413966), ('FMGGGMG', 2.485801), ('TTVEVAES', 2.4110012), ('GIIGWRES', 2.3958423), ('ISIEVAES', 2.393613), 
    #  ('VTTITVH', 2.3881288), ('TSVSTVH', 2.3856363), ('IDTITVH', 2.3842466), ('ATVITVH', 2.3812316), 
    #  ('DVVSTVH', 2.3747982), ('KPSPIDM', 2.3723977), ('DITTTVH', 2.3708809), ('TATTTVH', 2.3688111), 
    #  ('VTTVVSR', 2.3673486), ('IDTVVSR', 2.3634664), ('ATVVVSR', 2.3604514), ('DIVEVEAS', 2.3518113), ('VTTITAR', 2.3469099)]
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
    
