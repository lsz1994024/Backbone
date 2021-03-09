# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 09:27:56 2020

@author: slaiad
"""

import pyteomics.fasta as fas
from pyteomics import parser
import re
from Parameters import MAX_LENGTH_OF_PEP, MIN_LENGTH_OF_PEP, NUM_CLUSTERS_1ST_LAYER, MISSED_CLEAVAGES
from collections import Counter
import random
import numpy as np
from Utils.Consts import AA_RES_MASS, ATOM_MASS, H2O, NH3
from Utils.Funcs import divideInputList

def getAllProts(fastaDir):
    fastaData = fas.FASTA(fastaDir)
    
    strOfProts = fastaData.read()
    
    listOfProts = re.split(r'>sp\|.*\n', strOfProts) 
    listOfProts = [i for i in listOfProts if i != ""]
    
    for j in range(len(listOfProts)):
        listOfProts[j] = listOfProts[j].replace("\n","")
        
    return listOfProts
  
def writePdbPymolCmd(uniprotToPdbFile):
    pdbList = open(uniprotToPdbFile, 'r')
    pdbs = pdbList.readlines()
    pdbList.close()
    
    # if os.path.exists('TempData/pyMolCmd.txt'):
    #     os.remove('TempData/pyMolCmd.txt')
    pyMolCmd = open('TempData/pdbFilesOnlyName.txt', 'w')

    for pdb in pdbs:
        pyMolCmd.writelines([pdb.split()[1].lower(), '.pdb', '\n'])
        # pyMolCmd.writelines([pdb.split()[1], ','])

    pyMolCmd.close()    
    
def getInnerTags(cifFileDir):
    with open(cifFileDir, 'r') as cifFile:
        lines = cifFile.readlines()
    
    coord = []
    for line in lines:
        if line[0:4] == 'ATOM':
            line = line.split()
            coord.append([float(x) for x in line[10:13]])
    print(len(coord))
    return coord
    
def digest(prot):
    
    peps = list(parser.cleave(prot, parser.expasy_rules["trypsin"], missed_cleavages = MISSED_CLEAVAGES, min_length = MIN_LENGTH_OF_PEP))
        
    peps = [pep for pep in peps if (len(pep) <= MAX_LENGTH_OF_PEP 
                                    and 'B' not in pep
                                    and 'J' not in pep 
                                    and 'X' not in pep 
                                    and 'Z' not in pep 
                                    and 'O' not in pep 
                                    and 'U' not in pep)]
    additionalPeps = []
    
    for pep in peps:
        if pep[0] == 'M' and prot.find(pep) == 0:
            additionalPeps.append(pep[1:])
    peps.extend(additionalPeps)
    
    return peps
  
    
def getTheoPeaks(pep):
    dtype = [('mass', float),('peak',  np.unicode_, 10)]
    theoPeaks = np.array([(calcuSeqMass(pep), 'pc')], dtype)
    bPeaks = [0] # nterm peak
    yPeaks = [H2O]  # cterm peak
    
    bH2OPeaks = [0]
    yH2OPeaks = [H2O]
    
    bNH3Peaks = [0]
    yNH3Peaks = [H2O]
    # print('len', len(bPeaks))
    for i in range(1,len(pep)):
        left = pep[0:i]
        right = pep[i:]
        
        # print('b%d' % i, left)
        # print('y%d' % (len(pep)-i), right)
        #b
        b = calcuSeqMass(left) - H2O
        theoPeaks = np.concatenate((theoPeaks, np.array([(b, 'b%d' % i)], dtype)), axis = 0)
        bPeaks.append(b)
        
        #y
        y = calcuSeqMass(right)
        theoPeaks = np.concatenate((theoPeaks, np.array([(y, 'y%d' % (len(pep)-i))], dtype)), axis = 0)
        yPeaks.append(y)
        
        #b-H2O  b-NH3
        theoPeaks = np.concatenate((theoPeaks, np.array([(b - H2O, 'b%d-H2O' % i)], dtype)), axis = 0)
        theoPeaks = np.concatenate((theoPeaks, np.array([(b - NH3, 'b%d-NH3' % i)], dtype)), axis = 0)
        bH2OPeaks.append(b - H2O)
        bNH3Peaks.append(b - NH3)
        
        #y-H2O  y-NH3
        theoPeaks = np.concatenate((theoPeaks, np.array([(y - H2O, 'y%d-H2O' % (len(pep)-i))], dtype)), axis = 0)
        theoPeaks = np.concatenate((theoPeaks, np.array([(y - NH3, 'y%d-NH3' % (len(pep)-i))], dtype)), axis = 0)
        yH2OPeaks.append(y - H2O)
        yNH3Peaks.append(y - NH3)
        
        #pc
    pepMass = calcuSeqMass(pep)
    bPeaks.append(pepMass)
    yPeaks.append(pepMass)
    bH2OPeaks.append(pepMass - H2O)
    yH2OPeaks.append(pepMass - H2O)
    bNH3Peaks.append(pepMass - NH3)
    yNH3Peaks.append(pepMass - NH3)
    # print('lenend', len(bPeaks))
    return sorted(bPeaks), sorted(yPeaks), sorted(bH2OPeaks), sorted(yH2OPeaks), sorted(bNH3Peaks), sorted(yNH3Peaks), np.sort(theoPeaks, order = ['mass'])


def calcuSeqMass(seq):
    sumMass = 0
    
    for aa in seq:
        sumMass += AA_RES_MASS[aa]
    sumMass += ATOM_MASS['H']*2 + ATOM_MASS['O']
    
    return sumMass

    
def getAllPeps(protList):
    dtype = [('seq', np.unicode_, MAX_LENGTH_OF_PEP), ('pcMass', float)]
    allPeps = np.array([], dtype)
    
    # i  = 0
    for prot in protList:
        prot = prot.replace("L", "I")
        peps = digest(prot)
        # for pep in peps:
        allPeps = np.concatenate((allPeps, peps), axis = 0)
        # i+=1
        # if i%100 ==0:
        #     print(i)
    allPeps = (np.sort(allPeps, order = ['pcMass']))[::-1]
    return np.unique(allPeps)

def getAllPepsForParallel(protList):
    allPeps = []
    for prot in protList:
        prot = prot.replace("L", "I")
        peps = digest(prot)
        allPeps.extend(peps)
    return allPeps

def getTags(pep):
    tags = []
    for i in range(len(pep) - 2):
        tags.append(pep[i:i+3])
        
    tags = list(set(tags))
    for i in range(len(tags)):
        if tags[i][0] > tags[i][-1]:
            tags[i] = tags[i][::-1]
            
    return list(set(tags))

def getTagsList(pep):
    tags = []
    for i in range(len(pep) - 2):
        tags.append(pep[i:i+3])
        
    for i in range(len(tags)):
        if tags[i][0] > tags[i][-1]:
            tags[i] = tags[i][::-1]
            
    return tags

def getTagsWithCount(pep):
    tags = []
    for i in range(len(pep) - 2):
        tags.append(pep[i:i+3])
        
    for i in range(len(tags)):
        if tags[i][0] > tags[i][-1]:
            tags[i] = tags[i][::-1]
            
    tagsWithCount = Counter(tags)
    return tagsWithCount
    
def getTagsDict(allPeps):
    tagsOfPep = {}
    for i in range(len(allPeps)):
        tagsOfPep[allPeps[i]] = getTags(allPeps[i])
    # print(allTags)
    return tagsOfPep

def getTagsWithCountDict(allPeps):
    tagsOfPep = {}
    for i in range(len(allPeps)):
        tagsOfPep[i] = getTagsWithCount(allPeps[i])
#    print(tagsOfPep)
    return tagsOfPep
    
def getAllTags(tagsDict):
    allTags = []
    for pep in tagsDict:
        allTags.extend(tagsDict[pep])
        
    allTags = list(set(allTags))            

    return allTags

def getTagsInCommon(tags1,tags2):
    num = len([tag for tag in tags1 if tag in tags2])
    return num

def getSimilar(tagsWithCount1, tagsWithCount2):
    similarity = 0
    for tag in tagsWithCount1:
        if tag in tagsWithCount2:
            similarity += tagsWithCount1[tag] + tagsWithCount2[tag]
    return similarity
     
def getRandPepIndexForCluster(tagsDictWithCount, allPeps):
    randInt = random.randint(0, len(allPeps))
    pepIndex = [randInt]
    
    while len(pepIndex) < NUM_CLUSTERS_1ST_LAYER:
        similarity = 0
        randI = random.randint(0, len(allPeps) - 1)
        for j in range(len(pepIndex)):
            if j >= len(allPeps) or randI >= len(allPeps):
                print("lsz list out of range j  randI", j, randI)
            similarity += getSimilar(tagsDictWithCount[j], tagsDictWithCount[randI])
        if similarity == 0:
            pepIndex.append(randI)
            
    return pepIndex
   
def updateCluster(clusterOfTags, clusterOfPeps, tagsDictWithCount):
    # newClusterOfTags = clusterOfTags ##only to get the keys
    
    for cluster in clusterOfPeps:
        tagsDict = {}
        for pepIndex in clusterOfPeps[cluster]:
            pepTags = tagsDictWithCount[pepIndex]
            for tag in pepTags:
                if tag in tagsDict:
                   tagsDict[tag] += pepTags[tag]  ##need to rethink, how to evaluate the similarity of a pep and a cluster of peps
                else:    
                   tagsDict[tag] = pepTags[tag]
        clusterOfTags[cluster] = tagsDict
        # newClusterOfTags[cluster] = tagsDict
        
    for cluster in clusterOfPeps:
        clusterOfPeps[cluster] = []
    # return newClusterOfTags, clusterOfPeps
     
if __name__ == '__main__':
#    protList = getProtsList("test.fasta")
    peaks = getTheoPeaks('QTVAVGVIKAVDKK')
    #'VAVAHAGHR', 'TATAVAHCK', 'IKRAVAHK', 'AEDGHAVAK', 'QIAVAHEK']
    
    print(peaks[-1])

        

