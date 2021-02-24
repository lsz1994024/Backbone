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
from Utils.Consts import AA_RES_MASS, ATOM_MASS
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
    # print(prot)
    dtype = [('seq', np.unicode_, MAX_LENGTH_OF_PEP), ('pcMass', float)]
    pepsWithMass = np.array([], dtype)
    
    peps = list(parser.cleave(prot, parser.expasy_rules["trypsin"], missed_cleavages = MISSED_CLEAVAGES, min_length = MIN_LENGTH_OF_PEP))
    # if 'MESYHKPDQQK' in peps:
    #     print(peps)
    # if len(peps) != 0:
    #     firstPep = peps[0]
    #     if firstPep ==  'MESYHKPDQQK':
    #         print(firstPep)
    #     nTermPep = firstPep[0].lower() + firstPep[1:]
    #     peps.append(nTermPep)
        
    nTermPeps = []
    oxMPeps = []
    for pep in peps:
        # print(pep)
        # a = prot.index(pep)
        # if pep not in prot:
        #     print(pep)
        if prot.index(pep) == 0:
            nTermPep = pep[0].lower() + pep[1:]
            nTermPeps.append(nTermPep)
        elif 'M' in pep:
            indexs = [M.start() for M in re.finditer('M', pep)]
            
            for i in indexs:
                oxMPep = pep[0:i]+'m'+pep[i+1:]
                oxMPeps.append(oxMPep)
                
    peps.extend(nTermPeps)
    peps.extend(oxMPeps)
    peps = [pep for pep in peps if (len(pep) <= MAX_LENGTH_OF_PEP 
                                    and 'B' not in pep
                                    and 'J' not in pep 
                                    and 'X' not in pep 
                                    and 'Z' not in pep 
                                    and 'O' not in pep 
                                    and 'U' not in pep)]
    
    pepsWithMass = np.array([tuple([pep, calcuSeqMass(pep)]) for pep in peps], dtype)
    return pepsWithMass
  
    
def calcuSeqMass(seq):
    sumMass = 0
    
    for aa in seq:
        if aa.islower():
            if aa == seq[0]:
                sumMass += 42.011 #nterm acetly
            elif aa == 'm':
                sumMass += 15.99491463  #variable M oxi
            sumMass += AA_RES_MASS[aa.upper()]
        else:
            sumMass += AA_RES_MASS[aa]
    sumMass += ATOM_MASS['H']*2 + ATOM_MASS['O']
    
    return sumMass

def testGetPepsFromProt(prot):
    
    prot = prot.replace("L", "I")
    peps = digest(prot)
    return np.unique(peps)
    
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
    dtype = [('seq', np.unicode_, MAX_LENGTH_OF_PEP), ('pcMass', float)]
    allPeps = np.array([], dtype)
    
    for prot in protList:
        prot = prot.replace("L", "I")
        peps = digest(prot)
        allPeps = np.concatenate((allPeps, peps), axis = 0)
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

def plotCommonPepsFor(index):
    maxScore = 0
    maxIndex = -1
    scoreList = []
    for i in range(0, len(allPeps)):
        if i == index:
            continue
        score = getTagsInCommon(tagsDict[allPeps[index]], tagsDict[allPeps[i]])
        
        if score > 2:
            scoreList.append(score)
            
        if score > maxScore:
            maxScore = score
            maxIndex = i
    # return Counter(scoreList)
    print(Counter(scoreList))
    print(allPeps[index])
    print(allPeps[maxIndex])
    print(maxScore)
    # plt.hist(scoreList, edgecolor = 'k', alpha = 0.35, align ='left', log = True) # 设置直方边线颜色为黑色，不透明度为 0.35
    # plt.show()
    
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
    protList = getAllProts('G:/Database/fasta/uniprot_homo_sapiens.fasta')
    
    print("length of Prots = ", len(protList))
    
    allPeps = getAllPeps(protList)
    print("Number of all peps = ", len(allPeps))
    
    aas = list('ACDEFGHIKMNPQRSTVWY')
    pepWithoutAA = {}
    
    for aa in aas:
        pepWithoutAA[aa] = []
    for pep in allPeps:
        for aa in aas:
            if aa in pep:
                pepWithoutAA[aa].append(pep)
    
    tagsDict = getTagsDict(allPeps)
    
    tagsDictWithCount = getTagsWithCountDict(allPeps)
   
    allTags = getAllTags(tagsDict)
    print("Total tags number ", len(allTags))


        

