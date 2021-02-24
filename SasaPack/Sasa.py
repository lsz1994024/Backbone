# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 09:53:36 2021

@author: slaiad
"""

import os
import numpy as np
from ReadTagsPack.consts import ATOM_VAN_AREA, gaussian, AA_S_NAME, NOT_R_TOPO,\
                                MIN_LENGTH_OF_PEP, R_PART_NORM, FONT
import concurrent.futures as cf
from collections import deque, Counter
import re
import itertools as it
import matplotlib.pyplot as plt
import pandas as pd

def getChainsFromInfo(atomsInfo):
    chains = {}
    chainKeys = list(set(atomsInfo['aaChain']))
    dtype = [('atomName',  np.unicode_, 3), ('aaName', np.unicode_, 1), ('aaChain', np.unicode_, 1), ('aaIndex', int), ('saArea', float), ('msArea', float)]
    for key in chainKeys:
        chains[key] = np.array([atom for atom in atomsInfo if atom['aaChain'] == key], dtype)
        
    return chains
    
def getChainsWithAreaInfo(chains):
    
    dtype = [('aaName', np.unicode_, 1), ('saArea', float), ('msArea', float)]
    
    allChainsWithAreaW = []
        
    for key in chains:
        chain = chains[key]
        allaaIndex = list(set(chain['aaIndex']))
        allaaIndex = sorted(allaaIndex)
        
        chainWithAreaW = np.array([],dtype)
        for aaIndex in allaaIndex:
            aaWithAreaW = np.array([(  [atom['aaName'] for atom in chain if atom['aaIndex'] == aaIndex][0],\
                                      sum([gaussian(atom['saArea']/ATOM_VAN_AREA[atom['atomName']], 0.2) for atom in chain if atom['aaIndex'] == aaIndex]),\
                                      sum([gaussian(atom['msArea']/ATOM_VAN_AREA[atom['atomName']], 0.2) for atom in chain if atom['aaIndex'] == aaIndex])    )], dtype)
            # if aaWithAreaW['aaName'] == 'K':
            #     print('K')
            #     # print([gaussian(atom['saArea']/ATOM_VAN_AREA[atom['atomName']]) for atom in chain if atom['aaIndex'] == aaIndex])
            #     print([atom['saArea']/ATOM_VAN_AREA[atom['atomName']] for atom in chain if atom['aaIndex'] == aaIndex])

            chainWithAreaW = np.concatenate((chainWithAreaW, aaWithAreaW), axis = 0)
            
        allChainsWithAreaW.append(chainWithAreaW)
    return allChainsWithAreaW
        
def readValidAtomsFromPdb(fullfile):
    with open(fullfile, 'r') as fid:
        line = fid.readline()

        dtype = [('atomName',  np.unicode_, 3), ('aaName', np.unicode_, 1), ('aaChain', np.unicode_, 1), ('aaIndex', int), ('saArea', float), ('msArea', float)]
        atomsInfo = np.array([], dtype)          
        
        while line is not None and line != '':
            if line[0:4] != 'ATOM':
                break
            
            #for special cases
            #ATOM   1559  CG1AVAL B 255      -9.359 -37.294  -6.556  1.9    0.135  6.840 0
            lineL = line[0:16]
            lineM = line[17:22]
            lineR = line[22:]
            
            lineL = lineL.split()
            lineM = lineM.split()
            lineR = lineR.split()
            if (lineM[0] not in AA_S_NAME) or (lineL[2] in NOT_R_TOPO) or (not lineR[0].isnumeric()) or line[12] == 'H':
                #if aa name not right i.e. it is DNA name, skip and continue
                #or if atom name is in ~R part, i.e. on main part of aa, skip and continue. Only care about R part
                line = fid.readline()
                continue
            
                                 #lineL[2][0]
            newAtom = np.array([(line[13], AA_S_NAME[lineM[0]], lineM[1], lineR[0], lineR[5], lineR[6])], dtype) 
            atomsInfo = np.concatenate((atomsInfo, newAtom), axis = 0)
            line = fid.readline()
            
    return atomsInfo    
    
def extractTagsWithArea(allPeptides):
    dtype = [('tag3', np.unicode_, 3), ('saArea', float), ('msArea', float)]
    allTags = np.array([], dtype)
    for pep in allPeptides:
        
        for i in range(len(pep) - 2):
            tag   = ''.join(pep['aaName'][i:i+3])
            saArea = np.sum(pep['saArea'][i:i+3])
            msArea = np.sum(pep['msArea'][i:i+3])
            
            tagWithArea = np.array([(tag, saArea, msArea)], dtype)
            allTags = np.concatenate((allTags, tagWithArea), axis = 0)
    return allTags
        
    
def cleaveChain(chain, missed_cleavages = 0, min_length = None):
    """modified based on parser.cleave() from pyteomics
    """
    rule = r'([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))' #trpsin rule
    
    sequence = ''.join(chain['aaName'])
    
    peptidesWithArea = [] 

    ml = missed_cleavages + 2
    trange = range(ml)
    cleavage_sites = deque([0], maxlen=ml)
    if min_length is None:
        min_length = 1
    cl = 1
    
    for i in it.chain([x.end() for x in re.finditer(rule, sequence)], [None]):
        
        cleavage_sites.append(i)
        if cl < ml:
            cl += 1
        for j in trange[:cl-1]:
            # seq = sequence[cleavage_sites[j]:cleavage_sites[-1]]
            seq = chain[cleavage_sites[j]:cleavage_sites[-1]]
            if len(seq) >= min_length:
                peptidesWithArea.append(seq)          
    return peptidesWithArea

def cleaveChainsWithAreaInfo(chains):
    allPeptides = []
    for chain in chains:
        peptides = cleaveChain(chain, missed_cleavages = 0, min_length = MIN_LENGTH_OF_PEP)
        allPeptides.extend(peptides)
    return allPeptides 
    
def mergeRepeatedTags(allTags):
    monoTags = list(set(allTags['tag3']))

    dtype = [('tag3', np.unicode_, 3), ('saArea', float), ('msArea', float)]
    setAllTags = np.array([], dtype)
    
    for tag in monoTags:
        setTag = np.array([(tag, \
                            round(np.sum([oriTag['saArea'] for oriTag in allTags if oriTag['tag3'] == tag]), 3),\
                            round(np.sum([oriTag['msArea'] for oriTag in allTags if oriTag['tag3'] == tag]), 3) )], dtype)
        setAllTags = np.concatenate((setAllTags, setTag), axis = 0)
    return setAllTags


def getNormedTagWeight(allTags):
    test = []
    tagWeightDict = {}
    for tag in allTags:
        normTerm = R_PART_NORM[tag['tag3'][0]] + R_PART_NORM[tag['tag3'][1]] + R_PART_NORM[tag['tag3'][2]]
        weight = np.array([gaussian(tag['saArea']/normTerm), gaussian(tag['msArea']/normTerm)])
        tagWeightDict[tag['tag3']] = tagWeightDict.get(tag['tag3'], np.array([0, 0])) + weight
        
        test.append(tag['saArea']/normTerm)
    
    return tagWeightDict, test


def getSetAllTags(allTags):
    
    
    allTagsDict = {}
    count = 0
    for tag in allTags:
        areas = np.array([tag['saArea'], tag['msArea']])
        allTagsDict[tag['tag3']] = allTagsDict.get(tag['tag3'], np.array([0, 0])) + areas

        count += 1
        if count % 1000000 == 0:
            print('%d of 130,000,000 tags has been set.' % count)

    freqDict = Counter(allTags['tag3'])
    
    dtype = [('tag3', np.unicode_, 3), ('saArea', float), ('msArea', float), ('freq', int)]
    setAllTags = np.array([], dtype)
    for tag in allTagsDict:
        newTag = np.array([(tag, allTagsDict[tag][0], allTagsDict[tag][1], freqDict[tag])], dtype)
        setAllTags = np.concatenate((setAllTags, newTag), axis = 0)
    return setAllTags


def extractOnePdb(fullFilePath):
    
    atomsInfo = readValidAtomsFromPdb(fullFilePath)

    chains = getChainsFromInfo(atomsInfo)
    
    chainsWithAreaInfo = getChainsWithAreaInfo(chains)
    
    allPeptides = cleaveChainsWithAreaInfo(chainsWithAreaInfo)
    
    allTags = extractTagsWithArea(allPeptides)

    return allTags

    
def extractTag3FromPdbs(files):
    dtype = [('tag3', np.unicode_, 3), ('saArea', float), ('msArea', float)]
    totalTags = np.array([], dtype)
    errorFiles = []
    processedCount = 0
    failedCount = 0
    succCount = 0
    for file in files:
        # extractOnePdb(file)
        try:
            tags = extractOnePdb(file)
            processedCount += 1
            if processedCount % 300 == 0:
                print('Processed ', processedCount, 'files in ', file[1:3])
        except:
            failedCount += 1
            errorFiles.append(file)
            print('Failed %s' % file)
            if failedCount % 100 == 0:
                print('Failed ', failedCount, ' files ', file[1:3])
        else:
            totalTags = np.concatenate((totalTags, tags), axis = 0)
            succCount += 1
            if succCount % 100 == 0:
                print('Finished %d files' % succCount, ' in ', file[1:3])
            continue  
    return errorFiles, totalTags

def getTagsDict(tags):
    tagDict = {}
    for tag in tags:
        tagDict[tag['tag3']] = [tag['saArea'], tag['msArea']]
    return tagDict

def generateInputList(filesFolder):
    inputList = []
    
    for i in range(1,10):
        fullPath = filesFolder + '/F%d/'%i
        files = os.listdir(fullPath)
        inputFiles = [('/F%d/'%i + x) for x in files if x[-3:] == 'txt' and len(x) == 8 and os.path.exists(fullPath + x[0:4] + '_residue.txt')]
        inputList.append(inputFiles)
        
    return inputList
        
def divideInputInto36(inputList):
    inputList36 = []
    
    for aList in inputList:
        quater = int(np.floor(len(aList)/4))
        inputList36.append(aList[0:quater])
        inputList36.append(aList[quater:2*quater])
        inputList36.append(aList[2*quater:3*quater])
        inputList36.append(aList[3*quater:])
    
    return inputList36

#%%main
if __name__ == '__main__':
    
    #%%generateInputList
    filesFolder = '../DataBase/PDBsHomeSapiens'
    
    inputList = generateInputList(filesFolder)
    
    inputList36 = divideInputInto36(inputList)
    #%%single test
    #fullFilePath = r'D:/oneDrive/OneDrive - HKUST Connect/Database/BSA/%s' % file
    files = ['/F3/2vh6.txt']
    errorFiles, totalTags = extractTag3FromPdbs(files)
    
    #%% parallel compute
    dtype = [('tag3', np.unicode_, 3), ('saArea', float), ('msArea', float)]
    totalTagsList = np.array([], dtype)
    totalErrorList = []
    
    with cf.ProcessPoolExecutor(max_workers = 40) as executor:
        results = executor.map(extractTag3FromPdbs, inputList36)

    for result in results:
        totalErrorList.extend(result[0])
        totalTagsList = np.concatenate((totalTagsList, result[1]), axis = 0) 
    print('lsz finished computation')
    #%%revise for error files
    # newErrorList, newSuccList = extractTag3FromPdbs(totalErrorList)
    # newSucc = extractOnePdb(totalErrorList[0])
    #%% post processing
    # setAllTags = mergeRepeatedTags(totalTagsList)
    setAllTags = getSetAllTags(totalTagsList)
    # tagsDict = getTagsDict(setAllTags)
    
    saAveAreaRankedTags = (np.sort(setAllTags, order = ['saArea']))[::-1]
    msAveAreaRankedTags = (np.sort(setAllTags, order = ['msArea']))[::-1]
    
    # dataSaRank = pd.DataFrame({'tag3' : saAreaRankedTags['tag3'], 'SA area /${\AA}^2$' : saAreaRankedTags['saArea'], 'MS area /${\AA}^2$' : saAreaRankedTags['msArea']})
    # toFile = '../TempData/SASA/saAreaRankedTags.csv' 
    # dataSaRank.to_csv(toFile, index = True, header = True, sep = ',')
    
    # dataMsRank = pd.DataFrame({'tag3' : msAreaRankedTags['tag3'], 'SA area' : msAreaRankedTags['saArea'], 'MS area ' : msAreaRankedTags['msArea']})
    # toFile = '../TempData/SASA/msAreaRankedTags.csv'
    # dataMsRank.to_csv(toFile, index = True, header = True, sep = ',')
    
    # dataAllRanks = pd.DataFrame({'saRank' : saAreaRankedTags['tag3'], 'msRank' : msAreaRankedTags['tag3'], 'fqRank': tagFrqRankedTags['tag3']})
    # toFile = '../TempData/SASA/allRankedTags.csv'
    # dataAllRanks.to_csv(toFile, index = True, header = True, sep = ',')
    
    dataAve = pd.DataFrame({'tag3' : setAllTags['tag3'], 'saArea' : setAllTags['saArea'], 'msArea' : setAllTags['msArea'], 'freq' : setAllTags['freq'], 'aveSa' : setAllTags['saAveArea'], 'aveMs' : setAllTags['msAveArea']})
    toFile = '../TempData/SASA/dataAll.csv'
    dataAve.to_csv(toFile, index = True, header = True, sep = ',')
    #%%drawing
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111)
    
    ax.scatter(np.sqrt(setAllTags['saArea']), np.sqrt(setAllTags['msArea']), color = "black")
    ax.grid(True)
    ax.axis('equal')
    ax.plot([4, 2000], [4, 2000], color='green', linestyle='dashed')
    plt.xlabel('SAS area / ${\AA}^2$', FONT)
    plt.ylabel('MS area / ${\AA}^2$', FONT)
    plt.xlim([4, 2000])
    plt.ylim([4, 2000])
    
    plt.show()
    fig.savefig('../TempData/SASA/distributionTagsWithArea.eps', dpi = 200, format = 'eps')
    
