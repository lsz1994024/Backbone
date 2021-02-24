# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 15:57:51 2021

@author: slaiad
"""
from pyteomics import mzxml
import numpy as np
from ReadTagsPack.ReadMS import readTagsFromMS2
import pandas as pd
from ReadHKPack.ReadHK import readFromHK
from math import sqrt
from SasaPack.Sasa import extractTag3FromPdbs, getSetAllTags
from ReadTagsPack.consts import ATOM_MASS

def readXML(bsaMzXmlPath):
    
    dtype = [('scanNo', int), ('pcCharge', int), ('pcMass', float)]
    
    with mzxml.read(bsaMzXmlPath) as spectra:
        specInfo = np.array([], dtype)
        specDict = {}
        for spectrum in spectra:
            # print(spectrum)
            # break
            # precurMass.append(spectrum['precursorMz'][0]['precursorMz'])
            # precurChag.append(spectrum['precursorMz'][0]['precursorCharge'])
            # specDict[spectrum['num']] = np.insert(spectrum['m/z array'], 0, 42.011)
            spec = np.array([(int(spectrum['num']),\
                              int(spectrum['precursorMz'][0]['precursorCharge']),
                              spectrum['precursorMz'][0]['precursorCharge'])], dtype)
                
            specDict[int(spectrum['num'])] = spectrum['m/z array']
                                              
            # if int(spectrum['num']) == 36781 :
            #     print(spectrum)
    return specDict

def getValidScanNo(percolatorResPath):
    table = pd.read_excel(percolatorResPath)
    scanNums = []
    for rawScanNo in table['PSMId']:
        scanNums.append(rawScanNo.split('_')[0])
    return scanNums

def getExpTags(scanNums, specDict):
    totalTags = []
    count = 0
    for scanNo in scanNums:
        count += 1
        if count % 200 == 0:
            print('Finished %d spectra' % count)
        spec = specDict[scanNo]
        
        spec = preProcMzs(spec)
        # modSpec = [i + ]
        print(spec)
        
        tags = readTagsFromMS2(spec)
        totalTags.extend(tags)
        
    for tag in totalTags:
        # tag[1] = round(tag[1]/len(tag[0]), 7)
        tag[1] = round(tag[1]/sqrt(len(tag[0])), 7)
    totalTags = sorted(totalTags, key = (lambda x:x[1]), reverse = True)
    
    return totalTags


    
#%% main
if __name__ == '__main__':
    bsaMzXmlPath = 'G:/Dataset/PXD013040/072617_GPOSU_03F_DDA1.mzXML'
    humanXmlPath = 'G:/Dataset/PXD022999/191122_MK_SIO13_P2-GM1.mzXML'
    # newBsaPath = 'G:/Dataset/PXD018758/acetyl/20190709_AcBSA_45nce.mzXML'
    specDict = readXML(humanXmlPath)
    
    bsaHKPath = 'G:/Dataset/PXD013040/072617_GPOSU_03F_DDA1.txt'
    humanHKPath = 'G:/Dataset/PXD022999/191122_MK_SIO13_P2-GM1.txt'
    spectra = readFromHK(humanHKPath)
    
    percolatorResBSAPath = 'G:/Dataset//PXD013040/072617_GPOSU_03F_DDA1.xlsx'
    percolatorResHumanPath = 'G:/Dataset//PXD022999/MK_SIO13_P2_GM1.xlsx'
    scanNums = getValidScanNo(percolatorResBSAPath)
    #%%get tag from spectra
    # scanNums = [44323,44315,44886,44001,33067,44894,47724,27980,27632,28255,39135,29164,27791,14918,44882,47794,19009,23619,27742,30689,27034,31995,44324,39076,28766,43981,18857,34831,40641,27741,34663,45327,27456,23523,22330,42736,8813,47100,33839,42877,27923,2091,45467,30887,8679,1997,31330,28760,28926,40628,47799,41329,41167,28559,27588,28563,30791,2022]#strange 15123
    sortedTotalTags = []
    scanNums = [44001]
    sortedTotalTags = getExpTags(scanNums, spectra)
    
    sortedDTotalTags = []
    sortedDTotalTags = getExpTags(scanNums, specDict)
    
    #%%get tags from pdb
    files = ['6oa6.txt', '2r0o.txt', '1ydi.txt']
    preDir = 'G:/Dataset/PXD022999/pdb/'
    fullFiles =  [preDir + file for file in files]
    # tags = extractOnePdb(fullFiles[2])
    errorFiles, totalTagsFromPdb = extractTag3FromPdbs(fullFiles)
    setAllTags = getSetAllTags(totalTagsFromPdb)
    saAveAreaRankedTags = (np.sort(setAllTags, order = ['saArea']))[::-1]
    msAveAreaRankedTags = (np.sort(setAllTags, order = ['msArea']))[::-1]
    
    dataAve = pd.DataFrame({'tag3' : setAllTags['tag3'], 'saArea' : setAllTags['saArea'], 'msArea' : setAllTags['msArea'], 'freq' : setAllTags['freq']})
    toFile = 'TempData/SASA/verifyO43707.csv'
    dataAve.to_csv(toFile, index = True, header = True, sep = ',')
    
    
    # for i in range(len(b)):
        