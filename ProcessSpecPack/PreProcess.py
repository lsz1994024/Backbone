# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 21:36:29 2021

@author: slaiad
"""

import ms_deisotope
from Utils.Consts import ATOM_MASS
import numpy as np

def preProcMzs(mzs):
    
    
    # #for origin mz
    # newMzs = np.insert(mzs, 0, ATOM_MASS['O']*1 + ATOM_MASS['H']*2 + ATOM_MASS['Z']) #worked for y1
    # newMzs = np.insert(newMzs, 0, ATOM_MASS['Z']) #worked for b1
    
    #for neutral mass
    newMzs = np.insert(mzs, 0, ATOM_MASS['O']*1 + ATOM_MASS['H']*2 ) #worked for y1
    newMzs = np.insert(newMzs, 0, 0) #worked for b1
    return newMzs

def deIsotope(peaks):
    dcPeaks, _ = ms_deisotope.deconvolute_peaks(peaks,
                                            averagine=ms_deisotope.peptide,
                                            truncate_after = 0.999, ##truncate_after is very important, loose a little bit
                                            scorer=ms_deisotope.MSDeconVFitter(10.),
                                            incremental_truncation = 0.8)
    
    mzs = [peak.neutral_mass for peak in dcPeaks]
    
    mzs = preProcMzs(mzs)
    
    return mzs
#%%main
if __name__ == '__main__':
    
    #%%read
    # peptide_averagine = ms_deisotope.Averagine({"C": 4.9384, "H": 7.7583, "N": 1.3577, "O": 1.4773, "S": 0.0417})

    

    
    # humanmlPath = 'G:/Dataset/PXD022999/191122_MK_SIO13_P2-GM1.mzML'
    humanXmlPath = 'G:/Dataset/PXD022999/191122_MK_SIO13_P2-GM1.mzXML'
    # path = 'G:/Dataset/PXD022999/test/three_test_scans.mzML'
    # # reader = ms_deisotope.MSFileLoader(path)
    # specDict = readPeaksFromXML(humanXmlPath)
    
    #%%do
    # peaks = specDict[44001]
    # bunch = next(reader)
    # dcPeaks, _ = ms_deisotope.deconvolute_peaks(peaks,
    #                                             averagine=ms_deisotope.peptide,
    #                                             truncate_after = 0.999, ##truncate_after is very important, loose a little bit
    #                                             scorer=ms_deisotope.MSDeconVFitter(10.),
    #                                             incremental_truncation = 0.8)
    # dcPeaks = []  
                                                      
    # for peak in deconvoluted_peaks.peaks:
        
                        

    # dataPeaks = pd.DataFrame({'neutral mass' : [peak.neutral_mass for peak in dcPeaks],
    #                           'charge'       : [peak.charge for peak in dcPeaks],
    #                           'ori mass'     : [peak.mz for peak in dcPeaks],
    #                           'score'        : [peak.score for peak in dcPeaks],
    #                           'chosen'       : [peak.chosen_for_msms for peak in dcPeaks],
    #                           'len envolope' : [len(peak.envelope) for peak in dcPeaks]})
        
    # toFile = 'G:/Code/TagTree/DeIsotopePack/testData/dcPeaks.csv'
    # dataPeaks.to_csv(toFile, index = True, header = True, sep = ',')
    # bunch.precursor.deconvolute(averagine=ms_deisotope.peptide, scorer=ms_deisotope.MSDeconVFitter(10.))
    # product = bunch.products[0]
    # tags = []
    # spec = [peak.neutral_mass for peak in dcPeaks]
    
    # # spec = [128.0949727 ,141.522112  ,146.1049062 ,147.5873739 ,172.1213704 ,184.1577627 ,200.1161977 ,212.1521322 ,217.1413442 ,257.1362478 ,280.0413991 ,281.174944  ,299.1830922 ,303.1962453 ,366.2650319 ,388.205309  ,503.2326223 ,602.2973501 ,630.818041  ,715.3831655 ,756.4145986 ,757.0632192 ,844.4215566 ,973.4652578 ,998.6246206 ,1044.505907 ,1101.525927 ,1214.618822 ,1301.652147 ,]
    # spec = preProcMzs(spec)
    # spec = [peak[0] for peak in peaks]
    # spec = np.insert(spec, 61, 1415.7122060123704) 
    # spec = np.insert(spec, 57, 1172.613695270183) 
    
    
    # tags = readTagsFromMS2(spec)
    # # totalTags.extend(tags)
        
    # for tag in tags:
    #     # tag[1] = round(tag[1]/len(tag[0]), 7)
    #     tag[1] = round(tag[1]/sqrt(len(tag[0])), 7)
    # tags = sorted(tags, key = (lambda x:x[1]), reverse = True)
    
