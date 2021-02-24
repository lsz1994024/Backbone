# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 20:03:02 2020

@author: slaiad
"""
import numpy as np
from Utils.Consts import NUM_PERMUTATION
from datasketch import MinHash, MinHashLSH, WeightedMinHashGenerator, MinHashLSHEnsemble, MinHashLSHForest

def getHashSig(tagsListOfPep):
    minHash = MinHash(num_perm = NUM_PERMUTATION)
    for tag in tagsListOfPep:
        minHash.update(tag.encode('utf-8'))
    
    return minHash.digest()
        
def getHashMat(tagsDict, allPeps):
    hashMat = np.zeros([NUM_PERMUTATION, len(tagsDict)], dtype=np.uint64)
    
    # pepNum = 0
    for i in range(len(allPeps)):
        hashMat[:, i] = getHashSig(tagsDict[allPeps[i]])
    # for pep in tagsDict:
    #     hashMat[:, pepNum] = getHashSig(tagsDict[pep])
    #     pepNum += 1
        
        if i % 20000 == 0:
            print("Peptide hash calculating %.2f %%" % (i / len(tagsDict) * 100))
        
    return hashMat

def calcuJaccard(hashSig1, hashSig2):
    return np.float(np.count_nonzero(hashSig1 == hashSig2)) /\
                np.float(len(hashSig1))
                
if __name__ == '__main__':

    # data1 = 'VGFGEEWEDAAWCN'
    data1 = ['TAG', 'VGF', 'GTB','EEW']
    # data1 = ['TAG', 'RBT', 'WCDS']
    
    data2 = 'TAGDSAFDVGFGTEEWEQWWRFRSDAAWCDSNBH'
    data3 = 'TAGSSAFDDBFDTWEEWTDWWRFRSCASWCDSQBH'
    
    data2 = [data2[i:i+3] for i in range(len(data2)-2)]
    data3 = [data3[i:i+3] for i in range(len(data3)-2)]
    
    # for i in range(len(data2)-2):
    #     print(data2[i:i+3])
    
    
    m1, m2 = MinHash(), MinHash()
    # for d in data1:
    d = 'AVB'
    m1.update(d.encode('utf8'))
    # for d in data2:
    d = 'BVA'
    m2.update(d.encode('utf8'))
    print("Estimated Jaccard for data1 and data2 is", m2.jaccard(m1))
    
    # sumSimi = 0
    # m2 = MinHash()
    # m2.update(data2.encode('utf8'))
    # for d in data1:
    #     print(d)
    #     m1 = MinHash()
    #     m1.update(d.encode('utf8'))
    #     sumSimi += m2.jaccard(m1)
    
    # print('sum similarity ', sumSimi)
    
    set1 = set(data1)
    set2 = set(data2)
    set3 = set(data3)
    
    # Create MinHash objects
    m1 = MinHash(num_perm=128)
    m2 = MinHash(num_perm=128)
    m3 = MinHash(num_perm=128)
    for d in set1:
        m1.update(d.encode('utf8'))
    for d in set2:
        m2.update(d.encode('utf8'))
    for d in set3:
        m3.update(d.encode('utf8'))
    
    # Create an LSH Ensemble index with a threshold
    lshensemble = MinHashLSHEnsemble(threshold=0.594, num_perm=128)
    
    # Index takes an iterable of (key, minhash, size)
    lshensemble.index([("m2", m2, len(set2)), ("m3", m3, len(set3))])
    
    # Check for membership using the key
    print("m2" in lshensemble)
    print("m3" in lshensemble)
    
    # Using m1 as the query, get an result iterator
    print("Sets with containment > 0.8:")
    for key in lshensemble.query(m1, len(set1)):
        print(key)
        
    # Create a MinHash LSH Forest with the same num_perm parameter
    forest = MinHashLSHForest(num_perm=128)
    
    # Add m2 and m3 into the index
    forest.add("m2", m2)
    forest.add("m3", m3)
    
    # IMPORTANT: must call index() otherwise the keys won't be searchable
    forest.index()
    # Using m1 as the query, retrieve top 2 keys that have the higest Jaccard
    result = forest.query(m1, 2)
    print("Top 2 candidates", result)
