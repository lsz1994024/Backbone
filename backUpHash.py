# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 11:15:31 2021

@author: slaiad
"""

#%% Kmeans with similarity
allTags = getAllTags(tagsDict)
print("Total tags number ", len(allTags))

initRandPepIndex = getRandPepIndexForCluster(tagsDictWithCount, allPeps)
    
print("Finish generating initial randPep index")

clusterIndexOfAllPeps = []
clusterOfPeps = {}
clusterOfTags = {}

for i in range(len(initRandPepIndex)):
    clusterOfPeps[i] = [initRandPepIndex[i]]
    clusterOfTags[i] = tagsDictWithCount[initRandPepIndex[i]]
    
    
for iteration in range(1):
    print("iteration for " , iteration)
    for i in range(len(allPeps)):
        
        maxSimi = 0
        clusterPepIndex = -1
        clusterIndex = -1
        for j in range(NUM_CLUSTERS_1ST_LAYER):
            simi = getSimilar(tagsDictWithCount[i], clusterOfTags[j])
            
            if simi > maxSimi:
                maxSimi = simi
                clusterPepIndex = initRandPepIndex[j]##no use
                clusterIndex = j
        if clusterIndex != -1:
            clusterOfPeps[clusterIndex].append(i)

        clusterIndexOfAllPeps.append(clusterPepIndex)
    
    # updateCluster(clusterOfTags, clusterOfPeps, tagsDictWithCount)

# print("cluster125", clusterOfPeps[452632])

sumNot0 = 0
for i in clusterOfPeps:
    sumNot0 += len(clusterOfPeps[i])
print("sumNot0", sumNot0) ##513109

plotCommonPepsFor(1)

#%% MinHash
tagsDict = loadmat('TempData/tagsDict.mat')
del tagsDict['__header__']
del tagsDict['__version__']
del tagsDict['__globals__']

hashMat = getHashMat(tagsDict, allPeps)
np.save('TempData/hashMat.npy', hashMat)

#%% Load HashMat
hashMat = np.load('TempData/hashMat.npy')

#%%  simple searching test
#based on msh
# gg = np.load('TempData/tagsDict.npy').item()
maxSimi = 0
index = -1
testPep = allPeps[251145]
#反例
#251145，

tagsIndexList = sample(list(range(len(tagsDict[testPep]))), 3)
# randTagIndice = np.random.randint(0, len(tagsDict[testPep]), size = 4)
# testTags = np.array([tagsDict[testPep][i] for i in tagsIndexList])
testTags = np.array(tagsDict[testPep][tagsIndexList])


testHash = getHashSig(testTags)
for i in range(525597):
    if i % 299000 == 0:
        print("searching for testHash  %.2f %%" % (float(i / 525597) * 100))
    jac = calcuJaccard(testHash, hashMat[:, i])
    if jac > maxSimi:
        maxSimi = jac
        index = i
print(maxSimi)
print("testPep is ", testPep)
print("Observed tags is ", testTags)
print(tagsDict[testPep])

print(index)
print("Found pep is ", allPeps[index])
print(tagsDict[allPeps[index]])
#%% minhash lsh clustering

def main(hashMat):
    
    pepNum = len(hashMat[0, :])               
    
    totalRow = NUM_PERMUTATION

    bandNum = 32
    bandWidth = totalRow / bandNum
                        
    totalBucketSets = []
    
    for bandIndex in range(bandNum):
        bucket = []
        start_time = time.time()
        tupleList = []
        for pepIndex in range(pepNum):
            if pepIndex % 199999 == 0:
                print("searching for testHash  %.2f %%" % (float(pepIndex / pepNum) * 100))
            
            try:
                tupleIndex = tupleList.index(tuple(hashMat[int(bandIndex*bandWidth) : int((bandIndex + 1)*bandWidth), pepIndex]))
            except ValueError:
                tupleList.append(tuple(hashMat[int(bandIndex*bandWidth) : int((bandIndex + 1)*bandWidth), pepIndex]))
                bucket.append([pepIndex])
            else:
                bucket[tupleIndex].append(pepIndex)
        
        for subBuck in bucket:
            if len(subBuck) > 1:
                totalBucketSets.append(subBuck)
        
        
        print("Bucket length ", len(bucket))
        print("Current Band No. ", bandIndex)
        print('Time consuming ', int((time.time() - start_time) / 60 ), 'min')
        
    np.save('TempData/totalBucketSets.npy', totalBucketSets)
    print("Unique sets in totalBucketSets No.", len(set(totalBucketSets))) 
