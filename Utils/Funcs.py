# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 09:51:05 2021

@author: slaiad
"""
import numpy as np
from math import floor

def divideInputList(inputList, size):
    dividedList = []
    
    quantile = int(np.ceil(len(inputList)/size))
    
    realSize = int(np.ceil(len(inputList)/quantile))
    for i in range(realSize):
        dividedList.append(inputList[i*quantile : min((i+1)*quantile, len(inputList))])
    
    return dividedList

def binarySearch(arr, l, r, x): 
    if r >= l: 
        mid = int(l + (r - l)/2)
        # print(mid, x)
        if arr[mid] == x: 
            return mid 
        elif arr[mid] > x: 
            return binarySearch(arr, l, mid-1, x) 
        else: 
            return binarySearch(arr, mid+1, r, x) 
    else: 
        return max(l-1, 0)
    
def binarySearchList(arr, l,r, sortX):
    ids = []
    
    if len(sortX) == 0:
        return ids
    
    midX = floor(len(sortX)/2)
    index = binarySearch(arr, l, r, sortX[midX])
    ids.append(index)
    
    Lx = sortX[0:midX]
    Rx = sortX[midX + 1:]
    if len(sortX) == 1:
        return ids
    
    ids.extend(binarySearchList(arr, 0, index, Lx))
    # print('idsaf',ids)
    ids.extend(binarySearchList(arr, index + 1, len(arr), Rx))
    return ids
    
  
def sendEmail():
    import yagmail
    yag = yagmail.SMTP(user = '2578027596@qq.com', password = 'alexlptryglrdhii', host = 'smtp.qq.com')
    yag.send(to = ['laishengzhi1994@163.com'], subject = 'Finished code', contents = ['Finished code', 'TempData/ptmRevise.xlsx'])

MU = 0
def gaussian(x, SIGMA2):
    
    f = 1/sqrt(2*pi*SIGMA2)*exp(-(x-MU)**2/2/SIGMA2)
    return f

if __name__ == '__main__':
    ides = []
    arr = [0,100,200,300,600,700,800,900,1000,1100,1200,2000,2500,3500,3700,4200,4600,8500] 
    X = [1,300,700, 1200,2000,8499]
    i = 0
    ids = binarySearchList(arr, 0, len(arr), X)
    print(ids)