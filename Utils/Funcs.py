# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 09:51:05 2021

@author: slaiad
"""
import numpy as np

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
  
        if arr[mid] == x: 
            return mid 
        elif arr[mid] > x: 
            return binarySearch(arr, l, mid-1, x) 
        else: 
            return binarySearch(arr, mid+1, r, x) 
    else: 
        return max(l-1, 0)
  
def sendEmail():
    import yagmail
    yag = yagmail.SMTP(user = '2578027596@qq.com', password = 'alexlptryglrdhii', host = 'smtp.qq.com')
    yag.send(to = ['laishengzhi1994@163.com'], subject = 'Finished code', contents = ['Finished code'])

if __name__ == '__main__':
    
    arr = [1,2,3,4,5,6,7,8,9,10,11,12] 
    x = 1.9
    print(binarySearch(arr, 0, len(arr)-1, x))