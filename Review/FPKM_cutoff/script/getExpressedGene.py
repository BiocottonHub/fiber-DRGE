'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-05-09 10:14:53
LastEditors: zpliu
LastEditTime: 2023-05-09 10:35:56
@param: 
'''
from multiprocessing import Condition
import pandas as pd 
import sys 


stageData=sys.argv[1]
ExpressionData=pd.read_csv(
    stageData,
    header=0,index_col=0,sep="\t"
) 

# ConditionFPKM=0.1

ConditionFPKM=float(sys.argv[2])

def filterExpression(sampleArray,conditions,sampleCount):
    totalExpressedFPKM=0
    expressedCount=0
    for expression in sampleArray:
        if expression>=conditions:
            expressedCount+=1
            totalExpressedFPKM+=expression
    if expressedCount/sampleCount>=0.05 and totalExpressedFPKM>=10:
        return True 
    else:
        return False


#* 
filterExpressioned=ExpressionData.apply(
    lambda x:filterExpression(x,ConditionFPKM,ExpressionData.shape[1]),axis=1
)
#* 
with open(sys.argv[3],'w') as File:
    for i in ExpressionData.loc[filterExpressioned].index:
        File.write(i+"\n")