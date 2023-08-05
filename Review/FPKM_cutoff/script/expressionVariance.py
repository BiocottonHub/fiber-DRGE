'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-05-09 20:54:09
LastEditors: zpliu
LastEditTime: 2023-05-09 21:08:43
@param: 
'''
#TODO 分析每个时期的基因在指定cutoff下，基因expression variation的变化范围
from email import header
import pandas as pd 
import numpy as np 
import sys 

#TODO 分析每个时期的基因在指定cutoff下，基因expression variation的变化范围
cutoff=sys.argv[1]
expressedGeneArray=[]
expressionDataArray=[]
for stage in ["0DPA", "4DPA", "8DPA", "12DPA", "16DPA", "20DPA"]:
    expressedGeneArray.append(
        pd.read_csv("{}/{}_expressed_{}.txt".format(
            stage,stage,cutoff
            ),index_col=None,header=None,sep="\t"
        ),
    )
    expressionDataArray.append(
        pd.read_csv(
            "/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/MappingFPKM/TM1/AllGenes/StringTie_V1/{}_AllSamples_FPKM.txt".format(
                    stage
            ),
                header=0,index_col=0,sep="\t"
            )
    )



def fold_change(expression_quantile_25,expression_quantile_75):
    if expression_quantile_25 <0.5 and expression_quantile_75<0.5:
        return 1
    elif expression_quantile_25!=0 and expression_quantile_75>0.5:
        return expression_quantile_75/expression_quantile_25
    elif expression_quantile_25==0:
        #! 默认两倍差异
        return 2
expressionVariantData=pd.DataFrame()
for expressedGeneList,expressData in zip(expressedGeneArray,expressionDataArray):
    expressionQuantile=expressData.loc[expressedGeneList[0].values].quantile(
        [.05,.95],axis=1
    )
    expressionFoldChange=expressionQuantile.apply(
    lambda x:fold_change(x[.05],x[.95]),axis=0
    )        
    foldChangeGroup=pd.cut(
        expressionFoldChange,
        bins=[
            1,2,4,8,16,np.inf
        ],
        labels=[
            '1-2','2-4','4-8','8-16','>16'
        ],
        right=False
    )
    expressionVariantData=pd.concat([expressionVariantData,foldChangeGroup.value_counts()],axis=1)

expressionVariantData.columns=["0DPA", "4DPA", "8DPA", "12DPA", "16DPA", "20DPA"]
expressionVariantData.to_csv(sys.argv[2],header=True,index=True,sep="\t")