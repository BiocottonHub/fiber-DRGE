'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-05-09 15:13:44
LastEditors: zpliu
LastEditTime: 2023-05-16 19:53:15
@param: 
'''
import pandas as pd 
from statsmodels.stats.multitest import multipletests
from scipy.stats import wilcoxon
import sys 

#TODO 根据expressioned 基因分析同源基因是否已经表达
Homoeologous=pd.read_csv(
    "/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/Homoeolog_geneId.txt",
    header=None,index_col=None,sep="\t"
)
# expressionData=pd.read_csv(
#     "/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/MappingFPKM/TM1/AllGenes/0DPA_AllSamples_FPKM.txt",
#     header=0,index_col=0,sep="\t"4
# )
# condition=0.01 
expressionData=pd.read_csv(
    sys.argv[1],
    header=0,index_col=0,sep="\t"
)
condition=float(sys.argv[2])

# ExpressedGeneList=pd.read_csv(
#     "0DPA/0DPA_expressed_0.01.txt",header=None,index_col=None,sep="\t"
# )
ExpressedGeneList=pd.read_csv(
    sys.argv[3],header=None,index_col=None,sep="\t"
)




def BiasSampleCount(homoeologousExpression,condition):
    BiasACount=0
    BiasDCount=0
    for At,Dt in homoeologousExpression:
        if At>=condition and Dt>=condition:
            if (At-Dt)/(At+Dt) >=0.33:
                BiasACount+=1 
            elif (At-Dt)/(At+Dt) <=-0.33:
                BiasDCount+=1 
        elif At>=condition and Dt<condition:
            BiasACount+=1 
        elif At<condition and Dt>=condition:
            BiasDCount+=1 
    return BiasACount,BiasDCount

def homoeologousExpressionNot(genePair,expressedGeneList,conditions,totalSampleCount):
    At,Dt=genePair
    if At not in ExpressedGeneList[0].values and Dt not in ExpressedGeneList[0].values:
        return ("BiasN_noExpression",1,"BiasN")
    else:
        #? At和Dt Bias Sample count
        At_Expression=expressionData.loc[At]
        Dt_Expression=expressionData.loc[Dt]
        v_,pavl=wilcoxon(
            At_Expression.values,
            Dt_Expression.values,
            alternative='two-sided'
        )
        # BiasACount,BiasDCount=BiasSampleCount(zip(At_Expression,Dt_Expression),conditions)
        #! FPKM expression的标准是两个亚基因组至少一个表达量超过0.5 FPKM
        BiasACount,BiasDCount=BiasSampleCount(zip(At_Expression,Dt_Expression),0.5)
        if BiasACount/totalSampleCount>=0.05 and BiasDCount/totalSampleCount<0.05:
            return ("BiasA",pavl,"BiasA") 
        elif BiasACount/totalSampleCount<0.05 and BiasDCount/totalSampleCount>=0.05:
            return ("BiasD",pavl,"BiasD") 
        elif BiasACount/totalSampleCount>=0.05 and BiasDCount/totalSampleCount>=0.05:
            return ("BidirectBias",pavl,"BidirectBias") 
        else:
            return ("BiasN_Expression",1,"BiasN")  

#TODO 分析每个同源基因对是否发生Bias
out=[]
for genePair in Homoeologous.values:
    BiasTypeV1,pval,BiasType=homoeologousExpressionNot(genePair,ExpressedGeneList,0.01,expressionData.shape[1])
    out.append(
        ("-".join(genePair),BiasTypeV1,pval,BiasType)
    )
out=pd.DataFrame(out,columns=['genePair','BiasType','pval','BiasTypeV2'])
NoBiasData=out.loc[out['BiasType'].isin(['BiasN_noExpression','BiasN_Expression'])]
NoBiasData['fdr']=1
BiasData=out.loc[~out['BiasType'].isin(['BiasN_noExpression','BiasN_Expression'])]
qvalueList=multipletests(BiasData['pval'].values,method='fdr_bh')[1]
BiasData['fdr']=qvalueList
BiasData['BiasTypeV2']=BiasData.apply(lambda x: x[3] if x[4]<=0.05 else 'BiasN',axis=1)
mergeData=pd.concat(
    [BiasData,NoBiasData],axis=0
)
mergeData.to_csv(sys.argv[4],header=True,index=False,sep="\t")
