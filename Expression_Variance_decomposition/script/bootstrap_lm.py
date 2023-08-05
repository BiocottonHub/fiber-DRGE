'''
Descripttion: 对基因表达变异进行解构
version: 
Author: zpliu
Date: 2021-10-30 22:40:50
LastEditors: zpliu
LastEditTime: 2021-11-04 11:00:16
@param: 
'''
from numpy.lib.scimath import log
from rpy2.rinterface_lib.openrlib import LOGICAL
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from math import log2
from sklearn import preprocessing
import pandas as pd
import numpy as np
import sys
import logging 
logging.basicConfig(level=logging.INFO)



logging.info("loading expression matrix..")
stageArray = ['0DPA', '4DPA', '8DPA', '12DPA', '16DPA', '20DPA']
readCountData = []
AllsampleName = pd.read_csv("./All_sampleName.txt",
                            header=None, index_col=None, sep="\t")
# * get read count of each gene at times scheme
for stage in stageArray:
    read_count = pd.read_csv(
        "../Gene_align/read_count/"+stage+"/geneFPKM.txt", header=0, index_col=0, sep="\t")
    read_count.columns = AllsampleName.loc[AllsampleName[0].isin(
        read_count.columns)].index
    readCountData.append(read_count)


def get_lm_matrix(readCountArray, geneId):
    ''' format data to run lme4 regression model.
        + readCountArray: each time point read count.
    '''
    # * scale time to mean zeri variance 1.
    # ?  six time point
    time_array = preprocessing.scale(range(1, 7, 1))
    # *  expression  sampleId, scaled time
    readMatrix = np.empty((0, 3))
    for index, readData in enumerate(readCountArray):
        # * log2(FPKM+1) transform the FPKM
        timePoint = time_array[index]
        for sampleIndex, FPKM in readData.loc[geneId].iteritems():
            tmpData = [log2(FPKM+1), sampleIndex, timePoint]
            readMatrix = np.vstack([readMatrix, tmpData])
    readMatrix = pd.DataFrame(
        readMatrix, columns=['FPKM', 'sampleId', 'time'])
    readMatrix.to_csv("test.txt",header=True,index=False,sep="\t")
    return readMatrix


pandas2ri.activate()
robjects.r('''
        library(lme4)  
        library(aod)
        ''')

RunAnalyseMLMlogit = robjects.r(
    '''
       f <- function(readMatrix){
        #*---------------------------
        #* 1. log2(FPKM+1) 
        #* 2. sampleId 
        #* 3. time
        #*---------------------------
        # readMatrix$sampleId=as.factor(readMatrix$sampleId)
        # readMatrix$timepoint =as.factor(readMatrix$time)
        #*------------------------------
        #* Variance decomposition 
        #*------------------------------
        model=lmer('FPKM ~ 1 + (1 | time) +(1|sampleId)', data = readMatrix)
        
        timeVariant=as.data.frame(summary(model)$varcor)[2,4]
        sampleVariant=as.data.frame(summary(model)$varcor)[1,4]
        residualVariant=as.data.frame(summary(model)$varcor)[3,4]
        result <- c(timeVariant,sampleVariant,residualVariant)
        return (result)
       } 
    '''
)


def bootstrap(lm_matrix, sampleTimes=1000):
    VariantResult = []
    for item in range(0, sampleTimes):
        tmeData = lm_matrix.sample(frac=1, replace=True)
        timeVariant, sampleVariant, residualVariant = RunAnalyseMLMlogit(
            tmeData)
        VariantResult.append(
            (
                timeVariant/(timeVariant+sampleVariant+residualVariant),
                sampleVariant/(timeVariant+sampleVariant+residualVariant)
            )
        )
    #! 95% confidence interval
    # * 2.5 Quantile and 97.5 Quantile
    VariantResult = pd.DataFrame(VariantResult, columns=['tiems', 'samples'])
    confidenceInterVal = []
    for index, Variant in VariantResult.iteritems():
        meanData = np.mean(Variant)
        lowData = np.percentile(Variant, 2.5)
        highData = np.percentile(Variant, 97.5)
        confidenceInterVal.append((lowData, meanData, highData))
    # ----------------------------------
    # * confidence InterVal Array
    # [
    # [ low, mean,high]
    # ]
    # ----------------------------------
    return confidenceInterVal



if __name__ =="__main__":
    geneList=pd.read_csv(sys.argv[1],header=None,index_col=None,sep="\t")
    out=[]
    for geneId in geneList[0]:
        logging.info(f'deal with {geneId}')
        dataMatrix=get_lm_matrix(readCountData,geneId)
        try:
            timeInterVal,sampleVal=bootstrap(dataMatrix)
            out.append(
                [geneId]+list(timeInterVal)+list(sampleVal)
            )
        except:
            out.append(
                [geneId,'error','error','error','error','error','error']
            )
    out=pd.DataFrame(out,columns=['GeneId','time_low','time_mean','time_high','sample_low','sample_mean','sample_high'])
    out.to_csv(sys.argv[2],header=True,index=False,sep="\t")

