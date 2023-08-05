'''
Descripttion: 
version: 
Author: zpliu
Date: 2022-11-11 10:10:40
LastEditors: zpliu
LastEditTime: 2022-11-15 18:08:48
@param: 
'''
import pandas as pd 
import pysam
import sys 
from pybedtools import BedTool
ChromDict={
    'Ghir_A01':1,'Ghir_A02':2,'Ghir_A03':3,'Ghir_A04':4,'Ghir_A05':5,'Ghir_A06':6,
    'Ghir_A07':7,'Ghir_A08':8,'Ghir_A09':9,'Ghir_A10':10,'Ghir_A11':11,'Ghir_A12':12,
    'Ghir_A13':13,
    'Ghir_D01':14,'Ghir_D02':15,'Ghir_D03':16,'Ghir_D04':17,'Ghir_D05':18,'Ghir_D06':19,
    'Ghir_D07':20,'Ghir_D08':21,'Ghir_D09':22,'Ghir_D10':23,'Ghir_D11':24,'Ghir_D12':25,
    'Ghir_D13':26,
}
ChromDict1=dict(zip(
    ChromDict.values(),ChromDict.keys()
))
genomeObject=pysam.TabixFile("../impute_V2/All_chrom_fastlmm.vcf.gz")
#TODO 计算基因型占比，只打乱极端短池的数据
poolSample=pd.read_csv(sys.argv[1],header=None,index_col=None,sep="\t")

poolSample[0]=poolSample[0].apply(lambda x:str(x))

#* 提取对应样本在vcf文件中属于哪一列
SampleIndex=zip(range(0,len(genomeObject.header[-1].split("\t"))),genomeObject.header[-1].split("\t"))
SampleIndex=pd.DataFrame(list(SampleIndex))
#* 保持长池不变，只打乱短池
shortSample=SampleIndex.loc[SampleIndex[1].isin(poolSample.loc[poolSample[1]=='pool1'][0])]
longSample=SampleIndex.loc[SampleIndex[1].isin(poolSample.loc[poolSample[1]=='pool4'][0])]

# -----------------------------------------------------------
# * 统计指定窗口中基因型和其中一个亲本一致的占比
# -----------------------------------------------------------
genomeBed = pd.read_csv("genome_window.bed",
                        header=None, index_col=None, sep="\t")
# genomeBed[0]=genomeBed[0].map(ChromDict) 

geneBedCount=[]
for chrom,start,end in genomeBed.values:
    SNPData=list(genomeObject.fetch(
            reference=chrom,
            start=start,end=end
    ))
    #* 长根池中，与ZY043基因型相同的占比
    SNPData=[i.split("\t") for i in SNPData]
    longInterval=[]
    shortInterval=[]
    SNPcount=0
    for item in SNPData:
        if item[9] in ['0/0','1/1'] and item[10] in ['0/0','1/1']:
            SNPcount+=1
            longSNP=item[9]
            longpoolCount=0
            shortpoolCount=0
            for longIndex in longSample[0].values:
                if item[longIndex]==longSNP:
                    longpoolCount+=1
            for shortIndex in shortSample[0].values:
                if item[shortIndex]==longSNP:
                    shortpoolCount+=1
            longInterval.append(
                longpoolCount
            )
            shortInterval.append(
                shortpoolCount
            )
        else:
            continue
    if len(longInterval)==0:
        longRatio=0 
    else:
        #* 该区域内所有SNP与长纤维样本相同的平均比例
        longRatio=sum(longInterval)/(50*len(longInterval)) 
    if len(shortInterval)==0:
        shortRatio=0 
    else:
        #* 该区域内所有SNP与长纤维样本相同的平均比例
        shortRatio=sum(shortInterval)/(50*len(shortInterval)) 
    geneBedCount.append(
        ( chrom,start,end,SNPcount,longRatio,shortRatio)
    )   
geneBedCount=pd.DataFrame(geneBedCount,columns=['Chrom','start','end','SNPCount','longGenotypeRatio','shortGenotypeRatio'])  
# geneBedCount['Chrom']=geneBedCount['Chrom'].map(ChromDict1) 
geneBedCount.to_csv(sys.argv[2],header=True,index=False,sep="\t")