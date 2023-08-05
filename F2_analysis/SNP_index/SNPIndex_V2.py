'''
Descripttion: 
version: 
Author: zpliu
Date: 2022-11-11 10:10:40
LastEditors: zpliu
LastEditTime: 2022-11-13 11:21:36
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
SNP_position=BedTool("../impute_V2/SNP_impute_MAF_0.2.txt")
genomeObject=pysam.TabixFile("../F2_Q1000_short.recode.vcf.gz")
#TODO 计算基因型占比，只打乱极端短池的数据
poolSample=pd.read_csv(sys.argv[1],header=None,index_col=None,sep="\t")

poolSample[0]=poolSample[0].apply(lambda x:str(x))
genomeObject.header[-1].split("\t")
#* 提取对应样本在vcf文件中属于哪一列
SampleIndex=zip(range(0,len(genomeObject.header[-1].split("\t"))),genomeObject.header[-1].split("\t"))
SampleIndex=pd.DataFrame(list(SampleIndex))
#* 保持长池不变，只打乱短池
shortSample=SampleIndex.loc[SampleIndex[1].isin(poolSample.loc[poolSample[1]=='pool1'][0])]
longSample=SampleIndex.loc[SampleIndex[1].isin(poolSample.loc[poolSample[1]=='pool4'][0])]

# -----------------------------------------------------------
# * 统计指定窗口中基因型和其中一个亲本一致的占比
# -----------------------------------------------------------
genomeBed = pd.read_csv("genome_window_V2.bed",
                        header=None, index_col=None, sep="\t")
# genomeBed[0]=genomeBed[0].map(ChromDict) 

geneBedCount=[]
for chrom,start,end in genomeBed.values:
    SNPData=list(genomeObject.fetch(
            reference=chrom,
            start=start,end=end
    ))
    #* 构建SNP数据框, 只筛选MAF不超过0.2的
    windowSNP=SNP_position.filter(lambda x: x[0]==chrom and int(x[1])>=start and int(x[1])<=end) 
    # print(windowSNP)
    SNPsite=[ i[1] for i in windowSNP ]
    if len(SNPsite)==0:
        #* 该窗口不包含SNP
        continue
    SNPData=[i.split("\t") for i in SNPData if i.split("\t")[1] in SNPsite]
    longInterval=[]
    shortInterval=[]
    homozygousCount=0
    for item in SNPData:
        longSNP=item[9].split(":")[0]
        if longSNP!="0/0" and longSNP!="1/1":
            #* ZY043样本基因型数据不纯,所有的位点如果都不存？？
            longInterval.append((0,0))
            shortInterval.append((0,0))
            continue
        homozygousCount+=1
        longpool_longReads=0
        longpool_shortReads=0
        shortpool_longReads=0
        shortpool_shortReads=0
        # print(item[0],item[1],item[9],item[10])
        for longIndex in longSample[0].values:
            if item[longIndex].split(":")[0]=="./.":
                # 该位点被过滤，没有read map
                continue
            ref,alt=[i for i in item[longIndex].split(":")[1].split(",")]
            
            try:
                ref=int(ref)
            except ValueError:
                ref=0
            try:
                alt=int(alt)
            except ValueError:
                alt=0
            itemSNP=item[longIndex].split(":")[0]
            if longSNP=="0/0":
                longpool_longReads+=ref
                longpool_shortReads+=alt
            else:
                longpool_longReads+=alt
                longpool_shortReads+=ref

        for shortIndex in shortSample[0].values:
            if item[shortIndex].split(":")[0]=="./.":
                continue
            ref,alt=[i for i in item[shortIndex].split(":")[1].split(",")]
            try:
                ref=int(ref)
            except ValueError:
                ref=0
            try:
                alt=int(alt)
            except ValueError:
                alt=0
            itemSNP=item[shortIndex].split(":")[0]
            if longSNP=="0/0":
                shortpool_longReads+=ref
                shortpool_shortReads+=alt
            else:
                shortpool_longReads+=alt
                shortpool_shortReads+=ref

        longInterval.append(
            (longpool_longReads,longpool_shortReads)
        )
        shortInterval.append(
            (shortpool_longReads,shortpool_shortReads)
        )
    longInterval=pd.DataFrame(longInterval)
    longtotalReads=longInterval[0].sum()+longInterval[1].sum()
    shortInterval=pd.DataFrame(shortInterval)
    shorttotalReads=shortInterval[0].sum()+shortInterval[1].sum()
    if longtotalReads==0:
        longRatio=0 
    else:
        #* 该区域内所有SNP与长纤维样本相同的平均比例
        longRatio=longInterval[0].sum()/longtotalReads
    if shorttotalReads==0:
        shortRatio=0 
    else:
        #* 该区域内支持长纤维基因的read占总read的比例
        shortRatio=shortInterval[0].sum()/shorttotalReads
    geneBedCount.append(
        ( chrom,start,end,homozygousCount,longRatio,shortRatio)
    )   
geneBedCount=pd.DataFrame(geneBedCount,columns=['Chrom','start','end','SNPCount','longGenotypeRatio','shortGenotypeRatio'])  
# geneBedCount['Chrom']=geneBedCount['Chrom'].map(ChromDict1) 
geneBedCount.to_csv(sys.argv[2],header=True,index=False,sep="\t")