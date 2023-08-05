'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-02-23 10:42:43
LastEditors: zpliu
LastEditTime: 2023-02-23 10:47:04
@param: 
'''
#---------------------------------------------------------
#TODO 随机挑选1000个eQTL对应的lead SNP进行permutation 分析
#---------------------------------------------------------
from pysnptools.snpreader import Bed
import pandas as pd 
import numpy as np 
import sys 
import torch 
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from scipy.stats import norm
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

snp_on_disk=Bed("/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/Genotype/Trait_376.bed")
snpdata = snp_on_disk.read()
#------------------------------
#* SNP数据框
#------------------------------
snpdata=pd.DataFrame(snpdata.val,columns=snpdata.sid,index=snp_on_disk.iid[:,1]) 
#* eQTL数据
cis_eQTL=pd.read_csv("./All_stage_cis_eQTL.txt",header=0,index_col=None,sep="\t")
trans_eQTL=pd.read_csv("./All_stage_trans_eQTL.txt",header=0,index_col=None,sep="\t")
MergeData=pd.concat(
    [cis_eQTL,trans_eQTL],axis=0
)
MergeData=MergeData.loc[MergeData['rank']!="minor"]
MergeData.index=range(0,MergeData.shape[0])
#* 随机选取5000个eQTL进行permutate分析
np.random.seed(2023)
randomeQTLIndex=np.random.choice(MergeData.index,5000,replace=False) 

#* 随机选取相应时期的eQTL

stage=sys.argv[1]
geneExpression=pd.read_csv(
    "/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/MappingFPKM/TM1/{}_AllSamples_All_stage.txt".format(
        stage
    ),
    header=0,index_col=0,sep="\t"
    )
geneExpression.columns=[i.replace('Sample','S') for i in geneExpression.columns]    
randmoData=MergeData.loc[randomeQTLIndex]
randomData_0DPA=randmoData.loc[randmoData['stage']==stage] 

#-------------------------------------------------
#* 基因型映射
#? './.' as na 作为参考基因型  0
#? '0/1' as 1  作为非参考基因型  -1
#? '0/0' as 2  作为参考基因型  1
#? '1/1' as 0  作为参考基因型  -1
#-------------------------------------------------
genotypeMap={
    0:-1,
    2.0:1,
    1.0:-1,
}
out=[]
normScale=StandardScaler()
count=1
for geneId,SNPid in randomData_0DPA[['eGene','top_variant']].values:
    count+=1
    #* 执行10000次permutate
    sampleCount=geneExpression.shape[1]
    randomIndex=np.array(
        [np.random.permutation(range(0,sampleCount)) for i in range(10000)]
        )
    observedIndex=np.array(range(0,sampleCount)).reshape(1,-1)    
    permutate_inx=torch.LongTensor(
        np.vstack(
        [observedIndex,randomIndex]
        )
    ).to(device)
    #* 基因表型向量
    geneFPKM_t=torch.tensor(
        geneExpression.loc[geneId].values,
        dtype=torch.float
    ).to(device)
    #* 基因permutate后的表达量向量
    permutate_t=geneFPKM_t[permutate_inx] 
    #* 映射后的基因型信息
    leadSNPData=snpdata[[SNPid]].loc[
        geneExpression.columns
        ].applymap(lambda x:genotypeMap.get(x,0)).values

    leadSNPData_t=torch.tensor(
        leadSNPData,
        dtype=torch.float
    ).to(device) 
    #* 矩阵乘法
    FPKM_diff=np.array(
    torch.mm(permutate_t,leadSNPData_t
    )).reshape(1,-1)[0] 
    #* 将permutate数据进行正态转换
    normData=normScale.fit_transform(FPKM_diff.reshape(-1,1))
    print(normData)
    out.append(
        normData.reshape(1,-1)[0]
    )

pd.DataFrame(out).to_csv("./eQTL_noise_check/{}_permutate.txt".format(stage),header=None,index=False,sep="\t")