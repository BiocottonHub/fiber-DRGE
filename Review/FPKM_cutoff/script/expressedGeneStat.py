'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-05-09 16:03:26
LastEditors: zpliu
LastEditTime: 2023-05-09 16:05:03
@param: 
'''
import pandas as pd 
import sys 
import re 
#TODO 根据基因所在染色体来区分At、Dt以及Scaffold上
geneBed=pd.read_csv(
    "/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/03express_gene/eQTLgenes_Allstages_340samples/fastQTL_peer_20/Gene_information/All_gene_TSS.bed",
    header=0,index_col=4,sep="\t"
)
stage_0DPA=pd.read_csv(sys.argv[1],header=None,index_col=None,sep="\t")
stage_4DPA=pd.read_csv(sys.argv[2],header=None,index_col=None,sep="\t")
stage_8DPA=pd.read_csv(sys.argv[3],header=None,index_col=None,sep="\t")
stage_12DPA=pd.read_csv(sys.argv[4],header=None,index_col=None,sep="\t")
stage_16DPA=pd.read_csv(sys.argv[5],header=None,index_col=None,sep="\t")
stage_20DPA=pd.read_csv(sys.argv[6],header=None,index_col=None,sep="\t") 


mergeData=pd.concat([stage_0DPA,stage_4DPA],axis=0)
mergeData=pd.concat([mergeData,stage_8DPA],axis=0)
mergeData=pd.concat([mergeData,stage_12DPA],axis=0)
mergeData=pd.concat([mergeData,stage_16DPA],axis=0)
mergeData=pd.concat([mergeData,stage_20DPA],axis=0)

geneStat=[]
for gene in mergeData[0].unique():
    Chrom=geneBed.loc[gene][0]
    if re.match("Ghir_A",Chrom):
        geneStat.append('At')
    elif re.match("Ghir_D",Chrom):
        geneStat.append('Dt')
    else:
        geneStat.append('Scaffold')
pd.value_counts(geneStat).to_csv(sys.argv[7],header=False,index=True,sep="\t")