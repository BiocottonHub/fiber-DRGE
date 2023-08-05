'''
Descripttion: 
version: 
Author: zpliu
Date: 2022-01-14 16:16:32
LastEditors: zpliu
LastEditTime: 2023-03-06 20:29:37
@param: 
#! lastz换一种输出格式
'''
import pandas as pd 
import pysam
from lastz_v2 import lastz
import sys 
from extract_seq import extract_sequence
from tempfile import NamedTemporaryFile

def get_seq_fragment(site,fragmentLength):
    '''上下游指定范围的序列
    '''
    if site>fragmentLength:
        eStart=site-fragmentLength
    else:
        eStart=1
    eEnd=site+fragmentLength
    return (eStart,eEnd)

#* 临时文件



def getConservedRegion(QTLVariant):
    genomeFile = '/data/cotton/MaojunWang/WMJ_fiberFullPopulationRNAseq/MappingFPKM/Ghir_Genome_Index/Ghirsutum_genome.fasta'
    genomeObject = pysam.FastaFile(genomeFile)
    eQTL_seqFile = NamedTemporaryFile(mode='w+', encoding='utf-8')
    eGene_seqFile = NamedTemporaryFile(mode='w+', encoding='utf-8')
    generalFile = NamedTemporaryFile(mode='w+', encoding='utf-8')
    out=[]
    for value in QTLVariant.values:
        #* 清空临时文件
        eQTL_seqFile.truncate()
        eGene_seqFile.truncate()
        generalFile.truncate()
        eChrom=value[1]
        #* 左右各50bp
        eSite=value[2]
        eStart,eEnd=get_seq_fragment(eSite,200)
        eQTLid=value[3]
        subChrom=value[6]
        #* 上下游各1M范围
        subTSS=value[7]
        subStart,subEnd=get_seq_fragment(subTSS,2000000)
        subId=value[8]
        #* 序列长度信息
        QTLseq=extract_sequence(eChrom,eStart,eEnd,eQTLid,genomeObject)
        subseq=extract_sequence(subChrom,subStart,subEnd,subId,genomeObject)
        #print(subChrom,subStart,subEnd) 
        standSame_or_not=value[9]
        outDataFram=lastz(
            QTLseq,subseq,eQTLid,subStart,subChrom,eQTL_seqFile,
            eGene_seqFile,generalFile
        )
        # out.append("{}_{}_{}_{}_{}".format(alignChrom,alignSite,alignScore,alignStand,mappingRegionCount))
        # print(outDataFram.columns)
        for lastzItem in outDataFram.values:
            out.append(
                (
                    value[0],eQTLid,subId,lastzItem[0],lastzItem[12],lastzItem[14]
                )
            )
    eQTL_seqFile.close()
    eGene_seqFile.close()
    generalFile.close()
    return out
if __name__=="__main__":
    QTLVariant=pd.read_csv(sys.argv[1],header=0,index_col=None,sep="\t")

    #* 过滤掉error行，并且修改数据类型
    QTLVariant=QTLVariant.loc[QTLVariant['CausalSite']!="error"]
    QTLVariant['CausalSite']=QTLVariant['CausalSite'].astype(int)

    QTLVariant['conservedRegion']=getConservedRegion(QTLVariant)

    QTLVariant.to_csv(sys.argv[2],header=True,index=False,sep="\t")
