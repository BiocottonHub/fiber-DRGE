'''
Descripttion: 
version: 
Author: zpliu
Date: 2022-01-14 15:52:12
LastEditors: zpliu
LastEditTime: 2023-03-06 20:21:26
@param: 
'''
#! 获取lastz比对的指标
import os
import pysam
import pandas as pd 
from tempfile import NamedTemporaryFile
from pyliftover import LiftOver
from extract_seq import extract_sequence
lastZ = '/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/software/lastz-1.04.03/lastz-distrib/bin/lastz'
axtChain = '/public/home/software/opt/bio/software/ucsc_kentUtils/v389/bin/axtChain'
chainSwap = '/public/home/software/opt/bio/software/ucsc_kentUtils/v389/bin/chainSwap'
chainSort = '/public/home/software/opt/bio/software/ucsc_kentUtils/v389/bin/chainSort'


def lastz(seqrchSeq, targetSeq, searchId, targetSeqReal_start, targetChrom,
          eQTL_seqFile, eGene_seqFile, axtFile):
    '''map sequence in subgenome
        eQTL_seqFile: NamedTemporaryFile Object
    '''
    # * 将序列数据写入文件
    eQTL_seqFile.write(seqrchSeq)
    eGene_seqFile.write(targetSeq)
    eQTL_seqFile.seek(0)
    eGene_seqFile.seek(0)
    # * 执行比对
    # --filter=identity:85 --filter=coverage:75
    os.system(
        "{} {} {} K=3000 L=3000 H=2000 Y=5000 E=55 T=2 \
            O=600 --filter=identity:80 --filter=coverage:85 \
            --verbosity=0  --format=general  >{}".format(
            lastZ, eGene_seqFile.name, eQTL_seqFile.name, axtFile.name
        )
    )
    outData=pd.read_csv(axtFile.name,header=0,index_col=None,sep="\s+")
    return outData
    
  