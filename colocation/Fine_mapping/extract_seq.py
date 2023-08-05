'''
Descripttion: 
version: 
Author: zpliu
Date: 2022-01-14 14:14:58
LastEditors: zpliu
LastEditTime: 2022-01-14 14:16:12
@param: 
'''
import pandas as pd
import pysam


def extract_sequence(Chrom, start, end,seqId,genomeObject):
    '''extract QTL sequence
    '''
    sequence = genomeObject.fetch(
        start=start-1, end=end, region=Chrom)
    # print(">"+QTLid+"\n"+sequence+"\n")
    return ">"+seqId+"\n"+sequence 

if __name__ =="__main__":
    genomeFile = '/data/cotton/MaojunWang/WMJ_fiberFullPopulationRNAseq/MappingFPKM/Ghir_Genome_Index/Ghirsutum_genome.fasta'
    genomeObject = pysam.FastaFile(genomeFile)
    print(extract_sequence('Ghir_A01',1,10,'test',genomeObject))