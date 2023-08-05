'''
Descripttion: 
version: 
Author: zpliu
Date: 2022-02-11 20:31:48
LastEditors: zpliu
LastEditTime: 2022-02-12 10:43:19
@param: 
'''
import pandas as pd 
from dap_g import dap_g_formed
from conserved_causal import getConservedRegion
from scipy.stats import pearsonr
import sys 
import pysam
import re 
import logging
logging.basicConfig(level=logging.INFO)
#------------------------------------------
#* 把每个eVariant cluster和eGene进行匹配，第三列是连锁不平衡r2值
#* 得到Block的区间范围,以及对应的SNP id
#------------------------------------------
# BiasQTL=pd.read_csv("./0DPA_Bias_eQTL.txt",header=0,index_col=None,sep="\t")
# BiasQTL=BiasQTL.loc[(BiasQTL['eGene']!="-")&(BiasQTL['QTLlocal']!="trans")]
# BiasQTL['SNPcluster']=BiasQTL['SNPcluster'].apply(lambda x:x.split("-")[0])
# BiasQTL[['eGene','SNPcluster']].values


#! 对每个QTL进行分析
def QTL_sub_seq(eGene,eVariant,SNPCluster,VCFFile,HomoeologGene,GeneLocation,HomoeologExpression,GenomeObject):
    SNPBlock_QTL=[]
    tmpData=SNPCluster.loc[SNPCluster['SNP_A']==eVariant]
    if tmpData.empty:
            SNPBlock_QTL.append(
                (eGene,eVariant,1)
            )
            
    else:
        for SNP,r2 in  tmpData[['SNP_B','R2']].values:
                SNPBlock_QTL.append(
                    (eGene,SNP,r2)
                )
    SNPBlock_QTL=pd.DataFrame(SNPBlock_QTL) 
    logging.info("change QTL format...")  
    #* 将每个SNP个格式进行整理，获得进行lastz比对的数据框
    QTLVariant=dap_g_formed(SNPBlock_QTL,VCFFile,HomoeologGene,GeneLocation)
    #* 过滤掉一些SNP
    QTLVariant=QTLVariant.loc[QTLVariant['CausalSite']!="error"]
    QTLVariant['CausalSite']=QTLVariant['CausalSite'].astype(int)
    QTLVariant['sub_TSS']=QTLVariant['sub_TSS'].astype(int)
    #* 将SNP进行lastz比对，获取亚组间保守的区域
    logging.info("run lastz...")
    QTLVariant['conservedRegion']=getConservedRegion(QTLVariant)
    #* 过滤掉没有比对到保守区域的SNPs
    QTLVariant=QTLVariant.loc[QTLVariant['conservedRegion']!="-_-_-_-"]
    if QTLVariant.empty:
        #*过滤完之后；Block区域的SNP在另外一个亚组中都没有找到
        return ("-","-","-") 
    else:
        tabix_Object = pysam.TabixFile(VCFFile)
        cor,pval,SNPCount=seq_conserved_BiasRatio(QTLVariant,GenomeObject,tabix_Object,HomoeologExpression)
        return cor,pval,SNPCount


#------------------------------------------------------------
#* 对定位到的QTL区间，分析其在两个亚组中碱基相同的数目，与Bias程度的关系
#------------------------------------------------------------
def get_seq(x):
    ref,alt=x[3:5]
    seqDict={
        '0/0':ref,
        '1/1':alt,
        '0/1':'N',
        './.':'N'
    } 
    resultSeq=''
    for seqCode in x[9:]:
        base=seqDict.get(seqCode)
        resultSeq+=base+"_" 
    return resultSeq.strip("_")

def BiasType(x):
    '''
    只分析Bias程度的绝对度量
    '''
    At, Dt = x
    At, Dt = float(At), float(Dt)
    if At < 0.5 and Dt < 0.5:
        return 0
    elif (At >= 0.5 and Dt < 0.5) or (At < 0.5 and Dt > 0.5):
        return 1
    else:
        return abs((At-Dt)/(At+Dt))


def seq_conserved_BiasRatio(conservedSNps,genomeObject,tabix_Object,Expression):
    #*提取这些SNP在每个样本中的基因型
    vcfHeader=tabix_Object.header[0].split("\t") 
    eChrom=conservedSNps.iloc[0,1]
    BlockStart=conservedSNps.iloc[0,2]
    BlockEnd=conservedSNps.iloc[-1,2]
    eGene=conservedSNps.iloc[0,0]
    subGene=conservedSNps.iloc[0,8]

    block_SNPs=pd.DataFrame(tabix_Object.fetch(start=BlockStart-1,end=BlockEnd,reference=eChrom))[0].str.split("\t",expand=True)
    block_SNPs.columns=vcfHeader
    filterSNPs=block_SNPs.loc[block_SNPs['ID'].isin(conservedSNps['CausalVar'])]
    #* 转置后最终得到Block区域对应的序列
    Block_seq=filterSNPs.apply(get_seq,axis=1).str.split("_",expand=True)
    Block_seq.index=conservedSNps['CausalVar']
    Block_seq.columns=vcfHeader[9:]
    Block_seq=Block_seq.T.apply(lambda x:"".join(x),axis=1)
    #* 另外一个亚组对应区域的序列信息
    subSeq=''    
    for region in conservedSNps['conservedRegion']:
        tmp1,tmp2,site=region.split("_")[0:3]
        chrom=tmp1+"_"+tmp2
        site=int(site)
        base=genomeObject.fetch(start=site-1,end=site,reference=chrom)
        subSeq+=base
    #* 计算表达水平
    if re.match('Ghir_A', eGene):
        genePair = eGene+"-"+subGene
    else:
        genePair = subGene+"-"+eGene
    geneExpression = Expression.loc[genePair].str.split(
                    "-", expand=True)
    geneExpression.columns = ['At', 'Dt']
    geneExpression['Bias'] = geneExpression.apply(BiasType, axis=1) 
    #* 分析每个样本中差异SNP的数量
    DifferenceBaseCount=[]
    for sample in geneExpression.index:
        count=0
        for a,b in zip(subSeq,Block_seq[sample]):
            #*基因位于不同的链方向，这个没有考虑
            if a!=b:
                count+=1
        DifferenceBaseCount.append(count/conservedSNps.shape[0])
    geneExpression['DifferenceBaseCount']=DifferenceBaseCount
    # geneExpression.to_csv("./test_expression.txt",header=True,index=True,sep="\t")
    cor,pvalue=pearsonr(geneExpression['DifferenceBaseCount'],geneExpression['Bias'])
    # print(geneExpression)
    # print(cor,pvalue)
    #* 返回相关系数以及LD-Block区域内SNP的数目
    return cor,pvalue,conservedSNps.shape[0]
    





if __name__=="__main__":
    #* top-eVariant的连锁信息
    QTLData=pd.read_csv(sys.argv[1],header=None,index_col=None,sep="\s+")
    SNPCluster=pd.read_csv(sys.argv[2],header=0,index_col=None,sep="\s+")
    # eGene,eVariant='Ghir_D01G000900','SNP1596657'
    vcfFile=sys.argv[3]
    HomoeologGeneFile=sys.argv[4]
    GeneTSSFile=sys.argv[5]
    stage=sys.argv[6]
    HomoeologExpression = pd.read_csv("/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/Allhomolog_" +
                                  str(stage)+"_FPKM_V2.txt", header=0, index_col=0, sep="\t")
    HomoeologExpression.columns = [
        re.sub("[0-9]*D", "", i) for i in HomoeologExpression.columns]
    genomeFile = '/data/cotton/MaojunWang/WMJ_fiberFullPopulationRNAseq/MappingFPKM/Ghir_Genome_Index/Ghirsutum_genome.fasta'
    genomeObject = pysam.FastaFile(genomeFile)
    out=[]
    print(QTLData)
    for eGene,eVariant in QTLData.values:
        cor,pval,SNPCount=QTL_sub_seq(
        eGene,eVariant,SNPCluster,vcfFile,HomoeologGeneFile,GeneTSSFile,HomoeologExpression,genomeObject
        )
        out.append(
          (eGene,eVariant,cor,pval,SNPCount)
        )
    out=pd.DataFrame(out)
    out.to_csv(sys.argv[7],header=False,index=False,sep="\t")
    