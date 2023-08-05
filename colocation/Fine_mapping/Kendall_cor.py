'''
Descripttion: 
version: 
Author: zpliu
Date: 2022-02-09 20:06:45
LastEditors: zpliu
LastEditTime: 2022-02-09 20:51:18
@param: 
'''
import pandas as pd
import sys
from scipy.stats import kendalltau
import logging
import pysam
import re
logging.basicConfig(level=logging.INFO)
# lo = LiftOver('./../../../sub_variants/TM1_eQTL.chain')
QTLMap = {
    1: "Ghir_A01", 2: "Ghir_A02", 3: "Ghir_A03", 4: "Ghir_A04", 5: "Ghir_A05", 6: "Ghir_A06", 7: "Ghir_A07",
    8: "Ghir_A08", 9: "Ghir_A09", 10: "Ghir_A10", 11: "Ghir_A11", 12: "Ghir_A12", 13: "Ghir_A13", 14: "Ghir_D01",
    15: "Ghir_D02", 16: "Ghir_D03", 17: "Ghir_D04", 18: "Ghir_D05", 19: "Ghir_D06", 20: "Ghir_D07",
    21: "Ghir_D08", 22: "Ghir_D09", 23: "Ghir_D10", 24: "Ghir_D11", 25: "Ghir_D12", 26: "Ghir_D13"
}
genomeFile = '/data/cotton/MaojunWang/WMJ_fiberFullPopulationRNAseq/MappingFPKM/Ghir_Genome_Index/Ghirsutum_genome.fasta'
genomeObject = pysam.FastaFile(genomeFile)

# --------------------------------------------
# * 提取基因表达文件和基因型文件
# --------------------------------------------
stage = sys.argv[1]
HomoeologExpression = pd.read_csv("/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/Allhomolog_" +
                                  str(stage)+"_FPKM_V2.txt", header=0, index_col=0, sep="\t")
HomoeologExpression.columns = [
    re.sub("[0-9]*D", "", i) for i in HomoeologExpression.columns]
tabixfile = pysam.TabixFile(
    "/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/Genotype/"+stage+".filter2.vcf.gz")
vcfHeader = tabixfile.header[0].split("\t")


def QTLVariant_Map(AllelicCode, ref, Alt):
    # * 提取每个样本的基因型数据
    if AllelicCode == "0/0":
        return ref+ref
    elif AllelicCode == "1/1":
        return Alt+Alt
    elif AllelicCode == "0/1":
        return 'heterozgous'
    elif AllelicCode == "0/1":
        return 'heterozgous'
    else:
        return 'None'


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


if __name__ == "__main__":

    QTLData = pd.read_csv(sys.argv[2], header=0, index_col=None, sep="\t")
    Bias_subGenotype = []

    for value in QTLData.values:
        eGene = value[0]
        subGene = value[8]
        QTLchrom, QTLsite = value[1:3]
        subSite = value[10]
        SNPid = value[3]
        subStand = value[9]
        # * 拼凑基因对名称
        if re.match('Ghir_A', eGene):
            genePair = eGene+"-"+subGene
        else:
            genePair = subGene+"-"+eGene
        #! 分析亚组间的基因型；亚组间没有保守的位点
        if subSite == "-_-_-_-":
            Bias_subGenotype.append(("-", '-'))
        else:
            try:
                chrom1, chrom2, subBase, lastZscore, lastZStand = subSite.split(
                    "_")
                referenceChrom = chrom1+"_"+chrom2
            except ValueError:
                # * 亚组基因出现在scaffold上
                chrom1, subBase, lastZscore, lastZStand = subSite.split("_")
                referenceChrom = chrom1
            geneExpression = HomoeologExpression.loc[genePair].str.split(
                "-", expand=True)
            geneExpression.columns = ['At', 'Dt']
            geneExpression['Bias'] = geneExpression.apply(BiasType, axis=1)

            # * 提取QTL的基因型数据
            QTLVariant = tabixfile.fetch(
                start=QTLsite-1, end=QTLsite, reference=QTLchrom).next().split("\t")
            refer, Alt = QTLVariant[3:5]
            QTLVariant = [QTLVariant_Map(i, refer, Alt)
                          for i in QTLVariant[9:]]
            QTLVariantData = pd.DataFrame(
                QTLVariant, columns=[SNPid], index=vcfHeader[9:]
            )
            # * 提取亚组的基因型数据
            subBase = int(subBase)
            subGenomeSeq = genomeObject.fetch(
                start=subBase-1, end=subBase, reference=referenceChrom)
            QTLVariantData[SNPid+'_sub'] = subGenomeSeq*2
            mergeData = pd.concat([QTLVariantData, geneExpression], axis=1)
            # * 过滤掉QTLVariant杂合的和基因型未知的样本
            mergeData = mergeData[mergeData[SNPid].isin([refer*2, Alt*2])]
            #! 亚组间基因型相同时个体的genetic设置为0，不同设置为1；需要考虑同源基因链方向相反的情况
            # * 获取互补的基因型数据
            reversedDict = {
                'AA': 'TT',
                'TT': 'AA',
                'CC': "GG",
                'GG': 'CC',
                'NN': 'NN'
            }
            if subStand == "sameStand":
                mergeData['SNPCode']=mergeData[SNPid].apply(lambda x:0 if x==subGenomeSeq*2 else 1)
            else:
                #* 基因在相反的链因此对基因型进行互补
                subGenomeSeq2=reversedDict.get(subGenomeSeq*2)
                mergeData['SNPCode']=mergeData[SNPid].apply(lambda x:0 if x==subGenomeSeq2 else 1)
            cor,pvalue=kendalltau(mergeData['SNPCode'],mergeData['Bias'])
            Bias_subGenotype.append((cor,pvalue))

    Bias_subGenotype = pd.DataFrame(Bias_subGenotype, columns=['cor','pvalue'])
    outData = pd.concat([QTLData, Bias_subGenotype], axis=1) 
    outData.to_csv(sys.argv[3],header=True,index=False,sep="\t",na_rep='NA')           
