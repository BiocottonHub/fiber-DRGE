# -------------------------------------------------------------------------------
# TODO: 因果变异位点发生Bias的基因型；并与对应的亚组进行比较
# * 根据Bias-GWAS效应的方向，就可以判断是reference基因型导致Bias，还是Alt导致Bias ？
# * 用odd的值来确定QTL中发生Bias的基因型
# * reference中发生Bias样本的比例，Alter中发生Bias的比例；
#! 如果reference中Bias比例大于1则，表达reference是导致Bias的基因型
# -------------------------------------------------------------------------------
import pandas as pd
import sys
from scipy.stats import fisher_exact
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
    At, Dt = x
    At, Dt = float(At), float(Dt)
    if At < 0.5 and Dt < 0.5:
        return 'BiasN'
    elif (At >= 0.5 and Dt < 0.5) or (At < 0.5 and Dt > 0.5):
        return 'Bias'
    elif abs((At-Dt)/(At+Dt)) >= 0.33:
        return 'Bias'
    else:
        return 'BiasN'


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
        #! 分析亚组间的基因型
        if subSite == "-_-_-_-" :
            Bias_subGenotype.append(("-", "-", "-", '-'))
        else:
            try:
                chrom1, chrom2, subBase, lastZscore,lastZStand = subSite.split("_")
                referenceChrom = chrom1+"_"+chrom2
            except ValueError:
                # * 亚组基因出现在scaffold上
                chrom1, subBase, lastZscore,lastZStand = subSite.split("_")
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
            #mergeData.to_csv("./SNP1599991.txt",header=True,index=True,sep="\t")

            Bias1 = mergeData.loc[(mergeData['Bias'] == "Bias") & (
                mergeData[SNPid] == refer*2)].shape[0]
            BiasN1 = mergeData.loc[(mergeData['Bias'] == "BiasN") & (
                mergeData[SNPid] == refer*2)].shape[0]
            Bias2 = mergeData.loc[(mergeData['Bias'] == "Bias") & (
                mergeData[SNPid] == Alt*2)].shape[0]
            BiasN2 = mergeData.loc[(mergeData['Bias'] == "Bias") & (
                mergeData[SNPid] == Alt*2)].shape[0]
            # print(mergeData)
            # -----------------------------------
            # * odd >= 1: Bias allele 为reference
            # * odd <1: Bias allele  为Alt
            # -----------------------------------
            odd, pvalue = fisher_exact(
                [
                    (Bias1, BiasN1),
                    (Bias2, BiasN2)
                ],
                alternative="two-sided"
            )
            # ---------------------------------------------------
            # * 提取Bias的基因型，分析其是否与另外一个亚组的基因型相同
            # ? sameGenotype 基因型相同时发生Bias（Bias样本数目设置为负数）
            # ? sameGenotype 基因型相同时发生Bias（Bias样本数目设置为负数）
            # ----------------------------------------------------
            if odd >= 1:
                BiasAllele = refer*2
            else:
                BiasAllele = Alt*2
            # * 判断Bias的基因型是否与另外一个亚组的基因型相同
            if BiasAllele == subGenomeSeq*2 and subStand == "sameStand":
                # * 基因型相同时发生了Bias
                Bias_subGenotype.append(
                    (odd, pvalue, 'sameGenotype', -(Bias1+Bias2))
                )
            elif BiasAllele == subGenomeSeq*2 and subStand == "DiffStand":
                # * 可能正好是位于两条相反的链方向；这是基因型不同时发生了Bias
                Bias_subGenotype.append(
                    (odd, pvalue, 'DifferGenotype1', Bias1+Bias2)
                )
            elif BiasAllele != subGenomeSeq*2 and subStand == "sameStand":
                # * 可能正好是位于两条相同的链；这是基因型不同时发生了Bias
                Bias_subGenotype.append(
                    (odd, pvalue, 'DifferGenotype2', Bias1+Bias2)
                )
            else:
                # BiasAllele!=subGenomeSeq*2  and subStand=="DiffStand"
                # * 虽然位于不同的链的方向，就假定是它们的不同导致的Bias
                Bias_subGenotype.append(
                    (odd, pvalue, 'DifferGenotype3', (Bias1+Bias2))
                )
    Bias_subGenotype = pd.DataFrame(Bias_subGenotype, columns=[
                                    'odd', 'pvalue', 'BiasType', 'BiasSampleCount'])
    outData = pd.concat([QTLData, Bias_subGenotype], axis=1)
    
    GwasPath = '/public/home/yjwang/zpliu/'+stage+'/output/'

    def getQTL_effect(x):
        eGene = x[0]
        subGene = x[8]
        SNPid = x[3]
        if re.match("Ghir_A", eGene):
            GWASData = pd.read_csv(
                "{}{}-{}_out".format(GwasPath, eGene, subGene), header=0, index_col=None, sep="\t")
        else:
            GWASData = pd.read_csv(
                "{}{}-{}_out".format(GwasPath, subGene, eGene), header=0, index_col=None, sep="\t")
        SNPData = GWASData.loc[GWASData['SNP'] == SNPid]
        if SNPData.empty:
            # * 因果变异没有跑GWAS
            return ("-*-")
        else:
            return "{}*{}".format(SNPData.iloc[0, 9], SNPData.iloc[0, 10])

    QTLeffect = outData.apply(getQTL_effect, axis=1)
    tmpData = QTLeffect.str.split("*", expand=True)
    tmpData.columns = ['beta', 'se']
    outData = pd.concat([outData, tmpData], axis=1)
    outData.to_csv(sys.argv[3], header=True, index=False, sep="\t")
