'''
Descripttion: 
version: 
Author: zpliu
Date: 2022-01-16 20:44:18
LastEditors: zpliu
LastEditTime: 2022-02-11 21:44:49
@param: 
'''
import pandas as pd
import pysam
import sys 
QTLchromMap = {
    'Ghir_A01': 1, 'Ghir_A02': 2, 'Ghir_A03': 3, 'Ghir_A04': 4, 'Ghir_A05': 5, 'Ghir_A06': 6,
    'Ghir_A07': 7, 'Ghir_A08': 8, 'Ghir_A09': 9, 'Ghir_A10': 10, 'Ghir_A11': 11, 'Ghir_A12': 12,
    'Ghir_A13': 13,
    'Ghir_D01': 14, 'Ghir_D02': 15, 'Ghir_D03': 16, 'Ghir_D04': 17, 'Ghir_D05': 18, 'Ghir_D06': 19,
    'Ghir_D07': 20, 'Ghir_D08': 21, 'Ghir_D09': 22, 'Ghir_D10': 23, 'Ghir_D11': 24, 'Ghir_D12': 25,
    'Ghir_D13': 26
}
QTLchromMap2 = dict(zip(QTLchromMap.values(), QTLchromMap.keys()))

def dap_g_formed(CaVEMaNout, VCFfile, HomoeologFile, Gene_localFile):
    '''get CaVeMaN result 
      result_best: CaVeMaN predict result
      out: new Format
    '''
    tabix_Object = pysam.TabixFile(VCFfile)
    HomoeologGenePair = pd.read_csv(
        HomoeologFile, header=None, index_col=None, sep="\t")
    All_gene_bed = pd.read_csv(Gene_localFile, header=0, index_col=4, sep="\t")
    def itero_SNP_item(x):
        '''
            Ghir_A12G023660 SNP1424193      0.438276
        '''
        eGeneId = x[0]
        At, Dt = HomoeologGenePair.loc[
            HomoeologGenePair[0].isin([eGeneId]
                                      ) | HomoeologGenePair[1].isin([eGeneId])].iloc[0, :]
        # * 区分eGene与另外一个亚组的同源基因
        subGenome = list(set([At, Dt])-set([eGeneId]))[0]
        eGene_stand = All_gene_bed.loc[eGeneId][3]
        eGeneChrom,eGeneTSS,eGeneTTS,eGeneStand=All_gene_bed.loc[eGeneId]
        subGene_chrom, subGeneTSS, subGeneTTS, subGene_stand = All_gene_bed.loc[subGenome]
        stand = ''
        if eGene_stand == subGene_stand:
            stand = 'sameStand'
        else:
            stand = 'DiffStand'
        ChromId = eGeneChrom
        #* 搜索1M范围内的所有SNP
        if eGeneTSS>=2000000:
            starSite=eGeneTSS-2000000
        else:
            starSite=1
        snp_items = tabix_Object.fetch(
            start=starSite-1, end=eGeneTSS+2000000-1, reference=ChromId)
        #*根据SNP id信息获取SNP的坐标信息; 
        SMPlist=[i.split("\t")[0:3] for i in list(snp_items)]
        SNPid=x[1]
        filterData=[i for i in SMPlist if i[2]==SNPid]
        if filterData==[]:
            #* 这里有可能有的SNP超出了2M的范围
            return "{}*{}*{}*{}*{}*{}*{}*{}*{}*{}".format(
                eGeneId, ChromId, "error", SNPid, "-", x[2], subGene_chrom, subGeneTSS, subGenome, stand
            )
        else:
            return "{}*{}*{}*{}*{}*{}*{}*{}*{}*{}".format(
                eGeneId, ChromId, filterData[0][1], SNPid, "-", x[2], subGene_chrom, subGeneTSS, subGenome, stand
            )
    outData = CaVEMaNout.apply(itero_SNP_item, axis=1).str.split("*", expand=True)
    
    outData.columns = ['eGene', 'eChrom', 'CausalSite', 'CausalVar', 'dap_g_pval',
                       'dap_g_pip', 'sub_chrom', 'sub_TSS', 'sub_gene', 'Homoeolog_stand']
    return outData


if __name__ == "__main__":
    '''
    '../CaVEMaN/0DPA/results.best'
    '/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/Genotype/0DPA.filter2.vcf.gz'
    "/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/Homoeolog_geneId.txt"
    "/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/03express_gene/eQTLgenes_Allstages_340samples/fastQTL_peer_20/Gene_information/All_gene_TSS.bed"
    '''
    CaVEMaNout = pd.read_csv(sys.argv[1], header=None, index_col=None, sep="\t")
    out = dap_g_formed(CaVEMaNout, sys.argv[2],sys.argv[3],sys.argv[4])
    out.to_csv(sys.argv[5], header=True, index=False, sep="\t")