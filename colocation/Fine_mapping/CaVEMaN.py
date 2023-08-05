'''
Descripttion: 
version: 
Author: zpliu
Date: 2022-01-14 14:03:49
LastEditors: zpliu
LastEditTime: 2022-01-14 14:39:14
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


def CaVEMaN_formed(CaVeMaNFile, VCFfile, HomoeologFile, Gene_localFile):
    '''get CaVeMaN result 
      result_best: CaVeMaN predict result
      out: new Format
    '''
    CaVEMaNout = pd.read_csv(CaVeMaNFile, header=0, index_col=None, sep="\t")
    tabix_Object = pysam.TabixFile(VCFfile)
    HomoeologGenePair = pd.read_csv(
        HomoeologFile, header=None, index_col=None, sep="\t")
    All_gene_bed = pd.read_csv(Gene_localFile, header=0, index_col=4, sep="\t")

    def CavEMA_item(x):
        '''
        Ghir_A01G000740_1_572110_C_G	1	572110	C	G	-0.597719	3.292840e-37	0.157314	0.418123
        '''
        eGeneId = "{}_{}".format(x[0].split("_")[0], x[0].split("_")[1])
        At, Dt = HomoeologGenePair.loc[
            HomoeologGenePair[0].isin([eGeneId]
                                      ) | HomoeologGenePair[1].isin([eGeneId])].iloc[0, :]
        # * 区分eGene与另外一个亚组的同源基因
        subGenome = list(set([At, Dt])-set([eGeneId]))[0]
        eGene_stand = All_gene_bed.loc[eGeneId][3]
        subGene_chrom, subGeneTSS, subGeneTTS, subGene_stand = All_gene_bed.loc[subGenome]
        stand = ''
        if eGene_stand == subGene_stand:
            stand = 'sameStand'
        else:
            stand = 'DiffStand'
        ChromId = QTLchromMap2.get(x[1])
        snp_items = tabix_Object.fetch(
            start=x[2]-1, end=x[2], reference=ChromId)
        SNPid = list(snp_items)[0].split("\t")[2]
        return "{}*{}*{}*{}*{}*{}*{}*{}*{}*{}".format(
            eGeneId, ChromId, x[2], SNPid, x[6], x[8], subGene_chrom, subGeneTSS, subGenome, stand
        )
    outData = CaVEMaNout.apply(CavEMA_item, axis=1).str.split("*", expand=True)
    outData.columns = ['eGene', 'eChrom', 'CausalSite', 'CausalVar', 'CaVEMaN_pval',
                       'CaVEMaN_pip', 'sub_chrom', 'sub_TSS', 'sub_gene', 'Homoeolog_stand']
    return outData


if __name__ == "__main__":
    '''
    '../CaVEMaN/0DPA/results.best'
    '/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/Genotype/0DPA.filter2.vcf.gz'
    "/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/Homoeolog_geneId.txt"
    "/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/03express_gene/eQTLgenes_Allstages_340samples/fastQTL_peer_20/Gene_information/All_gene_TSS.bed"
    '''
    out = CaVEMaN_formed(sys.argv[1], sys.argv[2],sys.argv[3],sys.argv[4])
    out.to_csv(sys.argv[5], header=True, index=False, sep="\t")
