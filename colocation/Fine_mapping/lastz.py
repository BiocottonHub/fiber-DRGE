'''
Descripttion: 
version: 
Author: zpliu
Date: 2022-01-14 15:52:12
LastEditors: zpliu
LastEditTime: 2022-09-29 09:10:05
@param: 
'''
# *
import os
import pysam
from tempfile import NamedTemporaryFile
from pyliftover import LiftOver
from extract_seq import extract_sequence
lastZ = '/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/software/lastz-1.04.03/lastz-distrib/bin/lastz'
axtChain = '/public/home/software/opt/bio/software/ucsc_kentUtils/v389/bin/axtChain'
chainSwap = '/public/home/software/opt/bio/software/ucsc_kentUtils/v389/bin/chainSwap'
chainSort = '/public/home/software/opt/bio/software/ucsc_kentUtils/v389/bin/chainSort'


def lastz(seqrchSeq, targetSeq, searchId, targetSeqReal_start, targetChrom,
          eQTL_seqFile, eGene_seqFile, axtFile, chainFile, Chain_changeOrder, standSame_or_not):
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
        "{} {} {} K=3000 L=3000 H=2000 Y=5000 E=55 T=2 O=600 --filter=identity:80 --filter=coverage:85  --verbosity=10  --format=axt  >{}".format(
            lastZ, eGene_seqFile.name, eQTL_seqFile.name, axtFile.name
        )
    )
    #* 宽松的标准
    # os.system(
    #     "{} {} {} K=3000 L=3000 H=2000 Y=5000 E=55 T=2 O=600  --verbosity=10  --format=axt  >{}".format(
    #         lastZ, eGene_seqFile.name, eQTL_seqFile.name, axtFile.name
    #     )
    # )
    # * axt转换为chain
    eQTL_seqFile.seek(0)
    eGene_seqFile.seek(0)
    axtFile.seek(0)
    os.system(
        "{} -minScore=3000 -linearGap=medium {} -faT -faQ {} {} {}".format(
            axtChain, axtFile.name, eGene_seqFile.name, eQTL_seqFile.name, chainFile.name
        )
    )
    # * 交换方向
    chainFile.seek(0)
    os.system(
        "{} {} stdout | {} stdin  {}".format(
            chainSwap, chainFile.name, chainSort, Chain_changeOrder.name
        )
    )
    # ----------------------------------------------------
    # * 进行搜索
    # ----------------------------------------------------
    Chain_changeOrder.seek(0)
    lo = LiftOver(Chain_changeOrder.name)
    result = lo.convert_coordinate(searchId, 50)
    eQTL_seqFile.seek(0)
    # with open("qtl.fa",'w') as File:
    #     File.write(eQTL_seqFile.read())
    # eGene_seqFile.seek(0)
    # with open("gene.fa",'w') as File:
    #     File.write(eGene_seqFile.read())
    # print(searchId)
    # * 将得分最高的mapping到真实的坐标中
    # axtFile.seek(0)
    # print(axtFile.read())
    # print(targetSeqReal_start,result)
    # print(result,targetSeqReal_start)
    if result == None or result == []:
        # * 整段序列都没有mapping到
        # * 当前查询的位置没有mapping到
        return ('-', '-', '-', "-", 0)
    else:
        mappingRegionCount = len(result)
    if mappingRegionCount == 1:
        # * 只有一个保守片段
        site = result[0][1]+targetSeqReal_start
        Score = result[0][3]
        stand = result[0][2]
        return (targetChrom, site, Score, stand, mappingRegionCount)
    elif standSame_or_not == "DiffStand":
        filterData = [i for i in result if i[2] == "-"]
        if filterData == None or filterData == []:
            # * 整段序列都没有mapping到
            # * 当前查询的位置没有mapping到
            return ('-', '-', '-', "-", mappingRegionCount)
        else:
            # * 比对的时候是不同方向的位点; 同时获取得分比较高 mapping,index从0开始
            site = filterData[0][1]+targetSeqReal_start
            Score = filterData[0][3]
            stand = filterData[0][2]
            return (targetChrom, site, Score, stand, mappingRegionCount)
    else:
        filterData = [i for i in result if i[2] == "+"]
        if filterData == None or filterData == []:
            # * 整段序列都没有mapping到
            # * 当前查询的位置没有mapping到
            return ('-', '-', '-', "-", mappingRegionCount)
        else:
            # * 比对的时候是同方向的位点; 同时获取得分比较高 mapping,index从0开始

            site = filterData[0][1]+targetSeqReal_start
            Score = filterData[0][3]
            stand = filterData[0][2]
            return (targetChrom, site, Score, stand, mappingRegionCount)


if __name__ == "__main__":
    genomeFile = '/data/cotton/MaojunWang/WMJ_fiberFullPopulationRNAseq/MappingFPKM/Ghir_Genome_Index/Ghirsutum_genome.fasta'
    genomeObject = pysam.FastaFile(genomeFile)
    eQTL_seqFile = NamedTemporaryFile(mode='w+', encoding='utf-8')
    eGene_seqFile = NamedTemporaryFile(mode='w+', encoding='utf-8')
    axtFile = NamedTemporaryFile(mode='w+', encoding='utf-8')
    chainFile = NamedTemporaryFile(mode='w+', encoding='utf-8')
    Chain_changeOrder = NamedTemporaryFile(mode='w+', encoding='utf-8')
    # * 提取序列信息
    eChrom, eStart, eEnd, eQTLid = ('Ghir_A01', 571910, 572310, 'SNP172')
    subChrom, subStart, subEnd, subId = (
        'Ghir_D01', 1, 1572110, 'Ghir_D01G000710')
    searchSeq = extract_sequence(eChrom, eStart, eEnd, eQTLid, genomeObject)
    TargetSeq = extract_sequence(
        subChrom, subStart, subEnd, subId, genomeObject)
    data = lastz(
        searchSeq, TargetSeq, eQTLid, subStart, subChrom, eQTL_seqFile,
        eGene_seqFile, axtFile, chainFile, Chain_changeOrder
    )
    # print(data)
    eQTL_seqFile.close()
    eGene_seqFile.close()
    axtFile.close()
    chainFile.close()
    Chain_changeOrder.close()
