###
# @Descripttion:
# @version:
# @Author: zpliu
# @Date: 2023-02-03 21:22:49
 # @LastEditors: zpliu
 # @LastEditTime: 2023-08-05 15:15:57
# @@param:
###
geneBedRegion='/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/03express_gene/eQTLgenes_Allstages_340samples/fastQTL_peer_20/Gene_information/All_gene_bed.txt'

mosdepthPath='/cotton/Liuzhenping/Pan-genome/software/mosdepth'
BamPath='/data/cotton/MaojunWang/WMJ_fiberFullPopulationRNAseq/MappingFPKM/12DPA/Bamfiles/'
#* 检测read覆盖度,
for sample in $(ls ${BamPath} | grep bai | sed 's/_.*//g'); do
    bsub -q q2680v2 -n 1 -J ${sample} -e test.err -o test.out "
     ${mosdepthPath} -t 1 -b ${geneBedRegion} -n -x ${sample} ${BamPath}/${sample}_mapGhir_srt.bam
    "
done

#* 对同源基因进行二项分布检验
for sample in $(ls splitData); do
    bsub -q normal  -n 1 -J 8DPA_${sample} -e test.err -o test.out "
    python /public/home/zpliu/LZP_Fiber_GWAS/geneReadsCount/geneBiasTypeCount.py ./splitData/${sample} ./splitData/${sample}_contri ./splitData/${sample}_sample
    "
done


#TODO 统计mapping到基因外显子区域的read数目
# HTSeq: v2.0.2
BamPath='/data/cotton/MaojunWang/WMJ_fiberFullPopulationRNAseq/MappingFPKM/0DPA/Bamfiles/'
gftFile='/public/home/zpliu/work/Alternative_review/genomeData/Ghirsutum_genome_HAU_v1.1/Ghirsutum_gene_model.gtf'
for sample in $(ls ${BamPath} | grep bai | sed 's/_.*//g'); do
    bsub -q normal  -n 2 -J 4DPA_${sample} -e test.err -o test.out "
     python -m HTSeq.scripts.count -f bam -r pos -s no  -t exon \
                -i gene_id -m union --nonunique none \
                -n 2 -c ${sample}_htseq_readCout.tsv \
                ${BamPath}/${sample}_mapGhir_srt.bam  ${gftFile}
    "
done

#TODO: 将在5%样本中Bias的基因进行




