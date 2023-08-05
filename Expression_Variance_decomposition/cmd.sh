###
 # @Descripttion: 
 # @version: 
 # @Author: zpliu
 # @Date: 2021-10-30 21:33:02
 # @LastEditors: zpliu
 # @LastEditTime: 2021-12-10 22:37:56
 # @@param: 
### 
#--------------------------------------------------
#TODO： 对基因表达的变异进行结构
#* ref:https://pubmed.ncbi.nlm.nih.gov/33930287/
#? 探讨时间和个体间随机因素对基因表达的解释
#* 随机变量是对组内样本进行方差分析，计算各个随机变量的解释率
#---------------------------------------------------

#* 筛选需要进行变异分析的基因
#! 基因在至少一个时期的80%的样本中，基因表达水平大于5.
geneID=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/03express_gene/eQTLgenes_Allstages_340samples/fastQTL_peer_20/Gene_information/All_gene_TSS.bed


#* 对每个基因进行100次Bootstrap
#? 每次进行有放回取样，计算一次随机变量的方差解释率
for i in $(ls splitFile); do
    bsub -q normal -n 1 -e test.err -o test.out "python script/bootstrap_lm.py splitFile/${i} splitFile/${i}_out "
done
cat splitFile/*out |awk 'NR==1{print $0}NR>1&&$1~/Ghir/{print $0}' >Gene_expression_decomposition.txt



