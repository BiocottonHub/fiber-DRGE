###
 # @Descripttion: 
 # @version: 
 # @Author: zpliu
 # @Date: 2023-05-09 10:33:57
 # @LastEditors: zpliu
 # @LastEditTime: 2023-05-16 21:58:28
 # @@param: 
### 
#TODO 预览任务情况
snakemake --cores 1  -np -s  expression_conditions.smk 

#TODO 运行任务
snakemake --cores 6  -s  expression_conditions.smk 
#? 文章中采样的是perl计算cutoff为0.1时得到的结果,python计算cutoff 0.1时有一点点出入

#TODO 根据表达的基因，分析存在Expression bias的基因数目
snakemake --profile lsf  -j 100 -s  expression_conditions.smk 


#TODO 统计每个cutoff下在至少一个时期间存在表达的基因数目
snakemake --cores 7   -s  expression_conditions.smk 

#? 在6个时期都稳定Bias表达的同源基因数目和BiasN的同源基因数目
getData.ipynb



#TODO 统计expression variation的基因占比 
#* 对于每个cutoff统计unique 表达的基因以及单个时期分析时union数目;
#? 分析这些基因分别有多少比例存在相应的fold changes



#TODO 统计不同阈值条件下，发生Bias表达的基因对 
#* rule: expressionBiasConditions 
#! 提取每个时期中所有的eGene
#? eGene 目录： 03express_gene/eQTLgenes_Allstages_340samples/eQTL_effect_V2





















