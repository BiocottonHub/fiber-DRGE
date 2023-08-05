###
# @Descripttion:
# @version:
# @Author: zpliu
# @Date: 2022-03-20 18:07:58
 # @LastEditors: zpliu
 # @LastEditTime: 2022-03-29 17:02:51
# @@param:
###
#TODO
#* 提取SNP id
#? 1. 与lead SNP紧密连锁的SNP
#? 2. cis调控区域的SNP
#? 3. trans调控区域的SNP
#! 因此对于每个eGene只需要构造三个文件
stage=0DPA
plinkFile=/data/cotton/MaojunWang/WMJ_fiberFullPopulationRNAseq/fastlmm_eQTLs_V2/${stage}/${stage}_FastGENE
cis_eGeneFile=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/ExpressionBias_eQTL/Bias_Category/coloc_V3/Genetic_Parition/${stage}_cis_eGene
gcta=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/software/gcta_1.94.0beta/gcta64
rawDataFile=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/BayesS/${stage}/
phenotypeFile=${stage}_expression_normal.txt
covariat=${stage}_peer_covariat.txt
for eGene in $(cat ${cis_eGeneFile} | grep None -v | grep Gh | cut -f4); do
    grep ${eGene} ${cis_eGeneFile} | cut -f1-4 >${eGene}_cis_region.txt
    phenotypeId=$(awk '$4=="'${eGene}'"{print $7}' ${cis_eGeneFile})
    #* 有的eGene可能没有鉴定到cis-eQTL
    leadSNP=$(awk '$4=="'${eGene}'"{print $5}' ${cis_eGeneFile})
    plink --bfile ${plinkFile} --r2 --ld-snp ${leadSNP} --ld-window-kb 50 --ld-window 49999 --ld-window-r2 0.6 --out ${eGene}
    awk 'NR>1{print $6}' ${eGene}.ld >${eGene}_SNPCluster.txt
    bsub -q normal -n 1 -e test.err "
    plink --bfile ${plinkFile}  --extract range ${eGene}_cis_region.txt --write-snplist --out  ${eGene}
    #* 获取出SNP cluster外的其他cis SNPs
    cat ${eGene}_SNPCluster.txt ${eGene}.snplist|sort |uniq -d |cat - ${eGene}.snplist|sort|uniq -u >${eGene}_noSNPCluster.txt
    #* 制作grem文件
    ${gcta} --bfile ${plinkFile} --autosome-num  26 --extract ${eGene}_SNPCluster.txt --make-grm --out ${eGene}_SNPCluster
    ${gcta} --bfile ${plinkFile} --autosome-num  26 --extract ${eGene}_noSNPCluster.txt --make-grm --out ${eGene}_noSNPCluster
    ${gcta} --bfile ${plinkFile} --autosome-num  26 --exclude ${eGene}.snplist --make-grm --out ${eGene}_trans 
    #* 多个grem文件一同跑gcta
    echo "$(pwd)/${eGene}_SNPCluster" >${eGene}_mgrm.txt 
    echo "$(pwd)/${eGene}_noSNPCluster" >>${eGene}_mgrm.txt
    echo "$(pwd)/${eGene}_trans" >>${eGene}_mgrm.txt
    ${gcta}  --mgrm ${eGene}_mgrm.txt --pheno ${rawDataFile}/${phenotypeFile} --reml --qcovar ${rawDataFile}/${covariat} --thread-num 1  --mpheno  ${phenotypeId} --out ${eGene}_result
    rm ${eGene}_mgrm.txt ${eGene}_*.grm.bin  ${eGene}_*.grm.id  ${eGene}_*.grm.N.bin -rf 
    rm ${eGene}.snplist  ${eGene}_noSNPCluster.txt ${eGene}_SNPCluster.txt  ${eGene}*.log -rf 
    #rm ${eGene}.nosex  ${eGene}.ld ${eGene}_cis_region.txt -rf 
    rm ${eGene}.ld ${eGene}.log ${eGene}.nosex  ${eGene}_cis_region.txt  -rf 
    "
done

#-----------------------------------------------------
#* cis区域没有SNP的eGene
#-----------------------------------------------------
for stage in 0DPA 4DPA 8DPA 12DPA 16DPA 20DPA; do
    cd ${stage}
    plinkFile=/data/cotton/MaojunWang/WMJ_fiberFullPopulationRNAseq/fastlmm_eQTLs_V2/${stage}/${stage}_FastGENE
    cis_eGeneFile=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/ExpressionBias_eQTL/Bias_Category/coloc_V3/Genetic_Parition/${stage}_cis_eGene
    gcta=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/software/gcta_1.94.0beta/gcta64
    rawDataFile=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/BayesS/${stage}/
    phenotypeFile=${stage}_expression_normal.txt
    covariat=${stage}_peer_covariat.txt
    for eGene in $(cat ${cis_eGeneFile} | grep None | grep Gh | cut -f4); do
        grep ${eGene} ${cis_eGeneFile} | cut -f1-4 >${eGene}_cis_region.txt
        phenotypeId=$(awk '$4=="'${eGene}'"{print $7}' ${cis_eGeneFile})
        #* 有的eGene可能没有鉴定到cis-eQTL
        bsub -q normal -n 1 -e test.err "
        plink --bfile ${plinkFile}  --extract range ${eGene}_cis_region.txt --write-snplist --out  ${eGene}
        #* 制作grem文件
        ${gcta} --bfile ${plinkFile} --autosome-num  26 --extract ${eGene}.snplist --make-grm --out ${eGene}_noSNPCluster
        ${gcta} --bfile ${plinkFile} --autosome-num  26 --exclude ${eGene}.snplist --make-grm --out ${eGene}_trans 
        #* 多个grem文件一同跑gcta
        echo "$(pwd)/${eGene}_noSNPCluster" >>${eGene}_mgrm.txt
        echo "$(pwd)/${eGene}_trans" >>${eGene}_mgrm.txt
        ${gcta}  --mgrm ${eGene}_mgrm.txt --pheno ${rawDataFile}/${phenotypeFile} --reml --qcovar ${rawDataFile}/${covariat} --thread-num 1  --mpheno  ${phenotypeId} --out ${eGene}_result
        rm ${eGene}_mgrm.txt ${eGene}_*.grm.bin  ${eGene}_*.grm.id  ${eGene}_*.grm.N.bin -rf 
        rm ${eGene}.snplist  ${eGene}_noSNPCluster.txt   ${eGene}*.log -rf 
        #rm ${eGene}.nosex   ${eGene}_cis_region.txt -rf 
        rm  ${eGene}.log ${eGene}.nosex  ${eGene}_cis_region.txt  -rf 
        "
    done
    cd ..
done

#----------------------------------------------------------
#! lead SNP cluster 改成 ld friend
#----------------------------------------------------------
gcta64 --bfile /data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/Genotype/Trait_376
--ld All_leadSNP.txt --ld-wind 100 --ld-sig 0.05 --out ld_friend

for stage in 4DPA 8DPA 12DPA 16DPA 20DPA; do
    cd ${stage}
    plinkFile=/data/cotton/MaojunWang/WMJ_fiberFullPopulationRNAseq/fastlmm_eQTLs_V2/${stage}/${stage}_FastGENE
    cis_eGeneFile=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/ExpressionBias_eQTL/Bias_Category/coloc_V3/Genetic_Parition/${stage}_cis_eGene
    gcta=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/software/gcta_1.94.0beta/gcta64
    rawDataFile=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/BayesS/${stage}/
    phenotypeFile=${stage}_expression_normal.txt
    ldFriend=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/ExpressionBias_eQTL/Bias_Category/coloc_V3/Genetic_Parition/test/ld_friend.snp.ld2
    extractSNPCluster=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/ExpressionBias_eQTL/Bias_Category/coloc_V3/Genetic_Parition/test/extract_SNP.sh
    covariat=${stage}_peer_covariat.txt
    for eGene in $(cat ${cis_eGeneFile} | grep None -v | grep Gh | cut -f4); do
        grep ${eGene} ${cis_eGeneFile} | cut -f1-4 >${eGene}_cis_region.txt
        phenotypeId=$(awk '$4=="'${eGene}'"{print $7}' ${cis_eGeneFile})
        #* 有的eGene可能没有鉴定到cis-eQTL
        leadSNP=$(awk '$4=="'${eGene}'"{print $5}' ${cis_eGeneFile})
        bsub -q q2680v2 -n 1 -e test.err "
        #* 获取SNP friend
        bash ${extractSNPCluster} ${leadSNP} ${ldFriend} ${eGene}_SNPCluster.txt
        plink --bfile ${plinkFile}  --extract range ${eGene}_cis_region.txt --write-snplist --out  ${eGene}
        #* 获取出SNP cluster外的其他cis SNPs
        cat ${eGene}_SNPCluster.txt ${eGene}.snplist|sort |uniq -d |cat - ${eGene}.snplist|sort|uniq -u >${eGene}_noSNPCluster.txt
        #* 制作grem文件
        ${gcta} --bfile ${plinkFile} --autosome-num  26 --extract ${eGene}_SNPCluster.txt --make-grm --out ${eGene}_SNPCluster
        ${gcta} --bfile ${plinkFile} --autosome-num  26 --extract ${eGene}_noSNPCluster.txt --make-grm --out ${eGene}_noSNPCluster
        ${gcta} --bfile ${plinkFile} --autosome-num  26 --exclude ${eGene}.snplist --make-grm --out ${eGene}_trans 
        #* 多个grem文件一同跑gcta
        echo "$(pwd)/${eGene}_SNPCluster" >${eGene}_mgrm.txt 
        echo "$(pwd)/${eGene}_noSNPCluster" >>${eGene}_mgrm.txt
        echo "$(pwd)/${eGene}_trans" >>${eGene}_mgrm.txt
        ${gcta}  --mgrm ${eGene}_mgrm.txt --pheno ${rawDataFile}/${phenotypeFile} --reml --qcovar ${rawDataFile}/${covariat} --thread-num 1  --mpheno  ${phenotypeId} --out ${eGene}_result
        rm ${eGene}_mgrm.txt ${eGene}_*.grm.bin  ${eGene}_*.grm.id  ${eGene}_*.grm.N.bin -rf 
        rm ${eGene}.snplist  ${eGene}_noSNPCluster.txt ${eGene}_SNPCluster.txt  ${eGene}*.log -rf 
        #rm ${eGene}.nosex  ${eGene}.ld ${eGene}_cis_region.txt -rf 
        rm ${eGene}.ld ${eGene}.log ${eGene}.nosex  ${eGene}_cis_region.txt  -rf 
        "
    done
    cd ..
done
#----------------------------------------------------
#* 基于ld score对lead SNP进行分层
#----------------------------------------------------
for stage in 4DPA 8DPA 12DPA 16DPA 20DPA; do
    cd ${stage}
    plinkFile=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/Genotype/Trait_376
    cis_eGeneFile=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/ExpressionBias_eQTL/Bias_Category/coloc_V3/Genetic_Parition/${stage}_cis_eGene
    gcta=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/software/gcta_1.94.0beta/gcta64
    rawDataFile=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/BayesS/${stage}/
    phenotypeFile=${stage}_expression_normal.txt
    ldFriend=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/ExpressionBias_eQTL/Bias_Category/coloc_V3/Genetic_Parition/test/ld_friend.snp.ld2
    extractSNPCluster=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/ExpressionBias_eQTL/Bias_Category/coloc_V3/Genetic_Parition/test/extract_SNP.sh
    covariat=${stage}_peer_covariat.txt
    LeadSNPData=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/ExpressionBias_eQTL/Bias_Category/coloc_V3/Genetic_Parition/${stage}_leadSNP_conserved.txt
    stratifySNP=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/ExpressionBias_eQTL/Bias_Category/coloc_V3/Genetic_Parition/stratify_SNP.R
    extractBlockSNP=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/ExpressionBias_eQTL/Bias_Category/coloc_V3/Genetic_Parition/extract_geneticSNP.py
    #* 提取LD Block内的SNP
    for eGene in $(cat ${LeadSNPData} | grep None -v | cut -f2); do
        phenotypeId=$(awk '$4=="'${eGene}'"{print $7}' ${cis_eGeneFile})
        leadSNP=$(awk '$2=="'${eGene}'"{print $3}' ${LeadSNPData})
        leadSNPGenetic=$(awk '$2=="'${eGene}'"{print $4}' ${LeadSNPData})
        bsub -q normal -n 1 -e test.err "
        bash ${extractSNPCluster} ${leadSNP} ${ldFriend} ${eGene}_SNPCluster.txt
        #* 制作plink 文件
        plink --bfile ${plinkFile} --extract  ${eGene}_SNPCluster.txt --make-bed --out ${eGene}
        ${gcta} --bfile ${eGene}  --autosome-num 26 --ld-score-region 100 --out ${eGene}
        #* 进行分层
        Rscript  ${stratifySNP}  ${eGene}
        #* 判断文件是否存在
        if [ -f "${eGene}_snp_group1.txt" ]; then
            ${gcta} --bfile ${eGene} --autosome-num 26 --extract ${eGene}_snp_group1.txt --make-grm --out ${eGene}_snp_group1
            echo $(pwd)/${eGene}_snp_group1 >>${eGene}_mgrm.txt 
        fi
        if [ -f "${eGene}_snp_group2.txt" ]; then
            ${gcta} --bfile ${eGene} --autosome-num 26 --extract ${eGene}_snp_group2.txt --make-grm --out ${eGene}_snp_group2
            echo $(pwd)/${eGene}_snp_group2 >>${eGene}_mgrm.txt
        fi
        if [ -f "${eGene}_snp_group3.txt" ]; then
            ${gcta} --bfile ${eGene} --autosome-num 26 --extract ${eGene}_snp_group3.txt --make-grm --out ${eGene}_snp_group3
            echo $(pwd)/${eGene}_snp_group3 >>${eGene}_mgrm.txt
        fi
        if [ -f "${eGene}_snp_group4.txt" ]; then
            ${gcta} --bfile ${eGene} --autosome-num 26 --extract ${eGene}_snp_group4.txt --make-grm --out ${eGene}_snp_group4
            echo $(pwd)/${eGene}_snp_group4 >>${eGene}_mgrm.txt
        fi
        ${gcta} --mgrm ${eGene}_mgrm.txt --pheno ${rawDataFile}/${phenotypeFile} --reml --qcovar ${rawDataFile}/${covariat} --thread-num 1 --mpheno ${phenotypeId} --out ${eGene}_result
        #* 提取stratify 后的有效应的SNP集合,效应的定义为V(G)/V(P)>=0.05
        python ${extractBlockSNP} ${eGene}  ${leadSNPGenetic} ${leadSNP}
        rm ${eGene}_snp_group*  ${eGene}_SNPCluster.txt  ${eGene}.score.ld  -rf 
        rm ${eGene}_result*  ${eGene}_mgrm.txt ${eGene}.nosex  -rf 
        rm ${eGene}.log ${eGene}.fam ${eGene}.bim ${eGene}.bed -rf 
        "
    done
    cd ../
done

#-----------------------------------------
#* 从1000个eGene中随机抽取随机的SNP
#! 每个eGene抽取100个随机的SNP，总共抽取500个eGene
#-----------------------------------------
ciseGeneFile='/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/ExpressionBias_eQTL/Bias_Category/coloc_V3/Genetic_Parition/random_SNPs/random_eGene.txt'
plinkFile=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/Genotype/Trait_376
for eGene in $(cat ${ciseGeneFile} | grep Gh | cut -f4); do
    #* 临时cis region文件
    awk '$4=="'${eGene}'"{print $0}' ${ciseGeneFile} | cut -f1-4 | sort | uniq >${eGene}_cis_region.txt
    plink --bfile ${plinkFile} --extract range ${eGene}_cis_region.txt --write-snplist --out ${eGene}
    shuf -n 200 ${eGene}.snplist | awk '{print "'${eGene}'""\t"$0"\t0"}' >>random_snp.txt
    rm ${eGene}* -rf
    #* 提取200个SNP
done

#---------------------------------------------
#* 使用LDAK进行遗传力估计
#---------------------------------------------
#* 估计整个基因组中所有SNP的权重
plinkFile=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/Genotype/Trait_376
ldak=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/software/ldak/ldak5.2.linux
SNP_weight=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/ExpressionBias_eQTL/Bias_Category/coloc_V3/Genetic_Parition/LDAK/weights.all2
#${ldak} --thin thin --bfile ${plinkFile} --window-prune .98 --window-kb 100

for stage in 0DPA 4DPA 8DPA 12DPA 16DPA 20DPA; do
    cd $stage
    cis_eGeneFile=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/ExpressionBias_eQTL/Bias_Category/coloc_V3/Genetic_Parition/${stage}_cis_eGene
    # cis_eGeneFile=/public/home/zpliu/LZP_Fiber_GWAS/GCTA/${stage}_tmpData
    gcta=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/software/gcta_1.94.0beta/gcta64
    rawDataFile=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/BayesS/${stage}/
    phenotypeFile=${stage}_expression_normal.txt
    ldFriend=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/ExpressionBias_eQTL/Bias_Category/coloc_V3/Genetic_Parition/test/ld_friend.snp.ld2
    extractSNPCluster=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/ExpressionBias_eQTL/Bias_Category/coloc_V3/Genetic_Parition/test/extract_SNP.sh
    covariat=${stage}_peer_covariat.txt
    LeadSNPData=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/ExpressionBias_eQTL/Bias_Category/coloc_V3/Genetic_Parition/${stage}_leadSNP_conserved.txt
    stratifySNP=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/ExpressionBias_eQTL/Bias_Category/coloc_V3/Genetic_Parition/stratify_SNP.R
    extractBlockSNP=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/ExpressionBias_eQTL/Bias_Category/coloc_V3/Genetic_Parition/extract_geneticSNP.py
    for eGene in $(cat ${cis_eGeneFile} | grep None -v | grep Gh | cut -f4); do
        grep ${eGene} ${cis_eGeneFile} | cut -f1-4 >${eGene}_cis_region.txt
        phenotypeId=$(awk '$4=="'${eGene}'"{print $7}' ${cis_eGeneFile})
        #* 有的eGene可能没有鉴定到cis-eQTL
        leadSNP=$(awk '$4=="'${eGene}'"{print $5}' ${cis_eGeneFile})
        bsub -q normal -n 1 -e test.err "
            #* 获取SNP friend
            bash ${extractSNPCluster} ${leadSNP} ${ldFriend} ${eGene}_SNPCluster.txt
            plink --bfile ${plinkFile}  --extract range ${eGene}_cis_region.txt --write-snplist --out  ${eGene}
            cat ${eGene}_SNPCluster.txt ${eGene}.snplist|sort |uniq -d |cat - ${eGene}.snplist|sort|uniq -u >${eGene}_noSNPCluster.txt
            ${ldak} --calc-kins-direct  ${eGene}_noSNPCluster --bfile ${plinkFile}  --weights ${SNP_weight} --power -.25   --extract ${eGene}_noSNPCluster.txt
            ${ldak} --calc-kins-direct  ${eGene}_SNPCluster --bfile ${plinkFile}  --weights ${SNP_weight} --power -.25   --extract ${eGene}_SNPCluster.txt
            ${ldak} --calc-kins-direct  ${eGene}_trans --bfile ${plinkFile}  --weights  ${SNP_weight} --power -.25   --exclude  ${eGene}.snplist
            echo  ${eGene}_SNPCluster >${eGene}_mlist.txt
            echo  ${eGene}_noSNPCluster >>${eGene}_mlist.txt
            echo  ${eGene}_trans >>${eGene}_mlist.txt
            ${ldak} --reml ${eGene} --pheno ${rawDataFile}/${phenotypeFile} --mpheno ${phenotypeId} --mgrm ${eGene}_mlist.txt --covar ${rawDataFile}/${stage}_peer_covariat.txt --constrain YES --he-starts NO
            rm ${eGene}_noSNPCluster* ${eGene}_SNPCluster*  ${eGene}_trans* -rf 
            rm ${eGene}.vars ${eGene}.snplist ${eGene}.share ${eGene}.progress -rf 
            rm ${eGene}_mlist.txt ${eGene}_cis_region.txt ${eGene}.coeff ${eGene}.cross -rf 
            rm ${eGene}.indi.blp ${eGene}.indi.res ${eGene}.nosex ${eGene}.log -rf  
            "
    done
    cd ..
done

#--------------------------------------------------------
#* 没有cis-QTL的eGene
#--------------------------------------------------------
for stage in 0DPA 4DPA 8DPA 12DPA 16DPA 20DPA; do
    cd $stage
    cis_eGeneFile=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/ExpressionBias_eQTL/Bias_Category/coloc_V3/Genetic_Parition/${stage}_cis_eGene
    # cis_eGeneFile=/public/home/zpliu/LZP_Fiber_GWAS/GCTA/${stage}_tmpData
    gcta=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/software/gcta_1.94.0beta/gcta64
    rawDataFile=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/BayesS/${stage}/
    phenotypeFile=${stage}_expression_normal.txt
    ldFriend=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/ExpressionBias_eQTL/Bias_Category/coloc_V3/Genetic_Parition/test/ld_friend.snp.ld2
    extractSNPCluster=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/ExpressionBias_eQTL/Bias_Category/coloc_V3/Genetic_Parition/test/extract_SNP.sh
    covariat=${stage}_peer_covariat.txt
    LeadSNPData=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/ExpressionBias_eQTL/Bias_Category/coloc_V3/Genetic_Parition/${stage}_leadSNP_conserved.txt
    stratifySNP=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/ExpressionBias_eQTL/Bias_Category/coloc_V3/Genetic_Parition/stratify_SNP.R
    extractBlockSNP=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/ExpressionBias_eQTL/Bias_Category/coloc_V3/Genetic_Parition/extract_geneticSNP.py
    for eGene in $(cat ${cis_eGeneFile} | grep None | grep Gh | cut -f4); do
        grep ${eGene} ${cis_eGeneFile} | cut -f1-4 >${eGene}_cis_region.txt
        phenotypeId=$(awk '$4=="'${eGene}'"{print $7}' ${cis_eGeneFile})
        #* 有的eGene可能没有鉴定到cis-eQTL
        bsub -q normal -n 1 -e test.err "
            plink --bfile ${plinkFile}  --extract range ${eGene}_cis_region.txt --write-snplist --out  ${eGene}
            #* 制作grem文件
            ${ldak} --calc-kins-direct  ${eGene}_noSNPCluster --bfile ${plinkFile}  --weights ${SNP_weight} --power -.25   --extract ${eGene}.snplist
            ${ldak} --calc-kins-direct  ${eGene}_trans --bfile ${plinkFile}  --weights ${SNP_weight} --power -.25   --exclude  ${eGene}.snplist
            #* 多个grem文件一同跑gcta
            echo  ${eGene}_noSNPCluster >>${eGene}_mlist.txt
            echo  ${eGene}_trans >>${eGene}_mlist.txt
            ${ldak} --reml ${eGene} --pheno ${rawDataFile}/${phenotypeFile} --mpheno ${phenotypeId} --mgrm ${eGene}_mlist.txt --covar ${rawDataFile}/${stage}_peer_covariat.txt --constrain YES --he-starts NO
            rm ${eGene}_noSNPCluster*  ${eGene}_trans* -rf 
            rm ${eGene}.vars ${eGene}.snplist ${eGene}.share ${eGene}.progress -rf 
            rm ${eGene}_mlist.txt ${eGene}_cis_region.txt ${eGene}.coeff ${eGene}.cross -rf 
            rm ${eGene}.indi.blp ${eGene}.indi.res ${eGene}.nosex ${eGene}.log -rf  
            "
    done
    cd ..
done 
