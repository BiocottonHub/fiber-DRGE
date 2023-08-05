###
# @Descripttion:
# @version:
# @Author: zpliu
# @Date: 2022-11-12 23:17:33
 # @LastEditors: zpliu
 # @LastEditTime: 2022-11-15 19:55:38
# @@param:
###

#*V1 版本是使用双亲存在差异的位点，而无论该位点在F2中的MAF是否大于0.2 
#? 对于某一个区域，统计该区域内每个SNP在pool池中的ZY043基因型的占比，对所有SNP求平均值
python SNPIndex_V1.py ../pool_permutate/poolsample_${i}.txt splitData/random_${i}.txt






#* V2版本
#? 在F2和亲本材料中，有两种基因型频率并且MAF > 0.2
for i in $(seq 1 100); do
    bsub -q normal -n 1 -e test.err -o test.out "
    python SNPIndex_V2.py ../pool_permutate/poolsample_${i}.txt splitData/random_${i}.txt
"
done 



