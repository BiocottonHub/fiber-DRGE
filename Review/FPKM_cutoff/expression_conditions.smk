# cutoff = config["cutoff"]


rule all:
    input:
        filterExpressedData=expand(
            "{stage}/{stage}_expressed_{cutoff}.txt",
            stage=["0DPA", "4DPA", "8DPA", "12DPA", "16DPA", "20DPA"],
            cutoff=[0.01,0.05, 0.1, 0.5, 1, 5, 10, 20],
        ),
        stageMergeData=expand(
            "ExpressedGene/expressedGeneStats_{cutoff}.txt",
            cutoff=[0.01,0.05, 0.1, 0.5, 1, 5, 10, 20],
        ),
        homoeologousExpressionBias=expand(
            "{stage}/{stage}_expressionBias_{cutoff}.txt",
            stage=["0DPA", "4DPA", "8DPA", "12DPA", "16DPA", "20DPA"],
            cutoff=[0.01,0.05, 0.1, 0.5, 1, 5, 10, 20],
        ),
        expressVariantFile=expand(
            "ExpressedGene/expressionVariation_{cutoff}.txt",
            cutoff=[0.01,0.05, 0.1, 0.5, 1, 5, 10, 20],
        ),
        HomoeologExpressionBiasV2=expand(
            "{stage}/{stage}_expressionBiasV2_{fold_change}_accession_{percentage}.txt",
            stage=["0DPA", "4DPA", "8DPA", "12DPA", "16DPA", "20DPA"],
            fold_change=[2, 3, 4],
            percentage=[0.05, 0.1, 0.15],
        ),


rule filterExpressedGenes:
    input:
        rawExpression_stringTieV1="/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/MappingFPKM/TM1/AllGenes/StringTie_V1/{stage}_AllSamples_FPKM.txt",
    output:
        expressionedGeneId="{stage}/{stage}_expressed_{cutoff}.txt",
    shell:
        """
        python ./script/getExpressedGene.py {input.rawExpression_stringTieV1} {wildcards.cutoff} {output.expressionedGeneId}
        """


rule HomoeologousBias:
    input:
        expressionData="/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/MappingFPKM/TM1/AllGenes/{stage}_AllSamples_FPKM.txt",
        expressedGeneList="{stage}/{stage}_expressed_{cutoff}.txt",
    output:
        biasGeneStat="{stage}/{stage}_expressionBias_{cutoff}.txt",
    resources:
        mem_mb=10000,
    threads: 1
    shell:
        #统计每个标准下每个时期内发生Bias的基因数
        """
        python ./script/HomoeologousBiasType.py {input.expressionData} {wildcards.cutoff} {input.expressedGeneList} {output.biasGeneStat}
        """


rule stageExpressedMerge:
    input:
        stage_0DPA="0DPA/0DPA_expressed_{cutoff}.txt",
        stage_4DPA="4DPA/4DPA_expressed_{cutoff}.txt",
        stage_8DPA="8DPA/8DPA_expressed_{cutoff}.txt",
        stage_12DPA="12DPA/12DPA_expressed_{cutoff}.txt",
        stage_16DPA="16DPA/16DPA_expressed_{cutoff}.txt",
        stage_20DPA="20DPA/20DPA_expressed_{cutoff}.txt",
    output:
        statFile="ExpressedGene/expressedGeneStats_{cutoff}.txt",
    shell:
        """
        python ./script/expressedGeneStat.py {input.stage_0DPA} {input.stage_4DPA} {input.stage_8DPA} {input.stage_12DPA} {input.stage_16DPA} {input.stage_20DPA} {output.statFile}
        """


rule expressionVariation:
    input:
        expressionArray=expand(
            "/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/MappingFPKM/TM1/AllGenes/StringTie_V1/{stage}_AllSamples_FPKM.txt",
            stage=["0DPA", "4DPA", "8DPA", "12DPA", "16DPA", "20DPA"],
        ),
        expressedGeneList_0DPA="0DPA/0DPA_expressed_{cutoff}.txt",
        expressedGeneList_4DPA="4DPA/4DPA_expressed_{cutoff}.txt",
        expressedGeneList_8DPA="8DPA/8DPA_expressed_{cutoff}.txt",
        expressedGeneList_12DPA="12DPA/12DPA_expressed_{cutoff}.txt",
        expressedGeneList_16DPA="16DPA/16DPA_expressed_{cutoff}.txt",
        expressedGeneList_20DPA="20DPA/20DPA_expressed_{cutoff}.txt",
    output:
        expressVariantFile="ExpressedGene/expressionVariation_{cutoff}.txt",
    shell:
        """
        python ./script/expressionVariance.py {wildcards.cutoff} {output.expressVariantFile}
        """


rule expressionBiasConditions:
    input:
        expressionData="/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/MappingFPKM/TM1/AllGenes/{stage}_AllSamples_FPKM.txt",
        expressedGeneList="{stage}/{stage}_expressed_0.1.txt",
    output:
        expressionBias_V2="{stage}/{stage}_expressionBiasV2_{fold_change}_accession_{percentage}.txt",
    shell:
        """
        python  ./script/HomoeologousBiasType_V2.py  {input.expressionData} {wildcards.fold_change} {wildcards.percentage} {input.expressedGeneList} {output.expressionBias_V2}
        """
