# quickHapR

## 介绍

**通过vcf文件对基因进行单倍型分析**

1. 基于二代测序得到的大量群体的基因型信息进行单倍型分析，

2. 结合表型数据筛选出优异单倍型进行后续研究
   [quickHapR](https://gitee.com/zhangrenl/quickHapR)

## 基本逻辑

1. 导入数据（vcf数据，GFF注释，表型数据）
2. 通过vcf文件计算单倍型，导出单倍型数据
3. 单倍型数据结合注释信息将变异位点标注在基因示意图上；展示各单倍型之间的进化关系
4. 比较不同单倍型之间的表型差异，筛选优势单倍型

## 安装教程

**安装准备**

1. 安装Rtools软件
2. 安装git软件
3. 安装R packages: `devtools`，`BiocManager`
4. quickHapR 安装命令

```R
if(!require("quickHapR")) 
devtools::install_git("https://gitee.com/zhangrenl/quickhapr")
```

* 如果上述命令安装失败可前往Gitee下载最新预编译R包，选择本地安装
* 本地安装quickHapR前需手动安装依赖的R packages：`ggpubr`, `vcfR`, `tidyverse`, `stringr`, `reshape2`, `randomcoloR`,

  `rtracklayer`, `trackViewer`, `GenomicRanges`, `IRanges`

## 使用说明

#### 

```R
# 加载quickhapR
library(quickHapR)
data("quickHap_test") # 加载测试数据,处理自己的数据时不必执行该行

# 设定工作目录
setwd("/your/working/directory")


# 导入数据
vcf = import_vcf("Seita.1G001600_136756_144094_-_3k_final.vcf.gz")
gff = import_gff("Yugu1.gff3")
phenos = import_pheno("allPheno.txt")

# 计算并输出单倍型结果
# hap,data.frame:第一列与最后一列分别固定为HAP和Accession，中间列为位置及对应的基因型
# 前四行为注释信息分别是：CHR，POS，ALLELE,INFO
hap = get_hap(vcf,                 # import_vcf() 导入的vcfR
              filter_Chr = F,      # 筛选染色体选项
              Chr = "scaffold_1",  # 通过染色体对vcf信息进行筛选
              filter_POS = T,      # 通过位置进行筛选
              startPOS = 136756,   # Numeric, 起始位置，通过位置对vcf信息进行筛选
              endPOS = 144094)     # Numeric, 终止位置，通过位置对vcf信息进行筛选

# hapResult,data.frame: 第一列固定为HAP，最后两列分别固定为Accession和freq，中间列为位置及对应的基因型
# 前四行为注释信息分别是：CHR，POS，ALLELE,INFO
hapResult = hap_result(hap,        # hap 结果
                       out =T,     # 是否输出文件，如果为TRUE， 必须指定输出路径file
                       file = "results/Seita.1G001600_hapResult.txt")  # 输出文件路径（tab分隔的表格）


# 可视化单倍型结果
plotGeneStructure(gff,                # gff注释信息
                  hapResult,          # 单倍型结果
                  Chr = "scaffold_1", # 基因所在染色体
                  startPOS = 136756,  # 基因结构示意图的起始位点
                  endPOS = 144094,    # 基因结构示意图的终止位置
                  type = "pin",       # SNP类型
                  cex = 1,            # circle大小
                  CDS_h = 0.05,       # 不同基因结构的高度
                  fiveUTR_h = 0.02, 
                  threeUTR_h = 0.01) 
plotHapTable(hapResult,               # 单倍型结果
             hapPrefix = "H",         # 单倍型前缀（阿拉伯数字前的字母）
             geneID = "",             # 基因ID， 作为图表Title
             title.color = "grey90")  # 表头底色

# 单倍型与表型的关联分析
res = hapVsPheno(hap,        # data.frame:第一列与最后一列分别固定为HAP和Accession，中间列为位置及对应的基因型
                 phenos,      # data.frame: 第一列固定为Accession，随后各列为表型数据，phenoName作为colnames
                 phenoName = "yourPhenoName",  # 本次分析中使用的表型名称
                 hapPrefix = "H",            # 单倍型编号的前缀
                 geneID = "Seita.1G000000",  # 基因ID， 作为表头信息
                 mergeFigs = T,    # 是否将两图融合
                 minAcc = 5)       # 需要分析的单倍型包含的数据量最小值
                 
# plot(res$fig_pvalue)
# plot(res$fig_Violin)

plot(res$figs)
```


