# quickHapR

#### 介绍
{**通过vcf文件对基因进行单倍型分析**
基于二代测序得到的大量群体的基因型信息进行单倍型分析，结合表型数据筛选出优异单倍型进行后续研究
[quickHapR](https://gitee.com/zhangrenl/quickHapR)}

#### 软件架构


#### 安装教程
```R
if(!require("quickHapR")) 
devtools::install_git("https://gitee.com/zhangrenl/quickhapr")
```

#### 使用说明

```R
# 加载quickhapR
library(quickHapR)

# 设定工作目录
setwd("/your/working/directory")

# 导入数据
vcf = import_vcf("vcf/Seita.1G001600_136756_144094_-_3k_final.vcf.gz")
gff = import_gff("gff/Yugu1.gff3")
phenos = import_pheno("pheno/allPheno.txt")

# 计算并输出单倍型结果
# hap and 
hap = get_hap(vcf,                 # vcf imported by import_vcf()
              filter_Chr = F,      # Filter vcf by Chrom or not,defalt is false.
              Chr = "scaffold_1",  # Needed if filter_Chr is TRUE.
              filter_POS = T,      # Filter vcf by Position or not, defalt is TRUE
              startPOS = 136756,   # Numeric, start postion. Needed if filter_POS is TRUE.
              endPOS = 144094)     # Numeric, end position. Needed if filter_POS is TRUE.
hapResult = hap_result(hap,        # hap produced by get_hap()
                       out =T,     # Write the results or not, defalt is TRUE and file is needed
                       file = "results/Seita.1G001600_hapResult.txt") 
  
# 可视化单倍型结果
plotGeneStructure(gff, hapResult)
plotHapTable(hapResult = hapResult)

# 单倍型与表型的关联分析
hapVsPheno(hap, phenos, phenoName = colnames(phenos)[1], hapPrefix = "H",geneID = "Seita.1G001600")

```








#### 特技

1.  使用 Readme\_XXX.md 来支持不同的语言，例如 Readme\_en.md, Readme\_zh.md
2.  Gitee 官方博客 [blog.gitee.com](https://blog.gitee.com)
3.  你可以 [https://gitee.com/explore](https://gitee.com/explore) 这个地址来了解 Gitee 上的优秀开源项目
4.  [GVP](https://gitee.com/gvp) 全称是 Gitee 最有价值开源项目，是综合评定出的优秀开源项目
5.  Gitee 官方提供的使用手册 [https://gitee.com/help](https://gitee.com/help)
6.  Gitee 封面人物是一档用来展示 Gitee 会员风采的栏目 [https://gitee.com/gitee-stars/](https://gitee.com/gitee-stars/)
