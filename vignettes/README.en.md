# quickHapR

#### Description
{**When you're done, you can delete the content in this README and update the file with details for others getting started with your repository**}

#### Software Architecture
Software architecture description

#### Installation

1.  xxxx
2.  xxxx
3.  xxxx

#### Instructions

```R
# 加载quickhapR
library(quickHapR)

# 设定工作目录
setwd("/your_working_directory")

# daoru
vcf = import_vcf("vcf/cleanvcf/Seita.1G001600_136756_144094_-_3k_final.vcf.gz")
gff = import_gff("gff/Yugu1.gff3")
phenos = import_pheno("/data/zhangrenliang/GeneFamily/kinesin/Haptype/pheno/allPheno.txt")

# generate haps
hap = get_hap(vcf)
hapResult = hap_result(hap, out =T)

# plot and visualization your results
plotGeneStructure(gff, hapResult)
plotHapTable(hapResult = hapResult)
hapVsPhenos(hap, phenos[,1:2],hapPrefix = "H",geneID = "Seita.0G000000")
```

#### Contribution

1.  Fork the repository
2.  Create Feat_xxx branch
3.  Commit your code
4.  Create Pull Request


#### Gitee Feature

1.  You can use Readme\_XXX.md to support different languages, such as Readme\_en.md, Readme\_zh.md
2.  Gitee blog [blog.gitee.com](https://blog.gitee.com)
3.  Explore open source project [https://gitee.com/explore](https://gitee.com/explore)
4.  The most valuable open source project [GVP](https://gitee.com/gvp)
5.  The manual of Gitee [https://gitee.com/help](https://gitee.com/help)
6.  The most popular members  [https://gitee.com/gitee-stars/](https://gitee.com/gitee-stars/)
