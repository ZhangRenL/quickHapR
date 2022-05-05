# quickHapR

#### 介绍
{**以下是 Gitee 平台说明，您可以替换此简介**
Gitee 是 OSCHINA 推出的基于 Git 的代码托管平台（同时支持 SVN）。专为开发者提供稳定、高效、安全的云端软件开发协作平台
无论是个人、团队、或是企业，都能够用 Gitee 实现代码托管、项目管理、协作开发。企业项目请看 [https://gitee.com/enterprises](https://gitee.com/enterprises)}

#### 软件架构
软件架构说明


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
setwd("/your_working_directory")

# 导入数据
vcf = import_vcf("vcf/cleanvcf/Seita.1G001600_136756_144094_-_3k_final.vcf.gz")
gff = import_gff("gff/Yugu1.gff3")
phenos = import_pheno("/data/zhangrenliang/GeneFamily/kinesin/Haptype/pheno/allPheno.txt")

# 计算并输出单倍型结果
hap = get_hap(vcf)
hapResult = hap_result(hap, out =T, )

# 可视化单倍型结果
plotGeneStructure(gff, hapResult)
plotHapTable(hapResult = hapResult)

# 单倍型与表型的关联分析
hapVsPhenos(hap, phenos[,1:2],hapPrefix = "H",geneID = "Seita.0G000000")

```




#### 特技

1.  使用 Readme\_XXX.md 来支持不同的语言，例如 Readme\_en.md, Readme\_zh.md
2.  Gitee 官方博客 [blog.gitee.com](https://blog.gitee.com)
3.  你可以 [https://gitee.com/explore](https://gitee.com/explore) 这个地址来了解 Gitee 上的优秀开源项目
4.  [GVP](https://gitee.com/gvp) 全称是 Gitee 最有价值开源项目，是综合评定出的优秀开源项目
5.  Gitee 官方提供的使用手册 [https://gitee.com/help](https://gitee.com/help)
6.  Gitee 封面人物是一档用来展示 Gitee 会员风采的栏目 [https://gitee.com/gitee-stars/](https://gitee.com/gitee-stars/)
