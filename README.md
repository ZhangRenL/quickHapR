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
vcf = import_vcf("vcf/Seita.1G001600_136756_144094_-_3k_final.vcf.gz")
gff = import_gff("gff/Yugu1.gff3")
phenos = import_pheno("pheno/allPheno.txt")

# 计算并输出单倍型结果
hap = get_hap(vcf, filter_Chr = T, Chr = "scaffold_1", filter_POS = T, startPOS = 136756, endPOS = 144094)
hapResult = hap_result(hap, out =T, file = "results/Seita.1G001600_hapResult.txt")
  
# 可视化单倍型结果
plotGeneStructure(gff, hapResult)
plotHapTable(hapResult = hapResult)

# 单倍型与表型的关联分析
hapVsPheno(hap, phenos, phenoName = colnames(phenos)[1], hapPrefix = "H",geneID = "Seita.1G001600")

```



library(quickHapR)
setwd("/data/zhangrenliang/GeneFamily/kinesin/Haptype/")
# 自定义参数区 ####
vcfPath = "/data/zhangrenliang/GeneFamily/kinesin/Haptype/vcf/cleanvcf/"
hapPath = "/data/zhangrenliang/GeneFamily/kinesin/Haptype/haps/"
resultDir = "/data/zhangrenliang/GeneFamily/kinesin/Haptype/Results/Architecture"
phenoFile = "/data/zhangrenliang/Hap/pheno/S_Architecture_株型.txt"
gffFile = "/data/zhangrenliang/GeneFamily/kinesin/Haptype/gff/Sitalica_312_v2.2.gene_exons.gff3"
hapTypei = "Gene"
dir.exists(hapPath)
file.exists(gffFile)
file.exists(phenoFile)
if(!stringr::str_ends(resultDir,"/")) resultDir = paste0(resultDir,"/")
if(!stringr::str_ends(vcfPath,"/")) resultDir = paste0(resultDir,"/")
if(!dir.exists(resultDir)) dir.create(resultDir)

gff <- import_gff(gffFile = gffFile)
phenos = import_pheno(phenoFile = phenoFile, comment.char = "#")
for (file in dir(vcfPath)){
  GeneID <- stringr::str_extract(file,"Seita.[0-9]{0,1}[GJ][0-9]{6}")
  GeneID <- stringr::str_split(file,"_")[[1]][1]
  startPOS <- stringr::str_split(file,"_")[[1]][2]  %>% as.numeric()
  endPOS <- stringr::str_split(file,"_")[[1]][3] %>% as.numeric()
  strand <- stringr::str_split(file,"_")[[1]][4]
  message(GeneID)
  vcf = import_vcf(vcf_file = paste0(vcfPath, file))
  hap = get_hap(vcf, filter_Chr = F, Chr = "scaffold_1",filter_POS = T, startPOS = startPOS, endPOS = endPOS)
  hapResult = hap_result(hap, out =T, file = paste0(resultDir, GeneID, "_hapResult.txt"))
  #
  pdf(file = paste0(resultDir,"Architectures_",GeneID,".pdf"), width = 16, height = 8)
  plotGeneStructure(gff, hapResult, startPOS = startPOS, endPOS = endPOS)
  hapTable = plotHapTable(hapResult = hapResult,geneID = GeneID)
  plot(hapTable)
  for (pheno in colnames(phenos))
  {
    PlotResults = hapVsPhenos(hap, phenos[,1:2],hapPrefix = "H",geneID = GeneID)
    plot(PlotResults$figs)
    message(GeneID,": ", pheno)
  }
  dev.off()
}





#### 特技

1.  使用 Readme\_XXX.md 来支持不同的语言，例如 Readme\_en.md, Readme\_zh.md
2.  Gitee 官方博客 [blog.gitee.com](https://blog.gitee.com)
3.  你可以 [https://gitee.com/explore](https://gitee.com/explore) 这个地址来了解 Gitee 上的优秀开源项目
4.  [GVP](https://gitee.com/gvp) 全称是 Gitee 最有价值开源项目，是综合评定出的优秀开源项目
5.  Gitee 官方提供的使用手册 [https://gitee.com/help](https://gitee.com/help)
6.  Gitee 封面人物是一档用来展示 Gitee 会员风采的栏目 [https://gitee.com/gitee-stars/](https://gitee.com/gitee-stars/)
