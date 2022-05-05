#' @name plotHapTable
#' @title plotHapTable
#' @importFrom randomcoloR randomColor
#' @importFrom stringr str_starts
#' @importFrom stringr str_length
#' @import tidyr
#' @import ggplot2
#' @param hapResult hapResult
#' @param hapPrefix hapPrefix
#' @param geneID geneID
#' @export
plotHapTable = function(hapResult, hapPrefix = "H", geneID = ""){
  requireNamespace('tidyr')
  if("Accession" %in% colnames(hapResult)) hapResult = hapResult[,colnames(hapResult) != 'Accession']
  ALLELE = hapResult[hapResult[,1] == "ALLELE",]
  hps = hapResult[stringr::str_starts(hapResult[,1],hapPrefix),]
  hps = rbind(ALLELE, hps)
  foot = c()
  nfoot = 1
  for(i in 2:length(ALLELE)){
    if(stringr::str_length(ALLELE[i]) > 3){
      note = paste0("*",nfoot,collapse = '')
      nfoot = nfoot + 1
      foot = c(foot,paste0(note,": ", ALLELE[i]))
      ALLELE[i] = note
    }
  }

  meltHapRes = reshape2::melt(hps,1)
  colnames(meltHapRes) = c('Var1','Var2',"value")
  lab = meltHapRes
  lab$value = stringr::str_replace_all(lab$value, c("AA"="A", "TT"="T","CC"="C","GG"="G","[+]{2}"="+","--"="-"))
  meltHapRes$value[stringr::str_detect(meltHapRes$value,"[0-9]")] = NA
  levels = as.vector(unique(meltHapRes$Var1))
  levels = levels[order(levels, decreasing = T)]
  meltHapRes$Var1 = factor(meltHapRes$Var1, levels = levels)
  if(is.null(foot)) foot = " " else foot = paste(foot,collapse = ";")
  fig0 = ggplot2::ggplot(data = meltHapRes,
                         mapping = ggplot2::aes_(x=~Var2, y=~Var1, fill=~value)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes_(x=~Var2, y=~Var1,
                                    label = lab$value)) +
    ggplot2::scale_fill_discrete(na.value = "white") +
    ggplot2::labs(caption = foot) +
    ggplot2::ggtitle(label = geneID) +  ggplot2::scale_y_discrete() +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(position = "top")) +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x =  ggplot2::element_blank(),
      axis.title.y =  ggplot2::element_blank(),
      panel.grid.major =  ggplot2::element_blank(),
      panel.border =  ggplot2::element_blank(),
      panel.background =  ggplot2::element_blank(),
      axis.ticks =  ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle=ifelse(ncol(hapResult) >= 9, 45, 0),
                                          vjust = ifelse(ncol(hapResult) >= 9, 0.1, 0.5),
                                          hjust = ifelse(ncol(hapResult) >= 9, 0.1, 0.5)),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top", title.hjust = 0.5))
return(fig0)
}



#' @name import_gff
#' @title  import_gff
#' @importFrom rtracklayer import
#' @param gffFile gffFile
#' @export
import_gff = function(gffFile){
  gff = rtracklayer::import(gffFile)
  return(gff)
}



#' @name plotGeneStructure
#' @title plotGeneStructure
#' @importFrom trackViewer lolliplot
#' @importFrom  GenomicRanges GRanges
#' @importFrom  GenomicRanges strand
#' @importFrom IRanges IRanges %over%
#' @import tidyr
#' @param gff gff
#' @param hapResult hapResult
#' @param Chr Chr
#' @param startPOS startPOS
#' @param endPOS endPOS
#' @export
plotGeneStructure = function(gff, hapResult,
                             Chr,
                             startPOS, endPOS){
# lolliplot
  requireNamespace("trackViewer")
  requireNamespace("tidyr")
  if(missing(gff)) {
    message("gfdf");
    stop("missing gff")}
  if(missing('hapResult')) {
    message("gf2f");
    stop("missing hapResult")}

  geneElement = c("CDS","three_prime_UTR","five_prime_UTR")
  if("Accession" %in% colnames(hapResult)) hapResult = hapResult[,colnames(hapResult) != 'Accession']
  meta = hapResult[1:4,-1]
  POS = as.numeric(meta[2,])
  SNP = meta[4,]

  if(missing(Chr)) Chr = meta[1,1]
  if(missing(startPOS)) startPOS = min(POS)
  if(missing(endPOS)) endPOS = max(POS)

  SNP.gr <- GenomicRanges::GRanges(Chr, IRanges::IRanges(POS, width=1,
                                 names = paste0(POS,"(",SNP,")")),
                    color = sample.int(6, length(SNP), replace=TRUE),
                    score = sample.int(5, length(SNP), replace = TRUE),
                    angle = 45)



  gene = GenomicRanges::GRanges(Chr,
                 IRanges::IRanges(start = min(startPOS,endPOS),
                         end = max(startPOS,endPOS)))
  over = gff[gff %over% gene]
  over$height[over$type == "CDS"] = 0.05
  over$height[over$type == "three_prime_UTR"] = 0.01
  over$height[over$type == "five_prime_UTR"] = 0.015

  features = over[over$type %in% geneElement]
  strands = as.character(GenomicRanges::strand(features))
  layerID = unlist(features$Parent)
  layerID = stringr::str_remove_all(layerID,".v2.2")
  layerID = paste0(layerID,"(",ifelse(strands == "+", "5'->3'","3'<-5'"),")")
  features$featureLayerID = layerID
  names(features) = features$featureLayerID
  l = length(unique(names(features)))
  if (l < 6) fillc = c((1:l) + 1) else fillc = randomcoloR::randomColor(l)

  names(fillc) = unique(names(features))
  features$fill = fillc[names(features)]
  # features = c(gene, features)

  trackViewer::lolliplot(SNP.gr, features, type = "pin", jitter = NULL, ylab = "",cex = 0.5,yaxis = F)
}

