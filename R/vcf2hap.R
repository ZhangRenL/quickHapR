#' @name get_hap
#' @title generat haps from vcf
#' @description  generate hap format from vcf
#' @usage get_hap(vcf, hap_prefix = "H",
#' filter_Chr = F, Chr = "scaffold_1",
#' filter_POS = F, startPOS = as.numeric(), endPOS = as.numeric(),
#' hyb_remove = T, na.drop = T)
#' @author Zhangrenl
#' @description  import vcf file
#' @import tidyr
#' @import vcfR
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom  stats na.omit
#' @param vcf vcf
#' @param filter_Chr filter vcf by Chrome or not, defalt is FALSE. If TRUE, the Chr is needed.
#' @param filter_POS filter vcf by Position or not, defalt is FALSE. If TRUE, startPOS and endPOS are needed.
#' @param hap_prefix hap_prefix, defalt is "H"
#' @param hyb_remove Remove accessions contains hybrid site or not. Defalt is TRUE.
#' @param na.drop Drop Accessions contains unknown allele site or note. Defalt is TRUE
#' @param Chr Needed for filter vcf by Chrom
#' @param startPOS,endPOS Needed for filter vcf by position. startPOS must less than endPOS
#' @export
get_hap <- function(vcf,
                    hap_prefix = "H",
                    filter_Chr = F, Chr = "scaffold_1",
                    filter_POS = F, startPOS = as.numeric(), endPOS = as.numeric(),
                    hyb_remove = T,
                    na.drop = T) {
  # to do: add a filter_vcf function
  # for choose promoter, gene, downstream range for haplotyoe analysis
  requireNamespace('tidyr')
  requireNamespace('dplyr')
  if(filter_Chr){
    if(missing(Chr)) stop("Chr must be character")
    vcfChr <- vcf@fix[,1]
    probe <- vcfChr %in% Chr
    vcf@fix <- vcf@fix[probe,]
    vcf@gt <- vcf@gt[probe,]
  }

  # filter Postion according given range
  if(filter_POS){
    if(missing(startPOS)) stop("startPOS must be numeric")
    if(missing(endPOS)) stop("endPOS must be numeric")
    if(startPOS >= endPOS) stop("startPOS must less tan endPOS")
    vcfPOS <- as.numeric(vcf@fix[,2])
    probe <- c(vcfPOS >= startPOS & vcfPOS <= endPOS)

    #
    if(!(TRUE %in% probe)) {
      e = paste0("There is no variant in selected range. \nPlease check vcf file between ",
                 startPOS," and ", endPOS, ".")
      return(e)
    }
    vcf@fix <- vcf@fix[probe,]
    vcf@gt <- vcf@gt[probe,]
  }

  # vcf2data.frame for analysis
  gt <- vcfR::extract_gt_tidy(vcf)
  CHR <- vcfR::getCHROM(vcf)
  POS <- vcfR::getPOS(vcf)
  REF <- vcfR::getREF(vcf)
  ALT <- vcfR::getALT(vcf)
  ALLELE <- paste0(REF,"/",ALT)
  INFO <- vcfR::getINFO(vcf)
  hap <- tidyr::pivot_wider(data = gt,
                           id_cols = .data$Key,
                           names_from = .data$Indiv,
                           values_from = .data$gt_GT_alleles)
  hap <- dplyr::select(hap, -c(.data$Key))
  hap <- as.matrix(hap)
  rownames(hap) <- POS

  # convert Indel into +/-
  for(l in 1:nrow(hap)){
    if(stringr::str_length(ALLELE[l]) > 3){
      if(stringr::str_length(REF[l]) > stringr::str_length(ALT[l])){
        gety <- c("++", "+-","-+","--")
      } else {
        gety <- c("--", "-+","+-","++")
      }
      REFl <- REF[l]
      ALTl <- ALT[l]
      names(gety) <- paste(c(REFl,REFl,ALTl,ALTl),c(REFl,ALTl,REFl,ALTl),sep = "/")
      probe <- hap[l,]
      hap[l,] <- gety[probe]
    } else {
      hap[l,] <- stringr::str_remove_all(hap[l,], "[/]")
    }
  }

  hap <- t(hap)


  # deal with heterozygosis site set as "H"
  probe_hyb <- c("AA","CC","GG","TT","++","--")
  hap[!(hap %in% probe_hyb)] <- "H"
  if(hyb_remove) {
    hap[!(hap %in% probe_hyb)] <- NA
    hap <- na.omit(hap)
  }

  # drop na rows
  hap[hap == "."] <- NA
  if(na.drop) hap <- na.omit(hap)


  # reform the genotypes
  # A/A -> A; T/T ->T; C/C -> C; G/G ->G
  for(i in probe_hyb) {
    hap[hap == i] <- stringr::str_sub(i,1,1)
  }


  # name haps
  hap <- data.frame(hap, check.rows = F, check.names = F)

  HapID <- tidyr::unite(hap, dplyr::matches("[0-9]{1,}"),col = "IDs", sep = "")
  HapID <- HapID$IDs
  hap <- cbind(Hap = HapID, hap)
  hap$Accession <- row.names(hap)
  haps <- table(hap$Hap)
  haps = haps[order(haps,decreasing = T)]
  hapnms <- stringr::str_pad(c(1:length(haps)),3,"left","0")
  hapnms <- paste0(hap_prefix, hapnms)
  names(hapnms) <- names(haps)
  hap$Hap <- hapnms[hap$Hap]
  hap <- hap[order(hap$Hap),]

  # add infos
  meta <- rbind(c("CHR",CHR,""),
               c("POS",POS,""),
               c("INFO", INFO,""),
               c("ALLELE",ALLELE,""))
  colnames(meta) <- colnames(hap)
  hap <- rbind(meta, hap)
  rownames(hap) <- c(1:nrow(hap))


  # removed Redundancy cols
  removecols = c()
  for(c in 1:ncol(hap)){
    namec = colnames(hap)[c]
    if(!(namec %in% c("Hap", "Accession", "freq"))){
      gtc = unique(hap[-c(1:4),c])
      if(length(gtc) == 1) {
        removecols = c(removecols, c)
      }
    }
  }
  if(!is.null(removecols)) hap = hap[, -removecols]
  return(hap)
}


#' @name hap_result
#' @title generate hap results
#' @description summarize hap result and output a txt file
#' @import tidyr
#' @importFrom utils write.table
#' @param hap hap
#' @param out write hap results to a txt file
#' @param file file path
#' @export
hap_result <- function(hap, out = T, file = "hapResult.txt"){
  requireNamespace('tidyr')
  hapResults <- hap %>% data.frame(check.names = F)
  hapfre <- table(hapResults[,1])
  hapfre <- hapfre[stringr::str_starts(names(hapfre),"H")]
  hapResults <- hapResults %>% tidyr::chop(cols = "Accession")
  hapResults$freq[5:nrow(hapResults)] <- hapfre[hapResults[5:nrow(hapResults),1]]
  Acc <- c()
  for(i in 1:length(hapResults$Accession))  Acc[i] <- paste(hapResults$Accession[[i]],collapse = ";")
  hapResults$Accession <- Acc

  if(out)  utils::write.table(hapResults, file = file, sep = "\t",quote = F,row.names = F,col.names = F)
  return(hapResults)
}
