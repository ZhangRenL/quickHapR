#' @name import_vcf
#' @title import vcf from file
#' @description read and import vcf file
#' @importFrom vcfR read.vcfR
#' @author Zhangrenl
#' @description  import vcf file
#' @param vcf_file file path of vcf
#' @param ... pass to vcfR::read.vcfR
#' @usage import_vcf(vcf_file = vcf_file, ...)
#' @export
import_vcf <- function(vcf_file = vcf_file, ...) {
  vcf = vcfR::read.vcfR(vcf_file, ...)
  return(vcf)
}


#' @name get_hap
#' @title generat haps from vcf
#' @description  generate hap format from vcf
#' @author Zhangrenl
#' @description  import vcf file
#' @import tidyr
#' @import vcfR
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom  stats na.omit
#' @param vcf vcf
#' @param filter_Chr filter_Chr
#' @param filter_POS filter_POS
#' @param hap_prefix hap_prefix
#' @param  hyb_remove hyb_remove
#' @param na.drop na.drop
#' @export
get_hap <- function(vcf,
                    hap_prefix = "H",
                    filter_Chr = F,
                    filter_POS = F,
                    hyb_remove = T,
                    na.drop = T) {
  requireNamespace('tidyr')
  requireNamespace('dplyr')
  gt = vcfR::extract_gt_tidy(vcf)
  CHR = vcfR::getCHROM(vcf)
  POS = vcfR::getPOS(vcf)
  REF = vcfR::getREF(vcf)
  ALT = vcfR::getALT(vcf)
  ALLELE = paste0(REF,"/",ALT)
  INFO = vcfR::getINFO(vcf)
  hap = tidyr::pivot_wider(data = gt,
                           id_cols = .data$Key,
                           names_from = .data$Indiv,
                           values_from = .data$gt_GT_alleles)
  hap =  dplyr::select(hap, -c(.data$Key))
  hap = as.matrix(hap)
  rownames(hap) = POS
  # convert Indel into +/-
  for(l in 1:nrow(hap)){
    if(stringr::str_length(ALLELE[l]) > 3){
      if(stringr::str_length(REF[l]) > stringr::str_length(ALT[l])){
        gety = c("++", "+-","-+","--")
      } else {
        gety = c("--", "-+","+-","++")
      }
      REFl = REF[l]
      ALTl = ALT[l]
      names(gety) = paste(c(REFl,REFl,ALTl,ALTl),c(REFl,ALTl,REFl,ALTl),sep = "/")
      probe = hap[l,]
      hap[l,] = gety[probe]
    } else {
      hap[l,] = stringr::str_remove_all(hap[l,], "[/]")
    }
  }

  #
  hap = t(hap)
  hap[hap == "."] = NA
  if(na.drop) hap = na.omit(hap)

  # deal with heterozygosis &
  probe_hyb = c("AA","CC","GG","TT","++","--")
  hap[!(hap %in% probe_hyb)] = "H"
  if(hyb_remove) {
    hap[!(hap %in% probe_hyb)] = NA
    hap = na.omit(hap)
  }

  for(i in probe_hyb) {
    hap[hap == i] = stringr::str_sub(i,1,1)
  }

  # name haps
  hap = data.frame(hap, check.rows = F, check.names = F)
  HapID = tidyr::unite(hap, dplyr::matches("[0-9]{1,}"),col = "IDs", sep = "")
  HapID = HapID$IDs
  hap = cbind(Hap = HapID, hap)
  hap$Accession = row.names(hap)
  haps = table(hap$Hap)
  haps = haps[order(haps,decreasing = T)]
  hapnms = stringr::str_pad(c(1:length(haps)),3,"left","0")
  hapnms = paste0(hap_prefix, hapnms)
  names(hapnms) = names(haps)
  hap$Hap = hapnms[hap$Hap]
  hap = hap[order(hap$Hap),]

  # add infos
  meta = rbind(c("CHR",CHR,""),
               c("POS",POS,""),
               c("INFO", INFO,""),
               c("ALLELE",ALLELE,""))
  colnames(meta) = colnames(hap)
  hap = rbind(meta, hap)
  rownames(hap) = c(1:nrow(hap))
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
hap_result = function(hap, out = T, file = "hapResult.txt"){
  requireNamespace('tidyr')
  hapResults = hap %>% data.frame(check.names = F)
  hapfre = table(hapResults[,1])
  hapfre = hapfre[stringr::str_starts(names(hapfre),"H")]
  hapResults = hapResults %>% tidyr::chop(cols = "Accession")
  hapResults$freq[5:nrow(hapResults)] = hapfre[hapResults$POS[5:nrow(hapResults)]]
  Acc = c()
  for(i in 1:length(hapResults$Accession))  Acc[i] = paste(hapResults$Accession[[i]],collapse = ";")
  hapResults$Accession = Acc
  if(out)  utils::write.table(hapResults, file = file, sep = "\t",quote = F,row.names = F,col.names = F)
  return(hapResults)
}

#### Import the pipe operator from magrittr ####
#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL


#' @name %over%
#' @keywords internal
#' @export
#' @importFrom IRanges %over%
NULL


