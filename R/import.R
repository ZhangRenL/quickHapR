#' @name import_vcf
#' @title import vcf from file
#' @usage import_vcf(vcf_file = vcf_file, ...)
#' @description read and import vcf file
#' @importFrom vcfR read.vcfR
#' @author Zhangrenl
#' @description  import vcf file
#' @param vcf_file file path of vcf
#' @param ... pass to vcfR::read.vcfR
#' @usage import_vcf(vcf_file = vcf_file, ...)
#' @export
import_vcf <- function(vcf_file = vcf_file, ...) {
    vcf <- vcfR::read.vcfR(vcf_file, ...)
    return(vcf)
}


#' @name import_pheno
#' @title imports phenos from file
#' @usage import_pheno(phenoFile, comment.char = "#", ...)
#' @description 第一列为Accession
#' 第二列起为表型，
#' 第一行为表头（表型名称）：“表型名称（计量单位）.地点[年份]”
#' @importFrom utils read.delim
#' @param phenoFile pheno file path, should be a table separated by tab
#' @param comment.char comment.char, start with comment.char will be ignored
#' @param ... parameters will pass to read.delim
#' @export
import_pheno <- function(phenoFile, comment.char = "#", ...){
    phenos <- utils::read.delim(phenoFile,
                                check.names = FALSE,
                                row.names = 1,
                                comment.char = comment.char, ...)
    return(phenos)
}


#' @name import_gff
#' @title  import_gff
#' @usage import_gff(gffFile)
#' @importFrom rtracklayer import
#' @param gffFile gffFile
#' @export
import_gff <- function(gffFile){
    gff <- rtracklayer::import(gffFile)
    return(gff)
}



#### Import the pipe operator from magrittr ####
#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @importFrom IRanges %over%
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

