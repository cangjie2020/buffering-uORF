library(GenomicFeatures)
library(rtracklayer)
library(data.table)

# ensembl info ###############################################################################
extractGeneInfo <- function(gtf.path){
  txdb <- makeTxDbFromGFF(gtf.path)
  
  txlen <- transcriptLengths(txdb, with.cds_len = TRUE, with.utr5_len = TRUE, with.utr3_len = TRUE)
  x <- rtracklayer::import(gtf.path, 'gtf')
  x <- as.data.frame(x)
  scols <- c('gene_id', 'gene_name', 'gene_biotype', 'transcript_id',
             'transcript_biotype', 'protein_id', 'seqnames', 'strand')
  scols <- scols[scols %in% colnames(x)]
  txtype <- x[, scols]
  setDT(txtype)
  txtype <- unique(txtype)
  res <- merge(txlen, txtype, by.x = c('gene_id', 'tx_name'), by.y = c('gene_id', 'transcript_id'))
  res <- as.data.table(res)
  res <- res[order(gene_id, tx_name, protein_id)][!duplicated(tx_name)]
  setnames(res, 'seqnames', 'chrom')
  return(res)
}

args <- commandArgs(trailingOnly = TRUE)
res <- extractGeneInfo(args[1])
fwrite(res, file = args[2], sep = '\t')