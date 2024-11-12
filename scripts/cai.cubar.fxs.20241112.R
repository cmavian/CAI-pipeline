#
###this function counts the trna anticodon genome copies from gff:
FXprepGFF<- function(inpGFF) {
  #
  intGff<- read.table(
    file = inpGFF, header = FALSE, sep = '\t', comment.char = '#',
    quote = '', col.names = c('seqid', 'source', 'type', 'start', 'end', 'score',
      'strand', 'phase', 'attributes'), blank.lines.skip = TRUE,
    fill = TRUE)
  #
  #Remove mitochondrial genes first:
  mitChroms<- unique(
    intGff[grepl(pattern = 'mgt|mitochondrion', x = intGff$attributes, ignore.case = TRUE), "seqid"])
  intGff<- intGff[intGff$type=='tRNA' & !intGff$seqid %in% mitChroms,]
  #
  ##Now we gotta parse the anticodons:
  intGff$AntiCodon<- sapply(
    X = strsplit(x = intGff$attributes, split = ';'),
    FUN = "[", 2)
  ##Actually, there's already the numbering,
  ##But it probably easier to make the table
  intGff$AntiCodon<- gsub(
    pattern = "-[0-9]+$|[0-9]-[0-9]+$", replacement = '',
    x = intGff$AntiCodon)
  intGff$AntiCodon<- gsub(
    pattern = ".*-", replacement = '', x = intGff$AntiCodon)
  intGff$AA<- sapply(
    X = strsplit(x = intGff$attributes, split = 'product=tRNA-'),
    FUN = "[", 2)
  unique(intGff$AA) ##Good.
  ##Now we gotta make the counts:
  inttRNAcounts<- data.frame()
  #
  for(i in unique(intGff$AA)) {
    indDF<- intGff[intGff$AA==i,]
    indCounts<- as.data.frame(table(indDF$AntiCodon))
    colnames(indCounts)<- c('AntiCodon', "Abundance")
    indCounts$AA<- i
    inttRNAcounts<- rbind(inttRNAcounts, indCounts)
  }
  #
  #
  ##Remove Met and Trp, as they are only coded by one triplet
  #Met: anticodon = cau ; codon = ATG
  #Trp: anticodon = cca ; codon = TGG
  #
  #inttRNAcounts<- inttRNAcounts[!inttRNAcounts$AA %in% c('Met', 'Trp'),]
  #DO IT LATER, save the complete table for now
  #####ALSO ALSO: cubar takes anticodons with T, not U!!!
  inttRNAcounts$AntiCodon<- toupper(inttRNAcounts$AntiCodon)
  inttRNAcounts$AntiCodon<- gsub(
    pattern = 'U', replacement = 'T', x = inttRNAcounts$AntiCodon)
  #
  return(inttRNAcounts)
} ## End of FXprepGFF --> prepare tRNA genomic copies table
#
#
#
#
#############################
#
#
FXprepFA<- function(
    inpFAfile, inpCtable, inpCDSMinLen=cdsMinLen,
    strtChk = TRUE) {
  #
  inpCDS<- readDNAStringSet(filepath = inpFAfile)
  #
  ##Keep only seqss with at least 500bps
  inpCDS<- check_cds(
    seqs = inpCDS, codon_table = inpCtable,
    min_len = inpCDSMinLen,
    check_start = strtChk, ##individual viral pept might not have iMet
    check_stop = FALSE) ##by default it also removes iMet and stop codons,
  ##Also makes sure length is multiple of 3.
  # 
  return(inpCDS)
} ## End of FXprepFA --> quality check input fasta
#
#
#############################
#
#
#
#
#############################
#
#
FX_calcTAI<- function(inpMtx, inptRNAcounts, inpCodTable) { #(...,inpTopGenes)
  require(cubar)
  require(Biostrings)
  #
  #Step 1: 
  #Estimate tRNA weights per codon:
  #A) anticodon counts:
  tRNAvec<- inptRNAcounts$Abundance
  names(tRNAvec)<- toupper(inptRNAcounts$AntiCodon)
  #
  #Sij is a penalty on non-canonical codon anticodon pairings and
  ### differs among different species (Sabi and Tuller 2014). Cubar uses average Sij
  # values for eukaryotes (Sabi and Tuller 2014) by default.
  ##Absolute adaptiveness values are then normalized by maximum
  #
  #
  ## Get tRNA weight:
  W<- est_trna_weight(
    trna_level = tRNAvec, codon_table = inpCodTable)
  #
  ######Now get tRNA Adaptation Index (TAI):
  #
  TAI<- get_tai(cf = inpMtx, trna_w = W)
  ##It's a numeric vector, with names indicating the CDS name.
  ##Highly expressed genes present high tAI values
  ##(> 0.4 ?? -- >  https://github.com/mt1022/cubar ),
  ## which means that their codon usage resembles the genomic structure of tRNA genes.)
  #
  #
  return(TAI)
} #End of FX_calcTAI --> calculates host cds TAI
#
#############################
#
#
