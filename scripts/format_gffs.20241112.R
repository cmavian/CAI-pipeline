#
#GFFs have slightly different format, gotta fix that one by one:
#
dir.create(path = './trna_counts', showWarnings = FALSE)
#
host_path<- '' ###define path to tRNA count dir
#
#
myCounts<- FXprepGFF(
  inpGFF = paste(
    host_path,
    'hsap_grch38p14_refseq/ncbi_dataset/data/GCF_000001405.40/genomic.gff',
    sep = ''))
write.table(
  x = myCounts, file = 'trna_counts/hsapiens.anticodons.counts.tsv',
  sep = '\t', quote = FALSE, row.names = FALSE)
unique(myCounts$AA)
##Some of them have stop codon in them, but it's fine.
##Remember to remove them later though, together with trp and met
rm(myCounts)
#rm(FXprepGFF)
#
#
#
