#
#
rm(list = ls())
require(cubar)
require(Biostrings)
##
source('cai.cubar.fxs.R')
#
########################Set up paths and input files:
#
cds_path<- '' ##path to host genome CDS dir
tRNAcounts_path<- '' ##Path to tRNA counts dir
testGenes_path<- "" ##Path to dir containing viral genomes to test
#
myHosts<- data.frame(
  hosts = c(
    "aaegypti", "cbrevitarsis", "cjacchus",
    "cquinquefasciatus", "hsapiens","sapella"),
  trnaCountsFile = c(
    "aaegypti.anticodons.counts.tsv",
    "cbrevitarsis.anticodons.counts.tsv",
    "cjacchus.anticodons.counts.tsv",
    "cquinquefasciatus.anticodons.counts.tsv",
    "hsapiens.anticodons.counts.tsv",
    "sapella.anticodons.counts.tsv"),
  cdsFile = c(
    'aaegL5_0_refseq/cds_from_genomic.fna',
    'cbrev_refseq/data/GCF_036172545.1/cds_from_genomic.fna',
    'cjacchus_mCalJa1_2_pat_X/ncbi_dataset/data/GCF_011100555.1/cds_from_genomic.fna',
    'cquinquefasciatus_VPISU_Cqui_1_0_pri_paternal/ncbi_dataset/data/GCF_015732765.1/cds_from_genomic.fna',
    'hsap_grch38p14_refseq/ncbi_dataset/data/GCF_000001405.40/cds_from_genomic.fna',
    'sapella_gscmonkey1/ncbi_dataset/data/GCF_009761245.1/cds_from_genomic.fna'),
  geneCode = 1,
  stringsAsFactors = FALSE)
#
myViruses<- data.frame(
  gene = c(
    'segL','segM', "segS"),
  fileName = c(
    "OROV_segL_all_aln_edit.fa",
    "OROV_segM_all_aln_edit.fa",
    "OROV_segS_all_aln_edit.fa"),
  stringsAsFactors = FALSE)
#
#
#
########################Set up some parameters:
#
#
#inpTestAln --> fasta file with genes to test for CAI
#hostCDS --> host CDS fasta file  
#geneCode --> genetic code (most times you want 1 (standard code))
###geneCode<- 1
#cdsMinLen --> suggested value is 500
cdsMinLen<- 500
#Amino acids to keep. Feel free to discard stop codons + 
## Met and Trp, as they are only coded by one triplet and
#they can skew the analyses for proteins with high content of those
rmMetTrp<- TRUE
chkTestStart<- FALSE ## viral pept might not have iMet
cdsStop<- 'Sec'
#
#
#inpANN --> input host genome annotation (right now only works with gff)
#annFmt (gff, gbk) (right now only works with GFF)
#topNgenes (for now, use 500)
topNgenes<- 500
#outSfx --> for downstream convenience
#outDir --> name of output directory
host_hegs<- paste(getwd(), 'host_heg/', sep = '/')
out_dir<- paste(getwd(), 'cai_out/', sep = '/')
#
##############
#
##Full output table:
myCAIdf<- data.frame()
#
#
#
#########################Finally, let's run the analyses:
#
#
#since the finding of the host hge is the most time expensive process,
##Do that first for each host, so not to repeat it by set of test genes:
#
#
for(i in 1:nrow(myHosts)) {
  indCtable<- get_codon_table(gcid = as.character(myHosts[
    i, "geneCode"]))
  indHost<- myHosts[i, "hosts"]
  #
  #
  ###########infer host highly expressed genes:
  ##Load host cds:
  indCDS<- FXprepFA(
    inpFAfile = paste(cds_path, myHosts[i, "cdsFile"], sep = ''),
    inpCtable = indCtable,
    inpCDSMinLen=cdsMinLen)
  #
  ##Get host tRNA anticodon counts:
  indtRNACounts<- read.delim(
    file = paste(
      tRNAcounts_path, myHosts[i,"trnaCountsFile"], sep = ''),
    sep = '\t', header = TRUE,
    stringsAsFactors = FALSE)
  #
  #Remove met and Trp?
  if(rmMetTrp) {
    indtRNACounts<- indtRNACounts[!indtRNACounts$AA %in% c('Met', 'Trp'),]
  }
  #Remove stop codons:
  indtRNACounts<- indtRNACounts[indtRNACounts$AA!=cdsStop,]
  #
  #
  #########Get codon matrix for host CDS:
  CodMtx<- count_codons(indCDS)
  ##Remove the met and trp? ###sholud ot hard code it like this
  ##BUt for now it'll do. We should remove AAs with 1 triplet
  ##To be ore correct; fix it when you have the time
  ##Or if working on hosts using non-standard code
  #Met: anticodon = cau ; codon = ATG
  #Trp: anticodon = cca ; codon = TGG
  if(rmMetTrp) {
    #CodMtx<- CodMtx[, -which(colnames(CodMtx) %in% c('ATG', 'TGG'))]
    CodMtx[,colnames(CodMtx)=='ATG']<- 0
    CodMtx[,colnames(CodMtx)=='TGG']<- 0
  }
  #
  #
  ########infer host highly expressed genes:
  #
  ##Get host CDS TAI:
  hTAI<- FX_calcTAI(
    inpMtx = CodMtx, inptRNAcounts = indtRNACounts, inpCodTable = indCtable)
  #
  #
  #########Now use tai score to get host (presumably) highly expressed CDSs:
  hTAI<- hTAI[order(hTAI)]
  #
  ###Take into account rounding errors:
  myTol<- 1e-7
  TopNames<- names(hTAI)[
    hTAI>=(hTAI[(length(hTAI)-(topNgenes-1))] - myTol)]
  #
  #Now extract the top genes:
  TopCDS<- indCDS[names(indCDS) %in% TopNames]
  #
  ##And get the matrix for these ones:
  TopMtx<- CodMtx[rownames(CodMtx) %in% TopNames,]
  #
  #
  #########Now calculate rscu for these top genes:
  TopRscu<- est_rscu(cf = TopMtx, codon_table = indCtable) #, weight = 
  ##We could add the weights from TAI, should we?. Nah, leave it alone.
  #
  ########Save these results:
  dir.create(path = host_hegs, showWarnings = FALSE)
  write.table(
    x = TopRscu, sep = '\t', quote = FALSE,
    row.names = FALSE,
    file = paste(
      host_hegs, indHost, '.heg.rscu.tsv', sep = ''))
  #
  #####
  #
  ###############Now test CAI for each viral alignment:
  #
  #import viral alignment:
  for(z in 1:nrow(myViruses)) {
    #
    #Import alignment:
    testAln<- FXprepFA(
      inpFAfile = paste(
        testGenes_path,
        myViruses[z,"fileName"], sep = ''),
      inpCtable = indCtable,
      inpCDSMinLen=6, strtChk = chkTestStart) ##don't filter these by length
    #
    #Get codon matrix for test genes:
    testMtx<- count_codons(testAln)
    ##Now we gotta remove the met and trp:
    #Met: anticodon = cau ; codon = ATG
    #Trp: anticodon = cca ; codon = TGG
    if(rmMetTrp) {
      #testMtx<- testMtx[, -which(colnames(testMtx) %in% c('ATG', 'TGG'))]
      testMtx[,colnames(testMtx)=='ATG']<- 0
      testMtx[,colnames(testMtx)=='TGG']<- 0
    }
    #
    #
    #########Finally, calculate CAI:
    #
    indTestCAI<- get_cai(cf = testMtx, rscu = TopRscu)
    # 
    #Save as table:
    dir.create(path = out_dir, showWarnings = FALSE)
    #
    caiDF<- data.frame(
      seqName = names(indTestCAI),
      CAI = indTestCAI,
      Gene = myViruses[z, "gene"],
      #Year = myViruses[z,"Year"],
      host = indHost,
      stringsAsFactors = FALSE)
    #dir.create(path = out_dir, showWarnings = FALSE)
    write.table(
      x = caiDF, sep = '\t', row.names = FALSE, quote = FALSE,
      file = paste(
        out_dir, myViruses[z, "gene"], '.',
        #myViruses[z,"Year"], '.',
        indHost, '.cai.tsv', sep = ''))
    #
    #
    #Also add to main output:
    myCAIdf<- rbind(myCAIdf, caiDF)
  } ##End of z loop by viral alignment
} ##End of i loop by host.
#
#
#
write.table(
  x = myCAIdf, sep = '\t', row.names = FALSE, quote = FALSE,
  file = paste(
    out_dir, 'all.cai.tsv', sep = ''))
#
#Done :)





