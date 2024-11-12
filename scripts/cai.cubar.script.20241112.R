#
#
####################Required inputs:
#
#inpTestAln --> fasta file with genes to test for CAI
#hostCDS --> host CDS fasta file  
#geneCode --> genetic code (most times you want 1 (standard code))
#cdsMinLen --> suggested value is 500
#inpANN --> input host genome annotation (right now only works with gff)
#annFmt (gff, gbk) (right now only works with GFF)
#topNgenes (for now, use 500)
#outSfx --> for downstream convenience
#outDir --> name of output directory
#
##############
#
#
library(Biostrings)
library(cubar)
#
#
#
##############################Step00: prepare input:
#
#
#
##Define genetic code (default: standard) and prepare  codon table
#1: Standard
Ctable<- get_codon_table(gcid = as.character(geneCode))
#
#
#import alignment of genes to test and clean them up for downstream testing:
testAln<- FXprepFA(inpFAfile = inpTestAln, inpCtable = Ctable,
                   inpCDSMinLen=cdsMinLen)
#
#
##import and prepare Host CDS:
CDS<- FXprepFA(inpFAfile = hostCDS, inpCtable = Ctable, inpCDSMinLen=cdsMinLen)
#
# 
#
##Get host tRNA anticodon counts:
if(annFmt=='gff') {
  tRNACounts<- FXprepGFF(inpGFF = inpANN)
} ##Right now it only works with gff. Change this when you got time.
#
#
#
#
##Get codon matrix for host CDS:
CodMtx<- count_codons(CDS)
##Now we gotta remove the met and trp:
#Met: anticodon = cau ; codon = ATG
#Trp: anticodon = cca ; codon = TGG
CodMtx[,colnames(CodMtx)=='ATG']<- 0
CodMtx[,colnames(CodMtx)=='TGG']<- 0
#
##Get codon matrix for test genes:
testMtx<- count_codons(testAln)
##Now we gotta remove the met and trp:
#Met: anticodon = cau ; codon = ATG
#Trp: anticodon = cca ; codon = TGG
testMtx[,colnames(testMtx)=='ATG']<- 0
testMtx[,colnames(testMtx)=='TGG']<- 0
#
#
#
#
###Get host CDS TAI:
hTAI<- FX_calcTAI(
  inpMtx = CodMtx, inptRNAcounts = tRNACounts, inpCodTable = Ctable)
#
#
#
#
####Now use tai score to get host (presumably) highly expressed CDSs:
hTAI<- hTAI[order(hTAI)]
#
###Take into account rounding errors:
myTol<- 1e-7
TopNames<- names(hTAI)[
  hTAI>=(hTAI[(length(hTAI)-(topNgenes-1))] - myTol)]
#
#Now extract the top genes:
TopCDS<- CDS[names(CDS) %in% TopNames]
#
##And get the matrix for these ones:
myTopMtx<- CodMtx[rownames(CodMtx) %in% TopNames,]
#
#
#
#
####Now calculate rscu for these top genes:
myTopRscu<- est_rscu(cf = myTopMtx, codon_table = Ctable) #, weight = 
##We could add the weights from TAI, should we?. Nah, leave it alone.
#
#
#
#
#####
# Finally, calculate CAI for test genes:
#
myTestCAI<- get_cai(cf = testMtx, rscu = myTopRscu)
# 
#Save as table:
caiDF<- data.frame(
  seqName = names(myTestCAI),
  CAI = myTestCAI,
  Notes = outSfx,
  stringsAsFactors = FALSE)
dir.create(path = outDir, showWarnings = FALSE)
write.table(
  x = caiDF, sep = '\t', row.names = FALSE, quote = FALSE,
  file = paste(
    outDir, '/cai.', outSfx, '.tsv', sep = ''))
#
#
