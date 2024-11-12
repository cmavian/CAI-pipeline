#
######################Template to define your input files
#
#
#inpTestAln --> fasta file with genes to test for CAI
inpTestAln<- 'genes.test.2seqs.fasta'
#hostCDS --> host CDS fasta file  
hostCDS<- '../../host_genomes/aaegL5_0_refseq/cds_from_genomic.fna'
#genCode --> genetic code (most times you want 1 (standard code))
geneCode<- 1
#cdsMinLen --> suggested value is 500
cdsMinLen<- 500
#inpANN --> input host genome annotation (right now only works with gff)
inpANN<- '../../host_genomes/aaegL5_0_refseq/genomic.gff'
#annFmt (gff, gbk) (right now only works with GFF)
annFmt<- 'gff'
#topNgenes (for now, use 500)
topNgenes<- 500
#
#outSfx --> for downstream convenience
outSfx<- 'test'
#outDir --> name of output directory
outDir<- 'testRes'




