# CAI-pipeline
Codon adaptability index pipeline
#
#
The codon adaptation index (CAI) is a measure of the relative adaptability of the codon usage of a gene/viral gene towards that of highly expressed genes within a given host (Sharp et al, 1987). The maximum CAI is 1 and the minimum is 0. In general, the higher the CAI value, the more efficient the mRNA can be translated.

The CAI values of OROV genomes were analyzed in the context of the following: *Homo sapiens* as host (NCBI accession no. GCF_000001405.40), and several known reservoirs hosts species implicated in the transmission of OROV: nonhuman primates such as capuchin monkeys (*Sarajus spp.*; genome available for the Sapajus apella - GCF_009761245.1) and marmosets (genome available for the common marmoset *Callithrix jacchus* - GCF_011100555.1). The main vector responsible for the transmission of OROV in sylvatic and urban cycle are midges *Culicoides paraensis*; the mosquito species *Coquillettidia venezuelensis*, *Aedes serratus*, *Culex quinquefasciatus* have also been reported as potential secondary vectors in sylvatic and urban cycles (Cardoso et al., 2015; McGregor et al., 2021; Travassos et al., 2017). In absence of genomes for *C. paraensis*, *C. sonorensis*, *C. venezuelensis*, *Ae. serratus*, we have tested *Culicoides brevitarsis* (GCF_036172545.1) as potential sylvatic vector, *Cx. quinquefasciatus* (GCF_015732765.1) as cosmopolitan vector (rural, peri‑urban, and urban environments), and *Ae. aegypti* (GCF_002204515.2) as potential sylvatic and urban vector (de Mendonca et al., 2021).

Because eukaryotic viruses generally do not have their own translation machinery, they typically evolve codon usage in response to their host tRNA pools; thus, viruses evolved their increased translation efficiency to match the most abundant cognate tRNA. Gene expression is proportional to respective tRNA anticodon genome copy numbers (Chithambaram et al., 2014a; Chithambaram et al., 2014b; Prabhakaran et al., 2014; Prabhakaran et al., 2015; van Weringh et al., 2011). Therefore, we have identified the highly expressed genes for each host were inferred by calculating the tRNA adaptation index (dos Reis e al., 2003) of host CDS sequences at least 500 base pairs long: sequences with no start codon or having an internal stop codon were excluded; these filters were implemented in order to reduce the chance of pseudogenes being part of the dataset. Then, we have identified best tAI were considered highly expressed (Navon and Pilpel, 2011) from the top 500 genes. Single triplet amino acids (Met and Trp) were excluded to reduce bias analyses in favor of proteins rich in such amino acids, as they would get a high CAI (Xia, 2007). Afterwards, we have calculated CAI for each each target viral sequence.

Analyses were performed in R v.4.4 (R Core team) and RStudio (R Studio team, 2015). Fasta files were imported into R using the Biostring package (Pagès et al., 2020); tAI and CAI were calculated using the R package cubar (Zhang et al., 2024) and CAI plots were drawn with the ggplot2 R package (Wickham H., 2016).

#
## **References**
Cardoso BF, Serra OP, Heinen LB, Zuchi N, Souza VC, Naveca FG, Santos MA, Slhessarenko RD. Detection of Oropouche virus segment S in patients and inCulex quinquefasciatus in the state of Mato Grosso, Brazil. Mem Inst Oswaldo Cruz. 2015 Sep;110(6):745-54. doi: 10.1590/0074-02760150123. PMID: 26517653; PMCID: PMC4667577.

Chithambaram S, Prabhakaran R, Xia X. Differential codon adaptation between dsDNA and ssDNA phages in Escherichia coli. Mol Biol Evol. 2014 (a) Jun;31(6):1606-17. doi: 10.1093/molbev/msu087. Epub 2014 Feb 27. PMID: 24586046; PMCID: PMC4032129.

Chithambaram S, Prabhakaran R, Xia X. The effect of mutation and selection on codon adaptation in Escherichia coli bacteriophage. Genetics. 2014 (b) May;197(1):301-15. doi: 10.1534/genetics.114.162842. Epub 2014 Feb 28. PMID: 24583580; PMCID: PMC4012488.

de Mendonça SF, Rocha MN, Ferreira FV, Leite THJF, Amadou SCG, Sucupira PHF, Marques JT, Ferreira AGA, Moreira LA. Evaluation of Aedes aegypti, Aedes albopictus, and Culex quinquefasciatus Mosquitoes Competence to Oropouche virus Infection. Viruses. 2021 Apr 25;13(5):755. doi: 10.3390/v13050755. PMID: 33923055; PMCID: PMC8145018.

dos Reis M, Wernisch L, Savva R. Unexpected correlations between gene expression and codon usage bias from microarray data for the whole Escherichia coli K-12 genome. Nucleic Acids Res. 2003 Dec 1;31(23):6976-85. doi: 10.1093/nar/gkg897. PMID: 14627830; PMCID: PMC290265.

McGregor BL, Connelly CR, Kenney JL. Infection, Dissemination, and Transmission Potential of North American Culex quinquefasciatus, Culex tarsalis, and Culicoides sonorensis for Oropouche Virus. Viruses. 2021 Feb 2;13(2):226. doi: 10.3390/v13020226. PMID: 33540546; PMCID: PMC7912852.

Navon S, Pilpel Y. The role of codon selection in regulation of translation efficiency deduced from synthetic libraries. Genome Biol. 2011;12(2):R12. doi: 10.1186/gb-2011-12-2-r12. Epub 2011 Feb 1. PMID: 21284851; PMCID: PMC3188794.

Pagès H, Aboyoun P, Gentleman R, DebRoy S. Biostrings: Efficient manipulation of biological strings 2020 \[Available from: https://bioconductor.org/packages/Biostrings![image](https://github.com/user-attachments/assets/ed3b51d8-abc6-4d10-aac8-1a49ab2d2a7c) ]
Prabhakaran R, Chithambaram S, Xia X. Aeromonas phages encode tRNAs for their overused codons. Int J Comput Biol Drug Des. 2014;7(2-3):168-82. doi: 10.1504/IJCBDD.2014.061645. Epub 2014 May 28. PMID: 24878728.

Prabhakaran R, Chithambaram S, Xia X. Escherichia coli and Staphylococcus phages: effect of translation initiation efficiency on differential codon adaptation mediated by virulent and temperate lifestyles. J Gen Virol. 2015 May;96(Pt 5):1169-1179. doi: 10.1099/vir.0.000050. Epub 2015 Jan 22. PMID: 25614589; PMCID: PMC4631060.

RCoreTeam. R: A language and environment for statistical computing.: R Foundation for Statistical Computing, Vienna,
Austria. URL; 2016 \[Available from: https://www.R-project.org/ ]

RStudioTeam. RStudio: Integrated Development for R: RStudio, Inc., Boston, MA; 2015 \[Available from: http://www.rstudio.com/ ]

Sharp PM, Li WH. The codon Adaptation Index--a measure of directional synonymous codon usage bias, and its potential applications. Nucleic Acids Res. 1987 Feb 11;15(3):1281-95. doi: 10.1093/nar/15.3.1281. PMID: 3547335; PMCID: PMC340524.

Travassos da Rosa JF, de Souza WM, Pinheiro FP, Figueiredo ML, Cardoso JF, Acrani GO, Nunes MRT. Oropouche Virus: Clinical, Epidemiological, and Molecular Aspects of a Neglected Orthobunyavirus. Am J Trop Med Hyg. 2017 May;96(5):1019-1030. doi: 10.4269/ajtmh.16-0672. Epub 2017 Feb 6. PMID: 28167595; PMCID: PMC5417190.

Wickham H. ggplot2: Elegant Graphics for Data Analysis: Springer-Verlag New York; 201![image](https://github.com/user-attachments/assets/65aad21c-1af0-467d-9a58-f17b7d2d2ef9)

Xia X. An improved implementation of codon adaptation index. Evol Bioinform Online. 2007 May 17;3:53-8. PMID: 19461982; PMCID: PMC2684136.

van Weringh A, Ragonnet-Cronin M, Pranckeviciene E, Pavon-Eternod M, Kleiman L, Xia X. HIV-1 modulates the tRNA pool to improve translation efficiency. Mol Biol Evol. 2011 Jun;28(6):1827-34. doi: 10.1093/molbev/msr005. Epub 2011 Jan 7. PMID: 21216840; PMCID: PMC3098512.

Zhang H, Liu M, Zi B (2024). cubar: Codon Usage Bias Analysis. R package version 1.0.0.9000, https://mt1022.github.io/cubar/, https://github.com/mt1022/cubar 
