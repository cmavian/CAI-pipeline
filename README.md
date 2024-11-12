# CAI-pipeline
Codon adaptability index pipeline
#
#

#######
TO BE COMPLETED

Codon adaptability index 

The codon adaptation index (CAI) was determined in order to evaluate the relative adaptability of the codon usage of a gene towards that of highly expressed genes within a given host (Sharp et al, 1987). The maximum CAI is 1 and the minimum is 0. In general, the higher the CAI value, the more efficient the mRNA can be translated. The CAI values of OROV genomes were analyzed in the context of the following: Homo sapiens as host (NCBI accession no. GCF_000001405.40), and several known reservoirs hosts species implicated in the transmission of OROV: nonhuman primates such as capuchin monkeys (Sarajus spp.; genome available for the Sapajus apella - GCF_009761245.1) and marmosets (genome available for the common marmoset Callithrix jacchus - GCF_011100555.1). The main vector responsible for the transmission of OROV in sylvatic and urban cycle are midges Culicoides paraensis; the mosquito species Coquillettidia venezuelensis, Aedes serratus, Culex quinquefasciatus have also been reported as potential secondary vectors in sylvatic and urban cycles (Cardoso et al., 2015; McGregor et al., 2021; Travassos et al., 2017). In absence of genomes for Culicoides paraensis, Culicoides sonorensis, Coquillettidia venezuelensis, Aedes serratus, we have tested Culicoides brevitarsis (GCF_036172545.1) as potential sylvatic vector, Cx. quinquefasciatus (GCF_015732765.1) as cosmopolitan vector (rural, peri‑urban, and urban environments), and Ae. aegypti (GCF_002204515.2) as potential urban vector, and as potential sylvatic and urban vector (de Mendonca et al., 2021). Because eukaryotic viruses generally do not have their own translation machinery, they typically evolve codon usage in response to their host tRNA pools, thus viruses evolved their increased translation efficiency match the mosting abundant cognate tRNA. Gene expression is proportional to respective tRNA anticodon genome copy numbers (van Weringh et al., 2011;     REF XXXXXXXX. ###CONTINUE FROM HERE




Thus, we have identified the highly expressed genes for each host were inferred by calculating the tRNA adaptation index (tAI REF ) of host CDS sequences at least 500 base pairs long: sequences with no start codon or having an internal stop codon were excluded; these filters were implemented in order to reduce the chance of pseudogenes being part of the dataset. Then, we have identified best tAI were considered highly expressed REF  from the top 500 genes REF excluding single triplet amino acids (Met and Trp) to reduce bias analyses in favor of proteins rich in such amino acids, as they would get a high CAI REF . Afterwards, we have calculated CAI for each each target viral sequence. Analyses were performed in R v.4.4 and RStudio. Fasta files were imported into R using the Biostring package REF ; tAI and CAI were calculated using the R package cubar REF  and CAI plots were drawn with the ggplot2 package REF .


References

Cardoso BF, Serra OP, Heinen LB, Zuchi N, Souza VC, Naveca FG, Santos MA, Slhessarenko RD. Detection of Oropouche virus segment S in patients and inCulex quinquefasciatus in the state of Mato Grosso, Brazil. Mem Inst Oswaldo Cruz. 2015 Sep;110(6):745-54. doi: 10.1590/0074-02760150123. PMID: 26517653; PMCID: PMC4667577.

de Mendonça SF, Rocha MN, Ferreira FV, Leite THJF, Amadou SCG, Sucupira PHF, Marques JT, Ferreira AGA, Moreira LA. Evaluation of Aedes aegypti, Aedes albopictus, and Culex quinquefasciatus Mosquitoes Competence to Oropouche virus Infection. Viruses. 2021 Apr 25;13(5):755. doi: 10.3390/v13050755. PMID: 33923055; PMCID: PMC8145018.

McGregor BL, Connelly CR, Kenney JL. Infection, Dissemination, and Transmission Potential of North American Culex quinquefasciatus, Culex tarsalis, and Culicoides sonorensis for Oropouche Virus. Viruses. 2021 Feb 2;13(2):226. doi: 10.3390/v13020226. PMID: 33540546; PMCID: PMC7912852.

Sharp PM, Li WH. The codon Adaptation Index--a measure of directional synonymous codon usage bias, and its potential applications. Nucleic Acids Res. 1987 Feb 11;15(3):1281-95. doi: 10.1093/nar/15.3.1281. PMID: 3547335; PMCID: PMC340524.

Travassos da Rosa JF, de Souza WM, Pinheiro FP, Figueiredo ML, Cardoso JF, Acrani GO, Nunes MRT. Oropouche Virus: Clinical, Epidemiological, and Molecular Aspects of a Neglected Orthobunyavirus. Am J Trop Med Hyg. 2017 May;96(5):1019-1030. doi: 10.4269/ajtmh.16-0672. Epub 2017 Feb 6. PMID: 28167595; PMCID: PMC5417190.

van Weringh A, Ragonnet-Cronin M, Pranckeviciene E, Pavon-Eternod M, Kleiman L, Xia X. HIV-1 modulates the tRNA pool to improve translation efficiency. Mol Biol Evol. 2011 Jun;28(6):1827-34. doi: 10.1093/molbev/msr005. Epub 2011 Jan 7. PMID: 21216840; PMCID: PMC3098512.





