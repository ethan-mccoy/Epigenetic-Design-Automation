# Epigenetic-Design-Automation

Ethan Dale McCoy's BioE134 Final Project

- Gathers and processes GEO data to find genes with hyper/hypo methylation within a target sample
- Creates dCas9 constructs for the target sample to adjust its methylation towards a control epigenetic profile
- Testing these constructs is meant to be automated with the Opentrons Python Protocol API but I didn't get that far

EpigeneticDesignAutomation.py is the main module

cli.py is a command line interface to play around with

tests.py contains some tests on GEO DataSets of psoriasis lesions and cocaine addicted brains

# Tests output:

Test results for Psoriasis test:

ID4

plasmid: pdCas9-DNMT3A-EGFP

- restriction site: XbaI
- forward oligo: ataatTCTAGA{49: 'gaacagcaatgggctcagac', 256: 'gttcccagttcctaaattcc'}
- reverse oligo: ataatTCTAGA

plasmid: pcDNA3.1-dCas9-MQ1(Q147L)-EGFP

- restriction site: EcoRV
- forward oligo: ataatGATATC{49: 'gaacagcaatgggctcagac', 256: 'gttcccagttcctaaattcc'}gtttaagagctatgctggaaacagcatagcaagtttaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgcttttttt
- reverse oligo: ataatGATATC

CDKN2B

plasmid: pdCas9-Tet1-CD

- restriction site: Acc65I
- forward oligo: ataatGGTACC{75: 'ggccagagcggctttgagct', 272: 'gagagtgcgccggagcagcg', 304: 'gaagagtgtcgttaagttta', 314: 'gttaagtttacggccaacgg'}gtttaagagctatgctggaaacagcatagcaagtttaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgcttttttt
- reverse oligo: ataatGGTACC

plasmid: pcDNA-dCas9-p300 Core

- restriction site: ClaI
- forward oligo: ataatATCGAT{75: 'ggccagagcggctttgagct', 272: 'gagagtgcgccggagcagcg', 304: 'gaagagtgtcgttaagttta', 314: 'gttaagtttacggccaacgg'}gtttaagagctatgctggaaacagcatagcaagtttaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgcttttttt
- reverse oligo: ataatATCGAT

CDKN1A

plasmid: pdCas9-Tet1-CD

- restriction site: Acc65I
- forward oligo: ataatGGTACC{70: 'gcggattcgccgaggcaccg'}gtttaagagctatgctggaaacagcatagcaagtttaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgcttttttt
- reverse oligo: ataatGGTACC

plasmid: pcDNA-dCas9-p300 Core

- restriction site: ClaI
- forward oligo: ataatATCGAT{70: 'gcggattcgccgaggcaccg'}gtttaagagctatgctggaaacagcatagcaagtttaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgcttttttt
- reverse oligo: ataatATCGAT

CDKN2A

plasmid: pdCas9-Tet1-CD

- restriction site: Acc65I
- forward oligo: ataatGGTACC{20: 'gagagcaggcagcgggcggc'}gtttaagagctatgctggaaacagcatagcaagtttaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgcttttttt
- reverse oligo: ataatGGTACC

plasmid: pcDNA-dCas9-p300 Core

- restriction site: ClaI
- forward oligo: ataatATCGAT{20: 'gagagcaggcagcgggcggc'}gtttaagagctatgctggaaacagcatagcaagtttaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgcttttttt
- reverse oligo: ataatATCGAT

Test results for Cocaine abuse test:

CDKN1A

plasmid: pdCas9-DNMT3A-EGFP

- restriction site: XbaI
- forward oligo: ataatTCTAGA{74: 'gcggattcgccgaggcaccg'}
- reverse oligo: ataatTCTAGA

plasmid: pcDNA3.1-dCas9-MQ1(Q147L)-EGFP

- restriction site: EcoRV
- forward oligo: ataatGATATC{74: 'gcggattcgccgaggcaccg'}gtttaagagctatgctggaaacagcatagcaagtttaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgcttttttt
- reverse oligo: ataatGATATC

CCL2

plasmid: pdCas9-DNMT3A-EGFP

- restriction site: XbaI
- forward oligo: ataatTCTAGA{}
- reverse oligo: ataatTCTAGA

plasmid: pcDNA3.1-dCas9-MQ1(Q147L)-EGFP

- restriction site: EcoRV
- forward oligo: ataatGATATC{}gtttaagagctatgctggaaacagcatagcaagtttaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgcttttttt
- reverse oligo: ataatGATATC

PVALB

plasmid: pdCas9-Tet1-CD

- restriction site: Acc65I
- forward oligo: ataatGGTACC{47: 'tccacccccacccgagttgc'}gtttaagagctatgctggaaacagcatagcaagtttaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgcttttttt
- reverse oligo: ataatGGTACC

plasmid: pcDNA-dCas9-p300 Core

- restriction site: ClaI
- forward oligo: ataatATCGAT{47: 'tccacccccacccgagttgc'}gtttaagagctatgctggaaacagcatagcaagtttaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgcttttttt
- reverse oligo: ataatATCGAT

SLC6A3

plasmid: pdCas9-Tet1-CD

- restriction site: Acc65I
- forward oligo: ataatGGTACC{29: 'gagcgggaggggaggcttcg'}gtttaagagctatgctggaaacagcatagcaagtttaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgcttttttt
- reverse oligo: ataatGGTACC

plasmid: pcDNA-dCas9-p300 Core

- restriction site: ClaI
- forward oligo: ataatATCGAT{29: 'gagcgggaggggaggcttcg'}gtttaagagctatgctggaaacagcatagcaagtttaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgcttttttt
- reverse oligo: ataatATCGAT
