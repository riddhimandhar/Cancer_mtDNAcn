# Cancer_mtDNAcn
Code for analysis of mtDNA CN variation and differential expression in cancer


There are three folders. 

A> Folder '0_normal_data_mtDNA_CN': 
   - contains codes for processing and calculation of mtDNA copy number of normal 
     samples - either normal tissue, normal blood or normal other (such as lymph node). 


B> Folder '1_mutation_data_processing':
   - contains codes for processing and analyzing mutation data and SVs in all cancer types 
     individually as well as pan-cancer. 
   - contains codes for predictive regression and randon forest models 


C> Folder '2_gene_expression_data_processing':
   - contains codes for differential expression of genes and transcripts between high- and low- mtDNA samples 
   - contains codes for analysis of hallmark gene sets 
   - contains codes for predictive modeling of expression of genes in hallmark data 


D> The code '3_drug_response_resistance_risk.R' is used for calculation of likelihood score for 
   developing chemotherapy resistance in cancer samples based on their normalized mtDNA copy number 
