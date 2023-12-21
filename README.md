# LPSStudy
This is the code of manuscript 'A gene network driving hypertension associated dampened LPS response in monocytes from coronary artery disease patients'

How the cardiovascular disease (CVD) risk factors (e.g. a number of biological factors and unhealthy lifestyles) affect the activities of circulating monocytes, one of the most important cell types in CVD, remains unclear.

This study is to examine the impact of CVD risk factors on monocyte transcriptional responses to an infectious stimulus.

The datasets and codes generated during and/or analyzed during this study are available from the corresponding author on reasonable request.

Input data:
----------------------------------
data/LPSgeneexpr.RData : processed LPS gene expression, and baseline, and lps response.

data/geneinfo.rds : gene names and id, WGCNA clustering (modules)

data/riskfactor14.rds: risk factors;

data/network.txt :  aracne network generated using ARACNe-AP

data/MEs.rds: Module eigengenes

data/drugs.rds: drug info

Some file and datasets should be downloaded in advance:
----------------------------------
GO file from GSEA website http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

TFs file download from http://humantfs.ccbr.utoronto.ca/

(Optional) LINCS CMap L1000 data GSE 92742 and GSE70138 from :https://clue.io/connectopedia/lincs_cmap_data, 
including GSE70138_Broad_LINCS_inst_info，GSE70138_Broad_LINCS_pert_info，GSE70138_Broad_LINCS_sig_info，GSE70138_Broad_LINCS_sig_metrics, GSE92742_Broad_LINCS_cell_info, GSE92742_Broad_LINCS_gene_info_delta_landmark, GSE92742_Broad_LINCS_gene_info, GSE92742_Broad_LINCS_inst_info, GSE92742_Broad_LINCS_pert_info, GSE92742_Broad_LINCS_pert_metrics, GSE92742_Broad_LINCS_sig_info GSE92742_Broad_LINCS_sig_metrics, GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx, GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx

Softwares
----------------------------------

For network visualization, you need to install Cytoscape

For GRN generation, you need to download ARACNe-AP



Data analysis 
------------------------------------
(You can run all the steps below and then run "Figure Visualization" to reproduce all the figures in the manuscript. However, since all outputs from 'Data analysis' are already saved in the 'data' folder, you can also skip this and just run "Figure Visualization").

Step 1: 1_DFWGCNA.rmd ("lIMMA and WGCNA analysis for LPS response gene matrix" (WGCNA result has been pre-calculated to prevent the clustering algorithm from generating random results, which may result in the figures in the manuscript not being reproduced.) 
input: data/geneinfo.rds : gene names and id; 
Output: 'data/DFWGCNA.rds'

Step 2: 2_Correlation_all.rmd , calculating all correlations
input:'data/DFWGCNA.rds' from step 1, 'data/riskfactor14.rds' ,data/LPSgeneexpr.RData'
output: corrall.rds


Step 3: 3_ORA_analysis.Rmd, ORA analysis for each WGCNA module
input:'data/corrall.rds' from step 2
output: GOAlllist.rds


Step 3: 3_Bayesian_network.Rmd, Bayesian network construction based on 10 modules
input:'data/corrall.rds' from step 2，LPSgeneexpr.RData
output: bayesiannetwork.rds (the results may slightly different from the Figure 2E if you rerun this step because of the bootstrapping.)


Step 3: 3_TF_analysis.Rmd, TF analysis using Aracne and viper
input:'data/corrall.rds' from step 2， LPSgeneexpr.RData
output: vipertfinfoMS.rds， cortfBPDMS.rds  (if you recalculate VIPER, you may get slightly different results for Figure S4B)
output : ARACNEnet.rds， ARACNEvertices.rds，ARACNEregulons.rds


Figure visualization 
--------------------------------
1. Figures related to the whole 7933 genes (Fig 1 2, S1,  S3): Figure_whole_01.rmd


2. Figures related to module Salmon (Fig 3 4, S4): Figure_Salmon_02.rmd


3. Figures S2 and confunding analysis: Add_Analysis_03.rmd








