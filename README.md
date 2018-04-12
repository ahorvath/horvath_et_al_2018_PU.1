# Horvath et al. 2018 

This repository contains the most important scripts (written in R and bash) used to create data and figures for the Horvath et al. 2018 manuscript.

# Collaborative transcription factor cistromes predict chromatin accessibility and identify labeled enhancers controlling macrophage polarization

Attila Horvath1,3,8, Lajos Szeles1,8, Bence Daniel2,8, Laszlo Steiner5, Zsolt Czimmerer1, Ixchelt Cuaranta-Monroy1, Nikolas Giannakis5, Mate Kiss1, Lilla Ozgyin1, Zoltan Simandi2,  Szilard Poliska1,5, Emanuele Raineri6,7, Ivo G. Gut6,7, Benedek Nagy4 and Laszlo Nagy1,2,3,*

1Department of Biochemistry and Molecular Biology, Faculty of Medicine, University of Debrecen, Debrecen, Hungary

2Sanford-Burnham-Prebys Medical Discovery Institute, Orlando, Florida, USA

3MTA-DE ‘Lendulet’ Immunogenomics Research Group, University of Debrecen, 4Department of Mathematics, Eastern Mediterranean University, 
Famagusta, North Cyprus, Mersin 10, Turkey

5UD-GenoMed, Ltd. Nagyerdei krt. 98. Debrecen, Hungary H-4012

6Centro Nacional de Analisis Genomico (CNAG-CRG), Center for Genomic Regulation (CRG), Barcelona Institute for Science and Technology 
(BIST), C/Baldiri Reixac 4, 08028 Barcelona, Spain.

7Universitat Pompeu Fabra (UPF), Barcelona, Spain

8These authors contributed equally to this work

*Correspondence: lnagy@sbpdiscovery.org

ABSTRACT

The emerging concept of tissue-specific gene expression posits that lineage-determining transcription factors (LDTFs) determine the enhancer repertoire making these genomic loci accessible to other transcription factors (TFs). However, the interrelationship and the guiding principles of LDTFs’ occupancy, chromatin accessibility and enhancer activity has not been systematically evaluated. We addressed this issue by integrating genome-wide technologies (ChIP-Seq, ATAC-Seq, GRO-Seq and mRNA-Seq) and computational approaches using murine bone marrow-derived macrophages as a model system. We found that only one third of the PU.1 (also known as SPI1) cistrome, an essential LDTF for macrophages, is associated with highly accessible chromatin. Using the Random Forest classifier, a machine learning method, we demonstrated that PU.1-related chromatin accessibility could be accurately predicted using only occupancy values of the main TFs of macrophages (PU.1, IRF8, JUNB, CEBPA and RUNX1), but not from PU.1 alone. These TFs contribute to chromatin accessibility and enhancer activation in a collaborative and hierarchical manner. Finally, we identified a novel class of regulatory regions (termed “labeled enhancers”), characterized as PU.1-bound, low accessible genomic regions in unstimulated cells, providing molecular beacons for the IL4-induced alternative polarization program. Collectively, our study demonstrates that combining epigenomic approaches and computational tools is an efficient strategy to identify the main determinants of chromatin openness and the inducible enhancer repertoire. 




*ATAC_0h_1h.R

*ChIP-seq_anal-v1_9mm10_fast.sh

*H4ac_0h_1h.R

*STAT6_0h_1h.R

*boxplot_Irf8_down.R

*boxplot_Irf8_up.R

*boxplot_PolII_S2_any_up.R

*boxplot_PolII_S5_any_up.R

*calculate_FPKM.sh

*create_meta_hist_general.sh

*diffbind_five_factor.R

*donuts.R

*randomForest_PU1_cistrome_openness_add.R

*remap_motif_mscore.sh

*spec_pie.R

*tsv2venn.sh
