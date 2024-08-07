---
layout: epistasis
permalink: /epistatic-shifts/
---

---


### Overview 

You can use this tool to explore epistatic shifts in mutational effects on ACE2-binding affinity (-log10 $$K_D$$) between SARS-CoV-2 receptor-binding domain (RBD) variants. Once you've selected a comparison of interest, you can investigate the mutation-level epistatic shifts. 

#### Instructions

To use this tool, select two SARS-CoV-2 variants that you wish to compare between. To do this, simply select a *'comparator'* in the dropdown menu below the plot and select a *'variant'* by clicking on a variant name in the legend above the plot. Now you're comparing the epistatic shift (Jensen-Shannon divergence between ACE2 binding affinities) at each RBD site between the variant backgrounds you selected. 

Now that you've selected two variants to compare between, you can investigate site level differences by clicking on the points in the higlighted line plot. Simply click on a point – you should see the size of the point change indicating your selection – then you will see the differences in affinities of each of the 20 amino acids measured in each variant appear in the scatter plot on the right. Hover over individual mutations to see exact numerical details.

#### Technical Details

The epistatic shift is calculated as the Jensen-Shannon divergence in the set of Boltzmann-weighted affinities for all amino acids at each site. Mutation affinities were experimentally measured via high-throughput ACE2-binding titrations with yeast-displayed RBDs. The "Barcodes" quantity in per-mutation tooltips indicates the number of internally replicated barcodes with which a mutation was measured, where higher numbers indicate higher-confidence measurements. Mutations with fewer than 3 barcodes in either background are excluded from the Jensen-Shannon divergence calculation.

Data for variants from Wuhan-Hu-1 through Omicron XBB.1.5 are from previously published studies [here](https://www.science.org/doi/10.1126/science.abo7896), [here](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1010951), and [here](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1011901).


Raw data can be found [here](https://github.com/tstarrlab/SARS-CoV-2-RBD_DMS_Omicron-EG5-FLip-BA286/blob/main/results/epistatic_shifts/JSD_by_target.csv) for a table of all pairwise RBD epistatic shifts, and [here](https://github.com/tstarrlab/SARS-CoV-2-RBD_DMS_Omicron-EG5-FLip-BA286/blob/main/results/final_variant_scores/final_variant_scores.csv) for individual measurements of RBD mutant affinities. The code used to make these plots can be found [here](https://github.com/tstarrlab/SARS-CoV-2-RBD_DMS_Omicron-EG5-FLip-BA286/blob/main/Epistatic-Shifts-Interactive-Visualization.ipynb).
