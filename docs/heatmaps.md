---
layout: heatmaps
permalink: /RBD-heatmaps/
---

---

### Overview

You can use this tool to explore the experimentally determined impacts of amino acid mutations on ACE2-binding affinity and expression in SARS-CoV-2 receptor-binding domain (RBD) variants.

#### Instruction

To use this tool, select the RBD variants that you wish to display in each heatmap by selecting a variant from the drop down menu corresponding to each plot. Then, you can select the metric that you wish to display in the heatmap (either change in ACE2 binding affinity (-log10 $$K_D$$), or change in RBD expression (log(MFI)) by selecting that metric in the corresponding drop down menu. Hover over individual mutations to see exact numerical details.

#### Technical Details

The impact on ACE2 receptor-binding affinity ($$\Delta$$ log10 $$K_D$$) or RBD expression ($$\Delta$$ log(MFI)) of every single amino-acid mutation in SARS-CoV-2 RBDs, as determined by high-throughput titration assays. Wildtype amino acids are indicated by an 'x', and gray squares indicate missing mutations from each library. The number of internally replicated barcodes with which a mutation was measured is visible as `Barcode Count` in the tooltips, where higher numbers indicate higher-confidence measurements. The experiments underlying these data can be found in our preprint [here](to_add).

Data for variants from Wuhan-Hu-1 through Omicron XBB.1.5 are from previously published studies [here](https://www.science.org/doi/10.1126/science.abo7896), [here](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1010951), and [here](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1011901).

### Data

Raw data  can be found [here](https://github.com/tstarrlab/SARS-CoV-2-RBD_DMS_Omicron-EG5-FLip-BA286/blob/main/results/final_variant_scores/final_variant_scores.csv). The code used to make these plots can be found [here](https://github.com/tstarrlab/SARS-CoV-2-RBD_DMS_Omicron-EG5-FLip-BA286/blob/main/RBD-Heatmaps-Interactive-Visualization.ipynb).
