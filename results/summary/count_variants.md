# Count variants in each sample
This Python Jupyter notebook counts occurrences of each barcode in each sample from Illumina barcode sequencing, and adds these counts to the codon variant table.

## Set up analysis
### Import Python modules.
Use [plotnine](https://plotnine.readthedocs.io/en/stable/) for ggplot2-like plotting.

The analysis relies heavily on the Bloom lab's [dms_variants](https://jbloomlab.github.io/dms_variants) package:


```python
import itertools
import multiprocessing
import multiprocessing.pool
import os
import warnings

import alignparse
import alignparse.targets

import dms_variants.codonvarianttable
from dms_variants.constants import CBPALETTE
import dms_variants.illuminabarcodeparser
import dms_variants.utils
import dms_variants.plotnine_themes

from IPython.display import display, HTML

import pandas as pd

from plotnine import *

import yaml
```

Set [plotnine](https://plotnine.readthedocs.io/en/stable/) theme to the gray-grid one defined in `dms_variants`:


```python
theme_set(dms_variants.plotnine_themes.theme_graygrid())
```

Versions of key software:


```python
print(f"Using alignparse version {alignparse.__version__}")
print(f"Using dms_variants version {dms_variants.__version__}")
```

    Using alignparse version 0.6.0
    Using dms_variants version 1.4.3


Ignore warnings that clutter output:


```python
warnings.simplefilter('ignore')
```

### Parameters for notebook
Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Make output directory if needed:


```python
os.makedirs(config['counts_dir'], exist_ok=True)
os.makedirs(config['figs_dir'], exist_ok=True)
```

## Input variant tables
Initialize the table of barcode-variant pairs from the respective `process_ccs` notebooks for each background.


```python
variants = pd.read_csv(config['codon_variant_table_file_EG5'], na_filter=None)
variants = variants.append(pd.read_csv(config['codon_variant_table_file_FLip'], na_filter=None))
variants = variants.append(pd.read_csv(config['codon_variant_table_file_BA286'], na_filter=None))

variants = variants.reset_index(drop=True)

display(HTML(variants.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>target</th>
      <th>library</th>
      <th>barcode</th>
      <th>variant_call_support</th>
      <th>codon_substitutions</th>
      <th>aa_substitutions</th>
      <th>n_codon_substitutions</th>
      <th>n_aa_substitutions</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>EG5</td>
      <td>pool1</td>
      <td>AAAAAAAAAAGATTGT</td>
      <td>3</td>
      <td>AGA179GGT</td>
      <td>R179G</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>EG5</td>
      <td>pool1</td>
      <td>AAAAAAAAACAATTGT</td>
      <td>1</td>
      <td>AGA27ACT TTC134GGT</td>
      <td>R27T F134G</td>
      <td>2</td>
      <td>2</td>
    </tr>
    <tr>
      <td>EG5</td>
      <td>pool1</td>
      <td>AAAAAAAAACTATACA</td>
      <td>1</td>
      <td>CAA79TGG</td>
      <td>Q79W</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>EG5</td>
      <td>pool1</td>
      <td>AAAAAAAAAGAAAGTT</td>
      <td>1</td>
      <td>TTG57GAT</td>
      <td>L57D</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>EG5</td>
      <td>pool1</td>
      <td>AAAAAAAAAGACTTCT</td>
      <td>3</td>
      <td>GTT52GCT</td>
      <td>V52A</td>
      <td>1</td>
      <td>1</td>
    </tr>
  </tbody>
</table>


Are there any barcodes in the same library that are shared across targets?
If so, we need to get rid of those as they will be confounded in barcode parsing:


```python
dup_barcodes = (
    variants
    .groupby(['library', 'barcode'])
    .size()
    .rename('duplicate_count')
    .reset_index()
    .query('duplicate_count > 1')
    )

print('Here are duplicated barcodes:')
display(HTML(dup_barcodes.head().to_html(index=False)))

print(f"\nRemoving the {len(dup_barcodes)} duplicated barcodes."
      f"Started with {len(variants)} barcodes:")
variants = (
    variants
    .merge(dup_barcodes, on=['library', 'barcode'], how='outer')
    .query('duplicate_count.isnull()', engine='python')
    )
print(f"After removing duplicates, there are {len(variants)} barcodes.")
```

    Here are duplicated barcodes:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>duplicate_count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool1</td>
      <td>AAAAGTGAAAGAGTAC</td>
      <td>2</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AAACAATAAATTATGA</td>
      <td>2</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AAAGGACTTTTAGAAA</td>
      <td>2</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AACCAGCGAAATCTAA</td>
      <td>2</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>ACAATTGAGGACAATG</td>
      <td>2</td>
    </tr>
  </tbody>
</table>


    
    Removing the 34 duplicated barcodes.Started with 512628 barcodes:
    After removing duplicates, there are 512560 barcodes.


Pull out a target sequence for matching to the barcode and flanking sequence regions. Note, in this pipeline this is ok because our different backgrounds don't have differing flanks or other features within the actual N16 region covered in Illumina sequencing. If ever placing in-line barcodes here in the future, we would need to modify this.


```python
# get wildtype gene sequence for primary target
targets = alignparse.targets.Targets(seqsfile=config['amplicons_EG5'],
                                     feature_parse_specs=config['feature_parse_specs_EG5'])
```

## Setup to parse barcodes
Read data frame with list of all barcode runs.


```python
# barcode runs with R1 files by semicolon string split
barcode_runs = (pd.read_csv(config['barcode_runs'])
                .assign(R1=lambda x: x['R1'].str.split('; '))
                )

display(HTML(barcode_runs.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>sample</th>
      <th>sample_type</th>
      <th>sort_bin</th>
      <th>concentration</th>
      <th>date</th>
      <th>number_cells</th>
      <th>R1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>1.0</td>
      <td>240125</td>
      <td>384974</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s1-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_01_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>1.0</td>
      <td>240125</td>
      <td>860724</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s1-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_01_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>1.0</td>
      <td>240125</td>
      <td>2596538</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s1-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_01_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>1.0</td>
      <td>240125</td>
      <td>6350060</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s1-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_02_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>2.0</td>
      <td>240125</td>
      <td>630205</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s2-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_02_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>2.0</td>
      <td>240125</td>
      <td>837133</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s2-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_02_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>2.0</td>
      <td>240125</td>
      <td>2854102</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s2-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_02_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>2.0</td>
      <td>240125</td>
      <td>5742632</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s2-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_03_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>3.0</td>
      <td>240125</td>
      <td>982344</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s3-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_03_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>3.0</td>
      <td>240125</td>
      <td>960420</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s3-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_03_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>3.0</td>
      <td>240125</td>
      <td>3274652</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s3-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_03_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>3.0</td>
      <td>240125</td>
      <td>4475613</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s3-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_04_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>4.0</td>
      <td>240125</td>
      <td>1173631</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s4-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_04_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>4.0</td>
      <td>240125</td>
      <td>1414259</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s4-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_04_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>4.0</td>
      <td>240125</td>
      <td>4475535</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s4-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_04_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>4.0</td>
      <td>240125</td>
      <td>2540383</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s4-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_05_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>5.0</td>
      <td>240125</td>
      <td>2329147</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s5-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_05_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>5.0</td>
      <td>240125</td>
      <td>3486729</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s5-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_05_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>5.0</td>
      <td>240125</td>
      <td>4271970</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s5-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_05_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>5.0</td>
      <td>240125</td>
      <td>629341</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s5-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_06_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>6.0</td>
      <td>240125</td>
      <td>5231946</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s6-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_06_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>6.0</td>
      <td>240125</td>
      <td>3850987</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s6-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_06_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>6.0</td>
      <td>240125</td>
      <td>1008171</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s6-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_06_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>6.0</td>
      <td>240125</td>
      <td>7913</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s6-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_07_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>7.0</td>
      <td>240125</td>
      <td>7613259</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s7-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_07_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>7.0</td>
      <td>240125</td>
      <td>2394536</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s7-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_07_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>7.0</td>
      <td>240125</td>
      <td>187898</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s7-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_07_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>7.0</td>
      <td>240125</td>
      <td>6682</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s7-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_08_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>8.0</td>
      <td>240125</td>
      <td>7710354</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s8-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_08_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>8.0</td>
      <td>240125</td>
      <td>1991881</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s8-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_08_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>8.0</td>
      <td>240125</td>
      <td>101047</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s8-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_08_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>8.0</td>
      <td>240125</td>
      <td>7211</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s8-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_09_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>9.0</td>
      <td>240125</td>
      <td>8171174</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s9-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_09_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>9.0</td>
      <td>240125</td>
      <td>1853609</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s9-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_09_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>9.0</td>
      <td>240125</td>
      <td>71179</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s9-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_09_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>9.0</td>
      <td>240125</td>
      <td>5544</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-s9-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_01_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>1.0</td>
      <td>240215</td>
      <td>1707353</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s1-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_01_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>1.0</td>
      <td>240215</td>
      <td>824849</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s1-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_01_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>1.0</td>
      <td>240215</td>
      <td>1655178</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s1-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_01_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>1.0</td>
      <td>240215</td>
      <td>5955302</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s1-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_02_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>2.0</td>
      <td>240215</td>
      <td>1715895</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s2-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_02_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>2.0</td>
      <td>240215</td>
      <td>773681</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp2-s2-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_02_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>2.0</td>
      <td>240215</td>
      <td>1705442</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp2-s2-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_02_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>2.0</td>
      <td>240215</td>
      <td>5568746</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp2-s2-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_03_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>3.0</td>
      <td>240215</td>
      <td>2327674</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp2-s3-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_03_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>3.0</td>
      <td>240215</td>
      <td>933060</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp2-s3-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_03_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>3.0</td>
      <td>240215</td>
      <td>2028246</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s3-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_03_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>3.0</td>
      <td>240215</td>
      <td>4450399</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s3-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_04_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>4.0</td>
      <td>240215</td>
      <td>3474208</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s4-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_04_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>4.0</td>
      <td>240215</td>
      <td>1673592</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s4-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_04_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>4.0</td>
      <td>240215</td>
      <td>2862696</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s4-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_04_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>4.0</td>
      <td>240215</td>
      <td>1862063</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s4-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_05_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>5.0</td>
      <td>240215</td>
      <td>4981328</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s5-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_05_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>5.0</td>
      <td>240215</td>
      <td>2488502</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s5-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_05_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>5.0</td>
      <td>240215</td>
      <td>2089100</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s5-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_05_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>5.0</td>
      <td>240215</td>
      <td>535959</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s5-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_06_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>6.0</td>
      <td>240215</td>
      <td>7238765</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s6-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_06_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>6.0</td>
      <td>240215</td>
      <td>1981461</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s6-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_06_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>6.0</td>
      <td>240215</td>
      <td>639508</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s6-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_06_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>6.0</td>
      <td>240215</td>
      <td>13393</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s6-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_07_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>7.0</td>
      <td>240215</td>
      <td>7196977</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s7-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_07_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>7.0</td>
      <td>240215</td>
      <td>1669131</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s7-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_07_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>7.0</td>
      <td>240215</td>
      <td>303510</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s7-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_07_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>7.0</td>
      <td>240215</td>
      <td>9907</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s7-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_08_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>8.0</td>
      <td>240215</td>
      <td>7996419</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s8-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_08_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>8.0</td>
      <td>240215</td>
      <td>1905421</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s8-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_08_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>8.0</td>
      <td>240215</td>
      <td>320064</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s8-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_08_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>8.0</td>
      <td>240215</td>
      <td>12281</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s8-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_09_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>9.0</td>
      <td>240215</td>
      <td>7552576</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s9-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_09_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>9.0</td>
      <td>240215</td>
      <td>1475069</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s9-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_09_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>9.0</td>
      <td>240215</td>
      <td>196742</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s9-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_09_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>9.0</td>
      <td>240215</td>
      <td>8068</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240215-exp3-s9-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>SortSeq_bin1</td>
      <td>SortSeq</td>
      <td>1</td>
      <td>NaN</td>
      <td>240125</td>
      <td>3581850</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-RBD-bin1-1_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-RBD-bin1-2_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-RBD-bin1-3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>SortSeq_bin2</td>
      <td>SortSeq</td>
      <td>2</td>
      <td>NaN</td>
      <td>240125</td>
      <td>4649830</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-RBD-bin2-1_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-RBD-bin2-2_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-RBD-bin2-3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>SortSeq_bin3</td>
      <td>SortSeq</td>
      <td>3</td>
      <td>NaN</td>
      <td>240125</td>
      <td>8346717</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-RBD-bin3-1_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-RBD-bin3-2_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-RBD-bin3-3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>SortSeq_bin4</td>
      <td>SortSeq</td>
      <td>4</td>
      <td>NaN</td>
      <td>240125</td>
      <td>7841439</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-RBD-bin4-1_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-RBD-bin4-2_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240125-exp1-RBD-bin4-3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>SortSeq_bin1</td>
      <td>SortSeq</td>
      <td>1</td>
      <td>NaN</td>
      <td>240201</td>
      <td>2142177</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240201-exp2-RBD-bin1-1_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240201-exp2-RBD-bin1-2_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240201-exp2-RBD-bin1-3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>SortSeq_bin2</td>
      <td>SortSeq</td>
      <td>2</td>
      <td>NaN</td>
      <td>240201</td>
      <td>4158809</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240201-exp2-RBD-bin2-1_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240201-exp2-RBD-bin2-2_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240201-exp2-RBD-bin2-3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>SortSeq_bin3</td>
      <td>SortSeq</td>
      <td>3</td>
      <td>NaN</td>
      <td>240201</td>
      <td>5801807</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240201-exp2-RBD-bin3-1_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240201-exp2-RBD-bin3-2_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240201-exp2-RBD-bin3-3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>SortSeq_bin4</td>
      <td>SortSeq</td>
      <td>4</td>
      <td>NaN</td>
      <td>240201</td>
      <td>6361778</td>
      <td>[/uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240201-exp2-RBD-bin4-1_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240201-exp2-RBD-bin4-2_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/u6042467/starr-group1/sequencing/TNS/2024/240313_Illumina_OmiBA286-EG5-FLip_animal-polish_sarbeco-combo-polish/240201-exp2-RBD-bin4-3_R1_001.fastq.gz]</td>
    </tr>
  </tbody>
</table>


Make sure library / sample combinations are unique:


```python
assert len(barcode_runs) == len(barcode_runs.groupby(['library', 'sample']))
```

Make sure the the libraries for which we have barcode runs are all in our variant table:


```python
unknown_libs = set(barcode_runs['library']) - set(variants['library'])
if unknown_libs:
    raise ValueError(f"Libraries with barcode runs not in variant table: {unknown_libs}")
```

Now we initialize an [IlluminaBarcodeParser](https://jbloomlab.github.io/dms_variants/dms_variants.illuminabarcodeparser.html#dms_variants.illuminabarcodeparser.IlluminaBarcodeParser) for each library.

First, get the length of the barcode from the alignment target after making sure the same length for all targets:


```python
bclen = len(targets.targets[0].get_feature('barcode').seq)

assert (bclen == len(target.get_feature('barcode').seq) for target in targets.targets)

print(f"Barcodes of length {bclen}")
```

    Barcodes of length 16


The other barcode parsing params come from the config file:


```python
parser_params = config['illumina_barcode_parser_params']

display(HTML(
    pd.Series(parser_params, name='value')
    .rename_axis(index='parameter')
    .reset_index()
    .to_html(index=False)
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>parameter</th>
      <th>value</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>upstream</td>
      <td>GGCCGC</td>
    </tr>
    <tr>
      <td>downstream</td>
      <td></td>
    </tr>
    <tr>
      <td>minq</td>
      <td>20</td>
    </tr>
    <tr>
      <td>upstream_mismatch</td>
      <td>1</td>
    </tr>
    <tr>
      <td>downstream_mismatch</td>
      <td>0</td>
    </tr>
  </tbody>
</table>


The parser needs to know the set of valid barcodes, which are stored in the variant table and are different for each library.
So we create a different parser for each library using these valid barcode sets:


```python
# create dict keyed by library, value is parser for library
parsers = {lib: dms_variants.illuminabarcodeparser.IlluminaBarcodeParser(
                    bclen=bclen,
                    valid_barcodes=variants.loc[variants['library']==lib]['barcode'],
                    **parser_params)
           for lib in set(variants['library'])}

print('Number of valid barcodes searched for by each parser:')
display(HTML(
    pd.DataFrame([(lib, len(p.valid_barcodes)) for lib, p in parsers.items()],
                 columns=['library', 'number of valid barcodes'])
    .to_html(index=False)
    ))
```

    Number of valid barcodes searched for by each parser:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>number of valid barcodes</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool1</td>
      <td>259862</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>252698</td>
    </tr>
  </tbody>
</table>


## Parse barcodes
We now parse the barcodes.
Since this will take a while, we utilize multiple CPUs via the Python [multiprocessing](https://docs.python.org/3.6/library/multiprocessing.html) module.
First, determine how many CPUs to use.
We use the minimum of the user-specified number hardcoded below and the number actually available.
(If you are running *interactively* on the Hutch cluster, you may need to reduce the number below in order to avoid an error as there is an enforced CPU limit on the home `rhino` nodes):


```python
ncpus = min(config['max_cpus'], multiprocessing.cpu_count())
print(f"Using {ncpus} CPUs")
```

    Using 8 CPUs


Parse the barcodes in parallel via a [multiprocessing.Pool](https://docs.python.org/3.6/library/multiprocessing.html#multiprocessing.pool.Pool) using all the available CPUs to get a list of the data frames with barcode counts / fates for each sample:


```python
def process_func(parser, r1files, library, sample):
    """Convenience function to be starmapped to multiprocessing pool."""
    return parser.parse(r1files, add_cols={'library': library, 'sample': sample})

# parallel computation of list of data frames
with multiprocessing.pool.Pool(processes=ncpus) as pool:
    bclist = pool.starmap(
                process_func,
                [(parsers[run.library], run.R1, run.library, run.sample)
                  for run in barcode_runs.itertuples()],
                )
```

Now concatenate the list into data frames of barcode counts and barcode fates:


```python
counts = pd.concat([samplecounts for samplecounts, _ in bclist],
                   sort=False,
                   ignore_index=True)

print('First few lines of counts data frame:')
display(HTML(counts.head().to_html(index=False)))

fates = pd.concat([samplefates for _, samplefates in bclist],
                  sort=False,
                  ignore_index=True)

print('First few lines of fates data frame:')
display(HTML(fates.head().to_html(index=False)))
```

    First few lines of counts data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>barcode</th>
      <th>count</th>
      <th>library</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>CACTTATCACTAAAAC</td>
      <td>2524</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>CTCGTTTTACGCAATA</td>
      <td>2363</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>ACCGTTTATGCAGACA</td>
      <td>2270</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>TTTACATGACCATAAC</td>
      <td>2104</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>CGCGATTCTACCTGCA</td>
      <td>2032</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
  </tbody>
</table>


    First few lines of fates data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>fate</th>
      <th>count</th>
      <th>library</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>valid barcode</td>
      <td>3435720</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>invalid barcode</td>
      <td>932236</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>low quality barcode</td>
      <td>380724</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>unparseable barcode</td>
      <td>63304</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>failed chastity filter</td>
      <td>0</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
  </tbody>
</table>


## Examine fates of parsed barcodes
First, we'll analyze the "fates" of the parsed barcodes.
These fates represent what happened to each Illumina read we parsed:
 - Did the barcode read fail the Illumina chastity filter?
 - Was the barcode *unparseable* (i.e., the read didn't appear to be a valid barcode based on flanking regions)?
 - Was the barcode sequence too *low quality* based on the Illumina quality scores?
 - Was the barcode parseable but *invalid* (i.e., not in our list of variant-associated barcodes in the codon variant table)?
 - Was the barcode *valid*, and so will be added to variant counts.
 
First, we just write a CSV file with all the barcode fates:


```python
fatesfile = os.path.join(config['counts_dir'], 'barcode_fates.csv')
print(f"Writing barcode fates to {fatesfile}")
fates.to_csv(fatesfile, index=False)
```

    Writing barcode fates to results/counts/barcode_fates.csv


Next, we tabulate the barcode fates in wide format:


```python
display(HTML(fates
             .pivot_table(columns='fate',
                          values='count',
                          index=['library', 'sample'])
             .to_html()
             ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>fate</th>
      <th>failed chastity filter</th>
      <th>invalid barcode</th>
      <th>low quality barcode</th>
      <th>unparseable barcode</th>
      <th>valid barcode</th>
    </tr>
    <tr>
      <th>library</th>
      <th>sample</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="40" valign="top">pool1</th>
      <th>SortSeq_bin1</th>
      <td>0</td>
      <td>106369</td>
      <td>43208</td>
      <td>18984</td>
      <td>370773</td>
    </tr>
    <tr>
      <th>SortSeq_bin2</th>
      <td>0</td>
      <td>21480</td>
      <td>9804</td>
      <td>2039</td>
      <td>76914</td>
    </tr>
    <tr>
      <th>SortSeq_bin3</th>
      <td>0</td>
      <td>202116</td>
      <td>91533</td>
      <td>13769</td>
      <td>764498</td>
    </tr>
    <tr>
      <th>SortSeq_bin4</th>
      <td>0</td>
      <td>36203</td>
      <td>16526</td>
      <td>2518</td>
      <td>143736</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin1</th>
      <td>0</td>
      <td>932236</td>
      <td>380724</td>
      <td>63304</td>
      <td>3435720</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin2</th>
      <td>0</td>
      <td>258514</td>
      <td>118957</td>
      <td>18065</td>
      <td>979816</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin3</th>
      <td>0</td>
      <td>435455</td>
      <td>193086</td>
      <td>29438</td>
      <td>1634775</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin4</th>
      <td>0</td>
      <td>1729589</td>
      <td>866367</td>
      <td>112796</td>
      <td>6729150</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin1</th>
      <td>0</td>
      <td>155753</td>
      <td>64753</td>
      <td>10508</td>
      <td>581241</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin2</th>
      <td>0</td>
      <td>2212805</td>
      <td>997852</td>
      <td>153628</td>
      <td>8384239</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin3</th>
      <td>0</td>
      <td>536611</td>
      <td>203960</td>
      <td>34644</td>
      <td>1983587</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin4</th>
      <td>0</td>
      <td>1539434</td>
      <td>750389</td>
      <td>110386</td>
      <td>6010924</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin1</th>
      <td>0</td>
      <td>255350</td>
      <td>104328</td>
      <td>17418</td>
      <td>944405</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin2</th>
      <td>0</td>
      <td>255543</td>
      <td>115467</td>
      <td>17320</td>
      <td>939382</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin3</th>
      <td>0</td>
      <td>753982</td>
      <td>276825</td>
      <td>46398</td>
      <td>2754005</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin4</th>
      <td>0</td>
      <td>1249695</td>
      <td>620163</td>
      <td>83485</td>
      <td>4929456</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin1</th>
      <td>0</td>
      <td>2295307</td>
      <td>995472</td>
      <td>160113</td>
      <td>8667919</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin2</th>
      <td>0</td>
      <td>205515</td>
      <td>97441</td>
      <td>13644</td>
      <td>762249</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin3</th>
      <td>0</td>
      <td>1206838</td>
      <td>563330</td>
      <td>76811</td>
      <td>4465849</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin4</th>
      <td>0</td>
      <td>957788</td>
      <td>527557</td>
      <td>70215</td>
      <td>4060415</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin1</th>
      <td>0</td>
      <td>527818</td>
      <td>258782</td>
      <td>36528</td>
      <td>1958719</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin2</th>
      <td>0</td>
      <td>2855945</td>
      <td>1283748</td>
      <td>164682</td>
      <td>10055383</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin3</th>
      <td>0</td>
      <td>850329</td>
      <td>440048</td>
      <td>57870</td>
      <td>3281451</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin4</th>
      <td>0</td>
      <td>98355</td>
      <td>43805</td>
      <td>9907</td>
      <td>485327</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin1</th>
      <td>0</td>
      <td>1242148</td>
      <td>596066</td>
      <td>80669</td>
      <td>4514819</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin2</th>
      <td>0</td>
      <td>973709</td>
      <td>443879</td>
      <td>62286</td>
      <td>3587030</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin3</th>
      <td>0</td>
      <td>283564</td>
      <td>159000</td>
      <td>23788</td>
      <td>1319740</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin4</th>
      <td>0</td>
      <td>2880</td>
      <td>986</td>
      <td>167</td>
      <td>6984</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin1</th>
      <td>0</td>
      <td>2117651</td>
      <td>1078965</td>
      <td>135575</td>
      <td>7963230</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin2</th>
      <td>0</td>
      <td>6898425</td>
      <td>3680592</td>
      <td>480602</td>
      <td>27931300</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin3</th>
      <td>0</td>
      <td>55375</td>
      <td>27310</td>
      <td>5217</td>
      <td>279570</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin4</th>
      <td>0</td>
      <td>2040</td>
      <td>657</td>
      <td>127</td>
      <td>4840</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin1</th>
      <td>0</td>
      <td>2258903</td>
      <td>1055110</td>
      <td>139342</td>
      <td>8219675</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin2</th>
      <td>0</td>
      <td>1411193</td>
      <td>651434</td>
      <td>97979</td>
      <td>5179580</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin3</th>
      <td>0</td>
      <td>79986</td>
      <td>42341</td>
      <td>6068</td>
      <td>378872</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin4</th>
      <td>0</td>
      <td>2684</td>
      <td>753</td>
      <td>128</td>
      <td>5278</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin1</th>
      <td>0</td>
      <td>1581953</td>
      <td>737734</td>
      <td>109321</td>
      <td>5926033</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin2</th>
      <td>0</td>
      <td>224614</td>
      <td>109740</td>
      <td>16236</td>
      <td>935484</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin3</th>
      <td>0</td>
      <td>90735</td>
      <td>52599</td>
      <td>8892</td>
      <td>452880</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin4</th>
      <td>0</td>
      <td>2011</td>
      <td>631</td>
      <td>16279</td>
      <td>3914</td>
    </tr>
    <tr>
      <th rowspan="40" valign="top">pool2</th>
      <th>SortSeq_bin1</th>
      <td>0</td>
      <td>2178967</td>
      <td>962040</td>
      <td>135219</td>
      <td>7462507</td>
    </tr>
    <tr>
      <th>SortSeq_bin2</th>
      <td>0</td>
      <td>3831435</td>
      <td>2010386</td>
      <td>258920</td>
      <td>14958128</td>
    </tr>
    <tr>
      <th>SortSeq_bin3</th>
      <td>0</td>
      <td>1450711</td>
      <td>712898</td>
      <td>108031</td>
      <td>5588238</td>
    </tr>
    <tr>
      <th>SortSeq_bin4</th>
      <td>0</td>
      <td>1454523</td>
      <td>703851</td>
      <td>103342</td>
      <td>6381660</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin1</th>
      <td>0</td>
      <td>566744</td>
      <td>318702</td>
      <td>34919</td>
      <td>2179005</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin2</th>
      <td>0</td>
      <td>302890</td>
      <td>171084</td>
      <td>18819</td>
      <td>1176585</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin3</th>
      <td>0</td>
      <td>550151</td>
      <td>216921</td>
      <td>36521</td>
      <td>2111539</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin4</th>
      <td>0</td>
      <td>2149437</td>
      <td>1125487</td>
      <td>143371</td>
      <td>8933158</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin1</th>
      <td>0</td>
      <td>533530</td>
      <td>278551</td>
      <td>31043</td>
      <td>2011862</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin2</th>
      <td>0</td>
      <td>296052</td>
      <td>153517</td>
      <td>22751</td>
      <td>1139927</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin3</th>
      <td>0</td>
      <td>675647</td>
      <td>345938</td>
      <td>40043</td>
      <td>2554955</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin4</th>
      <td>0</td>
      <td>1981267</td>
      <td>1134244</td>
      <td>132746</td>
      <td>8306811</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin1</th>
      <td>0</td>
      <td>843752</td>
      <td>429173</td>
      <td>52605</td>
      <td>3228137</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin2</th>
      <td>0</td>
      <td>410576</td>
      <td>209283</td>
      <td>28470</td>
      <td>1581619</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin3</th>
      <td>0</td>
      <td>797528</td>
      <td>385553</td>
      <td>50637</td>
      <td>3072332</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin4</th>
      <td>0</td>
      <td>1600001</td>
      <td>885166</td>
      <td>108992</td>
      <td>6545689</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin1</th>
      <td>0</td>
      <td>1323216</td>
      <td>723717</td>
      <td>89837</td>
      <td>5100595</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin2</th>
      <td>0</td>
      <td>617699</td>
      <td>298692</td>
      <td>39575</td>
      <td>2368790</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin3</th>
      <td>0</td>
      <td>990071</td>
      <td>399833</td>
      <td>70820</td>
      <td>3895303</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin4</th>
      <td>0</td>
      <td>620106</td>
      <td>335851</td>
      <td>45148</td>
      <td>2737065</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin1</th>
      <td>0</td>
      <td>1641978</td>
      <td>974679</td>
      <td>103651</td>
      <td>6422605</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin2</th>
      <td>0</td>
      <td>790280</td>
      <td>435852</td>
      <td>58405</td>
      <td>3132331</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin3</th>
      <td>0</td>
      <td>788695</td>
      <td>414843</td>
      <td>51291</td>
      <td>3177216</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin4</th>
      <td>0</td>
      <td>185950</td>
      <td>117096</td>
      <td>13347</td>
      <td>896513</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin1</th>
      <td>0</td>
      <td>2773742</td>
      <td>1461010</td>
      <td>184701</td>
      <td>10913672</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin2</th>
      <td>0</td>
      <td>787205</td>
      <td>395784</td>
      <td>48923</td>
      <td>3281969</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin3</th>
      <td>0</td>
      <td>178959</td>
      <td>81451</td>
      <td>12271</td>
      <td>838993</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin4</th>
      <td>0</td>
      <td>28884</td>
      <td>8923</td>
      <td>7998</td>
      <td>46975</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin1</th>
      <td>0</td>
      <td>5356588</td>
      <td>3010060</td>
      <td>350364</td>
      <td>20966673</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin2</th>
      <td>0</td>
      <td>587004</td>
      <td>306564</td>
      <td>39204</td>
      <td>2520884</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin3</th>
      <td>0</td>
      <td>96521</td>
      <td>50338</td>
      <td>5879</td>
      <td>419086</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin4</th>
      <td>0</td>
      <td>9654</td>
      <td>1941</td>
      <td>165</td>
      <td>8625</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin1</th>
      <td>0</td>
      <td>3160370</td>
      <td>1692300</td>
      <td>207296</td>
      <td>12443737</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin2</th>
      <td>0</td>
      <td>700882</td>
      <td>431406</td>
      <td>45832</td>
      <td>2982071</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin3</th>
      <td>0</td>
      <td>105945</td>
      <td>65044</td>
      <td>6539</td>
      <td>442796</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin4</th>
      <td>0</td>
      <td>7738</td>
      <td>1620</td>
      <td>146</td>
      <td>7503</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin1</th>
      <td>0</td>
      <td>3004054</td>
      <td>1735398</td>
      <td>188001</td>
      <td>11562475</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin2</th>
      <td>0</td>
      <td>523773</td>
      <td>310533</td>
      <td>32730</td>
      <td>2167625</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin3</th>
      <td>0</td>
      <td>73273</td>
      <td>34267</td>
      <td>4236</td>
      <td>309089</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin4</th>
      <td>0</td>
      <td>12979</td>
      <td>2452</td>
      <td>236</td>
      <td>9697</td>
    </tr>
  </tbody>
</table>


Now we plot the barcode-read fates for each library / sample, showing the bars for valid barcodes in orange and the others in gray.
We see that the largest fraction of barcode reads correspond to valid barcodes, and most of the others are invalid barcodes (probably because the map to variants that aren't present in our variant table since we didn't associate all variants with barcodes). The exception to this is lib2 Titeseq_03_bin3; the PCR for this sample in the original sequencing run failed, so we followed it up with a single MiSeq lane. We did not filter out the PhiX reads from this data before parsing, so these PhiX reads will deflate the fraction of valid barcode reads as expected, but does not indicate any problems.


```python
barcode_fate_plot = (
    ggplot(
        fates
        .assign(sample=lambda x: pd.Categorical(x['sample'],
                                                x['sample'].unique(),
                                                ordered=True),
                fate=lambda x: pd.Categorical(x['fate'],
                                              x['fate'].unique(),
                                              ordered=True),
                is_valid=lambda x: x['fate'] == 'valid barcode'
                ), 
        aes('fate', 'count', fill='is_valid')) +
    geom_bar(stat='identity') +
    facet_grid('sample ~ library') +
    facet_grid('sample ~ library') +
    scale_fill_manual(CBPALETTE, guide=False) +
    theme(figure_size=(1.4 * (1 + fates['library'].nunique()),
                       1.7 * (1.2 + fates['sample'].nunique())),
          axis_text_x=element_text(angle=90),
          panel_grid_major_x=element_blank()
          ) +
    scale_y_continuous(labels=dms_variants.utils.latex_sci_not,
                       name='number of reads')
    )

_ = barcode_fate_plot.draw()
```


    
![png](count_variants_files/count_variants_42_0.png)
    


## Output csv of barcode counts in variant-barcode lookup table


```python
print(f"Writing variant counts to {config['variant_counts_file']}")
counts.to_csv(config['variant_counts_file'], index=False)
```

    Writing variant counts to results/counts/variant_counts.csv


The [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) has lots of nice functions that can be used to analyze the counts it contains.
However, we do that in the next notebook so we don't have to re-run this entire (rather computationally intensive) notebook every time we want to analyze a new aspect of the counts.


```python

```