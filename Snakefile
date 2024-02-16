"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------
import glob
import itertools
import os.path
import os
import textwrap
import urllib.request

import Bio.SeqIO

import dms_variants.codonvarianttable
import dms_variants.illuminabarcodeparser

import pandas as pd

# Configuration  --------------------------------------------------------------
configfile: 'config.yaml'

# run "quick" rules locally:
localrules: make_dag,
            make_summary,
            save_pinned_env,
            get_ccs,
            get_mut_antibody_escape,
            get_mut_bind_expr,
            get_SARSr_DMS

# Functions -------------------------------------------------------------------
def nb_markdown(nb):
    """Return path to Markdown results of notebook `nb`."""
    return os.path.join(config['summary_dir'],
                        os.path.basename(os.path.splitext(nb)[0]) + '.md')

# Global variables extracted from config --------------------------------------
pacbio_runs = (pd.read_csv(config['pacbio_runs'], dtype = str)
               .assign(pacbioRun=lambda x: x['library'] + '_' + x['bg'] + '_' + x['run'])
               )
assert len(pacbio_runs['pacbioRun'].unique()) == len(pacbio_runs['pacbioRun'])

# Information on samples and barcode runs -------------------------------------
barcode_runs = pd.read_csv(config['barcode_runs'])

# Rules -----------------------------------------------------------------------

# making this summary is the target rule (in place of `all`) since it
# is first rule listed.
rule make_summary:
    """Create Markdown summary of analysis."""
    input:
        dag=os.path.join(config['summary_dir'], 'dag.svg'),
        env='environment_pinned.yml',
        get_mut_bind_expr=config['mut_bind_expr'],
        get_mut_antibody_escape=config['mut_antibody_escape'],
        get_SARSr_DMS=config['SARSr_DMS'],
        process_ccs_EG5=nb_markdown('process_ccs_EG5.ipynb'),
        process_ccs_FLip=nb_markdown('process_ccs_FLip.ipynb'),
        process_ccs_BA286=nb_markdown('process_ccs_BA286.ipynb'),
        barcode_variant_table_EG5=config['codon_variant_table_file_EG5'],
        barcode_variant_table_FLip=config['codon_variant_table_file_FLip'],
        barcode_variant_table_BA286=config['codon_variant_table_file_BA286'],
        variant_counts_file=config['variant_counts_file'],
        count_variants=nb_markdown('count_variants.ipynb'),
        fit_titrations='results/summary/compute_binding_Kd.md',
        variant_Kds_file=config['Titeseq_Kds_file'],
        calculate_expression='results/summary/compute_expression_meanF.md',
        variant_expression_file=config['expression_sortseq_file'],
        collapse_scores='results/summary/collapse_scores.md',
        mut_phenos_file=config['final_variant_scores_mut_file'],
        epistatic_shifts='results/summary/epistatic_shifts.md',
        epistasis_viz=os.path.join(config['visualization_dir'], "epistasis.html"),
        heatmap_viz=os.path.join(config['visualization_dir'], "heatmap.html"),
        gisaid_rbd_mutations=nb_markdown('gisaid_rbd_mutations.ipynb'),
        gisaid_mutation_counts=config['gisaid_mutation_counts'],
    output:
        summary = os.path.join(config['summary_dir'], 'summary.md')
    run:
        def path(f):
            """Get path relative to `summary_dir`."""
            return os.path.relpath(f, config['summary_dir'])
        with open(output.summary, 'w') as f:
            f.write(textwrap.dedent(f"""
            # Summary

            Analysis run by [Snakefile]({path(workflow.snakefile)})
            using [this config file]({path(workflow.configfiles[0])}).
            See the [README in the top directory]({path('README.md')})
            for details.

            Here is the DAG of the computational workflow:
            ![{path(input.dag)}]({path(input.dag)})

            Here is the Markdown output of each Jupyter notebook in the
            workflow:

            1. Get prior RBD DMS mutation-level binding and expression data from [prior DMS study](https://github.com/tstarrlab/SARS-CoV-2-RBD_DMS_Omicron-XBB-BQ), and SARS-CoV-1 ACE2-binding DMS data from [this repo](https://github.com/tstarrlab/SARSr-CoV-RBD_DMS).
            
            2. Process PacBio CCSs for 
            [Omicron EG.5]({path(input.process_ccs_EG5)}), 
            [Omicron HK.3 (i.e. FLip)]({path(input.process_ccs_FLip)}), and 
            [Omicron BA.2.86]({path(input.process_ccs_BA286)}). Creates barcode-variant lookup tables for each background: [EG.5]({path(input.barcode_variant_table_EG5)}), [FLip]({path(input.barcode_variant_table_FLip)}), and [BA.2.86]({path(input.barcode_variant_table_BA286)}).
            
            3. [Count variants by barcode]({path(input.count_variants)}).
               Creates a [variant counts file]({path(input.variant_counts_file)})
               giving counts of each barcoded variant in each condition.

            4. [Fit titration curves]({path(input.fit_titrations)}) to calculate per-barcode K<sub>D</sub>, recorded in [this file]({path(input.variant_Kds_file)}).
            
            5. [Analyze Sort-seq]({path(input.calculate_expression)}) to calculate per-barcode RBD expression, recorded in [this file]({path(input.variant_expression_file)}).
            
            6. [Derive final genotype-level phenotypes from replicate barcoded sequences]({path(input.collapse_scores)}).
               Generates final phenotypes, recorded in [this file]({path(input.mut_phenos_file)}).
            
            7. [Count mutations in GISAID RBD sequences]({path(input.gisaid_rbd_mutations)}) to create [this counts file]({path(input.gisaid_mutation_counts)}).
            
            8. [Analyze patterns of epistasis in the DMS data and in SARS-CoV-2 genomic data]({path(input.epistatic_shifts)}).
            
            9. Make interactive data visualizations, available [here](https://tstarrlab.github.io/SARS-CoV-2-RBD_DMS_Omicron-XBB-BQ/)

            """
            ).strip())

rule make_dag:
    # error message, but works: https://github.com/sequana/sequana/issues/115
    input:
        workflow.snakefile
    output:
        os.path.join(config['summary_dir'], 'dag.svg')
    shell:
        "snakemake --forceall --dag | dot -Tsvg > {output}"

rule save_pinned_env:
    input:
        dag=os.path.join(config['summary_dir'], 'dag.svg'),
    log:
    	"environment_pinned.yml"
    shell:
        """
        conda env export > {log}
        """

rule interactive_jsd_plot:
    """ Make the interactive plot for visualizing epistatic shifts.
    """
    input: 
        scores=config['final_variant_scores_mut_file'],
        jsd=config['JSD_file']
    output:
        html=os.path.join(config['visualization_dir'], "epistasis.html")
    notebook: "Epistatic-Shifts-Interactive-Visualization.ipynb"


rule interactive_heatmap_plot:
    """ Make the interactive heatmap for expression and binding.
    """
    input: 
        scores=config['final_variant_scores_mut_file']
    params:
        annotations=config['RBD_sites']
    output:
        html=os.path.join(config['visualization_dir'], "heatmap.html")
    notebook: "RBD-Heatmaps-Interactive-Visualization.ipynb"


rule epistatic_shifts:
    input:
        config['final_variant_scores_mut_file'],
        config['mut_antibody_escape'],
        config['SARSr_DMS']
    output:
        config['JSD_file'],
        config['JSD_expr_file'],
        md='results/summary/epistatic_shifts.md',
        md_files=directory('results/summary/epistatic_shifts_files')
    envmodules:
        'R/4.1.3'
    params:
        nb='epistatic_shifts.Rmd',
        md='epistatic_shifts.md',
        md_files='epistatic_shifts_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule gisaid_rbd_mutations:
    input:
        config['gisaid_spikes'],
        config['wildtype_sequence'],
    output:
        config['gisaid_mutation_counts'],
        nb_markdown=nb_markdown('gisaid_rbd_mutations.ipynb'),
    params:
        nb='gisaid_rbd_mutations.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule collapse_scores:
    input:
        config['Titeseq_Kds_file'],
        config['expression_sortseq_file'],
        config['mut_bind_expr'],
        config['gisaid_mutation_counts'],
    output:
        config['final_variant_scores_mut_file'],
        md='results/summary/collapse_scores.md',
        md_files=directory('results/summary/collapse_scores_files')
    envmodules:
        'R/4.1.3'
    params:
        nb='collapse_scores.Rmd',
        md='collapse_scores.md',
        md_files='collapse_scores_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule fit_titrations:
    input:
        config['codon_variant_table_file_EG5'],
        config['codon_variant_table_file_FLip'],
        config['codon_variant_table_file_BA286'],
        config['variant_counts_file']
    output:
        config['Titeseq_Kds_file'],
        md='results/summary/compute_binding_Kd.md',
        md_files=directory('results/summary/compute_binding_Kd_files')
    envmodules:
        'R/4.1.3'
    params:
        nb='compute_binding_Kd.Rmd',
        md='compute_binding_Kd.md',
        md_files='compute_binding_Kd_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule calculate_expression:
    input:
        config['codon_variant_table_file_EG5'],
        config['codon_variant_table_file_FLip'],
        config['codon_variant_table_file_BA286'],
        config['variant_counts_file']
    output:
        config['expression_sortseq_file'],
        md='results/summary/compute_expression_meanF.md',
        md_files=directory('results/summary/compute_expression_meanF_files')
    envmodules:
        'R/4.1.3'
    params:
        nb='compute_expression_meanF.Rmd',
        md='compute_expression_meanF.md',
        md_files='compute_expression_meanF_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule count_variants:
    """Count codon variants from Illumina barcode runs."""
    input:
        config['codon_variant_table_file_EG5'],
        config['codon_variant_table_file_FLip'],
        config['codon_variant_table_file_BA286'],
        config['barcode_runs']
    output:
        config['variant_counts_file'],
        nb_markdown=nb_markdown('count_variants.ipynb')
    params:
        nb='count_variants.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"


rule get_SARSr_DMS:
    """Download SARS-related CoV mutation ACE2-binding and expression from URL."""
    output:
        file=config['SARSr_DMS']
    run:
        urllib.request.urlretrieve(config['SARSr_DMS_url'], output.file)
        
rule get_mut_bind_expr:
    """Download SARS-CoV-2 mutation ACE2-binding and expression from URL."""
    output:
        file=config['mut_bind_expr']
    run:
        urllib.request.urlretrieve(config['mut_bind_expr_url'], output.file)
        
rule get_mut_antibody_escape:
    """Download SARS-CoV-2 mutation antibody-escape data from URL."""
    output:
        file=config['mut_antibody_escape']
    run:
        urllib.request.urlretrieve(config['mut_antibody_escape_url'], output.file)


rule process_ccs_BA286:
    """Process the PacBio CCSs for BA286 background and build variant table."""
    input:
        expand(os.path.join(config['ccs_dir'], "{pacbioRun}_ccs.fastq.gz"),
               pacbioRun=pacbio_runs['pacbioRun']),
    output:
        config['processed_ccs_file_BA286'],
    	config['codon_variant_table_file_BA286'],
        nb_markdown=nb_markdown('process_ccs_BA286.ipynb')
    params:
        nb='process_ccs_BA286.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
rule process_ccs_FLip:
    """Process the PacBio CCSs for FLip background and build variant table."""
    input:
        expand(os.path.join(config['ccs_dir'], "{pacbioRun}_ccs.fastq.gz"),
               pacbioRun=pacbio_runs['pacbioRun']),
    output:
        config['processed_ccs_file_FLip'],
    	config['codon_variant_table_file_FLip'],
        nb_markdown=nb_markdown('process_ccs_FLip.ipynb')
    params:
        nb='process_ccs_FLip.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule process_ccs_EG5:
    """Process the PacBio CCSs for EG5 background and build variant table."""
    input:
        expand(os.path.join(config['ccs_dir'], "{pacbioRun}_ccs.fastq.gz"),
               pacbioRun=pacbio_runs['pacbioRun']),
    output:
        config['processed_ccs_file_EG5'],
    	config['codon_variant_table_file_EG5'],
        nb_markdown=nb_markdown('process_ccs_EG5.ipynb')
    params:
        nb='process_ccs_EG5.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

if config['seqdata_source'] == 'local':

    rule get_ccs:
        """Symbolically link CCS files."""
        input:
            ccs_fastq=lambda wildcards: (pacbio_runs
                                        .set_index('pacbioRun')
                                        .at[wildcards.pacbioRun, 'ccs']
                                        )
        output:
            ccs_fastq=os.path.join(config['ccs_dir'], "{pacbioRun}_ccs.fastq.gz")
        run:
            os.symlink(input.ccs_fastq, output.ccs_fastq)

elif config['seqdata_source'] == 'SRA':
    raise RuntimeError('getting sequence data from SRA not yet implemented')

else:
    raise ValueError(f"invalid `seqdata_source` {config['seqdata_source']}")
