# Uploading FASTQ files to the SRA

Details of how the raw PacBio sequencing files were uploaded to the NIH [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra).
The submitted files are in BioProject [PRJNA770094](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA770094).

The Python Jupyter notebooks [upload_to_SRA_BA286.ipynb](upload_to_SRA_BA286.ipynb), [upload_to_SRA_EG5.ipynb](upload_to_SRA_EG5.ipynb)  and [upload_to_SRA_FLip.ipynb](upload_to_SRA_FLip.ipynb) has instructions and does the uploading.

Because the FTP upload takes a while, you may want to run the Jupyter notebook using `slurm` so there is no timeout with::

    sbatch --wrap="jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 upload_to_SRA_BA286.ipynb" --time 2-0

And equivalent for EG5 and FLip