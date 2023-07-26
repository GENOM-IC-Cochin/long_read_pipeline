# Long Read Pipeline
This pipeline performs Differential Transcript Usage (DTU) according to Love et al.'s [workflow](https://doi.org/10.12688/f1000research.15398.3), which uses [DEXSeq (Anders et al.)](https://doi.org/10.1101/gr.133744.111) and [stageR (Van den Verge et al.)](https://doi.org/10.1186/s13059-017-1277-0), following alignment with [minimap2 (Li)](https://doi.org/10.1093/bioinformatics/btab705), alignment quality control with [RNA-SeQC (Graubert et al.)](https://doi.org/10.1093/bioinformatics/btab135) and quantification with [Isoquant (Prjibelski et al.)](https://doi.org/10.1038/s41587-022-01565-y).

## Usage
- Conda is needed to set up the pipeline's environment, visit [its website](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) for installation. [Mamba](mamba.readthedocs.io/) can speed up solving the environment and is straightforward to install.
- Set up the environment : `conda|mamba env create -f ont_hiv.yaml`
- Activate the environment : `conda|mamba activate ont_hiv`
- Fill the `config.yaml` file. It requires paths to a fastq directory, a reference genome (in the case of dual RNA-Seq, host and pathogen, the pathogen reference needs to be concactenated as a chromosome in the host genome), to a gtf file (containing the host and pathogen genome if no filtering is chosen) to a samplesheet, to a contrast csv file, and various options and output paths.
- Install renv ```r install.packages("renv")```
- [Optional test step] Try to restore the renv lock file (in the R console ```r renv::restore()```). Some packages have system requirements that might prevent the statistical analysis from running if not installed.
- Execute the pipeline with ```snakemake -c [number_of_threads]```.

