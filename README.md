# FECBUD reference-free engraftment analysis

Computational analysis workflow and results summary of the reference-free engraftment analysis of the FECBUD study.

---

## Short description

As a way to estimate engraftment from the Netherlands Donor Feces Bank (NDFB) donors
to ulcerative colitis (UC) patients, we are calculating whole-metagenome-based
Jaccard distances between samples using [Mash](https://github.com/marbl/Mash).

This workflow consists of four main parts:

1. Human-derived read removal (by [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
 and [SAMtools](https://github.com/samtools/samtools))  
2. Quality trimming and adapter removal (by [fastp](https://github.com/OpenGene/fastp))  
3. K-mer based Jaccard distance calculations between metagenomes (by [Mash](https://github.com/marbl/Mash))  
4. Visualizations and statistical calculations (using [R](https://www.r-project.org/))

The results of the human read removal and quality processing (together 'pre-processing', steps 1 and 2)
are summarized and visualized using [MultiQC](https://multiqc.info/),
of which the report can be found in [`doc/Preprocessing_report.html`](doc/Preprocessing_report.html).

The results of step 4 can be found in [this document](results/Engraftment_analysis-Mash.html).

(_Please note that the above linked HTML reports usually work best_ 
_when downloaded to your computer and then opened in a webbrowser._)

All the fastq processing is done with a Snakemake workflow.
The resulting tables with distances are imported into an RMarkdown script
for further analyses.

For details on program versions and parameters used, see the configuration files
in `envs/` and `config/`, and the `Snakefile`.

---

## Quick-start

The workflow is based on [Snakemake](https://snakemake.readthedocs.io/en/stable/) 
and uses [conda](https://docs.conda.io/en/latest/miniconda.html) 
environments to automatically install the required software tools.
Please make sure you have both of these installed before running the workflow.

For example, after installing conda:

```bash
conda create -n snakemake snakemake=5.10.0
conda activate snakemake
```

When using the default `config/config.yaml` file, you need to first create the directory `log/SLURM`:

```bash
mkdir -p log/SLURM
```

And the file `config/parameters.yaml` needs to be adjusted to point to the right input files.
For example, this YAML file:

```yaml
samples:
  - sampleA
  - sampleB
  - sambleC

raw_directory: "data/raw/"
```

Makes it so that the pipeline assumes your input files are:  
  - `data/raw/sampleA-trimmed.fastq.gz`  
  - `data/raw/sampleB-trimmed.fastq.gz`  
  - `data/raw/sampleC-trimmed.fastq.gz`

Furthermore, the workflow requires a bowtie2-index of the human genome to map the raw reads to to exclude human-derived reads.

```yaml
host_database: "/path/to/human/genome/"
```

There should be bowtie2-index files in this directory with the file prefix `reference`.
(See the rule 'host_removal' in the [Snakefile](Snakefile).)

The default configuration file assumes the workflow is being run a an HPC cluster
using the SLURM job scheduler.

The `slurm-cluster-status` project needs to be cloned from GitHub into the `bin` directory:

```bash
cd bin
git clone https://github.com/LUMC/slurm-cluster-status
cd ..
```

If you are not running on a SLURM cluster, please edit these lines in
`config/config.yaml` accordingly (for instance by removing them if you run on your local machine):

```yaml
cluster: "sbatch --parsable -N 1 -n 1 -c {threads} --mem={resources.mem} -t {resources.time} -D . -e log/SLURM/{name}-{jobid}.err -o log/SLURM/{name}-{jobid}.out"
cluster-status: bin/slurm-cluster-status/slurm-cluster-status.py
```

With these preparations ready, you should be able to run the workflow.

```bash
snakemake --profile config -n # First try a dry-run to check if the configuration is right

snakemake --profile config # If that was successful, run the workflow
```


---

## Project organisation

```
.
├── .gitignore
├── CITATION.cff
├── LICENSE
├── README.md
├── bin                <- Code and programs used in the project
├── config             <- Configuration files necessary fow the workflow
├── data               <- All project data, divided in subfolders
│   ├── processed      <- Final data, used for visualisation (e.g. tables)
│   ├── raw            <- Raw data, original, should not be modified (e.g. fastq files)
│   └── tmp            <- Intermediate data, derived from the raw data, but not yet ready for visualisation
├── doc                <- Documents such as quality control report
├── envs               <- Conda environments necessary to run the project
├── log                <- Log files from programs
└── results            <- Figures or reports generated from processed data
```

---

## License

This project is licensed under the terms of the [fill in your licence here](LICENSE).

---

## Citation

Please cite this project as described [here](CITATION.cff).
