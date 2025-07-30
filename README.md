# peakflow
A simple nextflow pipeline for calling peaks and producing tracks from ChIP-seq data.

## Description
Given a single ChIP and a single Input BAM files, `peakflow` estimates fragment size and calls peaks with MACS2 and produces coverage tracks for ChIP, Input, and ChIP/Input ratio with deepTools. Coverage tracks are binned according to a specified bin size and have reads extended to an estimated fragment size. BAM files are treated as single-end ones, even if they are paired-end.

## Input
The pipeline requires two kinds of inputs: read alignment files (BAM/BAI) and parameters. Paths to read alignment files are specified in `samplesheet.csv`.

### Read alignment files and `samplesheet.csv`
`peakflow` requires indexed BAM files for ChIP and Input samples with BAI index files in the same directory. `samplesheet.csv` specifies paths only to BAM files in the following format:

```
type,path
chip,/path/to/chip.bam
input,/path/to/input.bam

```

`samplesheet.csv` is a comma-delimited CSV file with a header and two data rows. The header consists of two columns: type and path. Type can be either "chip" or "input". Rows specify paths to BAM files with the corresponding types.

### Parameters
The pipeline has the following parameters:

| Name        | Type            | Default value                      | Description                                        |
|-------------|-----------------|------------------------------------|----------------------------------------------------|
| samplesheet | file path (csv) | `"./samplesheet.csv"`              | Path to the samplesheet                            |
| blacklist   | file path (bed) | `"./assets/hg19-blacklist.v2.bed"` | Path to the blacklist file                         |
| outdir      | file path       | `null`                             | Path to the output directory                       |
| prefix      | string          | `"sample"`                         | Filename prefix to start output files with         |
| binsize     | integer         | `1000`                             | Size of bins for coverage tracks                   |
| callpeaks   | boolean         | `true`                             | Whether to call peaks or not                       |
| extreads    | boolean         | `true`                             | Whether to extend reads for coverage tracks or not |

Parameters are supplied as a json file. An example can be found at `params.json`.

## How to run

### Basic run
```{bash}
nextflow run dmitrymyl/peakflow -profile PROFILE -params-file params.json
```
Available `PROFILE` values are `local` and those available form [nf-core](https://nf-co.re/configs/).

### How to run on CLIP CBE
Activate nextflow:
```
ml build-env/f2022
module --ignore-cache load "nextflow/23.10.1"
export NXF_OPTS='-Xms1g -Xmx4g'
export NXF_TEMP="/scratch/nextflow/"$USER"/nxf_temp"
export NXF_WORK="/scratch/nextflow/"$USER"/nxf_work"
export NXF_ANSI_LOG=false
```
Run the pipeline with:
```
nextflow -bg run dmitrymyl/peakflow -r main -profile cbe -params-file params.json
```
Profile `cbe` makes the pipeline use slurm for job submission and resource management.

### Useful options
In command line, nextflow options are specified with the single hyphen (such as `-profile` or `-params`), while the workflow parameters are specified with the double hyphen (such as `--samplesheet` and all the rest available in `params.json`). For example, instead of supplying `params.json` you can specify necessary workflow parameters directly on the command line:
```
nextflow -bg run dmitrymyl/peakflow -r main -profile cbe --samplesheet samplesheet.csv --blacklist blacklist.bed --outdir results --prefix sample
```
`-bg` option allows to run nextflow in the background. `-resume` option allows to skip successful steps in case of a rerun.

## Get the latest version of the pipeline

```
nextflow pull dmitrymyl/peakflow
```

## Output
The output directory has the following files:

1. `prefix.chip_track.bw` ChIP CPM coverage track.
2. `prefix.input_track.bw` Input CPM coverage track.
3. `prefix.ratio_track.bw` log2 ChIP/Input CPM ratio track.
4. `prefix.peaks.narrowPeak` MACS2 peaks.
5. `prefix.model.pdf` peak model plot produced by MACS2. 
6. `prefix.fragment_size.txt` estimated fragment size as a single number.

All bigwig tracks are binned and have reads extended to the fragment size.

## Software dependencies

There are two ways to manage dependencies in this pipeline: a conda environment and a container. The container is the recommended option to run the pipeline.

### Apptainer/Singularity container

The Apptainer/Singularity container is available at `oras://docker.io/gerlichlab/peakflow-apptainer:latest`. See `Apptainer` or `Dockerfile` definition file to build a container yourself. When running, nextflow automatically downloads the container.

### Conda environment
conda environment specification file is `conda.yml`. When running, nextflow creates a conda environment based on this specification

## Development

### Hosting singularity containers on dockerhub

Two step process: login and push.
```
apptainer remote login --username dh-user oras://docker.io
apptainer push container.sif oras://docker.io/dh-user/container:tag
```