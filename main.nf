#!/usr/bin/env nextflow


// Runs MACS2 to build the peak model and predict fragment size
process makePeakModel {
//    container 'oras://community.wave.seqera.io/library/deeptools_macs2_python_r-base:8454a235341a6609'
    conda './conda.yml'
    
    input:
        path chip_bam
    output:
        path 'model.R'
    script:
    """
    macs2 predictd -i $chip_bam -f BAM -g hs --rfile model.R
    """
}

// Parses R script from MACS2 to obtain the fragment size, saves it to a file and emits to stdout
process getFragmentSize {
    publishDir "${params.outdir}", mode: 'copy'
    conda './conda.yml'

    input:
        path model_script
    output:
        path 'fragment_size.txt'
        stdout emit: fragment_size
    script:
    """
    parse_model_script.py $model_script > fragment_size.txt
    cat fragment_size.txt
    """
}

// Runs R script from MACS2 to create pdf with the peak model
process plotPeakModel {
    publishDir "${params.outdir}", mode: 'copy'
    conda './conda.yml'
    
    input:
        path model_script
    output:
        path 'model.pdf'
    script:
    """
    Rscript $model_script
    mv ${model_script}_model.pdf model.pdf
    """
}

// Runs MACS2 to call peaks
process callPeaks {
    publishDir "${params.outdir}", mode: 'copy'
    conda './conda.yml'

    input:
        path chip_bam
        path input_bam
    output:
        path 'peaks.narrowPeak'
    script:
    """
    macs2 callpeak -t $chip_bam -c $input_bam -g hs -f BAM -n chip
    mv chip_peaks.narrowPeak peaks.narrowPeak
    """
}

// Runs deeptools bamCoverage to get CPM bigwigs
process makeCpmTrack {
    publishDir "${params.outdir}", mode: 'copy'
    conda './conda.yml'

    input:
        path bam
        path bai
        val kind
        path blacklist
        val binsize
        val fragsize
    output:
        path "${kind}_track.bw"
    script:
    """
    bamCoverage -b $bam -o ${kind}_track.bw -p 2 --normalizeUsing CPM --blackListFileName $blacklist --binSize $binsize --extendReads $fragsize --ignoreDuplicates
    """
}

// Runs deeptools bamCompare to get log2 ratio of ChIP and Input CPMs
process makeRatioTrack {
    publishDir "${params.outdir}", mode: 'copy'
    conda './conda.yml'

    input:
        path bam_chip
        path bai_chip
        path bam_input
        path bai_input
        path blacklist
        val binsize
        val fragsize
    output:
        path "ratio_track.bw"
    script:
    """
    bamCompare -b1 $bam_chip -b2 $bam_input -o ratio_track.bw -p 2 --normalizeUsing CPM --operation log2 --scaleFactorsMethod None --blackListFileName $blacklist --binSize $binsize --extendReads $fragsize --ignoreDuplicates
    """
}

params.samplesheet = './samplesheet.csv'
params.blacklist = './assets/hg19-blacklist.v2.bed'
params.binsize = 1000
params.outdir = "./results"

workflow {
    ch_bams = Channel.fromPath('./samplesheet.csv').splitCsv(header:true)
    bam_chip = ch_bams.filter { it.type == 'chip' }.map { row -> file(row.path) }
    bam_input = ch_bams.filter { it.type == 'input' }.map { row -> file(row.path) }
    blacklist = Channel.fromPath(params.blacklist)
    
    makePeakModel(bam_chip)
    plotPeakModel(makePeakModel.out)
    getFragmentSize(makePeakModel.out)
    callPeaks(bam_chip, bam_input)
    makeCpmTrack(ch_bams.map {row -> file(row.path)}, ch_bams.map {row -> file("${row.path}.bai")}, ch_bams.map { row -> row.type }, blacklist.first(), params.binsize, getFragmentSize.out.fragment_size.first())
    makeRatioTrack(bam_chip, bam_chip.map{"${it}.bai"}, bam_input, bam_input.map{"${it}.bai"}, blacklist.first(), params.binsize, getFragmentSize.out.fragment_size.first())
}