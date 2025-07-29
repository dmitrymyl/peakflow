#!/usr/bin/env nextflow


// Runs MACS2 to build the peak model and predict fragment size
process makePeakModel {
    label 'quick'
    
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
    label 'quick'

    input:
        path model_script
        val prefix
    output:
        path "${prefix}.fragment_size.txt"
        stdout emit: fragment_size
    script:
    """
    parse_model_script.py $model_script > ${prefix}.fragment_size.txt
    cat ${prefix}.fragment_size.txt
    """
}

// Runs R script from MACS2 to create pdf with the peak model
process plotPeakModel {
    publishDir "${params.outdir}", mode: 'copy'
    label 'quick'
    
    input:
        path model_script
        val prefix
    output:
        path "${prefix}.model.pdf"
    script:
    """
    Rscript $model_script
    mv ${model_script}_model.pdf ${prefix}.model.pdf
    """
}

// Runs MACS2 to call peaks
process callPeaks {
    publishDir "${params.outdir}", mode: 'copy'
    label 'peaks'

    input:
        path chip_bam
        path input_bam
        val prefix
    output:
        path "${prefix}.peaks.narrowPeak"
    script:
    """
    macs2 callpeak -t $chip_bam -c $input_bam -g hs -f BAM -n chip
    mv chip_peaks.narrowPeak ${prefix}.peaks.narrowPeak
    """
}

// Runs deeptools bamCoverage to get CPM bigwigs
process makeCpmTrack {
    publishDir "${params.outdir}", mode: 'copy'
    label 'tracks'

    input:
        path bam
        path bai
        val kind
        val prefix
        path blacklist
        val binsize
        val fragsize
    output:
        path "${prefix}.${kind}.bw"
    script:
    """
    bamCoverage -b $bam -o ${prefix}.${kind}.bw -p $task.cpus --normalizeUsing CPM --blackListFileName $blacklist --binSize $binsize --extendReads $fragsize --ignoreDuplicates
    """
}

// Runs deeptools bamCompare to get log2 ratio of ChIP and Input CPMs
process makeRatioTrack {
    publishDir "${params.outdir}", mode: 'copy'
    label 'tracks'

    input:
        path bam_chip
        path bai_chip
        path bam_input
        path bai_input
        val prefix
        path blacklist
        val binsize
        val fragsize
    output:
        path "${prefix}.ratio.bw"
    script:
    """
    bamCompare -b1 $bam_chip -b2 $bam_input -o ${prefix}.ratio.bw -p $task.cpus --normalizeUsing CPM --operation log2 --scaleFactorsMethod None --blackListFileName $blacklist --binSize $binsize --extendReads $fragsize --ignoreDuplicates
    """
}

params.samplesheet = './samplesheet.csv'
params.blacklist = './assets/hg19-blacklist.v2.bed'
params.outdir = "./results"
params.prefix = 'sample'
params.binsize = 1000
params.callpeaks = true
params.extreads = true

workflow {
    ch_bams = Channel.fromPath(params.samplesheet).splitCsv(header:true)
    bam_chip = ch_bams.filter { it.type == 'chip' }.map { row -> file(row.path) }
    bam_input = ch_bams.filter { it.type == 'input' }.map { row -> file(row.path) }
    blacklist = Channel.fromPath(params.blacklist)
    
    if (params.extreads) {
        makePeakModel(bam_chip)
        plotPeakModel(makePeakModel.out, params.prefix)
        getFragmentSize(makePeakModel.out, params.prefix)
        fragment_size = getFragmentSize.out.fragment_size.first()
    }
    else {
        fragment_size = Channel.value(0) // No read extension
    }

    if (params.callpeaks) {
        callPeaks(bam_chip, bam_input, params.prefix)
    }

    makeCpmTrack(ch_bams.map {row -> file(row.path)}, ch_bams.map {row -> file("${row.path}.bai")}, ch_bams.map { row -> row.type }, params.prefix, blacklist.first(), params.binsize, fragment_size)
    makeRatioTrack(bam_chip, bam_chip.map{"${it}.bai"}, bam_input, bam_input.map{"${it}.bai"}, params.prefix, blacklist.first(), params.binsize, fragment_size)
}