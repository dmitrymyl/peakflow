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

// Parses R script from MACS2 to obtain the fragment size
process getPeakSize {
    input:
        path model_script
    output:
        stdout
    script:
    """
    parse_model_script.py $model_script
    """
}

// Runs R script from MACS2 to create pdf with the peak model
process plotPeakModel {
    publishDir 'results', mode: 'copy'
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
    publishDir 'results', mode: 'copy'
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
    publishDir 'results', mode: 'copy'
    conda './conda.yml'

    input:
        path bam
        path bai
        path blacklist
        val binsize
        val fragsize
        val kind
    output:
        path "${kind}_track.bw"
    script:
    """
    bamCoverage -b $bam -o ${kind}_track.bw -p 2 --normalizeUsing CPM --blackListFileName $blacklist --binSize $binsize --extendReads $fragsize --ignoreDuplicates
    """
}

// Runs deeptools bamCompare to get log2 ratio of ChIP and Input CPMs
process makeRatioTrack {
    publishDir 'results', mode: 'copy'
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

params.bamchip = './assets/sample.chip.chr20.bam'
params.baminput = './assets/sample.input.chr20.bam'
params.blacklist = './assets/hg19-blacklist.v2.bed'
params.binsize = 1000
params.callpeaks = true

workflow {
    bam_chip = Channel.fromPath(params.bamchip)
    bam_input = Channel.fromPath(params.baminput)
    blacklist = Channel.fromPath(params.blacklist)
    makePeakModel(bam_chip)
    plotPeakModel(makePeakModel.out)
    getPeakSize(makePeakModel.out)

    if (params.callpeaks) {
        callPeaks(bam_chip, bam_input)
    }

    makeCpmTrack(bam_chip, bam_chip.map{"${it}.bai"}, blacklist, params.binsize, getPeakSize.out, 'chip')
    // makeCpmTrack(bam_input, bam_input.map{"${it}.bai"}, blacklist, params.binsize, getPeakSize.out, 'input')
    makeRatioTrack(bam_chip, bam_chip.map{"${it}.bai"}, bam_input, bam_input.map{"${it}.bai"}, blacklist, params.binsize, getPeakSize.out)

/*
    makeCpmTrack(bam_files)
    makeRatioTrack(bam_files)
*/

}