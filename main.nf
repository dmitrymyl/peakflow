#!/usr/bin/env nextflow

// Runs MACS2 to estimate fragment size and build a peak model.
process getPeakModelFragmentSize {
    publishDir "${params.outdir}", mode: 'copy'
    label 'quick'  // This process is supposed to run very quickly

    input:
        path chip_bam   // ChIP BAM file
        val  effgsize   // Effective genome size for normalization
        val  format     // Format of the BAM files (i.e., 'BAM' for single-end or 'BAMPE' for paired-end)
        val  prefix     // Prefix for output files
    output:
        path   "${prefix}.fragment_size.txt"  // Output file with fragment size
        path   "${prefix}.model.pdf"          // Output PDF with the peak model plot
        stdout emit: fragment_size            // Fragment size value emitted to stdout
    script:
    """
    macs2 predictd -i $chip_bam -f $format -g $effgsize --rfile model.R
    parse_model_script.py model.R > ${prefix}.fragment_size.txt
    cat ${prefix}.fragment_size.txt
    Rscript model.R > /dev/null
    mv model.R_model.pdf ${prefix}.model.pdf
    """
}

// Runs MACS2 to call peaks
process callPeaks {
    publishDir "${params.outdir}", mode: 'copy'
    label 'peaks'  // This process is supposed to run for a longer time and consume more RAM

    input:
        path chip_bam   // ChIP BAM file
        path input_bam  // Input BAM file
        val  prefix     // Prefix for output files
        val  effgsize   // Effective genome size for normalization
        val  format     // Format of the BAM files (i.e., 'BAM' for single-end or 'BAMPE' for paired-end)
    output:
        path "${prefix}.peaks.narrowPeak"  // Output peaks file in narrowPeak format
    script:
    """
    macs2 callpeak -t $chip_bam -c $input_bam -g $effgsize -f $format -n chip
    mv chip_peaks.narrowPeak ${prefix}.peaks.narrowPeak
    """
}

// Filter peaks by blacklist
process filterPeaksByBlacklist {
    publishDir "${params.outdir}", mode: 'copy'
    label 'quick'  // This process is supposed to run very quickly

    input:
        path peaks       // Peaks file in narrowPeak format
        path blacklist   // Blacklist file to filter out regions
        val  prefix      // Prefix for output files
    output:
        path "${prefix}.peaks_noblacklist.narrowPeak"  // Output filtered peaks file
    script:
    """
    bedtools intersect -v -a $peaks -b $blacklist > ${prefix}.peaks_noblacklist.narrowPeak
    """
}

// Runs deeptools bamCoverage to get CPM bigwigs
process makeCpmTrack {
    publishDir "${params.outdir}", mode: 'copy'
    label 'tracks'  // This process is supposed to be parallelized

    input:
        path bam        // BAM file to be converted to bigwig
        path bai        // BAM index file
        val  kind       // Type of the BAM file (i.e., 'chip' or 'input')
        val  prefix     // Prefix for output files
        path blacklist  // Blacklist file to filter out regions
        val  binsize    // Bin size for bigwig files
        val  fragsize   // Fragment size for read extension
        val  effgsize   // Effective genome size for normalization
    output:
        path "${prefix}.${kind}.bw"  // Output bigwig file with CPM normalization
    script:
    """
    bamCoverage -b $bam \
                -o ${prefix}.${kind}.bw \
                -p $task.cpus \
                --normalizeUsing CPM \
                --blackListFileName $blacklist \
                --binSize $binsize \
                --extendReads $fragsize \
                --ignoreDuplicates \
                --effectiveGenomeSize $effgsize
    """
}

// Runs deeptools bamCompare to get log2 ratio of ChIP and Input CPMs
process makeRatioTrack {
    publishDir "${params.outdir}", mode: 'copy'
    label 'tracks'  // This process is supposed to be parallelized

    input:
        path bam_chip   // ChIP BAM file
        path bai_chip   // ChIP BAM index file
        path bam_input  // Input BAM file
        path bai_input  // Input BAM index file
        val  prefix     // Prefix for output files
        path blacklist  // Blacklist file
        val  binsize    // Bin size for bigwig files
        val  fragsize   // Fragment size for read extension
        val  effgsize   // Effective genome size for normalization

    output:
        path "${prefix}.ratio.bw"  // Output bigwig file with log2 ratio

    script:
    """
    bamCompare -b1 $bam_chip \
               -b2 $bam_input \
               -o ${prefix}.ratio.bw \
               -p $task.cpus \
               --normalizeUsing CPM \
               --operation log2 \
               --scaleFactorsMethod None \
               --blackListFileName $blacklist \
               --binSize $binsize \
               --extendReads $fragsize \
               --ignoreDuplicates \
               --effectiveGenomeSize $effgsize
    """
}

params.samplesheet = './samplesheet.csv'               // Path to the samplesheet
params.blacklist   = './assets/hg19-blacklist.v2.bed'  // Path to the blacklist file
params.outdir      = "./results"                       // Path to the output folder
params.prefix      = 'sample'                          // Prefix for output files
params.binsize     = 1000                              // Bin size for bigwig files
params.callpeaks   = true                              // Whether to call peaks or not
params.maketracks  = true                              // Whether to make bigwig tracks or not
params.extreads    = true                              // Whether to extend reads based on fragment size
params.effgsize    = 2736124898                        // Effective genome size for normalization (hg19 by default)
params.pairedend   = false                             // Whether the input BAM files are paired-end or single-end

workflow {
    ch_bams   = Channel.fromPath( params.samplesheet ).splitCsv( header:true )       // Parse samplesheet
    bam_chip  = ch_bams.filter { it.type == 'chip' }.map { row -> file(row.path) }   // Separate channel for ChIP BAM
    bam_input = ch_bams.filter { it.type == 'input' }.map { row -> file(row.path) }  // Separate channel for Input BAM
    blacklist = Channel.fromPath(params.blacklist)                                   // Blacklist file

    if (params.pairedend) {
        bam_format = "BAMPE"
    }
    else {
        bam_format = "BAM"
    }
    
    if (params.extreads) {
        getPeakModelFragmentSize(
                         bam_chip,
                         params.effgsize,
                         bam_format,
                         params.prefix
                         )
        fragment_size = getPeakModelFragmentSize.out.fragment_size.first()
    }
    else {
        fragment_size = Channel.value(0) // No read extension
    }

    if (params.callpeaks) {
        callPeaks(
                 bam_chip,
                 bam_input,
                 params.prefix,
                 params.effgsize,
                 bam_format
                 )
        filterPeaksByBlacklist(
                              callPeaks.out,
                              blacklist.first(),
                              params.prefix
                              )
    }

    if (params.maketracks) {

        makeCpmTrack(
                    ch_bams.map {row -> file(row.path)},
                    ch_bams.map {row -> file("${row.path}.bai")},
                    ch_bams.map { row -> row.type },
                    params.prefix,
                    blacklist.first(),
                    params.binsize,
                    fragment_size,
                    params.effgsize
                    )

        makeRatioTrack(
                      bam_chip,
                      bam_chip.map{"${it}.bai"},
                      bam_input,
                      bam_input.map{"${it}.bai"},
                      params.prefix,
                      blacklist.first(),
                      params.binsize,
                      fragment_size,
                      params.effgsize
                      )
    }
}