#!/usr/bin/env bash

# Script to run PeakFlow on the CBE cluster using Nextflow
# Load necessary modules and set environment variables
ml build-env/f2022
module --ignore-cache load "nextflow/23.10.1"
export NXF_OPTS='-Xms1g -Xmx4g'
export NXF_TEMP="/scratch/nextflow/"$USER"/nxf_temp"
export NXF_WORK="/scratch/nextflow/"$USER"/nxf_work"
export NXF_ANSI_LOG=false

# Set notification email for IMBA people
NOTIFICATION_EMAIL=${USER}@imba.oeaw.ac.at
# Path to the parameters file -- modify it as needed
PARAMS_FILE="./params.json"

# Run the peakflow pipeline
nextflow -bg run dmitrymyl/peakflow -r main -profile cbe -params-file $PARAMS_FILE -N $NOTIFICATION_EMAIL
