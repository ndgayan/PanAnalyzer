# PanAnalyzer - A tool to analyze and visualize pan-genome data
# Developed by Gayan Nagoda
# Copyright (C) 2025 Gayan Nagoda


import os

import utility as app_utility

########################################################################################################################
# Options
########################################################################################################################

CONDA_PATH = os.path.expanduser("~/miniconda3/bin/conda")
CONDA_ENV = "GEM"

SAMPLES_PREFIX = "Study-All.csv"  # Name of the study file

SPADES_OUTPUT = "./OUTPUT/SPAdes_Results"
SAMPLE_FORWARD_READS_POSTFIX = "R1_001.trim.fastq.gz"
SAMPLE_REVERSE_READS_POSTFIX = "R2_001.trim.fastq.gz"

REFERENCE_FILE_EXTENSION = ".fna"

PIPE_SPADES = False
PIPE_ANVI_O = True

########################################################################################################################

if __name__ == "__main__":
    # Load sample information from Study-All.csv
    samples = app_utility.read_study_file(f"./DATA/{SAMPLES_PREFIX}")
    print(f"\nTotal samples: {len(samples)}")
    ####################################################################################################################
    # STEP 1 : SPeads Assembly
    ####################################################################################################################
    if PIPE_SPADES:
        # Clean the output directory before starting
        app_utility.clean_output_directory(SPADES_OUTPUT)

        # Validate Sample files (paired-end reads)
        try:
            # Get all files in the Samples directory
            sample_files = app_utility.get_files_from_directory("./DATA/Samples")
            print(f"\nTotal files in Samples directory: {len(sample_files)}")

            # Validate that all samples have required forward and reverse read files
            validated_samples = app_utility.validate_sample_files(
                samples_dict=samples,
                sample_files=sample_files,
                forward_postfix=SAMPLE_FORWARD_READS_POSTFIX,
                reverse_postfix=SAMPLE_REVERSE_READS_POSTFIX,
                filter_type="Sample",
            )
            print(f"\nValidated samples with required files: {len(validated_samples)}")

        except FileNotFoundError:
            print(
                "\n⚠️  PROCESS STOPPED: Cannot proceed without all required sample files"
            )
            exit(1)

        for idx, (sample_id, sample_data) in enumerate(validated_samples.items(), 1):
            r1_file = sample_data["forward"]
            r2_file = sample_data["reverse"]
            output_dir = os.path.join(SPADES_OUTPUT, sample_id)
            print(
                f"Sample {idx}: {sample_id} - Running SPAdes assembly using:\n{r1_file}\n{r2_file}"
            )

            command = [
                CONDA_PATH,
                "run",
                "-n",
                CONDA_ENV,
                "spades.py",
                "--isolate",
                "-o",
                output_dir,
                "-1",
                r1_file,
                "-2",
                r2_file,
            ]

            result = app_utility.bash_execute(command)
            print(f"Success: {result['success']}")
            if not result["success"]:
                print(f"Error running SPAdes for sample {sample_id}: {result['error']}")
                exit(1)

        print("\n✅ SPAdes assembly completed for all samples.")

    ####################################################################################################################
    # STEP 2 : Anvi'o Database Creation
    ####################################################################################################################

    if PIPE_ANVI_O:
        # Load reference genome information to the pan genome Analysis
        # Validate Reference files
        try:
            # Get all files in the Reference directory
            reference_files = app_utility.get_files_from_directory("./DATA/Ref")
            print(f"\nTotal files in Reference directory: {len(reference_files)}")

            # Validate that all reference samples have required files
            validated_references = app_utility.validate_reference_files(
                samples_dict=samples,
                reference_files=reference_files,
                file_extension=REFERENCE_FILE_EXTENSION,
                filter_type="Reference",
            )
            print(
                f"\nValidated references with required files: {len(validated_references)}"
            )

        except FileNotFoundError:
            print(
                "\n⚠️  PROCESS STOPPED: Cannot proceed without all required reference files"
            )
            exit(1)
