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
REFERENCE_FILE_EXTENSION = ".fna"

TEMP_OUTPUT = "./OUTPUT/TEMP"

SPADES_OUTPUT = "./OUTPUT/SPAdes_Results"
SAMPLE_FORWARD_READS_POSTFIX = "R1_001.trim.fastq.gz"
SAMPLE_REVERSE_READS_POSTFIX = "R2_001.trim.fastq.gz"

ANVIO_PROJECT_NAME = "PanAnalyzer"
ANVIO_GENOMES_DB = f"{ANVIO_PROJECT_NAME}-Genomes.db"
ANVIO_OUTPUT = "./OUTPUT/Anvio_Results"

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
        app_utility.clean_output_directory(ANVIO_OUTPUT)
        app_utility.clean_output_directory(TEMP_OUTPUT)
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

        try:
            validated_contigs = app_utility.validate_spades_output(
                samples_dict=samples,
                spades_output_dir=SPADES_OUTPUT,
                filter_type="Sample",
            )
        except FileNotFoundError:
            print(
                "\n⚠️  PROCESS STOPPED: Cannot proceed without all contigs.fasta files"
            )
            exit(1)

        combined_genomes = {}
        for sample_id, contig_path in validated_contigs.items():
            combined_genomes[sample_id] = contig_path
        for ref_id, ref_data in validated_references.items():
            combined_genomes[ref_id] = ref_data["file"]

        for genome_id, genome_file in combined_genomes.items():
            db_file_name = ANVIO_OUTPUT + f"/{genome_id}.db"

            # NCBI Reference genomes need to be reformatted
            temp_genome_file = None
            if genome_file.endswith(REFERENCE_FILE_EXTENSION):
                temp_genome_file = TEMP_OUTPUT + f"/{genome_id}.fna"
                command = [
                    CONDA_PATH,
                    "run",
                    "-n",
                    CONDA_ENV,
                    "anvi-script-reformat-fasta",
                    genome_file,
                    "-o",
                    temp_genome_file,
                    "-l",
                    "0",
                    "--simplify-names",
                ]
                result = app_utility.bash_execute(command)
                if not result["success"]:
                    print(
                        f"Error creating contigs DB for sample {sample_id}: {result['error']}"
                    )
                    exit(1)
                genome_file = temp_genome_file

            # Create contigs database for each genome
            command = [
                CONDA_PATH,
                "run",
                "-n",
                CONDA_ENV,
                "anvi-gen-contigs-database",
                "-f",
                genome_file,
                "-o",
                db_file_name,
                "-n",
                genome_id,
                "--project-name",
                ANVIO_PROJECT_NAME,
            ]
            result = app_utility.bash_execute(command)
            if not result["success"]:
                print(
                    f"Error creating contigs DB for sample {sample_id}: {result['error']}"
                )
                exit(1)
            if temp_genome_file and os.path.exists(temp_genome_file):
                os.remove(temp_genome_file)

            # Run HMMs on the contigs database
            command = [
                CONDA_PATH,
                "run",
                "-n",
                CONDA_ENV,
                "anvi-run-hmms",
                "-c",
                db_file_name,
                "--num-threads",
                "8",
            ]
            result = app_utility.bash_execute(command)
            if not result["success"]:
                print(f"Error running HMMs for sample {sample_id}: {result['error']}")
                exit(1)

            # Annotate genes with COGs
            command = [
                CONDA_PATH,
                "run",
                "-n",
                CONDA_ENV,
                "anvi-run-ncbi-cogs",
                "-c",
                db_file_name,
                "--num-threads",
                "8",
            ]
            result = app_utility.bash_execute(command)
            if not result["success"]:
                print(
                    f"Error running Annotations for sample {sample_id}: {result['error']}"
                )
                exit(1)

        print("\n✅ Anvi'o contigs databases created for all genomes.")
