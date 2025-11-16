########################################################################################################################
# PanAnalyzer - A tool to analyze and visualize pan-genome data
# Author: Gayan Nagoda
# Date: 2025-11-08
########################################################################################################################

import os

import utility as app_utility

########################################################################################################################
########################################################################################################################
# Options
########################################################################################################################
########################################################################################################################

THREADS = "8"
MEMORY = "8192"  # in MB

CONDA_PATH = os.path.expanduser("~/miniconda3/bin/conda")
CONDA_ENV = "GEM"

SAMPLES_PREFIX = "Study-All.csv"  # Name of the study file
REFERENCE_FILE_EXTENSION = ".fna"

SAMPLE_FORWARD_READS_POSTFIX = "R1_001.trim.fastq.gz"
SAMPLE_REVERSE_READS_POSTFIX = "R2_001.trim.fastq.gz"

ANVIO_PROJECT_NAME = "PanAnalyzer"
ANVIO_GENOMES_DB = f"{ANVIO_PROJECT_NAME}-GENOMES.db"

# What part of the study to run.
PIPE_FASTQC = False
PIPE_SPADES = False

PIPE_ANVI_O = False
PIPE_RELATIONSHIPS = False

PIPE_PROKKA = False
PIPE_ROARY = True

########################################################################################################################
########################################################################################################################
########################################################################################################################

# STOP: Do not change these output directory paths. Create them if they do not exist.
TEMP_OUTPUT = "./OUTPUT/TEMP"
FASTQC_OUTPUT = "./OUTPUT/FastQC_Results"
SPADES_OUTPUT = "./OUTPUT/SPAdes_Results"
ANVIO_OUTPUT = "./OUTPUT/Anvio_Results"
PROKKA_OUTPUT = "./OUTPUT/Prokka_Results"
ROARY_OUTPUT = "./OUTPUT/Roary_Results"


if __name__ == "__main__":
    # Load sample information from Study-All.csv
    samples = app_utility.read_study_file(f"./DATA/{SAMPLES_PREFIX}")
    print(f"\nTotal samples: {len(samples)}")

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
        print("\n⚠️  PROCESS STOPPED: Cannot proceed without all required sample files")
        exit(1)

    ####################################################################################################################
    # STEP 0 : FastQC Quality Check (Optional)
    ####################################################################################################################
    if PIPE_FASTQC:
        # Clean the output directory before starting
        app_utility.clean_output_directory(FASTQC_OUTPUT)

        for idx, (sample_id, sample_data) in enumerate(validated_samples.items(), 1):
            r1_file = sample_data["forward"]
            r2_file = sample_data["reverse"]

            for sample_file in [r1_file, r2_file]:
                # Run FastQC for the sample -> https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
                command = [
                    CONDA_PATH,
                    "run",
                    "-n",
                    CONDA_ENV,
                    "fastqc",
                    sample_file,
                    "-o",
                    FASTQC_OUTPUT,
                    "--threads",
                    THREADS,
                    "--memory",
                    MEMORY,
                ]
                result = app_utility.bash_execute(command)
                print(f"Success: {result['success']}")
                if not result["success"]:
                    print(
                        f"Error running FastQC for sample {sample_file}: {result['error'].output}"
                    )
                    exit(1)
                else:
                    print(
                        f"Sample {idx}: {sample_id} - FastQC completed for {sample_file}"
                    )

        print("\n✅ FastQC quality check completed for all samples.")

        # Open-source tool to aggregate bioinformatic analyses results into a single report.
        # Run MultiQC for the sample -> https://seqera.io/multiqc/
        command = [
            CONDA_PATH,
            "run",
            "-n",
            CONDA_ENV,
            "multiqc",
            FASTQC_OUTPUT,
            "-o",
            os.path.join(FASTQC_OUTPUT, "MultiQC_Results"),
            "--clean-up",
        ]
        result = app_utility.bash_execute(command)
        print(f"Success: {result['success']}")
        if not result["success"]:
            print(f"Error running MultiQC for FastQC Samples: {result['error'].output}")
            exit(1)

        print("\n✅ MultiQC report generated for FastQC results.")

    ####################################################################################################################
    # STEP 1 : SPAdes Assembly
    ####################################################################################################################
    if PIPE_SPADES:
        # Clean the output directory before starting
        app_utility.clean_output_directory(SPADES_OUTPUT)

        for idx, (sample_id, sample_data) in enumerate(validated_samples.items(), 1):
            r1_file = sample_data["forward"]
            r2_file = sample_data["reverse"]
            output_dir = os.path.join(SPADES_OUTPUT, sample_id)
            print(
                f"Sample {idx}: {sample_id} - Running SPAdes assembly using:\n{r1_file}\n{r2_file}"
            )

            # Run SPAdes assembly for the sample -> https://ablab.github.io/spades/getting-started.html
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
                print(
                    f"Error running SPAdes for sample {sample_id}: {result['error'].output}"
                )
                exit(1)
            else:
                print(f"Sample {idx}: {sample_id} - SPAdes assembly completed.")

        print("\n✅ SPAdes assembly completed for all samples.")

    ####################################################################################################################
    # STEP 2 : Anvi'o Pipeline
    ####################################################################################################################

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

    # Validate SPAdes output contig files
    try:
        validated_contigs = app_utility.validate_spades_output(
            samples_dict=samples,
            spades_output_dir=SPADES_OUTPUT,
            filter_type="Sample",
        )
    except FileNotFoundError:
        print("\n⚠️  PROCESS STOPPED: Cannot proceed without all contigs.fasta files")
        exit(1)

    # Combine sample contigs and reference genomes
    combined_genomes = {}
    for sample_id, contig_path in validated_contigs.items():
        combined_genomes[sample_id] = contig_path
    for ref_id, ref_data in validated_references.items():
        combined_genomes[ref_id] = ref_data["file"]

    if PIPE_ANVI_O:
        # Clean the output directory before starting
        app_utility.clean_output_directory(ANVIO_OUTPUT)
        app_utility.clean_output_directory(TEMP_OUTPUT)

        # Create Anvi'o contigs databases for all genomes (samples + references)
        for genome_id, genome_file in combined_genomes.items():
            db_file_name = os.path.join(ANVIO_OUTPUT, f"{genome_id}.db")

            # NCBI Reference genomes need to be reformatted -> https://anvio.org/help/main/programs/anvi-script-reformat-fasta/
            temp_genome_file = None
            if genome_file.endswith(REFERENCE_FILE_EXTENSION):
                temp_genome_file = os.path.join(TEMP_OUTPUT, f"{genome_id}.fna")
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
                        f"Error creating contigs DB for sample {sample_id}: {result['error'].output}"
                    )
                    exit(1)
                genome_file = temp_genome_file

            # Create contigs database for the genome -> https://anvio.org/help/main/programs/anvi-gen-contigs-database/
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
                    f"Error creating contigs DB for sample {sample_id}: {result['error'].output}"
                )
                exit(1)
            if temp_genome_file and os.path.exists(temp_genome_file):
                os.remove(temp_genome_file)

            # Run HMMs on the contigs database -> https://anvio.org/help/7/programs/anvi-run-hmms/
            command = [
                CONDA_PATH,
                "run",
                "-n",
                CONDA_ENV,
                "anvi-run-hmms",
                "-c",
                db_file_name,
                "--num-threads",
                THREADS,
            ]
            result = app_utility.bash_execute(command)
            if not result["success"]:
                print(
                    f"Error running HMMs for sample {sample_id}: {result['error'].output}"
                )
                exit(1)

            # Annotate genes with COGs -> https://anvio.org/help/main/programs/anvi-run-ncbi-cogs/
            command = [
                CONDA_PATH,
                "run",
                "-n",
                CONDA_ENV,
                "anvi-run-ncbi-cogs",
                "-c",
                db_file_name,
                "--num-threads",
                THREADS,
            ]
            result = app_utility.bash_execute(command)
            if not result["success"]:
                print(
                    f"Error running Annotations for sample {sample_id}: {result['error'].output}"
                )
                exit(1)

            print(f"Contigs database created for genome: {genome_id}")

        print("\n✅ Anvi'o contigs databases created for all genomes.")

        # Create external-genomes.txt file.
        try:
            external_genomes_file = app_utility.generate_external_genomes_file(
                combined_genomes=combined_genomes,
                anvio_output_dir=ANVIO_OUTPUT,
                output_file="external-genomes.txt",
            )
            print(f"\nExternal genomes file created: {external_genomes_file}")

        except FileNotFoundError:
            print("\n⚠️  PROCESS STOPPED: Cannot proceed without all database files")
            exit(1)

        # Create genome storage -> https://anvio.org/help/7/programs/anvi-gen-genomes-storage/
        command = [
            CONDA_PATH,
            "run",
            "-n",
            CONDA_ENV,
            "anvi-gen-genomes-storage",
            "-e",
            external_genomes_file,
            "-o",
            os.path.join(ANVIO_OUTPUT, ANVIO_GENOMES_DB),
        ]
        result = app_utility.bash_execute(command)
        if not result["success"]:
            print(f"Error creating genomes storage: {result['error'].output}")
            exit(1)

        print(f"\n✅ Genomes storage database created: {ANVIO_GENOMES_DB}")

        # Run Pan Genome Analysis -> https://anvio.org/help/main/programs/anvi-pan-genome/
        command = [
            CONDA_PATH,
            "run",
            "-n",
            CONDA_ENV,
            "anvi-pan-genome",
            "-g",
            os.path.join(ANVIO_OUTPUT, ANVIO_GENOMES_DB),
            "-n",
            ANVIO_PROJECT_NAME,
            "-o",
            os.path.join(ANVIO_OUTPUT, "PAN"),
            "--num-threads",
            THREADS,
            "--minbit",
            "0.5",  # Minimum alignment coverage (0.5 = 50%, good default)
            "--mcl-inflation",
            "10",
            "--use-ncbi-blast",  # Use NCBI BLAST instead of DIAMOND (more sensitive, slower)
        ]
        result = app_utility.bash_execute(command)
        if not result["success"]:
            print(f"Error running pangenome analysis: {result['error'].output}")
            exit(1)

        print("\n✅ Pangenome analysis completed successfully.")

        # Compute ANI (Average Nucleotide Identity) -> https://anvio.org/help/8/programs/anvi-compute-genome-similarity/
        command = [
            CONDA_PATH,
            "run",
            "-n",
            CONDA_ENV,
            "anvi-compute-genome-similarity",
            "--external-genomes",
            external_genomes_file,
            "--program",
            "pyANI",
            "--output-dir",
            os.path.join(ANVIO_OUTPUT, "ANI"),
            "--num-threads",
            THREADS,
            "--pan-db",
            os.path.join(ANVIO_OUTPUT, "PAN", f"{ANVIO_PROJECT_NAME}-PAN.db"),
        ]
        result = app_utility.bash_execute(command)
        if not result["success"]:
            print(f"Error computing genome similarity: {result['error'].output}")
            exit(1)

        print(
            "\n✅ Genome similarity analysis completed successfully. Please find GenomeDB and PanDB in the Anvio_Results directory. Please use server.py script to run the pan genome web server and visualize the pan."
        )

    if PIPE_RELATIONSHIPS:
        # Compute phylogenomic tree from core genes -> https://anvio.org/help/main/programs/anvi-get-sequences-for-gene-clusters/

        # ANIb_percentage_identity.newick: Best for strain typing and identifying closely related isolates
        # phylogenomic_tree.nwk: Best for evolutionary analysis and understanding gene-level relationships
        # Both trees may show similar clustering patterns for closely related genomes, but can differ significantly when comparing more divergent strains or when horizontal gene transfer is involved.
        print("\nExtracting core gene sequences for phylogenomic tree...")
        command = [
            CONDA_PATH,
            "run",
            "-n",
            CONDA_ENV,
            "anvi-get-sequences-for-gene-clusters",
            "-p",
            os.path.join(ANVIO_OUTPUT, "PAN", f"{ANVIO_PROJECT_NAME}-PAN.db"),
            "-g",
            os.path.join(ANVIO_OUTPUT, ANVIO_GENOMES_DB),
            # "--min-num-genomes-gene-cluster-occurs",
            str(len(combined_genomes)),
            "--concatenate-gene-clusters",
            "--output-file",
            os.path.join(ANVIO_OUTPUT, "PAN", "gene_clusters_aligned.faa"),
            "--force-overwrite",
            "-C",
            "DEFAULT",
        ]
        result = app_utility.bash_execute(command)
        if not result["success"]:
            print(f"Error extracting core gene sequences: {result['error'].output}")
            exit(1)

        print("\n✅ Core gene sequences extracted.")

        # Generate phylogenomic tree -> https://anvio.org/help/main/programs/anvi-gen-phylogenomic-tree/
        command = [
            CONDA_PATH,
            "run",
            "-n",
            CONDA_ENV,
            "anvi-gen-phylogenomic-tree",
            "-f",
            os.path.join(ANVIO_OUTPUT, "PAN", "gene_clusters_aligned.faa"),
            "-o",
            os.path.join(ANVIO_OUTPUT, "PAN", "phylogenomic_tree.nwk"),
        ]
        result = app_utility.bash_execute(command)
        if not result["success"]:
            print(f"Error generating phylogenomic tree: {result['error'].output}")
            exit(1)

        print("\n✅ Phylogenomic tree generated successfully.")

    ####################################################################################################################
    # STEP 3 : Prokka Annotation
    ####################################################################################################################

    # Collect all genomes to annotate (both contigs and references)
    genomes_to_annotate = {}

    # Add reference genomes from DATA/Ref
    for ref_id, ref_data in validated_references.items():
        genomes_to_annotate[ref_id] = ref_data["file"]

    # Add sample contigs from SPAdes
    for sample_id, contig_path in validated_contigs.items():
        genomes_to_annotate[sample_id] = contig_path

    if PIPE_PROKKA:
        app_utility.clean_output_directory(PROKKA_OUTPUT)

        print(f"\n✅ Total genomes to annotate: {len(genomes_to_annotate)}")

        # Run Prokka for each genome
        for genome_id, genome_path in genomes_to_annotate.items():
            print(f"\nAnnotating {genome_id}...")

            output_dir = os.path.join(PROKKA_OUTPUT, genome_id)

            command = [
                CONDA_PATH,
                "run",
                "-n",
                "GEM-PROKKA",
                "prokka",
                "--outdir",
                output_dir,
                "--prefix",
                genome_id,
                "--cpus",
                THREADS,
                "--force",  # Overwrite existing output directory
                genome_path,
            ]

            result = app_utility.bash_execute(command)
            if not result["success"]:
                print(f"Error running Prokka for {genome_id}: {result['error']}")
                continue

            print(f"✅ Prokka annotation completed for {genome_id}")

        print(
            f"\n✅ All Prokka annotations completed. Results stored in {PROKKA_OUTPUT}"
        )

    ####################################################################################################################
    # STEP 4 : ROARY Pangenome Analysis
    ####################################################################################################################

    if PIPE_ROARY:
        app_utility.clean_output_directory(ROARY_OUTPUT)

        # Collect all GFF files from Prokka output directory
        gff_files = []
        seen_labels = set()

        if not os.path.exists(PROKKA_OUTPUT):
            print(f"ERROR: Prokka output directory not found: {PROKKA_OUTPUT}")
            exit(1)

        # Iterate through sample directories in Prokka output
        for sample_id in validated_contigs.keys():
            # Look for exact directory match
            sample_dir = os.path.join(PROKKA_OUTPUT, sample_id)

            if not os.path.isdir(sample_dir):
                print(
                    f"WARNING: No Prokka directory found for {sample_id}; skipping..."
                )
                continue

            # Check for duplicate labels (shouldn't happen with exact matching)
            if sample_id in seen_labels:
                print(
                    f"WARNING: Duplicate sample label {sample_id} encountered; skipping duplicate entry."
                )
                continue

            # Check if GFF file exists
            gff_path = os.path.join(sample_dir, f"{sample_id}.gff")

            if not os.path.exists(gff_path):
                print(f"WARNING: GFF file missing for {sample_id}; skipping...")
                continue

            # Add full path to list
            gff_files.append(os.path.realpath(gff_path))
            seen_labels.add(sample_id)
            print(f"✓ Found {os.path.realpath(gff_path)}")

        if len(gff_files) == 0:
            print("ERROR: No GFF files collected for Roary.")
            exit(1)

        print(f"\n✅ Collected {len(gff_files)} GFF files for Roary analysis")

        # Run Roary core genome alignment -> https://github.com/sanger-pathogens/Roary

        print("\nRunning Roary core genome alignment...")

        # Change to ROARY_OUTPUT directory to run Roary there
        original_dir = os.getcwd()
        os.chdir(ROARY_OUTPUT)

        command = [
            CONDA_PATH,
            "run",
            "-n",
            "GEM-ROARY",
            "roary",
            "-e",  # Create a multiFASTA alignment of core genes using PRANK
            "-n",  # Fast core gene alignment with MAFFT
            "-p",
            THREADS,  # Number of threads
            "-f",
            ".",  # Output directory
        ] + gff_files  # Add all GFF file paths
        result = app_utility.bash_execute(command)
        if not result["success"]:
            print(f"Error running Roary: {result['error'].output}")
            exit(1)

        print(f"\n✅ Roary analysis completed. Results stored in {ROARY_OUTPUT}")

        # IQ TREE Phylogenetic Analysis

        # print("\nRunning IQ-TREE phylogenetic analysis...")

        # # Run IQ-TREE -> http://www.iqtree.org/
        # command = [
        #     CONDA_PATH,
        #     "run",
        #     "-n",
        #     "GEM-ROARY",  # or your IQ-TREE conda environment
        #     "iqtree",
        #     "-s",
        #     "core_gene_alignment.aln",  # Input alignment
        #     "-T",
        #     THREADS,  # Number of threads
        #     "-m",
        #     "GTR+G",  # GTR model with Gamma rate heterogeneity
        #     "-bb",
        #     "1000",  # Bootstrap replicates (1000 is standard)
        # ]
        # result = app_utility.bash_execute(command)
        # if not result["success"]:
        #     print(f"Error running IQ-TREE: {result['error'].output}")
        #     exit(1)

        # print("\n✅ IQ-TREE phylogenetic analysis completed.")
