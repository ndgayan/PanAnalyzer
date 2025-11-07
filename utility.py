import csv
import os
import shutil
import subprocess


def clean_output_directory(output_dir):
    """
    Delete all files and subdirectories in the output directory.

    Args:
        output_dir: Path to the output directory to clean
    """
    if os.path.exists(output_dir):
        print(f"\n🗑️  Cleaning output directory: {output_dir}")
        try:
            shutil.rmtree(output_dir)
            print(f"✓ Successfully removed: {output_dir}")
        except Exception as e:
            print(f"✗ Error removing directory: {e}")
            return False

    # Recreate the directory
    os.makedirs(output_dir, exist_ok=True)
    print(f"✓ Created clean directory: {output_dir}\n")
    return True


def get_files_from_directory(directory_path):
    """
    Get all files from a directory and return their absolute paths.

    Args:
        directory_path (str): Path to the directory

    Returns:
        list: Sorted list of absolute file paths

    Raises:
        FileNotFoundError: If the directory doesn't exist
        ValueError: If the path is not a directory
    """
    if not os.path.exists(directory_path):
        raise FileNotFoundError(f"Directory not found: {directory_path}")

    if not os.path.isdir(directory_path):
        raise ValueError(f"Path is not a directory: {directory_path}")

    file_paths = [
        os.path.abspath(os.path.join(directory_path, f))
        for f in os.listdir(directory_path)
        if os.path.isfile(os.path.join(directory_path, f))
    ]

    file_paths.sort()
    print(f"Found {len(file_paths)} files in {directory_path}")
    return file_paths


def read_study_file(file_path):
    """
    Read the study CSV file and return a dictionary.

    Args:
        file_path (str): Path to the CSV file

    Returns:
        dict: Dictionary with sample IDs as keys and [type, species_name] as values

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    sample_dict = {}
    with open(file_path, "r", encoding="utf-8") as file:
        for row in csv.reader(file):
            if len(row) >= 3:
                sample_dict[row[0].strip()] = [row[1].strip(), row[2].strip()]

    print(f"Successfully loaded {len(sample_dict)} samples from {file_path}")
    return sample_dict


def _filter_samples_by_type(samples_dict, filter_type):
    """Filter samples dictionary by type."""
    return {
        sample_id: info
        for sample_id, info in samples_dict.items()
        if info[0] == filter_type
    }


def _create_file_lookup(file_paths):
    """Create a filename to filepath mapping for quick lookup."""
    return {os.path.basename(f): f for f in file_paths}


def _find_paired_read_files(
    sample_id, available_files, forward_postfix, reverse_postfix
):
    """Find forward and reverse read files for a sample."""
    forward_file = reverse_file = None

    for filename, filepath in available_files.items():
        if filename.startswith(f"{sample_id}_"):
            if filename.endswith(forward_postfix):
                forward_file = filepath
            elif filename.endswith(reverse_postfix):
                reverse_file = filepath

    return forward_file, reverse_file


def _generate_missing_file_errors(
    sample_id, forward_file, reverse_file, forward_postfix, reverse_postfix
):
    """Generate error messages for missing files."""
    errors = []
    if not forward_file:
        errors.append(
            f"  - {sample_id}: Missing forward read (expected: {sample_id}_*{forward_postfix})"
        )
    if not reverse_file:
        errors.append(
            f"  - {sample_id}: Missing reverse read (expected: {sample_id}_*{reverse_postfix})"
        )
    return errors


def _raise_validation_error(missing_files, error_type="files"):
    """Raise a FileNotFoundError with formatted error message."""
    error_msg = (
        f"\n❌ VALIDATION FAILED: Missing {len(missing_files)} required {error_type}:\n"
        + "\n".join(missing_files)
    )
    print(error_msg)
    raise FileNotFoundError(error_msg)


def validate_sample_files(
    samples_dict, sample_files, forward_postfix, reverse_postfix, filter_type="Sample"
):
    """
    Validate that each sample has both forward and reverse read files.

    Args:
        samples_dict (dict): Dictionary with sample IDs as keys and [type, species_name] as values
        sample_files (list): List of file paths in the Samples directory
        forward_postfix (str): Expected postfix for forward reads
        reverse_postfix (str): Expected postfix for reverse reads
        filter_type (str): Filter samples by type (default: "Sample")

    Returns:
        dict: Validated samples with forward/reverse file paths

    Raises:
        FileNotFoundError: If required files are missing
    """
    filtered_samples = _filter_samples_by_type(samples_dict, filter_type)
    available_files = _create_file_lookup(sample_files)

    print(f"\nValidating {len(filtered_samples)} samples of type '{filter_type}'...")

    validated_samples = {}
    missing_files = []

    for sample_id, sample_info in filtered_samples.items():
        forward_file, reverse_file = _find_paired_read_files(
            sample_id, available_files, forward_postfix, reverse_postfix
        )

        if forward_file and reverse_file:
            validated_samples[sample_id] = {
                "forward": forward_file,
                "reverse": reverse_file,
                "type": sample_info[0],
                "species": sample_info[1],
            }
        else:
            missing_files.extend(
                _generate_missing_file_errors(
                    sample_id,
                    forward_file,
                    reverse_file,
                    forward_postfix,
                    reverse_postfix,
                )
            )

    if missing_files:
        _raise_validation_error(missing_files, "sample file(s)")

    print(
        f"✅ Validation successful: All {len(validated_samples)} samples have required paired-end reads"
    )
    return validated_samples


def validate_reference_files(
    samples_dict, reference_files, file_extension=".fna", filter_type="Reference"
):
    """
    Validate that each reference sample has its corresponding reference file.

    Args:
        samples_dict (dict): Dictionary with sample IDs as keys and [type, species_name] as values
        reference_files (list): List of file paths in the Reference directory
        file_extension (str): Expected file extension (default: ".fna")
        filter_type (str): Filter samples by type (default: "Reference")

    Returns:
        dict: Validated references with file paths

    Raises:
        FileNotFoundError: If required files are missing
    """
    filtered_samples = _filter_samples_by_type(samples_dict, filter_type)
    available_files = _create_file_lookup(reference_files)

    print(
        f"\nValidating {len(filtered_samples)} reference samples of type '{filter_type}'..."
    )

    validated_references = {}
    missing_files = []

    for sample_id, sample_info in filtered_samples.items():
        expected_filename = f"{sample_id}{file_extension}"

        # Try exact match, then case-insensitive match
        reference_file = available_files.get(expected_filename) or next(
            (
                filepath
                for filename, filepath in available_files.items()
                if filename.lower() == expected_filename.lower()
            ),
            None,
        )

        if reference_file:
            validated_references[sample_id] = {
                "file": reference_file,
                "type": sample_info[0],
                "species": sample_info[1],
            }
        else:
            missing_files.append(
                f"  - {sample_id}: Missing reference file (expected: {expected_filename})"
            )

    if missing_files:
        _raise_validation_error(missing_files, "reference file(s)")

    print(
        f"✅ Validation successful: All {len(validated_references)} reference samples have required files"
    )
    return validated_references


def bash_execute(command):
    """
    Execute a bash command and return the result.

    Args:
        command (list): Command and arguments as a list
    """
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        return {"success": True, "result": result}
    except subprocess.CalledProcessError as e:
        return {"success": False, "error": e}


def run_spades_for_sample(
    sample_id, r1_file, r2_file, output_dir, conda_path, conda_env
):
    """
    Run SPAdes for a single sample using conda run.

    Args:
        sample_id: Sample identifier
        r1_file: Path to forward reads file
        r2_file: Path to reverse reads file
        output_dir: Output directory for SPAdes results

    Returns:
        bool: True if successful, False otherwise
    """
    os.makedirs(output_dir, exist_ok=True)

    print(f"Running SPAdes for sample: {sample_id}")
    # Use 'conda run' to execute spades in the specified environment
    command = [
        conda_path,
        "run",
        "-n",
        conda_env,
        "spades.py",
        "--isolate",
        "-o",
        output_dir,
        "-1",
        r1_file,
        "-2",
        r2_file,
    ]

    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)

        print(f"✓ Completed: {sample_id}")
        print("-" * 60)
        print("")
        return True

    except subprocess.CalledProcessError as e:
        print(f"✗ Error running SPAdes for {sample_id}:")
        print(f"  {e.stderr}")
        print("-" * 60)
        print("")
        return False
