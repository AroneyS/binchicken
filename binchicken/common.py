import os
from pathlib import Path
import re
import ast
import time

def fasta_to_name(filename):
    sample_name = os.path.basename(filename)
    # Put .gz first so it is stripped off first.
    for extension in (".gz",".fna",".fq",".fastq",".fasta",".fa"):
        if sample_name.endswith(extension):
            sample_name = sample_name[0:(len(sample_name)-len(extension))]
            if extension != ".gz":
                break
    return sample_name

def parse_snake_dict(s):
    # Remove leading/trailing whitespace
    s = s.strip()
    # If it doesn't start with '{', it's not a dict
    if not s.startswith("{"):
        raise ValueError("Input does not look like a dict: " + s)
    # Add quotes around keys and values if missing
    # This regex finds keys and values that are not quoted and adds quotes
    s = re.sub(r'([{,]\s*)([^:,}{]+)\s*:', r'\1"\2":', s)
    s = re.sub(r':\s*([^,}{]+)', lambda m: ': "' + m.group(1).strip().strip('"\'') + '"', s)
    # Now safe to use ast.literal_eval
    return ast.literal_eval(s)

def pixi_run_func():
    # pixi.toml ships inside this package directory (see MANIFEST.in), so
    # under a normal pip/conda install that's the only copy and the only
    # sensible pixi project root. In a git checkout, though, a convenience
    # symlink at the repo root also points at this same file, and the
    # project's own pixi tasks (postinstall's `cwd = "binchicken"` + `-e ..`)
    # are written assuming *that* repo-root location is the pixi project
    # root -- dev tooling that runs `pixi run ...` from the repo root
    # creates environments on that basis. If this function instead pointed
    # at the real nested file, pixi would see the same project referenced
    # via two different literal manifest paths and, since 0.82.1, refuse to
    # continue with "the environment directory has moved". Prefer the
    # repo-root symlink when it's actually an alias for our own pixi.toml,
    # to stay consistent with whichever path created the environment; fall
    # back to this package's own location otherwise (the only copy that
    # exists in a real install).
    package_dir = Path(__file__).resolve().parent
    local_manifest = package_dir / "pixi.toml"
    repo_root_manifest = package_dir.parent / "pixi.toml"
    if repo_root_manifest.exists() and repo_root_manifest.resolve() == local_manifest.resolve():
        manifest_path = repo_root_manifest
    else:
        manifest_path = local_manifest
    return f"pixi run --frozen --manifest-path {manifest_path}"

pixi_run = pixi_run_func()


# A unique identifier for the current workflow invocation. This is used
# when creating log directories so that retries from the same workflow do
# not overwrite previous log files.
workflow_identifier = time.strftime("%Y%m%d_%H%M%S")


def setup_log(log_dir_base: str, attempt: int) -> str:
    """Return a unique log file path for a given rule attempt.

    Parameters
    ----------
    log_dir_base: str
        Directory in which log files for a rule should be stored.  This
        function will create a subdirectory named with the workflow
        identifier and place attempt specific log files inside it.
    attempt: int
        The Snakemake retry attempt number.

    Returns
    -------
    str
        Path to a log file unique to the workflow invocation and attempt.
    """

    log_dir = os.path.join(log_dir_base, workflow_identifier)
    os.makedirs(log_dir, exist_ok=True)
    return os.path.join(log_dir, f"attempt{attempt}.log")

BIN_INFO_COLUMNS = {
    "Name": str,
    "Completeness": str,
    "Contamination": str,
    "Completeness_Model_Used": str,
    "Translation_Table_Used": str,
    "Coding_Density": str,
    "Contig_N50": str,
    "Average_Gene_Length": str,
    "Genome_Size": str,
    "GC_Content": str,
    "Total_Coding_Sequences": str,
    "Total_Contigs": str,
    "Max_Contig_Length": str,
    "Additional_Notes": str,
}

CHECKM2_QUALITY_COLUMNS = {
    "Name": str,
    "Completeness": str,
    "Contamination": str,
    "Completeness_Model_Used": str,
    "Translation_Table_Used": str,
    "Coding_Density": str,
    "Contig_N50": str,
    "Average_Gene_Length": str,
    "Genome_Size": str,
    "GC_Content": str,
    "Total_Coding_Sequences": str,
    "Total_Contigs": str,
    "Max_Contig_Length": str,
    "Additional_Notes": str,
}

ELUSIVE_EDGES_COLUMNS = {
    "style": str,
    "cluster_size": int,
    "samples": str,
    "target_ids": str,
}

ELUSIVE_CLUSTERS_COLUMNS = {
    "samples": str,
    "length": int,
    "total_targets": float,
    "total_size": int,
    "recover_samples": str,
    "coassembly": str,
}
