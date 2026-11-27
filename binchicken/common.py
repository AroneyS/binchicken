import os
from importlib.resources import files, as_file
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

RECOVERED_BIN_NAME_REGEX = re.compile(r"^binchicken_(.+)\.\d+$")
ITERATION_PREFIX_REGEX = re.compile(r"^it\d+_")

def bin_coassembly(bin_name, coassemblies=None):
    """Return the coassembly (or single-assembly sample) a bin key belongs to.

    Bin keys are either recovered-dir names, binchicken_[it<N>_]<coassembly>.<n>,
    or <coassembly>-<genome> (from --aviary-outputs / --new-genomes).

    Coassembly names may themselves contain "-" or start with it<N>_, so when the known
    coassembly names are given, keys are resolved against them: recovered names are tried
    as-is before stripping an iteration prefix, and <coassembly>-<genome> keys take the
    longest matching name. Without known names (or no match), keys are parsed directly.
    """
    known = set(coassemblies) if coassemblies is not None else None
    match = RECOVERED_BIN_NAME_REGEX.match(bin_name)
    if match:
        name = match.group(1)
        stripped = ITERATION_PREFIX_REGEX.sub("", name, count=1)
        if known is not None:
            for candidate in [name, stripped]:
                if candidate in known:
                    return candidate
        return stripped

    if known is not None:
        matches = [c for c in known if bin_name.startswith(c + "-")]
        if matches:
            return max(matches, key=len)
    return bin_name.split("-")[0]

def split_bin_name(name, coassemblies=None):
    """Split a bin name, or a SingleM OTU table sample name (<bin>_transcripts), into (coassembly, genome).

    Recovered-dir names keep the whole name as the genome; <coassembly>-<genome> names keep
    the genome suffix (None if there is none).
    """
    base = name[:-len("_transcripts")] if name.endswith("_transcripts") else name
    if RECOVERED_BIN_NAME_REGEX.match(base):
        return bin_coassembly(base, coassemblies), name

    coassembly = bin_coassembly(name, coassemblies)
    genome = name[len(coassembly) + 1:] if name.startswith(coassembly + "-") else None
    return coassembly, genome

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
    pixi_source = files("binchicken").joinpath("pixi.toml")

    with as_file(pixi_source) as manifest_path:
        return f"pixi run --as-is --manifest-path {manifest_path}"

pixi_run = pixi_run_func()
# --as-is is only valid for `pixi run`, so build the explicit install separately
pixi_install = pixi_run.replace("run --as-is", "install -a", 1)


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
