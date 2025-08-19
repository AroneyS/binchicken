import os
import importlib.resources
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
    with importlib.resources.path("binchicken", "pixi.toml") as manifest_path:
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
