import os
import importlib.resources
import re
import ast

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
