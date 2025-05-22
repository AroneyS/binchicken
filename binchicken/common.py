import os
import re
import ast

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
    binchicken_basedir = os.path.join(os.path.dirname(__file__), "..")
    return f"pixi run --frozen --manifest-path {binchicken_basedir}/pixi.toml"

pixi_run = pixi_run_func()
