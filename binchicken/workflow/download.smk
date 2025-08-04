#############
### Setup ###
#############
import os
import polars as pl
from binchicken.common import pixi_run
os.umask(0o002)

output_dir = os.path.abspath("download")
logs_dir = output_dir + "/logs"
benchmarks_dir = output_dir + "/benchmarks"

tmpdir = config["tmpdir"]
singlem_metapackage = config["singlem_metapackage"]
checkm2_db = config["checkm2_db"]
gtdbtk_db = config["gtdbtk_db"]

download_arg = "--download "
if singlem_metapackage:
    download_arg += "singlem "
if checkm2_db:
    download_arg += "checkm2 "
if gtdbtk_db:
    download_arg += "gtdb "

def get_mem_mb(wildcards, threads, attempt):
    return 8 * 1000 * threads * attempt

def get_runtime(base_hours):
    def runtime_func(wildcards, attempt):
        return f"{attempt * base_hours}h"
    return runtime_func

#############
### Rules ###
#############
rule all:
    input:
        output_dir + "/aviary_downloads.done",
    localrule: True

rule aviary_download:
    output:
        output_dir + "/aviary_downloads.done",
    params:
        tmpdir = f"export TMPDIR={tmpdir} && " if tmpdir else "",
        singlem_metapackage_env = "export SINGLEM_METAPACKAGE_PATH=. && " if not singlem_metapackage else "",
        checkm2_db_env = "export CHECKM2DB=. && " if not checkm2_db else "",
        gtdbtk_db_env = "export GTDBTK_DATA_PATH=. && " if not gtdbtk_db else "",
        singlem_metapackage = "--singlem-metapackage-path " + singlem_metapackage if singlem_metapackage else "",
        checkm2_db = "--checkm2-db-path " + checkm2_db if checkm2_db else "",
        gtdbtk_db = "--gtdb-path " + gtdbtk_db if gtdbtk_db else "",
        download_arg = download_arg,
    threads: 16
    resources:
        mem_mb=get_mem_mb,
        runtime = get_runtime(base_hours = 16),
    log:
        logs_dir + "/aviary_downloads.log"
    shell:
        "{params.tmpdir} "
        "{params.singlem_metapackage_env} "
        "{params.checkm2_db_env} "
        "{params.gtdbtk_db_env} "
        "export EGGNOG_DATA_DIR=. && "
        f"{pixi_run} -e aviary "
        "aviary configure "
        "{params.singlem_metapackage} "
        "{params.checkm2_db} "
        "{params.gtdbtk_db} "
        "{params.download_arg} "
        "&> {log} "
        "&& touch {output} "
