#############
### Setup ###
#############
import os

output_dir = os.path.abspath("coassemble")
logs_dir = output_dir + "/logs"

singlem_bin = "/home/aroneys/bin/singlem-dev"

rule all:
    input:
        output_dir + "/target/coassemble_commands.sh",
        output_dir + "/target/recover_commands.sh"

##################################
### Map reads to matching bins ###
##################################
rule map_reads:


##############################
### Create Aviary commands ###
##############################
rule coassemble_commands:
    input:
        elusive_clusters=output_dir + "/target/elusive_clusters.tsv"
    output:
        coassemble_commands=output_dir + "/target/coassemble_commands.sh",
        recover_commands=output_dir + "/target/recover_commands.sh"
    log:
        logs_dir + "/coassemble_commands.log"
    script:
        "scripts/coassemble_commands.py"
