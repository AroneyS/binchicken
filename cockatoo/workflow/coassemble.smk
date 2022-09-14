#############
### Setup ###
#############
import os

output_dir = os.path.abspath("coassemble")
logs_dir = output_dir + "/logs"

reads_1 = config["reads_1"]
reads_2 = config["reads_2"]

rule all:
    input:
        output_dir + "/commands/coassemble_commands.sh",
        output_dir + "/commands/recover_commands.sh"

##################################
### Map reads to matching bins ###
##################################


##############################
### Create Aviary commands ###
##############################
rule coassemble_commands:
    input:
        elusive_clusters=config["elusive_clusters"],
    output:
        coassemble_commands=output_dir + "/commands/coassemble_commands.sh",
        recover_commands=output_dir + "/commands/recover_commands.sh"
    params:
        reads_1=reads_1,
        reads_2=reads_2,
    log:
        logs_dir + "/coassemble_commands.log"
    script:
        "scripts/coassemble_commands.py"
