"""
Author: Y. Ahmed-Braimah
--- Workflow for Trans Proteomic Pipeline (TPP)
--- using a TPP docker container (spctools/tpp)

"""

import json
import os
import re
from os.path import join, basename, dirname
from os import getcwd
from subprocess import check_output

##--------------------------------------------------------------------------------------##
## Functions
##--------------------------------------------------------------------------------------##

# To print process messages
def message(x):
  print()

# To remove suffix from a string
def rstrip(text, suffix):
    if not text.endswith(suffix):
        return text
    return text[:len(text)-len(suffix)]

## define environment variables

##--------------------------------------------------------------------------------------##
## Global config files:
##--------------------------------------------------------------------------------------##

configfile: 'config.yml'

PEP = config['PEP']

# Full path to a folder where final output files will be deposited.
OUT_DIR = config['OUT_DIR']
HOME_DIR = config['HOME_DIR']
LOGS_DIR = config['LOGS_DIR']

FILES = json.load(open(config['MZ_SAMPLES_JSON']))
SAMPLES = sorted(FILES.keys())


## Create the final output directory if it doesn't already exist
if not os.path.exists(OUT_DIR):
            os.makedirs(OUT_DIR)

if not os.path.exists(LOGS_DIR):
            os.makedirs(LOGS_DIR)


ENGINES = ['comet']

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule all:
    input:
        join(OUT_DIR, 'iProph', 'xinteract.iProph.prot.xml')
##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule decoy:
    input:
        pep = PEP
    output:
        decoy = join(OUT_DIR, 'database_plus_decoy.fa')
    params:
        i = join('/data', basename(PEP)),
        o = join('/data', basename(OUT_DIR), 'database_plus_decoy.fa')
    singularity:
        "docker://spctools/tpp"
    shell:
        "decoyFastaGenerator.pl {params.i} DECOY {params.o}"

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule tandem:
    input:
        mzml = lambda wildcards: FILES[wildcards.sample],
        tPARAM = join(HOME_DIR, 'params', 'tandem_params.xml'),
        decoy = join(OUT_DIR, 'database_plus_decoy.fa')
    output:
        tandem = join(OUT_DIR, 'tandem', '{sample}.tandem')
    params:
        i = "/data/mzML_files/{sample}.mzML",
        o = join('/data', basename(OUT_DIR), 'tandem', '{sample}.tandem'),
        t = "/data/params/tandem_params.xml",
        l = join('/data', basename(OUT_DIR), basename(LOGS_DIR), '{sample}_tandem.log')
    threads:
        8
    resources:
        mem_mb=16000
    singularity:
        "docker://spctools/tpp"
    shell:
        "sed 's/SAMPLE_NAME/{wildcards.sample}/g' {params.t} > {wildcards.sample}_tandem.params && "
        "tandem {wildcards.sample}_tandem.params > {params.l} 2>&1 && "
        "rm {wildcards.sample}_tandem.params &&"

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule Tandem2XML:
    input:
        tandem = join(OUT_DIR, 'tandem', '{sample}.tandem')
    output:
        xml = join(OUT_DIR, 'tandem', '{sample}.pep.xml')
    params:
        i = "/data/OUTPUT/tandem/{sample}.tandem",
        o = "/data/OUTPUT/tandem/{sample}.pep.xml",
        l = "/data/OUTPUT/LOGS/{sample}_T2XML.log"
    threads:
        8
    resources:
        mem_mb=16000
    log:
        "/home/yahmed/Proteomics/Aedes_aegypti/SM_OUTPUT/LOGS/{sample}_T2XML.log"
    singularity:
        "docker://spctools/tpp"
    shell:
        "Tandem2XML {params.i} {params.o} > {params.l} 2>&1"

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule comet:
    input:
        mzml = lambda wildcards: FILES[wildcards.sample],
        cPARAM = join(HOME_DIR, 'params', 'comet.params.high-high'),
        decoy = join(OUT_DIR, 'database_plus_decoy.fa')
    output:
        comet = join(OUT_DIR, 'comet', '{sample}.pep.xml')
    params:
        i = "/data/mzML_files/{sample}.mzML",
        o = join('/data', basename(OUT_DIR), 'comet', '{sample}.pep.xml'),
        t = "/data/params/comet.params.high-high",
        l = join('/data', basename(OUT_DIR), basename(LOGS_DIR), '{sample}_comet.log'),
        d = join('/data', basename(OUT_DIR), 'database_plus_decoy.fa')
    singularity:
        "docker://spctools/tpp"
    threads:
        8
    resources:
        mem_mb=16000
    shell:
        'comet'
        ' -P{params.t}'
        ' -D{params.d}'
        ' -N{wildcards.sample}'
        ' {params.i}'
        ' > {params.l} 2>&1'
        ' && mv {wildcards.sample}.pep.xml ' + join('/data', basename(OUT_DIR), 'comet')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule PeptideProphet:
    input:
        pepXML = expand(join(OUT_DIR, '{engine}', '{sample}.pep.xml'), sample = SAMPLES, engine = ENGINES)
    output:
        PepProph = join(OUT_DIR, 'PepProph', 'xinteract.pep.xml')
    params:
        # i = join('/data', basename(OUT_DIR), '{engine}', '{sample}.pep.xml'),
        o = join('/data', basename(OUT_DIR), 'PepProph', 'xinteract.pep.xml'),
        l = join('/data', basename(OUT_DIR), basename(LOGS_DIR), 'combined_PepProph.log'),
        d = join('/data', basename(OUT_DIR), 'database_plus_decoy.fa'),
        c = "/data/params/condition.xml"
    singularity:
        "docker://spctools/tpp"
    threads:
        8
    resources:
        mem_mb=32000
    shell:
        'xinteract'
        ' -N{params.o}'
        ' -p0.05 -l5 -PPM -OAP -D{params.d}'
        ' -L{params.c}'
        ' -THREADS=8'
        ' ' + join('/data', basename(OUT_DIR), '*', '*.pep.xml') +
        ' > {params.l} 2>&1'

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule InterProphet:
    input:
        # PepProph = join(OUT_DIR, 'PepProph', '{engine}', '{sample}.pep.xml')
        PepProph = join(OUT_DIR, 'PepProph', 'xinteract.pep.xml')
    output:
        # iProph = join(OUT_DIR, 'iProph', '{engine}', '{sample}.iProph.pep.xml')
        iProph = join(OUT_DIR, 'iProph', 'xinteract.iProph.pep.xml')
    params:
        i = join('/data', basename(OUT_DIR), 'PepProph', 'xinteract.pep.xml'),
        o = join('/data', basename(OUT_DIR), 'iProph', 'xinteract.iProph.pep.xml'),
        l = join('/data', basename(OUT_DIR), basename(LOGS_DIR), 'combined_iProph.log'),
        d = join('/data', basename(OUT_DIR), 'database_plus_decoy.fa')
    singularity:
        "docker://spctools/tpp"
    threads:
        8
    resources:
        mem_mb=32000
    shell:
        'InterProphetParser'
        ' THREADS=8'
        ' DECOY={params.d}'
        ' {params.i}'
        ' {params.o}'
        ' > {params.l} 2>&1'

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule ProteinProphet:
    input:
        # iProph = join(OUT_DIR, 'iProph', '{engine}', '{sample}.iProph.pep.xml')
        iProph = join(OUT_DIR, 'iProph', 'xinteract.iProph.pep.xml')
    output:
        # ProtProph = join(OUT_DIR, 'iProph', '{engine}', '{sample}.iProph.prot.xml')
        ProtProph = join(OUT_DIR, 'iProph', 'xinteract.iProph.prot.xml')
    params:
        i = join('/data', basename(OUT_DIR), 'iProph', 'xinteract.iProph.pep.xml'),
        o = join('/data', basename(OUT_DIR), 'iProph', 'xinteract.iProph.prot.xml'),
        l = join('/data', basename(OUT_DIR), basename(LOGS_DIR), 'combined_ProtProph.log'),
        c = "/data/params/condition.xml"
    singularity:
        "docker://spctools/tpp"
    shell:
        'ProteinProphet'
        ' {params.i}'
        ' {params.o}'
        ' IPROPHET NOGROUPWTS PLOTPNG EXCELPEPS EXCEL.90'
        ' LIBRA{params.c}'
        ' > {params.l} 2>&1'

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule InterProphet_comb:
    input:
        pepXML = expand("/home/yahmed/Proteomics/Aedes_aegypti/SM_OUTPUT/PepProph/{engine}/{sample}.pep.xml", sample = SAMPLES, engine = ENGINES)
    output:
        iProph = "/home/yahmed/Proteomics/Aedes_aegypti/SM_OUTPUT/iProph/{sample}.iProph_comb.pep.xml"
    params:
        ic = "/data/SM_OUTPUT/PepProph/comet/{sample}.pep.xml",
        it = "/data/SM_OUTPUT/PepProph/tandem/{sample}.pep.xml",
        o = "/data/SM_OUTPUT/iProph/{sample}.iProph_comb.pep.xml",
        l = "/data/SM_OUTPUT/LOGS/{sample}_iProph_comb_pep.log",
        d = "/data/database_plus_decoy.fa"
    log:
        "/home/yahmed/Proteomics/Aedes_aegypti/SM_OUTPUT/LOGS/{sample}_iProph_comb_pep.log"
    singularity:
        "docker://spctools/tpp"
    shell:
        'InterProphetParser'
        ' THREADS=8'
        ' DECOY={params.d}'
        ' {params.ic}'
        ' {params.it}'
        ' {params.o}'
        ' > {params.l} 2>&1'

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule ProteinProphet_comb:
    input:
        iProph = "/home/yahmed/Proteomics/Aedes_aegypti/SM_OUTPUT/iProph/{sample}.iProph_comb.pep.xml"
    output:
        iProph = "/home/yahmed/Proteomics/Aedes_aegypti/SM_OUTPUT/iProph/{sample}.iProph_comb_prot.prot.xml"
    params:
        i = "/data/SM_OUTPUT/iProph/{sample}.iProph_comb.pep.xml",
        o = "/data/SM_OUTPUT/iProph/{sample}.iProph_comb_prot.prot.xml",
        l = "/data/SM_OUTPUT/LOGS/{sample}_iProph_comb_prot.log",
        d = "/data/database_plus_decoy.fa"
    log:
        "/home/yahmed/Proteomics/Aedes_aegypti/SM_OUTPUT/LOGS/{sample}_iProph_comb_prot.log"
    singularity:
        "docker://spctools/tpp"
    shell:
        'ProteinProphet'
        ' {params.i}'
        ' {params.o}'
        ' IPROPHET NOGROUPWTS PLOTPNG EXCELPEPS EXCEL.90'
        ' > {params.l} 2>&1'
