"""
Author: Y. Ahmed-Braimah
--- Workflow for Trans Proteomic Pipeline (TPP)
--- using a TPP docker container (spctools/tpp)
--- (For TMT MS/MS)

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

ENGINES = ['comet', 'tandem']

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule TMT_quants:
    input:
        join(OUT_DIR, 'combined', 'iProph', 'xinteract.iProph_comb.prot.xml'),
        expand(join(OUT_DIR, '{engine}', 'TMT_quant.tsv'), engine = ENGINES)
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
    threads:
        8
    resources:
        mem_mb=16000
    message:
        """  ------  Generating protein decoy database  ------  """
    shell:
        "decoyFastaGenerator.pl {params.i} DECOY {params.o}"

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule tandem:
    input:
        mzml = join(HOME_DIR, 'mzML_files', '{sample}.mzML'),
        tPARAM = join(HOME_DIR, 'params', 'tandem_params.xml'),
        taxono = join(HOME_DIR, 'params', 'taxonomy.xml'),
        decoy = join(OUT_DIR, 'database_plus_decoy.fa')
    output:
        tandem = join(OUT_DIR, 'tandem', '{sample}.tandem')
    params:
        i = "/data/mzML_files/{sample}.mzML",
        o = join('/data', basename(OUT_DIR), 'tandem', '{sample}.tandem'),
        t = "/data/params/tandem_params.xml",
        x = "/data/params/taxonomy.xml",
        l = join('/data', basename(OUT_DIR), basename(LOGS_DIR), '{sample}_tandem.log')
    threads:
        16
    resources:
        mem_mb=32000
    singularity:
        "docker://spctools/tpp"
    message:
        """  ------  Runing X!Tandem search for mzML file {wildcards.sample}  ------  """
    shell:
        "sed 's/SAMPLE_NAME/{wildcards.sample}/g' {params.t} > {wildcards.sample}_tandem.params && "
        "sed 's/OUTPUT_DIR/" + basename(OUT_DIR) + "/g' {params.x} > taxonomy_{wildcards.sample}.xml && "
        "tandem {wildcards.sample}_tandem.params > {params.l} 2>&1 && "
        "mv {wildcards.sample}.tandem " + basename(OUT_DIR) + "/tandem/ &&"
        "rm {wildcards.sample}_tandem.params taxonomy_{wildcards.sample}.xml"

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule Tandem2XML:
    input:
        tandem = join(OUT_DIR, 'tandem', '{sample}.tandem')
    output:
        xml = join(OUT_DIR, 'tandem', '{sample}.pep.xml')
    params:
        i = join('/data', basename(OUT_DIR), 'tandem', '{sample}.tandem'),
        o = join('/data', basename(OUT_DIR), 'tandem', '{sample}.pep.xml'),
        l = join('/data', basename(OUT_DIR), basename(LOGS_DIR), '{sample}_T2XML.log')
    threads:
        8
    resources:
        mem_mb=16000
    singularity:
        "docker://spctools/tpp"
    message:
        """  ------  Converting X!Tandem search results for mzML file {wildcards.sample} to pep.xml ------  """
    shell:
        "Tandem2XML {params.i} {params.o} > {params.l} 2>&1"

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule comet:
    input:
        mzml = join(HOME_DIR, 'mzML_files', '{sample}.mzML'),
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
    message:
        """  ------  Runing Comet search for mzML file {wildcards.sample}  ------  """
    threads:
        16
    resources:
        mem_mb=32000
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
        PepProph = join(OUT_DIR, '{engine}', 'PepProph', 'xinteract.pep.xml')
    params:
        o = join('/data', basename(OUT_DIR), '{engine}', 'PepProph', 'xinteract.pep.xml'),
        l = join('/data', basename(OUT_DIR), basename(LOGS_DIR), '{engine}_PepProph.log'),
        d = join('/data', basename(OUT_DIR), 'database_plus_decoy.fa'),
        c = "/data/params/condition.xml"
    singularity:
        "docker://spctools/tpp"
    message:
        """  ------  Runing PeptideProphet search for combined mzMLs from {wildcards.engine} ------  """
    threads:
        16
    resources:
        mem_mb=32000
    shell:
        'xinteract'
        ' -N{params.o}'
        ' -p0.05'
        ' -l5'
        ' -PPM'
        ' -OAP'
        ' -D{params.d}'
        ' -dDECOY'
        ' -L{params.c}'
        ' -THREADS=8'
        ' ' + join('/data', basename(OUT_DIR), '{wildcards.engine}', '*.pep.xml') +
        ' > {params.l} 2>&1'

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule InterProphet:
    input:
        PepProph = join(OUT_DIR, '{engine}', 'PepProph', 'xinteract.pep.xml')
    output:
        iProph = join(OUT_DIR, '{engine}', 'iProph', 'xinteract.iProph.pep.xml')
    params:
        i = join('/data', basename(OUT_DIR), '{engine}', 'PepProph', 'xinteract.pep.xml'),
        o = join('/data', basename(OUT_DIR), '{engine}', 'iProph', 'xinteract.iProph.pep.xml'),
        l = join('/data', basename(OUT_DIR), basename(LOGS_DIR), '{engine}_iProph.log'),
        d = join('/data', basename(OUT_DIR), 'database_plus_decoy.fa')
    singularity:
        "docker://spctools/tpp"
    message:
        """  ------  Runing InterProphet search for combined mzMLs from {wildcards.engine} ------  """
    threads:
        16
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
        iProph = join(OUT_DIR, '{engine}', 'iProph', 'xinteract.iProph.pep.xml')
    output:
        ProtProph = join(OUT_DIR, '{engine}', 'iProph', 'xinteract.iProph.prot.xml'),
        quants = join(OUT_DIR, '{engine}', 'TMT_quant.tsv')
    params:
        i = join('/data', basename(OUT_DIR), '{engine}', 'iProph', 'xinteract.iProph.pep.xml'),
        o = join('/data', basename(OUT_DIR), '{engine}', 'iProph', 'xinteract.iProph.prot.xml'),
        l = join('/data', basename(OUT_DIR), basename(LOGS_DIR), '{engine}_ProtProph.log'),
        q = join('/data', basename(OUT_DIR), '{engine}', 'TMT_quant.tsv'),
        c = "/data/params/condition.xml"
    singularity:
        "docker://spctools/tpp"
    threads:
        16
    resources:
        mem_mb=32000
    message:
        """  ------  Runing ProteinProphet search for combined mzMLs from {wildcards.engine} ------  """
    shell:
        'ProteinProphet'
        ' {params.i}'
        ' {params.o}'
        ' IPROPHET'
        ' NOGROUPWTS'
        ' PLOTPNG'
        ' EXCELPEPS'
        ' EXCEL.90'
        ' LIBRA{params.c}'
        ' > {params.l} 2>&1'
        ' && mv quantitation.tsv {params.q}'

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule InterProphet_combined:
    input:
        expand(join(OUT_DIR, '{engine}', 'PepProph', 'xinteract.pep.xml'), engine = ENGINES)
    output:
        iProph = join(OUT_DIR, 'combined', 'iProph', 'xinteract.iProph_comb.pep.xml')
    params:
        ic = join('/data', basename(OUT_DIR), 'comet', 'PepProph', 'xinteract.pep.xml'),
        it = join('/data', basename(OUT_DIR), 'tandem', 'PepProph', 'xinteract.pep.xml'),
        o = join('/data', basename(OUT_DIR), 'combined', 'iProph', 'xinteract.iProph_comb.pep.xml'),
        l = join('/data', basename(OUT_DIR), basename(LOGS_DIR), 'combined_iProph.log'),
        d = join('/data', basename(OUT_DIR), 'database_plus_decoy.fa')
    singularity:
        "docker://spctools/tpp"
    threads:
        16
    resources:
        mem_mb=32000
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

rule ProteinProphet_combined:
    input:
        iProph = join(OUT_DIR, 'combined', 'iProph', 'xinteract.iProph_comb.pep.xml')
    output:
        iProph = join(OUT_DIR, 'combined', 'iProph', 'xinteract.iProph_comb.prot.xml')
    params:
        i = join('/data', basename(OUT_DIR), 'combined', 'iProph', 'xinteract.iProph_comb.pep.xml'),
        o = join('/data', basename(OUT_DIR), 'combined', 'iProph', 'xinteract.iProph_comb.prot.xml'),
        l = join('/data', basename(OUT_DIR), basename(LOGS_DIR), 'combined_ProtProph.log'),
        q = join('/data', basename(OUT_DIR), 'combined', 'TMT_quant.tsv'),
        c = "/data/params/condition.xml",
    singularity:
        "docker://spctools/tpp"
    threads:
        16
    resources:
        mem_mb=32000
    shell:
        'ProteinProphet'
        ' {params.i}'
        ' {params.o}'
        ' IPROPHET'
        ' NOGROUPWTS'
        ' PLOTPNG'
        ' EXCELPEPS'
        ' EXCEL.90'
        ' LIBRA{params.c}'
        ' > {params.l} 2>&1'
        ' && mv quantitation.tsv {params.q}'
