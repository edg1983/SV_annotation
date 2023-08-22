#!/usr/bin/env python3

'''
SV annotation
Take a SV VCF and annotate it using fields from other VCF or BED files
Can add gene based annotations as well using a GFF file

Author: Edoardo Giacopuzzi
Email: edoardo.giacopuzzi@well.ox.ac.uk

Each variant type is annotated only from the corresponding type
INS:ME / DEL:ME are annotated with data from INS / DEL
BND annotated by making a flanking region defined in config
'''

import json
import sys, os
import gzip
import argparse
import subprocess
import tempfile
import re
import pandas as pd
from datetime import datetime
from collections import OrderedDict
from cyvcf2 import VCF, Writer
#sys.path.append("/home/edg1983/Servers/well/gel/HICF2/software/BRC_tools")
#from pyCore.Core import get_stdout, tokenize
#from pyCore.VCFtools import INFO


VERSION = "0.2"
COL_OPERATIONS = {
    'String': 'unique',
    'Float': 'max',
    'Integer': 'max'
}
SV_TYPES = ['DEL','DUP','INS','INV','BND']
bedtools = ['bedtools','intersect','-wa','-wb']
coverage = ['bedtools','coverage']

class INFO():
	def __init__(self, info_field):
		exp = re.compile('(.+)=(.+)')
		self.infos = OrderedDict()
		info_tags = info_field.split(";")
		for value in info_tags:
			m = exp.match(value)
			if m:
				self.infos[m.group(1)] = m.group(2)
			else:
				self.infos[value] = True
	
	def selectINFO(self, info_keys):
		if isinstance(info_keys, list):
			self.infos = { key: self.infos[key] for key in info_keys }
		else:
			self.infos = { info_keys : self.infos[info_keys]}
	
	def dropINFO(self, info_keys):
		if isinstance(info_keys, list):
			for key in info_keys: self.infos.pop(key)
		else:
			self.infos.pop(key)
    
	def addINFO(self, info_key, value):
		self.infos[info_key] = value
    
	def updateINFO(self, info_key, value):
		if info_key in self.infos.keys():
			self.infos[info_key] = value
		else:
			self.addINFO(info_key, value)
	
	def getINFO(self,format):
		if format == "dict":
			return self.infos
		elif format == "string":
			output = [] 
			for key, value in self.infos.items():
				output.append(key + "=" + str(value))
		return ';'.join(output)

def log(message, level = "INFO"):
    timestring = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestring}] - {level} - {message}", file=sys.stderr)

def get_stdout(cmd, shell=False):
	proc = subprocess.Popen(cmd, shell=shell,
		stdout = subprocess.PIPE,
		stderr = subprocess.DEVNULL)
    
	for line in iter(proc.stdout.readline, b''):
		yield line.decode('utf-8')

def tokenize(line,sep):
	line = line.rstrip('\n')
	line = line.split(sep)
	return line

def reciprocalOpt(option):
    if option == "TRUE": 
        reciprocal = ['-r'] 
    else:
        reciprocal = []
    return reciprocal

def computeOverlap(dataset, group, a_file, b_file, overlap, reciprocal, cols_idx, fields, tmp_dir):
    filename = tmp_dir + "/SVannot.tmp_" + dataset + "_" + group
    out_file = open(filename, "w+")
    header = ['ID'] + [dataset + "_" + x for x in fields.split(",")]
    out_file.write("\t".join(header) + "\n")
    cols_idx = [3] + [int(x) + 3 for x in cols_idx.split(",")]

    for line in get_stdout(bedtools + overlap + reciprocal + ['-a',a_file,'-b',b_file]):
        line = tokenize(line, "\t")
        out_file.write("\t".join([line[i] for i in cols_idx]) + "\n")  
    out_file.close()    

    mydf = pd.read_csv(filename, sep="\t", index_col="ID")
    return mydf

def computeCoverage(dataset, group, a_file, b_file, tmp_dir):
    filename = tmp_dir + "/SVannot.tmp_cov_" + dataset + "_" + group
    out_file = open(filename, "w+")
    header = ['ID'] + [dataset + "_overlap"]
    out_file.write("\t".join(header) + "\n")
    cols_idx = [3, 7]
    
    log(f"Coverage command: {coverage + ['-a',a_file,'-b',b_file]}")

    for line in get_stdout(coverage + ['-a',a_file,'-b',b_file]):
        line = tokenize(line, "\t")
        out_file.write("\t".join([line[i] for i in cols_idx]) + "\n")  
    out_file.close()    

    mydf = pd.read_csv(filename, sep="\t", index_col="ID")
    return mydf

def setDataType(dataset, fields, data_types):
    tags = [dataset + "_" + x for x in fields.split(",")]
    data_types = data_types.split(",")
    if len(data_types) == 1:
        tags_dict = {k:data_types[0] for k in tags}
    else:
        tags_dict = {k:v for k,v in zip(tags, data_types)}
    
    return tags_dict

def readFileConfig(file_config):
    group = file_config[0]
    b_file = file_config[1]
    cols_idx = file_config[2] if len(file_config) > 2 else None
    fields = file_config[3] if len(file_config) > 3 else None
    data_types = file_config[4] if len(file_config) > 4 else None
    return group,b_file,cols_idx,fields,data_types

def prepare_tmp_files(vcf_file,sv_types,cipos,default_pad,tmp_folder):
    tmp_files = {}
    out = {}
    pads = {}
    for t in sv_types:
        tmp_files[t] = tmp_folder + "/SVannot.tmp_" + t + ".bed"
        out[t] = open(tmp_files[t], "a")
        pads[t] = {'minus': 0, 'plus': 0}
    
    skipped = 0
    counts = {t:0 for t in sv_types}

    vcf = VCF(vcf_file)
    for v in vcf:
        svtype = v.INFO.get(config['SVTYPE'])
        alt = v.ALT[0]
        if cipos:
            cipos_value = v.INFO.get(config['CIPOS'])
            if cipos_value:
                if isinstance(cipos_value, str): 
                    cipos_value = cipos_value.split(",")
                if len(cipos_value) == 2:
                    for t in sv_types: pads[t]['minus'] = abs(int(cipos_value[0]))
                    for t in sv_types: pads[t]['plus'] == abs(int(cipos_value[1]))
                else:
                    pads['INS']['minus'] = pads['INS']['plus'] = pads['BND']['minus'] = pads['BND']['minus'] = default_pad
            else:
                pads['INS']['minus'] = pads['INS']['plus'] = pads['BND']['minus'] = pads['BND']['minus'] = default_pad
        else:
            pads['INS']['minus'] = pads['INS']['plus'] = pads['BND']['minus'] = pads['BND']['minus'] = default_pad
        
        for t in sv_types:
            if svtype == t or re.match(t,alt): 
                svtype = t

        start = int(v.POS) - pads[svtype]['minus']
        stop = int(v.POS) + pads[svtype]['plus']

        try:
            out[svtype].write("{chrom}\t{start}\t{stop}\t{ID}\n".format(chrom=v.CHROM,start=start,stop=stop,ID=v.ID))
            counts[svtype] += 1
        except:
            skipped += 1

    return(tmp_files, counts, skipped)

def checkFiles(files):
    not_found = []
    success = True
    for f in files:
        if not os.path.isfile(f):
            not_found.append(f)
            success = False
    return (success, not_found)

#arguments
parser = argparse.ArgumentParser(description='Script to add annotSV annotations to VCF file')
parser.add_argument("-i", "--inputvcf", help="Input VCF file to be annotated", action="store", required=True)
parser.add_argument("-o", "--out", help="Output VCF file", action="store", required=False, default="stdout")
parser.add_argument("-t", "--tmpdir", help="Folder to store temp files", action="store", required=False, default=tempfile.gettempdir())
parser.add_argument("-b", "--build", help="Genome build", action="store", choices=["GRCh37", "GRCh38"], required=True)
parser.add_argument("-s", "--config", help="Config file (json)", action="store", required=True)
parser.add_argument("-r", "--remove_tmp", help="Set to remove tmp files", action="store_true", required=False)
args = parser.parse_args()

#Startup message
log(f"### SV annotation {VERSION} ###")

#Initial configuration
input_file = args.inputvcf
build = args.build
json_file = args.config
tmp_folder = args.tmpdir
os.makedirs(tmp_folder, exist_ok=True)

if args.out == "-":
    log("Writing output to stdout")
    outfile = sys.stdout
else:
    log(f"Writing output to {args.out}")
    outfile = open(args.out, "w")

with open(json_file) as json_config:
    config = json.load(json_config)

#If resource_folder is defined, this path is added to all files
resource_folder = config.get('resource_folder', False)
if resource_folder:
    for section in ['AF_datasets','custom_datasets','genes','coverage_datasets']:
        for dataset in config[section][build].keys():
            for i in range(len(config[section][build][dataset])): 
                config[section][build][dataset][i][1] = resource_folder + "/" + config[section][build][dataset][i][1]

#Check all files from config exists
input_files = [input_file]
for g in ['AF_datasets','custom_datasets','genes','coverage_datasets']:
    for v in config[g][build].values():
        for f in v: input_files.append(f[1])
success, not_found = checkFiles(input_files)
if not success:
    sys.exit("ERROR. Following input files not found: " + ";".join(not_found))

#Read input VCF and create temp bed files
log("Prepare intervals from input VCF")
tmp_files, counts, skipped = prepare_tmp_files(input_file, SV_TYPES, config.get('CIPOS', None),config.get('padding',10),tmp_folder)
for key, value in counts.items(): log(f"{key}: {value}")
log(f"Skipped vars: {skipped}")

#Run bedtools intersect on the single datasets
annotations = []
tags_dataTypes = {}

log("Processing AF datasets")
log(f"Overlap: {config['overlap']['AF_datasets'][0]}")
log(f"Reciprocal: {config['overlap']['AF_datasets'][1]}")

for dataset, items in config['AF_datasets'][build].items():
    log(dataset)
    for file_config in items:
        group,b_file,cols_idx,fields,data_types = readFileConfig(file_config)
        a_file = tmp_files[group]
        if group in ['BND', 'INS']:
            overlap = []
            reciprocal = []
        else:      
            overlap = ['-f', config['overlap']['AF_datasets'][0]]
            reciprocal = reciprocalOpt(config['overlap']['AF_datasets'][1]) 

        log (f"{group}: {b_file}")

        annotations.append(computeOverlap(dataset,group,a_file,b_file,overlap,reciprocal,cols_idx,fields,tmp_folder))
        tags_dataTypes.update(setDataType(dataset,fields,data_types))

log("Processing coverage datasets")
for dataset, items in config['coverage_datasets'][build].items():
    log(dataset)
    for file_config in items:
        group,b_file,cols_idx,fields,data_types = readFileConfig(file_config)
        a_file = tmp_files[group]
        log (f"{group}: {b_file}")

        annotations.append(computeCoverage(dataset,group,a_file,b_file,tmp_folder))
        tags_dataTypes.update(setDataType(dataset,"overlap","Float"))

log("Processing custom datasets")
log(f"Overlap: {config['overlap']['custom_datasets'][0]}")
log(f"Reciprocal: {config['overlap']['custom_datasets'][1]}")
 
for dataset, items in config['custom_datasets'][build].items():
    log(dataset)
    for file_config in items:
        group,b_file,cols_idx,fields,data_types = readFileConfig(file_config)
        a_file = tmp_files[group]
        if group in ['BND', 'INS']:
            overlap = []
            reciprocal = []
        else:   
            overlap = ['-F', config['overlap']['custom_datasets'][0]]
            reciprocal = reciprocalOpt(config['overlap']['custom_datasets'][1])

        log (f"{group}: {b_file}")

        annotations.append(computeOverlap(dataset,group,a_file,b_file,overlap,reciprocal,cols_idx,fields,tmp_folder))
        tags_dataTypes.update(setDataType(dataset,fields,data_types))

log("Processing genes datasets")
log(f"Overlap: {config['overlap']['genes'][0]}")
log(f"Reciprocal: {config['overlap']['genes'][1]}")


for dataset, items in config['genes'][build].items():
    log(dataset)
    for file_config in items:
        group,b_file,cols_idx,fields,data_types = readFileConfig(file_config)
        a_file = tmp_files[group]
        if group in ['BND', 'INS']:
            overlap = []
            reciprocal = []
        else:  
            overlap = ['-F', config['overlap']['genes'][0]]
            reciprocal = reciprocalOpt(config['overlap']['genes'][1]) 
        
        log (f"{group}: {b_file}")

        annotations.append(computeOverlap(dataset,group,a_file,b_file,overlap,reciprocal,cols_idx,fields,tmp_folder))
        tags_dataTypes.update(setDataType(dataset,fields,data_types))

#Concatenate all annotations in a single df and replace NA with 0
log("Concatenating annotations...")
full_annots = pd.concat(annotations, axis=0, join='outer', sort=False)
full_annots.fillna(0, inplace=True)

#Perform operations on columns to generate a unique annotation per variant
#Currently this annotate max AF and unique distinct genes
log("Collapsing annotations...")
col_operations = {}
for col in full_annots.columns:
    col_operations[col] = COL_OPERATIONS[tags_dataTypes[col]]

log(f"Columns operations: {col_operations}")
full_annots = full_annots.groupby('ID').agg(col_operations)
full_annots.index = full_annots.index.map(str)

log("==== Annotations extract ====")
print(full_annots.head(), file=sys.stderr)
log("=============================")

#Add the annotations and write new VCF
log("Adding annotations to VCF...")
if args.inputvcf.endswith("vcf.gz"):
    vcf = gzip.open(args.inputvcf,"rt")
elif args.inputvcf.endswith("vcf"):
    vcf = open(args.inputvcf,"r")

line = vcf.readline()
while line.startswith('##'):
    outfile.write(line)
    line = vcf.readline()

#Add INFO tags descriptions to header
for tag, data_type in tags_dataTypes.items():
    header_line = '##INFO=<ID={tag},Number=1,Type={datatype},Description="SV annotation">'.format(tag=tag,datatype=data_type)
    outfile.write(header_line + "\n")

while line.startswith('#'):
    outfile.write(line)        
    line = vcf.readline()

#Annotate vars
nline = 0
annotated = 0
skipped = 0
while line:
    new_annotation_added = False
    nline += 1
    if nline % 10000 == 0: 
        log(f"{nline} variants annotated")
    line = tokenize(line,"\t")
    ID = line[2]
    infos = INFO(line[7])
    for col in full_annots.columns:
        try:
            newannot = full_annots.loc[ID,col]
            if isinstance(newannot, (float, int)):
                newannot = str(round(newannot,4))
                infos.addINFO(col, newannot)
                new_annotation_added = True
            else:
                newannot = [str(x) for x in newannot if x != 0]
                if len(newannot) > 0:
                    newannot = ",".join(newannot)
                    infos.addINFO(col, newannot)
                    new_annotation_added = True
            newline = line[0:7]
            newline.append(infos.getINFO("string"))
            newline.extend(line[8:])   
        except:
            newline = line
    if new_annotation_added:
        annotated += 1
    else:
        skipped += 1
    outfile.write("\t".join(newline) + "\n")
    line = vcf.readline()

outfile.close()
print("", file=sys.stderr)
log(f"Variants annotated: {annotated}; not annotated: {skipped}")

if args.remove_tmp:
    log("Removing temp files")
    tmp_files = os.listdir(tmp_folder)
    for myfile in tmp_files:
        if myfile.startswith('SVannot.tmp_'):
            os.remove(tmp_folder + "/" + myfile)

log("Annotations done!")