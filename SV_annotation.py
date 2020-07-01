'''
SV annotation
Take a SV VCF and annotate it using fields from other VCF or BED files
Can add gene based annotations as well using a GFF file

Author: Edoardo Giacopuzzi
Email: edoardo.giacopuzzi@well.ox.ac.uk
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

COL_OPERATIONS = {
    'String': 'unique',
    'Float': 'max',
    'Integer': 'max'
}

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
    cols_idx = file_config[2]
    fields = file_config[3]
    data_types = file_config[4]
    return group,b_file,cols_idx,fields,data_types

#arguments
parser = argparse.ArgumentParser(description='Script to add annotSV annotations to VCF file')
parser.add_argument("-i", "--inputvcf", help="Input VCF file to be annotated", action="store", required=True)
parser.add_argument("-o", "--out", help="Output VCF file", action="store", required=True)
parser.add_argument("-t", "--tmpdir", help="Folder to store temp files", action="store", required=False, default=tempfile.gettempdir())
parser.add_argument("-b", "--build", help="Genome build", action="store", choices=["GRCh37", "GRCh38"], required=True)
parser.add_argument("-s", "--config", help="Config file (json)", action="store", required=True)
parser.add_argument("-r", "--remove_tmp", help="Set to remove tmp files", action="store_true", required=False)
args = parser.parse_args()

#Initial configuration
input_file = args.inputvcf
build = args.build
json_file = args.config
tmp_folder = args.tmpdir

bedtools = ['bedtools','intersect','-wa','-wb']

with open(json_file) as json_config:
    config = json.load(json_config)

tmp_files = {}
tmp_files['INS'] = tmp_folder + "/SVannot.tmp_INS.bed"
tmp_files['DEL'] = tmp_folder + "/SVannot.tmp_DEL.bed"
tmp_files['DUP'] = tmp_folder + "/SVannot.tmp_DUP.bed"
tmp_files['INV'] = tmp_folder + "/SVannot.tmp_INV.bed"

#Read input VCF and create temp bed files
vcf = VCF(input_file)
skipped = 0
out = {}
counts = {'INS': 0, 'DEL': 0, 'DUP': 0, 'INV': 0}

out['INS'] = open(tmp_files['INS'], "w+")
out['DEL'] = open(tmp_files['DEL'], "w+")
out['DUP'] = open(tmp_files['DUP'], "w+")
out['INV'] = open(tmp_files['INV'], "w+")
for v in vcf:
    try:
        out[v.INFO.get(config['SVTYPE'])].write("{chrom}\t{start}\t{stop}\t{ID}\n".format(chrom=v.CHROM,start=v.POS,stop=v.INFO.get(config['END']),ID=v.ID))
        counts[v.INFO.get(config['SVTYPE'])] += 1
    except:
        skipped += 1
out['INS'].close()
out['DEL'].close()
out['DUP'].close()
out['INV'].close()
for key, value in counts.items(): print(key, ": ", value)
print("Skipped vars: ", skipped)

#Run bedtools intersect on the single datasets
annotations = []
tags_dataTypes = {}

print("Processing AF datasets")
print("Overlap: ", config['overlap']['AF_datasets'][0])
print("Reciprocal: ", config['overlap']['AF_datasets'][1])
overlap = ['-f', config['overlap']['AF_datasets'][0]]
reciprocal = reciprocalOpt(config['overlap']['AF_datasets'][1]) 

for dataset, items in config['AF_datasets'][build].items():
    print(dataset)
    for file_config in items:
        group,b_file,cols_idx,fields,data_types = readFileConfig(file_config)
        a_file = tmp_files[group]
        print ("\t", group, ": ", b_file)

        annotations.append(computeOverlap(dataset,group,a_file,b_file,overlap,reciprocal,cols_idx,fields,tmp_folder))
        tags_dataTypes.update(setDataType(dataset,fields,data_types))

print("Processing custom datasets")
print("Overlap: ", config['overlap']['custom_datasets'][0])
print("Reciprocal: ", config['overlap']['custom_datasets'][1])
overlap = ['-F', config['overlap']['custom_datasets'][0]]
reciprocal = reciprocalOpt(config['overlap']['custom_datasets'][1]) 

for dataset, items in config['custom_datasets'][build].items():
    print(dataset)
    for file_config in items:
        group,b_file,cols_idx,fields,data_types = readFileConfig(file_config)
        a_file = tmp_files[group]
        print ("\t", group, ": ", b_file)

        annotations.append(computeOverlap(dataset,group,a_file,b_file,overlap,reciprocal,cols_idx,fields,tmp_folder))
        tags_dataTypes.update(setDataType(dataset,fields,data_types))

full_annots = pd.concat(annotations, axis=0, join='outer', sort=False)
full_annots.fillna(0, inplace=True)

col_operations = {}
for col in full_annots.columns:
    col_operations[col] = COL_OPERATIONS[tags_dataTypes[col]]

full_annots = full_annots.groupby('ID').agg(col_operations)
full_annots.index = full_annots.index.map(str)

print("Processing genes datasets")
print("Overlap: ", config['overlap']['genes'][0])
print("Reciprocal: ", config['overlap']['genes'][1])
overlap = ['-F', config['overlap']['genes'][0]]
reciprocal = reciprocalOpt(config['overlap']['genes'][1]) 

for dataset, items in config['genes'][build].items():
    print(dataset)
    for file_config in items:
        group,b_file,cols_idx,fields,data_types = readFileConfig(file_config)
        a_file = tmp_files[group]
        print ("\t", group, ": ", b_file)

        annotations.append(computeOverlap(dataset,group,a_file,b_file,overlap,reciprocal,cols_idx,fields,tmp_folder))
        tags_dataTypes.update(setDataType(dataset,fields,data_types))

full_annots = pd.concat(annotations, axis=0, join='outer', sort=False)
full_annots.fillna(0, inplace=True)

col_operations = {}
for col in full_annots.columns:
    col_operations[col] = COL_OPERATIONS[tags_dataTypes[col]]

full_annots = full_annots.groupby('ID').agg(col_operations)
full_annots.index = full_annots.index.map(str)

#Add the annotations and write new VCF
print("Adding annotations to VCF...")
outfile = open(args.out, "w")
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
    nline += 1
    print("Variant ",nline, end="\r")
    line = tokenize(line,"\t")
    ID = line[2]
    infos = INFO(line[7])
    for col in full_annots.columns:
        try:
            newannot = full_annots.loc[ID,col]
            if isinstance(newannot, (float, int)):
                newannot = str(round(newannot,4))
                infos.addINFO(col, newannot)
                annotated += 1
            else:
                newannot = [str(x) for x in newannot if x != 0]
                if len(newannot) > 0:
                    newannot = ",".join(newannot)
                    infos.addINFO(col, newannot)
                    annotated += 1
            newline = line[0:7]
            newline.append(infos.getINFO("string"))
            newline.extend(line[8:])   
        except:
            skipped += 1
            newline = line
            
    outfile.write("\t".join(newline) + "\n")
    line = vcf.readline()

outfile.close()
print("\nVariants annotated: ",annotated, "; skipped: ", skipped)

if args.remove_tmp:
    print("Removing temp files")
    tmp_files = os.listdir(tmp_folder)
    for myfile in tmp_files:
        if myfile.startswith('SVannot.tmp_'):
            os.remove(tmp_folder + "/" + myfile)

print("Annotations done!")