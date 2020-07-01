# SV annotation

Tool to annotate standard SV VCF based on region overlap
Configuration of regions is provided with a json file.

Works with VCF files for popular caller like Manta, Lumpy, Delly.
It should work with any SV VCF as long as it contains
- unique var ID in the ID colum
- an SV type and end coordinate annotation in INFO field (the exact tag can be configured)

Files usedo for annotation are standard BED. 
You can configure location of files and column to annotate from in the json
