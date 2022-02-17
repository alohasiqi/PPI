# network-ppi

This program is a pipeline to generate JSON and XML files for use with Cytoscape for visualizing a protein-protein interaction network.

## Installation
Python NetworkX

Cytoscape and the enhancedGraphics Package
The enhancedGraphics Package is used for multi-color nodes and can be installed from the link:
https://apps.cytoscape.org/apps/enhancedgraphics


## Usage
```
Usage:
       NetworkPPI.py -c <input candidate file> -a < annotations files> -n <network edge list> [ -o <output> ]

Required arguments:
  -c <input candidate file>	path to candidate input file. 
  -a <annotations files>		path to directory with annotation files. 
  -n <network edge list>		path to network edge list. 
                                                    
Options:
  -o <output>     name of output .json file. Default network.json
```

## Input Parameter File

### Candidate Input File

A list of genes separated by newlines. Genes in this list always appear in the network.

### Annotation Input
**Tab-delimited** file with annotations with the following columns. Nodes will be assigned annotation_name as its label.
 - Gene name
 - annotation_name (values must be binary 1 or 0 to be graphed)
 
Gene&nbsp;&nbsp;annotation_name

ABAT&nbsp;&nbsp;1

ABCA10&nbsp;&nbsp;1
 

### Network File Input
**Tab-delimited** edge list file. No header required.
- Gene Node 1
- Gene Node 2
- Database Edge (1 - Edge exists, 0 - Edge does not)
- Database Edge
- Database Edge
- DBCount (integer)


Gene_1&nbsp;&nbsp;Gene_2&nbsp;&nbsp;consensusPathDB&nbsp;&nbsp;SFARI&nbsp;&nbsp;GIANT

CACNG6&nbsp;&nbsp;RYR1&nbsp;&nbsp;0&nbsp;&nbsp;1&nbsp;&nbsp;0&nbsp;1

H4C13&nbsp;&nbsp;TP53BP1&nbsp;&nbsp;1&nbsp;&nbsp;1&nbsp;&nbsp;0&nbsp;2



