import csv
import networkx as nx
import argparse, sys
import os
import json

# Usage:
#        NetworkPPI.py -c <input candidate file> -a < annotations files> -n <network edge list> [ -o <output> -d <degree> -e <edge_dbcount>]

# Required arguments:
#   -c <input candidate file>     path to candidate input file. 
#   -a <annotations files>        path to directory with annotation files. 
#   -n <network edge list>        path to network edge list.

                                                    
# Options:
#   -d <degree>     Remove nodes with degree less than the specified value
#   -e <edge_dbcount>     Remove edges with a DBCount less than the specified value 
#   -o <output>     name of output .json file. Default network.json

parser = argparse.ArgumentParser()

parser.add_argument('--candidate','-c', dest="candidate")
parser.add_argument('--annotations','-a', dest="annotations")
parser.add_argument('--network','-n', default ='20211215.combinedDB.ConsensusPathDB_v35_StringDB_v11.5.GIANTv2.simplified.tsv', dest="network")
parser.add_argument('--degree','-d', dest="degree")
parser.add_argument('--edge_dbcount','-e', dest="dbcount")
parser.add_argument('--output','-o', default='network_test.json', dest="output")


args = parser.parse_args()

# Load in Edge List.
g = nx.read_edgelist(args.network,create_using=nx.Graph(), nodetype = str,data=(("consensus", bool),("STRING", bool),("GIANT", bool),("DBCount", int),))



#Function to read in a directory of annotation files.
#Inputs g - network, path - path to directory
def dir_annotation(g,path):
    attrs = dict()
    annotation_names = []
    for annotation_file in os.listdir(path):
        attrs.clear()
        with open(path + "/"+ annotation_file, "r") as file:
            reader = csv.reader(file, delimiter='\t')
            header = next(reader)
            for row in reader:
                gene = row[0]
                if gene in g:
                    if len(row) == 1:
                        attrs[gene] = 1
                    else:
                        attrs[gene] = row[1]
            annotation_names.append(header[0])
            nx.set_node_attributes(g,attrs, header[0])
    return annotation_names

#Function to read in a directory of candidate files
#Inputs g - network, path - path to directory
def file_annotation(g, path):
    attrs = dict()
    with open(path, "r") as file:
        reader = csv.reader(file, delimiter = '\t')
        header = next(reader)
        for row in reader:
            gene = row[0]
            if gene in g:
                if len(row) == 1:
                    attrs[gene] = 1
                else:
                    attrs[gene] = row[1]
        nx.set_node_attributes(g,attrs,header[0])
    return header[0]


candidate_name = 'candidate'
#Read in candidate files if directory
if os.path.isdir(args.candidate):
    candidate_name = dir_annotation(g,args.candidate)
#Read in candidate file
elif os.path.isfile(args.candidate):
    candidate_name = file_annotation(g,args.candidate)

#Assign annotation values.
if args.annotations is not None:
    annotation_names = dir_annotation(g,args.annotations)

#Search for nodes with candidate annotations.
if isinstance(candidate_name,str):
    #If the candidate set is from a single file.
    candidate_dict = nx.get_node_attributes(g, candidate_name)
    candidates = list(candidate_dict.keys())
else:
    #If the candidate set is from a list of files; get all the annotations.
    candidate_set = set()
    for i in candidate_name:
        a = nx.get_node_attributes(g,i)
        candidate_set.update(list(a.keys()))
    candidates = list(candidate_set)



#Select nodes with specific annotation name.
node_annotations = set()
for i in annotation_names:
    #One can do a specific check here for the desired annotations
    a = nx.get_node_attributes(g,i)
    node_annotations.update(list(a.keys()))


#Generate subgraph with nodes from candidates and desired annotation nodes
node_list = list(node_annotations) + candidates
h = g.subgraph(node_list)
subgraph = nx.Graph(h)
#Optional remove isolate nodes
subgraph.remove_nodes_from(list(nx.isolates(subgraph)))

#Optional pruning of edges by database count
if args.dbcount is not None:
    selected_edges = [(u,v) for u,v,e in subgraph.edges(data='DBCount') if e < int(args.dbcount)]
    subgraph.remove_edges_from(selected_edges)
    subgraph.remove_nodes_from(list(nx.isolates(subgraph)))

#Optional pruning of nodes by degree
if args.degree is not None:
    selected_nodes = [node for node,degree in dict(subgraph.degree()).items() if degree < int(args.degree)]
    subgraph.remove_nodes_from(selected_nodes)
    subgraph.remove_nodes_from(list(nx.isolates(subgraph)))


#Convert NetworkX format to JSON
cy = nx.readwrite.json_graph.cytoscape_data(subgraph)

#Set positional data. Required for Cytoscape.
pos = nx.spring_layout(subgraph)
for n,p in zip(cy['elements']['nodes'],pos.values()):
     n['position'] = {'x':2000 * p[0],'y':2000 * p[1]}

#Write to output file.
with open(args.output, "w") as outfile:
     json.dump(cy, outfile)

