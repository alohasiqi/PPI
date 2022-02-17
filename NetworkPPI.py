import csv
import networkx as nx
import argparse, sys
import os


# Usage:
#        NetworkPPI.py -c <input candidate file> -a < annotations files> -n <network edge list> [ -o <output> ]

# Required arguments:
#   -c <input candidate file> path to candidate input file. 
#   -a <annotations files>        path to directory with annotation files. 
#   -n <network edge list>        path to network edge list. 
                                                    
# Options:
#   -o <output>     name of output .json file. Default network.json

parser = argparse.ArgumentParser()

parser.add_argument('--candidate','-c', dest="candidate")
parser.add_argument('--annotations','-a', dest="annotations")
parser.add_argument('--network','-n', dest="network")
parser.add_argument('--output','-o', default='network.json', dest="output")

args = parser.parse_args()

# Load in Edge List.
g = nx.read_edgelist(args.network,create_using=nx.Graph(), nodetype = str,data=(("consensus", bool),("STRING", bool),("GIANT", bool),("DBCount", int),))


# Assign edge colors based on number of databases that indicate a PPI connection
edges = g.edges()
for u,v in edges:
    dbcount = g[u][v]['DBCount']
    if dbcount == 1:
        g[u][v]["color"] = "red"
    elif dbcount == 2:
        g[u][v]["color"] = "green"
    elif dbcount == 3:
        g[u][v]["color"] = "blue"


#Assign node label with candidate
attrs = dict()
with open(args.candidate, "r") as file:
    for line in file:
        gene = line.rstrip()
        if gene in g:
            attrs[gene] = 1
    nx.set_node_attributes(g,attrs,"candidate")


#Assign annotation values for the remaining nodes. Second column name becomes name of the annotation
for annotation_file in os.listdir(args.annotations):
    attrs.clear()
    with open(annotation_file, "r") as file:
        reader = csv.reader(file)
        header = next(reader)
        for row in reader:
            gene = row[0]
                if gene in g:
                attrs[gene] = row[1]

        nx.set_node_attributes(g,attrs, header[1])


#Select all candidate nodes
candidates = nx.get_node_attributes(g,'candidate')

#Select nodes with annotation name 'utr'
utr = nx.get_node_attributes(g,'utr')


node_list = list(utr.keys()) + list(candidates.keys())


#Generate subgraph with only nodes from node_list
h = g.subgraph(node_list)
a = nx.Graph(h)
#Optional remove isolate nodes
a.remove_nodes_from(list(nx.isolates(a)))

pos = nx.spring_layout(a)


#Convert NetworkX format to JSON
cy = nx.readwrite.json_graph.cytoscape_data(a)

#Set positional data. Required for Cytoscape.
for n,p in zip(cy['elements']['nodes'],pos.values()):
     n['position'] = {'x':2000 * p[0],'y':2000 * p[1]}

#Write to output file.
with open(args.output, "w") as outfile:
     json.dump(cy, outfile)

