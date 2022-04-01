import csv
import networkx as nx
import argparse, sys
import os
import json

# Usage:
#        NetworkPPI.py -c <input candidate file> -a < annotations files> -n <network edge list> [ -o <output> -d <degree> -e <edge_dbcount> -i <int_nodes>]

def arg_parser(argv=None):
    """
    Required arguments:
  -c <input candidate file>     path to candidate input file. 
  -a <annotations files>        path to directory with annotation files. 
  -n <network edge list>        path to network edge list.
  

                                                    
    Options:
  -d <degree>     Remove nodes with degree less than the specified value
  -e <edge_dbcount>     Remove edges with a DBCount less than the specified value 
  -i <int_nodes>   Flag to remove annotation nodes that do not connect to at least 2 candidate nodes
  -o <output>     name of output .json file. Default network.json
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('--candidates','-c', dest="candidates")
    parser.add_argument('--annotations','-a', dest="annotations")
    parser.add_argument('--network','-n', default ='20211215.combinedDB.ConsensusPathDB_v35_StringDB_v11.5.GIANTv2.simplified_filtered.tsv', dest="network")
    parser.add_argument('--degree','-d', dest="degree")
    parser.add_argument('--edge_dbcount','-e', dest="dbcount")
    parser.add_argument('--int_nodes','-i', default = False,dest="intermediate_nodes")
    parser.add_argument('--output','-o', default='network_out.json', dest="output")
    return parser.parse_args(argv)


def dir_annotation(g,path):
    """
    Read in a directory of files.
    
    Parameters
    g: networkX Graph Object
    path: Path to Directory

    Returns annotation name
    """
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

def file_annotation(g, path):
    """
    Read in a single file

    Parameters
    g: networkX Graph Object
    path: Path to Directory

    Returns annotation name
    """
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

def clean_intermediate_nodes(g,candidate_names,annotation_names):
    """
    Remove annotation nodes if they do not connect to at least two candidate nodes.

    Parameters
    g:networkX Graph Object
    candidate_nodes: Set of Candidates nodes
    annotation_nodes: Set of Annotation nodes 

    Returns networkX object
    """

    candidate_nodes = get_named_nodes(subgraph,candidate_names)
    annotation_nodes = get_named_nodes(subgraph,annotation_names)

    nodes_to_delete = set()

    #Nodes that are both part of the candidate and annotation set are not removed
    intermediates = annotation_nodes.difference(candidate_nodes)
    
    for i in intermediates:
        counter = 0
        for n in g.neighbors(i):
            if n in candidate_nodes:
                counter += 1
            if counter == 2:
                break
        if counter < 2:
            nodes_to_delete.add(i)
    g.remove_nodes_from(nodes_to_delete)
    return g

def get_named_nodes(g,name):
    """
    Get Nodes with annotation label

    Parameters
    g: networkX Graph Object
    name: list of names with desired labels

    Returns list of nodes
    """
    node_set = set()

    for i in name:
        a = nx.get_node_attributes(g,i)
        node_set.update(list(a.keys()))
    return node_set

def prune_network(g,degree_filter,edge_filter):
    """
    Prune network by node degree or edge value

    Parameters
    g: networkX Graph Object
    degree_filter(int): Nodes with degree less than set value are removed
    edge_filter(int): Edges with DBCount less than set value are removed

    Returns pruned networkX object
    """

    #Optional pruning of edges by database count
    if edge_filter is not None:
        selected_edges = [(u,v) for u,v,e in subgraph.edges(data='DBCount') if e < int(edge_filter)]
        subgraph.remove_edges_from(selected_edges)
        subgraph.remove_nodes_from(list(nx.isolates(subgraph)))

    #Optional pruning of nodes by degree
    if degree_filter is not None:
        selected_nodes = [node for node,degree in dict(subgraph.degree()).items() if degree < int(degree_filter)]
        subgraph.remove_nodes_from(selected_nodes)
        subgraph.remove_nodes_from(list(nx.isolates(subgraph)))
    return subgraph


def generate_network(args):
    """
    Generate NetworkX file

    Parameters:
    args: arg_parser arguments (See arg_parser for details)
    
    Returns networkX Object
    """
    node_list = set()
    if os.path.isdir(args.candidates):
        candidate_names = dir_annotation(g,args.candidates)
        candidates = get_named_nodes(g,candidate_names)

    #Read in candidate file
    elif os.path.isfile(args.candidates):
        candidate_names = []
        candidate_names.append(file_annotation(g,args.candidates))

        candidates = get_named_nodes(g,candidate_names)


  
    node_list = node_list.union(candidates)

    #Assign annotation values.
    if args.annotations is not None:
        annotation_names = dir_annotation(g,args.annotations)
        annotation_nodes = get_named_nodes(g, annotation_names)

        node_list = node_list.union(annotation_nodes)

    if not len(node_list):
        sys.exit("Candidate and/or annotations could not be read")

    #Generate subgraph with nodes from candidates and desired annotation nodes

    h = g.subgraph(node_list)
    subgraph = nx.Graph(h)

    #Optionally Remove Isolated Nodes from file
    #subgraph.remove_nodes_from(list(nx.isolates(subgraph)))

    return subgraph




if __name__ == "__main__":
    args = arg_parser()

    # Load in Edge List.
    g = nx.read_edgelist(args.network,create_using=nx.Graph(), nodetype = str,data=(("consensus", bool),("STRING", bool),("GIANT", bool),("DBCount", int),))
    
    subgraph = generate_network(args)
    subgraph = prune_network(subgraph, args.degree, args.dbcount)

    if args.intermediate_nodes:
        subgraph = clean_intermediate_nodes(subgraph,args.candidates,args.annotations)

    #Convert NetworkX format to JSON
    cy = nx.readwrite.json_graph.cytoscape_data(subgraph)
    #Set positional data. Required for Cytoscape.
    pos = nx.spring_layout(subgraph)
    for n,p in zip(cy['elements']['nodes'],pos.values()):
         n['position'] = {'x':2000 * p[0],'y':2000 * p[1]}

    #Write to output file.
    with open(args.output, "w") as outfile:
         json.dump(cy, outfile)
