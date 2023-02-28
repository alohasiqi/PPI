import pandas as pd
from collections import Counter
import csv
import networkx as nx
import argparse, sys
import os
import json
import warnings
from pathlib import Path


# Usage:
#        NetworkPPI.py -c <input candidate file> -a < annotations files> -n <network edge list> [ -o <output> -d <degree> -i <int_nodes> -e <edge_dbcount> -m <connection_number> -r <rm_dupintermediates>]

def arg_parser(argv=None):
    """
    Required arguments:
  -c <input candidate file>     path to candidate input file. 
  -a <annotations files>        path to directory with annotation files. 
  -n <network edge list>        path to network edge list.
  

                                                    
    Options:
  -d <degree>     Remove nodes with degree less than the specified value
  -e <edge_dbcount>     Remove edges with a DBCount less than the specified value
  -dg <degree_greater> Remove edges with a degree greater than the specified value 
  -i <int_nodes>   Flag to remove annotation nodes that do not connect to at least 2 candidate nodes
  -s <subgraph_id> Specific specific candidate gene lists to graph
  -o <output>     name of output .json file. Default network.json
  -m <connection_number>    Add intermediate nodes to the graph that connect to at least(>) a number of candidate genes (connection_count)
  -r <rm_dupintermediates>    Remove intermediate nodes that have the same connections to the candidate/desired nodes and randomly keep one representative intermediate node, valid only when -m <connection_number>  is not None
  """
    parser = argparse.ArgumentParser()

    parser.add_argument("--candidates","-c", dest="candidates", type=Path)
    parser.add_argument("--annotations","-a", dest="annotations", type=Path)
    parser.add_argument("--network","-n", dest="network", type=Path)
    parser.add_argument("--degree","-d", dest="degree", type=int, default=None)
    parser.add_argument("--edge_dbcount","-e", dest="dbcount", type=int, default=None)
    parser.add_argument("--connection_number","-m", dest="connection_number", type=int, default=None)
    parser.add_argument("--rm_dupintermediates","-r", dest="dedup", type=bool, default=False)
    parser.add_argument("--degree_greater","-dg", dest="degree_greater", type=int, default=None)
    parser.add_argument("--int_nodes","-i", action="store_true",dest="intermediate_nodes", default = False)
    parser.add_argument("--subgraph_id","-s",dest="subname", default=None)
    parser.add_argument("--output","-o",dest="output", type=Path, default="network_out.json")
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
        file = path / f'{annotation_file}'
        with open(file, "r") as file:
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

    candidate_nodes = get_named_nodes(g,candidate_names)
    annotation_nodes = get_named_nodes(g,annotation_names)

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
    g.remove_nodes_from(list(nx.isolates(g)))



    annotation_nodes = get_named_nodes(g,annotation_names)
    intermediates = annotation_nodes.difference(candidate_nodes)

    edges = g.edges(intermediates)

    selected_edges = set()
    for edge in edges:
        if edge[0] in intermediates and edge[1] in intermediates:
            selected_edges.add(edge)
    g.remove_edges_from(selected_edges)
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

def prune_network(g,degree_filter,edge_filter,degree_greater):
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
        selected_edges = [(u,v) for u,v,e in g.edges(data='DataBaseCount') if int(e) < int(edge_filter)]
        g.remove_edges_from(selected_edges)
        g.remove_nodes_from(list(nx.isolates(g)))
    
    if degree_greater is not None:
        selected_nodes = [node for node,degree in dict(g.degree()).items() if degree > int(degree_greater)]
        g.remove_edges_from(selected_edges)
        g.remove_nodes_from(list(nx.isolates(g)))

    #Optional pruning of nodes by degree
    if degree_filter is not None:
        selected_nodes = [node for node,degree in dict(g.degree()).items() if degree < int(degree_filter)]
        g.remove_nodes_from(selected_nodes)
        g.remove_nodes_from(list(nx.isolates(g)))
    return g

def allow_intermediate(df_edge_all,node_list,connection_number,dedup):
    """
    Optional: add intermediate nodes to the graph that connect to a number of candidate genes (connection_count)
    None means no intermediate nodes

    Parameters:
    df_edges_all: the interaction dataframe generated from the network
    node_list: the candidate and desired nodes generated
    connection_number: args.connection_number
    dedup:args.dedup
    
    Returns 
    df_edge_filt: the interaction dataframe that contains intermediated nodes
    intermediate_nodes: the selected intermediate nodes
    """
    #get edges with more than connection_number candiates or desired annotation nodes
    df_edge = df_edge_all[df_edge_all['node1'].isin(node_list) | df_edge_all['node2'].isin(node_list)]
    #generate a new dataframe with node1 (candiate genes, genes in the list), node2 (intermediate genes)
    df_edge = df_edge.copy()
    df_edge["genenode1"] = df_edge['node1'].isin(node_list)
    df_edge_sort = df_edge[df_edge["genenode1"]==True]
    df_edge_modify = df_edge[df_edge["genenode1"]==False]
    df_edge_modify = df_edge_modify.copy()
    df_edge_modify.rename(columns={"node1":"node2", "node2":"node1"}, inplace=True)
    network_clean = pd.concat([df_edge_sort, df_edge_modify])
    network_clean.drop(columns=["genenode1"], inplace=True)
    
    #count neighbors of intermediate nodes and select intermediate nodes with greater than connection_count interactions with candidates/desired annotation nodes
    nodes_counter = Counter(network_clean['node2'])
    intermediate_nodes = [i for i in nodes_counter if (nodes_counter[i] > int(connection_number))] 
    
    df_edge_filt = network_clean[network_clean['node2'].isin(intermediate_nodes)]
    
    if dedup:
        #Remove intermediate nodes with the same connections (i.e.connect with the same targeted genes) and randomly choose one representative intermediate node
        node2_connect = df_edge_filt.groupby(["node2"])["node1"].agg(list).to_frame()
        node2_connect["node1"] = node2_connect["node1"].apply(lambda x: ",".join(e for e in set(x)))
        node2_connect['Duplicate'] = node2_connect.duplicated(keep=False).map({True:'Yes', False:'No'}) 
        node2_connect_filt = node2_connect.drop_duplicates()
        
        intermediate_nodes = list(node2_connect_filt.index)
        df_edge_filt = df_edge_filt[df_edge_filt['node2'].isin(intermediate_nodes)]
    return df_edge_filt, intermediate_nodes

def generate_network(g,network,candidate_path,annotation_path,subname,connection_number,dedup):
    """
    Generate NetworkX file

    Parameters:
    args: arg_parser arguments (See arg_parser for details)
    
    Returns networkX Object, candidate_node_headers, annotation_name_headers
    """
    node_list = set()
    if os.path.isdir(candidate_path):
        candidate_names = dir_annotation(g,candidate_path)
        if subname is not None:
            candidate_names = subname
        candidates = get_named_nodes(g,candidate_names)

    #Read in candidate file
    elif os.path.isfile(candidate_path):
        candidate_names = []
        candidate_names.append(file_annotation(g,candidate_path))

        candidates = get_named_nodes(g,candidate_names)

    node_list = node_list.union(candidates)

    #Assign annotation values.
    if annotation_path is not None:
        annotation_names = dir_annotation(g, annotation_path)
        annotation_nodes = get_named_nodes(g, annotation_names)

        node_list = node_list.union(annotation_nodes)

    df_edge_all = pd.read_csv(network, sep="\t")
    network_edges = df_edge_all[(df_edge_all['node1'].isin(node_list)) & (df_edge_all['node2'].isin(node_list))]
    
    if len(network_edges)==0:
        warnings.formatwarning = lambda msg, *args, **kwargs: f'{msg}\n'
        warnings.warn("Warning: Input genes can not connect with each other, plase consider adding intermediate nodes by option -m", stacklevel=2)
    
    if connection_number is not None:
        df_edge_filt, intermediate_nodes = allow_intermediate(df_edge_all,node_list,connection_number,dedup)
        if len(intermediate_nodes)==0:
            warnings.warn("No intermediate nodes found")
            if len(network_edges)==0:
                sys.exit("No network connections found")
        if len(network_edges)==0:
                network_edges = df_edge_filt
        else:
            network_edges = pd.concat(network_edges, df_edge_filt)
    
        node_list = node_list.union(intermediate_nodes)
    #generate graph
    attr = list(network_edges.columns[2:])
    subgraph = nx.from_pandas_edgelist(network_edges, "node1", "node2", attr)
    #Remove Isolated Nodes from file
    subgraph.remove_nodes_from(list(nx.isolates(subgraph)))

    return subgraph,candidate_names,annotation_names
        

if __name__ == "__main__":
    args = arg_parser()

    # Load in Edge List.
    #node1  node2   consensus_interactionType   consensus_score STRING_neighborhood STRING_fusion   STRING_cooccurence  STRING_homology STRING_coexpression STRING_experiments  STRING_database STRING_textmining   STRING_combined_score   GIANT_score DataBaseCount   
    gmax = nx.read_edgelist(args.network,create_using=nx.Graph(), nodetype = str,data=(("consensus_interactionType", str),("consensus_score", str),("STRING_neighborhood", str),("STRING_fusion", str),("STRING_cooccurence", str),("STRING_homology", str),("STRING_coexpression", str),("STRING_experiments", str),("STRING_database", str),("STRING_textmining", str),("STRING_combined_score", str),("GIANT_score", str),("DataBaseCount",str),))
    
    subgraph,candidate_names,annotation_names = generate_network(gmax,args.network,args.candidates,args.annotations,args.subname,args.connection_number,args.dedup)

    subgraph = prune_network(subgraph, args.degree, args.dbcount,args.degree_greater)

    # with open('mirna_target_pair_network.txt','r') as file:
    #     reader = csv.DictReader(file, delimiter='\t')

    #     for row in reader:
    #         if row['Ortholog'] in list(subgraph):
    #             subgraph.add_edge(row['Ortholog'],row['miRNA'],mirna_target = 1) 

    if args.intermediate_nodes:
        subgraph = clean_intermediate_nodes(subgraph,candidate_names,annotation_names)

    #Convert NetworkX format to JSON
    cy = nx.readwrite.json_graph.cytoscape_data(subgraph)
    #Set positional data. Required for Cytoscape.
    pos = nx.spring_layout(subgraph)
    for n,p in zip(cy['elements']['nodes'],pos.values()):
         n['position'] = {'x':2000 * p[0],'y':2000 * p[1]}

    #Write to output file.
    with open(args.output, "w") as outfile:
         json.dump(cy, outfile)
         
         
