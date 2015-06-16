'''
Created on Aug 12, 2013

@author: ksahlin
'''

import os
import networkx as nx
import math
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError, RuntimeError:
    pass
from pylab import *
from networkx import *
 
def histogram(x_, param, bins=20, x_label='x', y_label='y', title='Histogram',nr_obs = 10000):
    x = x_[:nr_obs] # If many contigs/edges we only plot 10 000 (default) for time and visuability purposes
    dest = os.path.join(param.output_directory, 'plots', title + '.png')
    plt.hist(x, bins)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)

    plt.savefig(dest)
    plt.clf()
    return()

def dot_plot(x_, y_, param, x_label='x', y_label='y', title='Dotplot'):
    x = x_[:10000] # If many contigs/edges we only plot 10 000 for time and visuability purposes
    y = y_[:10000] # If many contigs/edges we only plot 10 000 for time and visuability purposes
    dest = os.path.join(param.output_directory, 'plots', title + '.png')
    plt.scatter(x, y, marker='o')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.grid(True)

    plt.savefig(dest)
    plt.clf()
    return


def multiple_histogram(list_of_datasets_, param, x_label='x', y_label='y', title='Stacked_histogram'):
    list_of_datasets = [list[:10000] for list in list_of_datasets_]
    dest = os.path.join(param.output_directory, 'plots', title + '.png')
    # filter out if any of the lists contains 0 elemnets
    list_of_datasets = filter(lambda x: len(x) > 0, list_of_datasets)
    for dataset in list_of_datasets:
        plt.hist(dataset, alpha=0.5)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)

    plt.savefig(dest)
    plt.clf()

    return()

def draw_graph(graph, Contigs, small_contigs, param, graph_layout='spectral',
               node_size=1600, node_color='blue', node_alpha=0.6,
               node_text_size=12,
               edge_color='blue', edge_alpha=0.3, edge_tickness=1,
               edge_text_pos=0.3,
               text_font='sans-serif'):

    G=nx.Graph()
    # make a list of node_names for each contig
    node_names = [Contigs[ctg_name].scaffold for ctg_name in Contigs.keys()] + [small_contigs[ctg_name].scaffold for ctg_name in small_contigs.keys()]
    # nodes = []
    # added_nodes = set()
    # for node in G.nodes():
    #     if node[0] not in added_nodes:
    #         nodes.append(node[0])
    #         added_nodes.add(node[0])

    # make a list of node sizes based on length of contigs
    node_sizes = [Contigs[ctg_name].length  for ctg_name in Contigs.keys()] + [small_contigs[ctg_name].length for ctg_name in small_contigs.keys()]
    smallest_node = min(node_sizes)

    print [int(size) for size in node_sizes]
    # make a list of node colors based on coverage of contigs
    # numeric values will convert according to 
    # http://matplotlib.org/api/colors_api.html#matplotlib.colors.Colormap
    # with cmap="jet"
    node_colors = [Contigs[ctg_name].coverage  for ctg_name in Contigs.keys()] + [small_contigs[ctg_name].coverage for ctg_name in small_contigs.keys()]

    node_dict = dict(zip(node_names, zip(node_colors, node_sizes)))
    # create a one node per contig graph based on the contig graph
    for node in node_names: 
        if len(graph.neighbors((node,'L')) + graph.neighbors((node,'R'))) == 2:
            G.add_node(node)
            continue
        for nbr in graph.neighbors((node,'L')):
            if node != nbr[0]:
                weight = graph[(node,'L')][nbr]["nr_links"]
                if weight >= 3:
                    G.add_edge(node, nbr[0] , w=weight)

        for nbr in graph.neighbors((node,'R')):
            if node != nbr[0]:
                weight = graph[(node,'R')][nbr]["nr_links"]
                if weight >= 3:
                    G.add_edge(node, nbr[0] , w=weight)

    # sort nodes in node_dict to get same order as in G.nodes()
    for node in G.nodes():
        if node_dict[node][1] < 200:
             G.remove_node(node)
    node_names = []
    node_colors = []
    node_sizes = []
    print len(G.nodes())
    for node in G.nodes():        
        node_names.append(node)
        node_colors.append(node_dict[node][0])
        #natural logarithmic scale on sizes
        node_sizes.append(math.log(node_dict[node][1]/float(smallest_node) +1) * 10)

           

    print [int(color) for color in node_colors]
    # get edge weights
    edge_labels = {} #dict(zip(graph, labels))
    for edge in G.edges():
        edge_labels[edge] = G[edge[0]][edge[1]]['w'] # dict(zip(graph, labels))



    # these are different layouts for the network you may try
    # shell seems to work best
    if graph_layout == 'spring':
        graph_pos=nx.spring_layout(G)
    elif graph_layout == 'spectral':
        graph_pos=nx.spectral_layout(G)
    elif graph_layout == 'random':
        graph_pos=nx.random_layout(G)
    else:
        graph_pos=nx.shell_layout(G)

    # draw graph
    nx.draw_graphviz(G,
         node_size = node_sizes,
         node_color = node_colors,
         with_labels = False)
    colorbar()
    # nx.draw_networkx_nodes(G,graph_pos,node_size=node_sizes, 
    #                        alpha=node_alpha, node_color=node_colors)
    # nx.draw_networkx_edges(G,graph_pos,width=edge_tickness,
    #                        alpha=edge_alpha,edge_color=edge_color)
    # nx.draw_networkx_labels(G, graph_pos,font_size=node_text_size,
    #                         font_family=text_font)


    # nx.draw_networkx_edge_labels(G, graph_pos, edge_labels=edge_labels, 
    #                              label_pos=edge_text_pos)

    # show graph
    matplotlib.pyplot.savefig(os.path.join(param.output_directory, 'graph.eps'))
    matplotlib.pyplot.clf()



# def VizualizeGraph(G, param, Information):
#    import os
#    try:
#        import matplotlib
#        matplotlib.use('Agg')

#        try:
#            os.mkdir(param.output_directory + '/graph_regions' + str(int(param.mean_ins_size)))
#        except OSError:
#            #directory is already created
#            pass
#        counter = 1
#        import copy

#        G_copy = copy.deepcopy(G)
#        RemoveIsolatedContigs(G_copy, Information)
#        CB = nx.connected_component_subgraphs(G_copy)
#        for cycle in CB:
#            nx.draw(cycle)
#            matplotlib.pyplot.savefig(param.output_directory + 'graph_regions' + str(int(param.mean_ins_size)) + '/' + str(counter) + '.png')
#            matplotlib.pyplot.clf()
#            counter += 1
#    except ImportError:
#        pass
#    return()
