'''
Created on Aug 12, 2013

@author: ksahlin
'''

import os
import networkx as nx
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError, RuntimeError:
    pass


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

def draw_graph(graph, Contigs, param, labels=None, graph_layout='shell',
               node_size=1600, node_color='blue', node_alpha=0.3,
               node_text_size=12,
               edge_color='blue', edge_alpha=0.3, edge_tickness=1,
               edge_text_pos=0.3,
               text_font='sans-serif'):

    G=nx.Graph()

    # make a list of node_names for each contig
    node_names = [Contigs[ctg_name].scaffold  for ctg_name in Contigs.keys()]

    # nodes = []
    # added_nodes = set()
    # for node in G.nodes():
    #     if node[0] not in added_nodes:
    #         nodes.append(node[0])
    #         added_nodes.add(node[0])

    # make a list of node sizes based on length of contigs
    node_sizes = [Contigs[ctg_name].length  for ctg_name in Contigs.keys()]
    # make a list of node colors based on coverage of contigs
    # numeric values will convert according to 
    # http://matplotlib.org/api/colors_api.html#matplotlib.colors.Colormap
    # with cmap="jet"
    node_colors = [Contigs[ctg_name].coverage  for ctg_name in Contigs.keys()]

    # create a one node per contig graph based on the contig graph
    for node in node_names:
        for nbr in graph.neighbors( (node,'L')):
            weight = graph[(node,'L')][(nbr, 'R')]["nr_links"]
            if node != nbr[0]:
                G.add_edge(, , w=weight)
        for nbr in graph.neighbors( (node,'R')):
            weight = graph[(node,'L')][(nbr, 'R')]["nr_links"]
            if node != nbr[0]:
                G.add_edge(, , w=weight)

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
    nx.draw_networkx_nodes(G,graph_pos,node_size=node_size, 
                           alpha=node_alpha, node_color=node_color)
    nx.draw_networkx_edges(G,graph_pos,width=edge_tickness,
                           alpha=edge_alpha,edge_color=edge_color)
    nx.draw_networkx_labels(G, graph_pos,font_size=node_text_size,
                            font_family=text_font)

    if labels is None:
        labels = range(len(graph))

    edge_labels = dict(zip(graph, labels))
    nx.draw_networkx_edge_labels(G, graph_pos, edge_labels=edge_labels, 
                                 label_pos=edge_text_pos)

    # show graph
    #plt.show()
    # nx.draw(graph)
    # matplotlib.pyplot.savefig(os.path.join(param.output_directory, + 'graph.png'))
    # matplotlib.pyplot.clf()



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
