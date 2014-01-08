'''
Created on Oct 12, 2013

@author: ksahlin
'''

import os

import argparse
from collections import defaultdict

from BESST import plots, Parameter

class ScoreType(object):
    def __init__(self, name):
        self.name = name
        self.true_positives = []
        self.false_positives = []

    def add_true_positive(self, score):
        self.true_positives.append(score)
    def add_false_positive(self, score):
        self.false_positives.append(score)

def read_in_files(args):
    true_graph = defaultdict(dict)
    for line in open(args.truth, 'r'):
        contig1, contig2, true_gap = line.split()
        true_graph[contig1][contig2] = int(true_gap)
        #true_graph[contig2] = (contig1, int(true_gap))

    total_edges_in_file = 0
    besst_graph = defaultdict(dict)
    for line in open(args.besst, 'r'):
        #TODO: Create ScoreType objects here, one for each additional column after gap column.
        contig1, contig2, besst_gap, nr_links, besst_score, weighted = line.split()
        #if contig1
        besst_graph[contig1][contig2] = (int(besst_gap), int(nr_links), float(besst_score), float(weighted))
        besst_graph[contig2][contig1] = (int(besst_gap), int(nr_links), float(besst_score), float(weighted))
        total_edges_in_file +=1


    return(true_graph, besst_graph, total_edges_in_file)


def true_positives(true_graph, besst_graph):
    true_positives_score = []
    true_positives_nr_links = []
    true_positives_weighted = []

    ##
    # true_negatives is not so important for now at least
    #true_negatives_score = []
    #true_negatives_nr_links = []


    for contig1 in true_graph:
        if contig1 in besst_graph:
            for contig2 in true_graph[contig1]:

                if contig2 in besst_graph[contig1]:
                    true_positives_score.append(besst_graph[contig1][contig2][2])
                    true_positives_nr_links.append(besst_graph[contig1][contig2][1])
                    true_positives_weighted.append(besst_graph[contig1][contig2][3])
                else:
                    ##
                    # true negatives here: in true graph but not in any besst linking
                    # Can bee a covarage issue, it is not necessraily bad to end up here

                    #true_negatives_score.append(besst_graph[contig1][contig2][2])
                    #true_negatives_score.append(besst_graph[contig1][contig2][1])
                    pass

        else:
            ##
            # true negatives here: in true graph but not in any besst linking
            # Can bee a covarage issue, it is not necessraily bad to end up here

            #for contig2 in true_graph[contig1]:
            #    true_negatives_score.append(besst_graph[contig1][contig2][2])
            #    true_negatives_nr_links.append(besst_graph[contig1][contig2][1])
            pass

    return(true_positives_score, true_positives_nr_links, true_positives_weighted)

def false_positives(true_graph, besst_graph):
    false_positives_score = []
    false_positives_nr_links = []
    false_positives_weighted = []
    #besst_edges_set = set()
    true_edges_set = set()
    for contig1 in true_graph:
        for contig2 in true_graph[contig1]:
            true_edges_set.add(tuple(sorted((contig1, contig2))))

    already_visited = set()
    for contig1 in besst_graph:
        for contig2 in besst_graph[contig1]:
            if contig2 not in already_visited and tuple(sorted((contig1, contig2))) not in true_edges_set:
                false_positives_score.append(besst_graph[contig1][contig2][2])
                false_positives_nr_links.append(besst_graph[contig1][contig2][1])
                # if besst_graph[contig1][contig2][1] > 1:
                #     print besst_graph[contig1][contig2][1]
                #     print tuple(sorted((contig1, contig2))) , besst_graph[contig1]
                false_positives_weighted.append(besst_graph[contig1][contig2][3])
        already_visited.add(contig1)




    return(false_positives_score, false_positives_nr_links, false_positives_weighted)

def plot(out, true_positives_score, true_positives_nr_links, true_positives_weighted, false_positives_score, false_positives_nr_links, false_positives_weighted):
    param = Parameter.parameter()
    param.output_directory = out
    try:
        os.mkdir(os.path.join(out, 'plots'))
    except OSError:
        pass
    plots.histogram(true_positives_score, param, 20, 'score', 'freq', "true_positives_score")
    plots.histogram(true_positives_nr_links, param, 20, 'nr_links', 'freq', "true_positives_nr_links")
    #plots.histogram(false_positives_score, param, 10, 'score', 'freq', "false_positives_score")
    #plots.histogram(false_positives_nr_links, param, 10, 'score', 'freq', "false_positives_nr_links")

    ##
    # do a dotplot here
    plots.dot_plot(true_positives_score, true_positives_nr_links, param, 'score', 'nr_links', "true_positives_dist")
    plots.dot_plot(false_positives_score, false_positives_nr_links, param, 'score', 'nr_links', "false_positives_dist")

    plots.dot_plot(true_positives_weighted, true_positives_nr_links, param, 'weighted', 'nr_links', "true_positives_weight")
    plots.dot_plot(false_positives_weighted, false_positives_nr_links, param, 'weighted', 'nr_links', "false_positives_weight")


if __name__ == '__main__':
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Evaluates different scoring functions by looking at scores \
    for true and spurious edges.")
    parser.add_argument('truth', type=str, help='Path to the file with true linkings. ')
    parser.add_argument('besst', type=str, help='Path to the file with BESST linkings and scores.')
    parser.add_argument('out', type=str, help='Path to output location for files. ')
    parser.add_argument("-n", action="store_true",
                  help="File with edges that are not chosen.")


    args = parser.parse_args()

    true_graph, besst_graph, total_edges_in_file = read_in_files(args)
    true_positives_score, true_positives_nr_links, true_positives_weighted = true_positives(true_graph, besst_graph)
    false_positives_score, false_positives_nr_links, false_positives_weighted = false_positives(true_graph, besst_graph)
    if not args.n:
        print 'TP:', len(true_positives_score), 'FP:', total_edges_in_file - len(true_positives_score)
    else:
        print 'TN:', total_edges_in_file - len(true_positives_score), 'FN', len(true_positives_score) 
    #plot(args.out, true_positives_score, true_positives_nr_links, true_positives_weighted, false_positives_score, false_positives_nr_links, false_positives_weighted)
