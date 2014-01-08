'''
Created on Oct 13, 2013

@author: ksahlin
'''
import sys
import os
from scipy.stats import ks_2samp

from mathstats.normaldist.truncatedskewed import param_est
from BESST import plots

def GiveScoreOnEdges(G, Scaffolds, small_scaffolds, Contigs, param, Information, plot):

    span_score_obs = []
    std_dev_score_obs = []
    gap_obs = []
    nr_link_obs = []
    cnt_sign = 0

    for edge in G.edges():
        mean_ = 0
        std_dev = 0
        if G[edge[0]][edge[1]]['nr_links'] != None:
            n = G[edge[0]][edge[1]]['nr_links']
            obs_squ = G[edge[0]][edge[1]]['obs_sq']
            mean_ = G[edge[0]][edge[1]]['obs'] / float(n)
            data_observation = (n * param.mean_ins_size - G[edge[0]][edge[1]]['obs']) / float(n)
            try:
                len1 = Scaffolds[ edge[0][0] ].s_length
            except KeyError:
                len1 = small_scaffolds[ edge[0][0] ].s_length
            try:
                len2 = Scaffolds[ edge[1][0] ].s_length
            except KeyError:
                len2 = small_scaffolds[ edge[1][0] ].s_length
            if 2 * param.std_dev_ins_size < len1 and 2 * param.std_dev_ins_size < len2:
                gap = param_est.GapEstimator(param.mean_ins_size, param.std_dev_ins_size, param.read_len, mean_, len1, len2)
            else:
                gap = data_observation

            G[edge[0]][edge[1]]['gap'] = int(gap)
            if -gap > len1 or -gap > len2:
                G[edge[0]][edge[1]]['score'] = 0
                #print gap, len1, len2
                continue

            #std_dev_d_eq_0 = param_est.tr_sk_std_dev(param.mean_ins_size, param.std_dev_ins_size, param.read_len, len1, len2, gap)

            if 2 * param.std_dev_ins_size < len1 and 2 * param.std_dev_ins_size < len2:
                std_dev_d_eq_0 = param_est.tr_sk_std_dev(param.mean_ins_size, param.std_dev_ins_size, param.read_len, len1, len2, gap)
            else:
                std_dev_d_eq_0 = 2 ** 32

            try:
                std_dev = ((obs_squ - n * mean_ ** 2) / (n - 1)) ** 0.5
                #chi_sq = (n - 1) * (std_dev ** 2 / std_dev_d_eq_0 ** 2)
            except ZeroDivisionError:
                std_dev = 2 ** 32
                #chi_sq = 0


            try:
                l1 = G[edge[0]][edge[1]][Scaffolds[edge[0][0]].name]
            except KeyError:
                l1 = G[edge[0]][edge[1]][small_scaffolds[edge[0][0]].name]
            try:
                l2 = G[edge[0]][edge[1]][Scaffolds[edge[1][0]].name]
            except KeyError:
                l2 = G[edge[0]][edge[1]][small_scaffolds[edge[1][0]].name]

            l1.sort()
            n_obs = len(l1)
            l1_mean = sum(l1) / float(n_obs)
            #l1_median = l1[len(l1) / 2]
            l1 = map(lambda x: x - l1_mean, l1)
            #l1 = map(lambda x: x - l1_median, l1)
            max_obs2 = max(l2)
            l2.sort(reverse=True)
            l2 = map(lambda x: abs(x - max_obs2), l2)
            l2_mean = sum(l2) / float(n_obs)
            #l2_median = l2[len(l2) / 2]
            l2 = map(lambda x: x - l2_mean, l2)
            #l2 = map(lambda x: x - l2_median, l2)
            KS_statistic, p_value = ks_2samp(l1, l2)


            #M_W_statistic, p_val = mannwhitneyu(l1, l2)

            #diff = map(lambda x: abs(abs(x[1]) - abs(x[0])), zip(l1, l2))
            #sc = sum(diff) / len(diff)

            if len(l1) < 3:
                span_score = 0
            else:
                span_score = 1 - KS_statistic


#            try:
#                span_score = 1 - sc / float(min((max_obs1 - min_obs1), (max_obs2 - min_obs2)))
#            except ZeroDivisionError:
#                span_score = 0
#            if span_score < 0:
#                span_score = 0
#                print  'ZEEERO', max_obs1 - min_obs1, max_obs2 - min_obs2, gap, sc, Scaffolds[edge[0][0]].contigs[0].name, Scaffolds[edge[1][0]].contigs[0].name

# if len(l1) > 3:
#     print >> Information , 'avg_diff: ', sc, 'span1: ', (max_obs1 - min_obs1), 'span2: ', (max_obs2 - min_obs2), 'Span score: ', span_score, 'pval: ', p_value, 'Est gap: ', gap, 'Nr_links: ', len(l1) #, Scaffolds[edge[0][0]].contigs[0].name, Scaffolds[edge[1][0]].contigs[0].name, len(diff)


                #print >> Information , l1
                #print >> Information , l2
                #print >> Information , diff
            #print  span_score

#
#            k = normal.MaxObsDistr(n, 0.95)
#            if 2 * param.read_len < len1 and 2 * param.read_len < len2:
#                #span_max1 = min(param.mean_ins_size + k * param.std_dev_ins_size - 2 * param.read_len, len1 - param.read_len + max(0, gap))
#                #span_max2 = min(param.mean_ins_size + k * param.std_dev_ins_size - 2 * param.read_len, len2 - param.read_len + max(0, gap))
#
#                span_max1 = min(param.mean_ins_size + k * param.std_dev_ins_size - param.read_len, len1 + max(0, gap))
#                span_max2 = min(param.mean_ins_size + k * param.std_dev_ins_size - param.read_len, len2 + max(0, gap))
#
#                try:
#                    span_obs1 = Scaffolds[ edge[0][0] ].upper_right_nbrs_obs[edge[1]] - Scaffolds[ edge[0][0] ].lower_right_nbrs_obs[edge[1]] if edge[0][1] == 'R' else Scaffolds[ edge[0][0] ].upper_left_nbrs_obs[edge[1]] - Scaffolds[ edge[0][0] ].lower_left_nbrs_obs[edge[1]]
#                except KeyError:
#                    span_obs1 = small_scaffolds[ edge[0][0] ].upper_right_nbrs_obs[edge[1]] - small_scaffolds[ edge[0][0] ].lower_right_nbrs_obs[edge[1]] if edge[0][1] == 'R' else small_scaffolds[ edge[0][0] ].upper_left_nbrs_obs[edge[1]] - small_scaffolds[ edge[0][0] ].lower_left_nbrs_obs[edge[1]]
#                try:
#                    span_obs2 = Scaffolds[ edge[1][0] ].upper_right_nbrs_obs[edge[0]] - Scaffolds[ edge[1][0] ].lower_right_nbrs_obs[edge[0]] if edge[1][1] == 'R' else Scaffolds[ edge[1][0] ].upper_left_nbrs_obs[edge[0]] - Scaffolds[ edge[1][0] ].lower_left_nbrs_obs[edge[0]]
#                except KeyError:
#                    span_obs2 = small_scaffolds[ edge[1][0] ].upper_right_nbrs_obs[edge[0]] - small_scaffolds[ edge[1][0] ].lower_right_nbrs_obs[edge[0]] if edge[1][1] == 'R' else small_scaffolds[ edge[1][0] ].upper_left_nbrs_obs[edge[0]] - small_scaffolds[ edge[1][0] ].lower_left_nbrs_obs[edge[0]]
#
#
#                #span_score1 = min((max(0, gap) + 2 * param.read_len + span_obs1) / float(span_max1) , float(span_max1) / (max(0, gap) + 2 * param.read_len + span_obs1)) if span_obs1 > 0 else 0
#                #span_score2 = min((max(0, gap) + 2 * param.read_len + span_obs2) / float(span_max2) , float(span_max2) / (max(0, gap) + 2 * param.read_len + span_obs2)) if span_obs2 > 0 else 0
#
#                span_score1 = min((max(0, gap) + param.read_len + span_obs1) / float(span_max1) , float(span_max1) / (max(0, gap) + param.read_len + span_obs1)) if span_obs1 > 0 else 0
#                span_score2 = min((max(0, gap) + param.read_len + span_obs2) / float(span_max2) , float(span_max2) / (max(0, gap) + param.read_len + span_obs2)) if span_obs2 > 0 else 0
#
#                span_score = min(span_score1, span_score2)
#
#                #span_score = (max(0, gap) + param.read_len + span_obs1) / float(span_max1)
#            else:
#                span_score = 0


            try:
                std_dev_score = min(std_dev / std_dev_d_eq_0, std_dev_d_eq_0 / std_dev) #+ span_score #+ min(n/E_links, E_links/float(n))
            except ZeroDivisionError:
                std_dev_score = 0
                # sys.stderr.write(str(std_dev) + ' ' + str(std_dev_d_eq_0) + ' ' + str(span_score) + ' ' + str(n) + '\n')

            if std_dev_score > param.scorethreshold and span_score > param.scorethreshold:
            # if span_score > param.scorethreshold:
                G[edge[0]][edge[1]]['score'] = std_dev_score + span_score
                # G[edge[0]][edge[1]]['score'] = span_score
            else:
                G[edge[0]][edge[1]]['score'] = 0
                print std_dev_score, span_score, n
            if param.plots:
                span_score_obs.append(span_score)
                std_dev_score_obs.append(std_dev_score)
                gap_obs.append(gap)
                nr_link_obs.append(n_obs)


    if param.plots:
        plots.histogram(span_score_obs, param, bins=20, x_label='score', y_label='frequency', title='Dispersity_score_distribuion' + plot + '.' + param.bamfile.split('/')[-1])
        plots.histogram(std_dev_score_obs, param, bins=20, x_label='score', y_label='frequency', title='Standard_deviation_score_distribuion' + plot + '.' + param.bamfile.split('/')[-1])
        plots.dot_plot(std_dev_score_obs, span_score_obs, param, x_label='std_dev_score_obs', y_label='span_score_obs', title='Score_correlation' + plot + '.' + param.bamfile.split('/')[-1])
        plots.dot_plot(std_dev_score_obs, gap_obs, param, x_label='std_dev_score_obs', y_label='estimated gap size', title='Gap_to_sigma' + plot + '.' + param.bamfile.split('/')[-1])
        plots.dot_plot(span_score_obs, gap_obs, param, x_label='span_score_obs', y_label='estimated gap size', title='Gap_to_span' + plot + '.' + param.bamfile.split('/')[-1])
        plots.dot_plot(span_score_obs, nr_link_obs, param, x_label='span_score_obs', y_label='Number links', title='Obs_to_span' + plot + '.' + param.bamfile.split('/')[-1])

    for edge in G.edges():
        if G[edge[0]][edge[1]]['nr_links'] != None:
            try:
                G[edge[0]][edge[1]]['score']
            except KeyError:
                pass
                # sys.stderr.write(str(G[edge[0]][edge[1]]) + ' ' + str(Scaffolds[edge[0][0]].s_length) + ' ' + str(Scaffolds[edge[1][0]].s_length))
    print >> Information, 'Number of significantly spurious edges:', cnt_sign

    # if param.development:
    #     tempout = open(os.path.join(param.output_directory, 'besst_edges.txt'), 'w')
    #     for edge in G.edges():
    #         if G[edge[0]][edge[1]]['nr_links'] != None:
    #             try:
    #                 c1 = Scaffolds[edge[0][0]].contigs[0].name
    #             except KeyError:
    #                 c1 = small_scaffolds[edge[0][0]].contigs[0].name
    #             try:
    #                 c2 = Scaffolds[edge[1][0]].contigs[0].name
    #             except KeyError:
    #                 c2 = small_scaffolds[edge[1][0]].contigs[0].name

    #             print >> tempout, c1 + '\t' + c2 + '\t' + str(G[edge[0]][edge[1]]['gap']) + '\t' + str(G[edge[0]][edge[1]]['nr_links']) + '\t' + str(G[edge[0]][edge[1]]['score']) + '\t' + str(G[edge[0]][edge[1]]['nr_links'] * G[edge[0]][edge[1]]['score'])

    return()
