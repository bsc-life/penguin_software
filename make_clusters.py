#!/usr/bin/env python
# coding: utf-8



import gc
import os

from ete3 import Tree, TreeStyle, CircleFace, faces, TextFace

import pandas as pd
import pyranges as pr
from IPython.display import display
import IPython

import numpy as np

from scipy.stats import fisher_exact
from sklearn.metrics import v_measure_score, homogeneity_score, completeness_score


from gmpy2 import popcount

from matplotlib import pyplot as plt
from matplotlib.colors import to_hex
from pathlib import Path
import pickle


import scipy.cluster.hierarchy as sch
from sklearn.metrics import calinski_harabasz_score, silhouette_score, davies_bouldin_score


from itertools import product
from collections import defaultdict

t = Tree( "((a,b),c);" )
circular_style = TreeStyle()
circular_style.mode = "c" # draw tree in circular mode
circular_style.scale = 20
#display(t.render('%%inline'))





data_path = './src/EPIN_reconstruction_data'
cell_line = "LNCaP"
cell_type_data_path = os.path.join(data_path, cell_line)
data_type = 'all-nodes' # *all-nodes* (or linker-nodes or bound-nodes)
metric    = 'overlap' # overlap or jaccard
#metric    = 'jaccard' # overlap or jaccard


#tf_path = os.path.join(data_path, 'LNCaP/ChIP-seq/')

type_filtered="filtered"
type_net = "Mega"

type_enha = "normalised"
loops = "P2PBckgr_1_AnnotatedData_LNCaP_pad5000_EP.txt" if type_enha == "unnormalised" else "E-P_loops_filtered-regions.tsv"
######

######
ctcf_bed_file = 'CTCF_lncap_ENCFF155SPQ.bed'

ctcf_file = os.path.basename(ctcf_bed_file).replace(".bed", "")

#######
paintor = "unfiltered_paintor"



#####

net_path = os.path.join(f'outputs_{cell_line}_EPINS', "tables" ,type_filtered , 'EP_graph_edges')
print("net_path: ",net_path)

result_path = os.path.join(f'outputs_{cell_line}_EPINS', "tables" ,type_filtered , "cluster_results")
if not os.path.exists(result_path):
    os.makedirs(result_path)
    

print("result_path: ",result_path)

## clustering functions
def iter_leaves(sons, c):
    for cc in sons[c]:
        if cc in sons:
            yield from iter_leaves(sons, cc)
        else:
            yield cc

def iter_childs(sons, c):
    yield c
    for cc in sons[c]:
        if cc in sons:
            yield from iter_childs(sons, cc)
        else:
            yield cc

def get_labels(Y, cut, sons=None):
    """
    Given a Scipy linkage table, and a distance cutoff, gives the cluster
    labelling of each input row.

    :params Y: Scipy linkage table
    :params cut: distance cutoff
    :params None sons: dctionary containing the relation between
       nodes: {mother: [daugther1, daugther2], ...}. By default it is calculated
       from Y, but it is slower.

    :returns: a list of labels of each label, and also each internal node of the clustering
    """
    size = len(Y) + 1
    label = 1
    labels = [0] * (size * 2 - 1)
    if sons is None:
        sons = dict((float(c), (i, j)) for c, (i, j, d, l) in enumerate(Y, size))
    dist = dict((float(c), d) for c, (i, j, d, l) in enumerate(Y, size))
    for c in sons:
        if dist[c] > cut:
            for cc in sons[c]:
                if not cc in sons:
                    labels[int(cc)] = label
                    label += 1
                elif dist[cc] < cut:
                    for ccc in iter_childs(sons, cc):
                        labels[int(ccc)] = label
                    label += 1
    return labels

def get_labels_fast(Y, cut, sons=None, dist=None, max_clusters=100):
    """
    Given a Scipy linkage table, and a distance cutoff, gives the cluster
    labelling of each input row.

    :params Y: Scipy linkage table
    :params cut: distance cutoff
    :params None sons: dctionary containing the relation between
       nodes: {mother: [daugther1, daugther2], ...}. By default it is calculated
       from Y, but it is slower.

    :returns: a list of labels of each label, and also each internal node of the clustering
    """
    size = len(Y) + 1
    label = 1
    labels = [0] * (size * 2 - 1)
    if sons is None:
        sons = dict((float(c), (i, j)) for c, (i, j, d, l) in enumerate(Y, size))
    if dist is None:
        dist = dict((float(c), d) for c, (i, j, d, l) in enumerate(Y, size))
    for cc in (cc for c in sons if dist[c] > cut for cc in sons[c]):
        if not cc in sons:
            labels[int(cc)] = label
            label += 1
            if label > max_clusters:
                break
        elif dist[cc] < cut:
            for ccc in iter_childs(sons, cc):
                labels[int(ccc)] = label
            label += 1
            if label > max_clusters:
                break
    return labels


def best_cut(Y, sons=None, method='CH', max_clusters=100):
    if method == 'CH':
        method = calinski_harabasz_score
    elif method == 'silhouette':
        method = silhouette_score
    elif method == 'DB':
        method = davies_bouldin_score
    else:
        raise NotImplementedError(method)
    size = len(Y) + 1
    if sons is None:
        sons = dict((float(c), (i, j)) for c, (i, j, d, l) in enumerate(Y, size))
    results = []
    dist = dict((float(c), d) for c, (i, j, d, l) in enumerate(Y, size))
    for n, c in enumerate([i[2] * 1.0000001 for i in Y if i[3] != 2][-max_clusters * 2:-1]):
        labels = get_labels_fast(Y, c, sons, dist, max_clusters=max_clusters)[:size - 1]
        # we want on average at least 4 nodes per cluster
        if len(set(labels)) >= size * 0.5:
            continue
        results.append((method(Y, labels), c))
    cut = max(results)[1] if len(results) > 0 else 0
    return cut

def myresample(m, size):
    x1 = np.random.randint(0, size, size**2 // 8)
    y1 = np.random.randint(0, size, size**2 // 8)
    x2 = np.random.randint(0, size, size**2 // 8)
    y2 = np.random.randint(0, size, size**2 // 8)
    m2 = m.copy()
    m2[x2, y2] = m2[x1, y1]
    return np.maximum(m2, m2.transpose() )

def bootstrap(m, ori_Y, cl_method='ward', metric='euclidean', n=1000, p=0.9):
    size = len(m)
    ori_sons = dict((float(c), (i, j)) for c, (i, j, d, l) in enumerate(ori_Y, size))
    ori_dad = {}
    for c, (i, j) in ori_sons.items():
        ori_dad[i] = c
        ori_dad[j] = c
    bootstrap = defaultdict(int)
    dady = {}
    for a in ori_sons:
        dady[tuple(sorted(iter_leaves(ori_sons, a)))] = a
    for _ in range(n):
        Y = sch.linkage(myresample(m, size), method=cl_method, metric=metric)
        sons = dict((float(c), (i, j)) for c, (i, j, d, l) in enumerate(Y, size))
        for n in sons:
            bootstrap[tuple(sorted(iter_leaves(sons, n)))] += 1
    labels = [0] * (size * 2)
    c = 1
    p = n * p
    for nodes in sorted(bootstrap.keys(), key=lambda x: len(x), reverse=True):
        try:
            ori_node = int(max(iter_childs(ori_sons, dady[nodes])))
        except KeyError:
            continue
        if bootstrap[nodes] >= p:
            try:
                if not labels[int(ori_dad[ori_node])]:
                    continue
            except KeyError:  # root of the tree
                pass
            labels[ori_node] = c
            c += 1
        else:
            if labels[int(ori_dad[ori_node])]:
                cc = labels[int(ori_dad[ori_node])]
                for n in iter_leaves(ori_sons, dady[nodes]):
                    labels[int(n)] = cc
#     for n, c in enumerate(labels):
#         try:
#             d = int(ori_dad[n])
#         except KeyError:
#             labels[n] = 0
#             continue
#         if c and not labels[d]:
#             labels[d] = c
#         if labels[d] != c:
#             labels[d] = 0
    labels[size * 2 - 2] = 0
    nums = set(labels)
    nums = dict((c, n) for n, c in enumerate(nums))
    labels = [nums[l] for l in labels]
    return labels

def dendrogram(Y, axe=None, cut=None, cut_score='CH', matrix=None, method=None, metric=None,
               max_clusters=100):
    if not axe:
        plt.figure(figsize=(12, 3))
        axe = plt.subplot(111)
    dmax = max(i[2] for i in Y)
    size = len(Y) + 1
    sons = dict((float(c), (i, j)) for c, (i, j, d, l) in enumerate(Y, size))
    pos = dict((c, (p, 0)) for p, c in enumerate(iter_leaves(sons, size * 2 - 2)))
    if cut_score == 'bootstrap':
        labels = bootstrap(matrix, Y, cl_method=method, metric=metric)
    else:
        if cut_score is not None:
            cut = best_cut(Y, sons=sons, method=cut_score, max_clusters=max_clusters)
        labels = get_labels(Y, cut, sons)
    color_keys = ['grey'] + [c['color'] for c in plt.rcParams['axes.prop_cycle']]
    for c, (i, j, d, l) in enumerate(Y, size):
        xi, yi = pos.get(i, (0, 0))
        xj, yj = pos.get(j, (0, 0))
        color = 'grey' if labels[c] == 0 else color_keys[1 + labels[c] % 10]
        axe.plot([xi, xi, xj, xj], [yi, d, d, yj], color=color, lw=1 if color == 'grey' else 2)
        pos[c] = (xi + xj) / 2, d
        if i < size:
            axe.text(xi, -dmax * 0.01, int(i), va='top', ha='center', size=8, color=color, rotation=90)
        if j < size:
            axe.text(xj, -dmax * 0.01, int(j), va='top', ha='center', size=8, color=color, rotation=90)
    plt.xlim(-0.5, size - 0.5)
    axe.axis('off')
    return [int(i) for i in iter_leaves(sons, size * 2 - 2)], labels[:size]

def plot_corr_mat(D, Y, labels=None, fig=None, vmin=None, vmax=None, title='distance',
                  normed=True, method='ward', metric='euclidean', color_threshold=0.7,
                  cut_score='CH', max_clusters=100):

    # Compute and plot first dendrogram.
    if not fig:
        fig = plt.figure(figsize=(12, 14))

    ax1 = fig.add_axes([0.2, 0.71, 0.7, 0.1])
    idx1, clusters = dendrogram(Y, cut=color_threshold, axe=ax1, cut_score=cut_score,
                                matrix=D, method=method, metric=metric, max_clusters=max_clusters)

    color_keys = ['grey'] + [c['color'] for c in plt.rcParams['axes.prop_cycle']]

    # plot matrix
    axmatrix = fig.add_axes([0.2, 0.1, 0.7, 0.6])
    D = D[idx1,:]
    D = D[:,idx1]
    im = axmatrix.imshow(D, vmin=vmin, vmax=vmax, cmap='viridis')
    axmatrix.set_xticks([])
    axmatrix.set_yticks(range(len(D)))
    labels = range(len(D)) if labels is None else labels
    # labels in the Y
    axmatrix.set_yticklabels(list('{}      '.format(l) for l in np.array(labels)[idx1]), ha='right', size=6)
    # numbers in the Y
    for i, l in enumerate(idx1):
        plt.text(-0.5, i, '{} '.format(l), ha='right', va='center', size=6,
                 color=color_keys[1 + clusters[l] % 10])
    axmatrix.tick_params(axis=u'both', which=u'both',length=0)
    vmin, vmax = im.get_clim()

    # colorbar
    axcolor = fig.add_axes([0.05, 0.75, 0.13, 0.05])
    plt.colorbar(im, cax=axcolor, orientation='horizontal')
    plt.xlabel(title, rotation=0, size=10)
    height, pos = np.histogram([v for l in D for v in l], bins=20, range=(vmin, vmax))
    height = (vmin + height / height.max()) * (vmax - vmin) * 0.95
    pos    = (vmin + pos    / pos.max()   ) * (vmax - vmin)
    axcolor.step(pos[:-1], height, color='white', linewidth=2)
    return clusters, labels



#for ctcf_bed_file in ['CTCF_lncap_ENCFF155SPQ.bed', 'TFs/CTCF_bert_39737_peaks_hg19.bed']:




# add the features here to test for enrichments
features = {}
features["CTCF"] = {'path': os.path.join(cell_type_data_path, ctcf_bed_file), 'skip': 0}
features['GWAS'] = {'path': os.path.join(cell_type_data_path, 'paintor_1causals.txt'), 'skip': 1}
features['GWAS_Cat_prostate.carcinoma'] = {'path': os.path.join(cell_type_data_path, 'GWAS_Cat_prostate.carcinoma'), 'skip': 1}
#features['H3K27ac'] = {'path': os.path.join(data_path, 'LNCaP/HiChIP/broad_narrow.bed'), 'skip': 0}


#here the gwas catalog SNPs.
#for gwas_file in os.listdir(os.path.join(data_path, 'GWAS_Catalog')):
#    features[gwas_file] = {'path': os.path.join(data_path, 'GWAS_Catalog', gwas_file), 'skip': 1}



feature_promoters = {}
feature_coordinates = {}

for feature in features:
    print("create_feature promotters: ", feature)
    if feature == "GWAS":
        tmp_feature_coordinates = pd.read_csv(
        features[feature]['path'], skiprows= features[feature]['skip'],
                sep='\t', names=['Chromosome', 'Start', 'End', 'Posterior_Prob', 'Pvalue'], usecols=[0, 1, 2, 3, 7])
        if paintor == "filtered_paintor":
            tmp_feature_coordinates = tmp_feature_coordinates[tmp_feature_coordinates.Pvalue <= 5e-08][['Chromosome' , 'Start' , 'End', 'Posterior_Prob', 'Pvalue']]
        else:
            tmp_feature_coordinates = tmp_feature_coordinates[['Chromosome' , 'Start' , 'End', 'Posterior_Prob', 'Pvalue']]
    else:
        tmp_feature_coordinates = pd.read_csv(
            features[feature]['path'], skiprows= features[feature]['skip'],
                    sep='\t', names=['Chromosome', 'Start', 'End'], usecols=[0, 1, 2, ])
    tmp_feature_coordinates['Chromosome'] = tmp_feature_coordinates['Chromosome'].replace(
        dict((i, i[3:]) for i in set(tmp_feature_coordinates['Chromosome'])))
    promoter_coordinates = pd.read_csv(os.path.join(cell_type_data_path, "genepos.txt"),
                sep='\t', names=["Chromosome", "Start", "End", "promoter"], usecols=[0,1,2,3],
            skiprows=0)
    promoter_coordinates = promoter_coordinates.drop_duplicates()
    promoter_coordinates = pr.PyRanges(promoter_coordinates)
    feature_coordinates[feature] = pr.PyRanges(tmp_feature_coordinates)
    promoter_coordinates = promoter_coordinates.nearest(feature_coordinates[feature])
    try:
        maxgap = 1 if feature != "CTCF" else 10000
        feature_promoters[feature] = set(promoter_coordinates[promoter_coordinates.Distance < maxgap].promoter)
    except Exception as e:
        print(e)
# ## load graph to get mark presence in enhancer/promoters# In[101]:



#metaloop_features = {}
#edges_by_metaloop = {}
feature_loop = {}
all_edges = {}

try:
    for enum_feature , feature in enumerate(features):
        print(feature)

        metaloop_features = {}
       
            
        if data_type == 'all-nodes':  # this is the good one
            filter_pairs = lambda pairs: False
            update_pairs = lambda pairs, vals: pairs.append(tuple(sorted(vals)))
        elif data_type == 'linker-nodes':
            filter_pairs = lambda pairs: any(p[:5] in ['ENHA_', 'PROM_'] for p in pairs)
            update_pairs = lambda pairs, vals: pairs.append(tuple(sorted(vals)))
        elif data_type == 'bound-nodes':  # in this case dear reader we only add one single protein
            filter_pairs = lambda pairs: not all(p[:5] in ['ENHA_', 'PROM_'] for p in pairs)
            update_pairs = lambda pairs, vals: pairs.extend(vals)

        metaloop_features[feature] = {}
        
        edges_by_metaloop = {}
        for fnam in os.listdir(net_path):
            if not len(metaloop_features[feature]) % 1000:
                print(' ', len(metaloop_features[feature]))
            metaloop = fnam[:-4]
            chromosomes = []
            starts = []
            ends = []
            pairs = []
            # check if metaloop promoter is near feature
            prom_feature = int(metaloop in feature_promoters[feature])
            for line in open(os.path.join(net_path, fnam)):
                prot1, prot2 = this_prot1, this_prot2 = line.strip().split("\t")
                if prot1.startswith('ENHA_'):
                    _, c, b, e, p = prot1.split('_', 4)
                    chromosomes.append(c)
                    starts.append(int(b))
                    ends.append(int(e))
                    this_prot1 = f'ENHA_{p}'
                if prot2.startswith('ENHA_'):
                    _, c, b, e, p = prot2.split('_', 4)

                    chromosomes.append(c)
                    starts.append(int(b))
                    ends.append(int(e))
                    this_prot2 = f'ENHA_{p}'
                if filter_pairs([prot1, prot2]):  # never go there when all-nodes
                    continue
                
                update_pairs(pairs, [this_prot1, this_prot2])
            tmp = pr.PyRanges(chromosomes=chromosomes, starts=starts, ends=ends)
            tmp = tmp.drop_duplicate_positions()
            tmp = tmp.nearest(feature_coordinates[feature])
            if len(tmp) > 0:
                maxgap = 1 if feature != "CTCF" else 10000
                enha_feature = sum(tmp.Distance < maxgap)
                metaloop_features[feature][metaloop] = {'promoter': prom_feature, 'enhancer': enha_feature}
            
            #save only if it is not available
            
            edges_by_metaloop[metaloop] = set(pairs)

        


        #for feature in features:
        where = 'both' if feature == 'CTCF' else 'enhancer'
        where = 'both' if feature == 'CTCF' else 'promoter'
        
        if where == 'both':
            feature_loop[feature] = dict([(k, metaloop_features[feature][k]['enhancer']) for k in metaloop_features[feature]
                                            if (metaloop_features[feature][k]['promoter'] and metaloop_features[feature][k]['enhancer']) !=0])
        elif where == 'any':
            feature_loop[feature] = dict([(k, metaloop_features[feature][k]['enhancer']) for k in metaloop_features[feature]
                                            if (metaloop_features[feature][k]['promoter'] or metaloop_features[feature][k]['enhancer']) !=0])
        elif where == 'promoter':
            feature_loop[feature] = dict([(k, metaloop_features[feature][k]['enhancer']) for k in metaloop_features[feature]
                                if (metaloop_features[feature][k]['promoter']) !=0])
        elif where == 'enhancer':
            feature_loop[feature] = dict([(k, metaloop_features[feature][k]['enhancer']) for k in metaloop_features[feature]
                                if (metaloop_features[feature][k]['enhancer']) !=0])
        else:
            raise Exception('escribe bien!!')   


        all_edges[feature] = set(k for metaloop in edges_by_metaloop
                                for k in edges_by_metaloop[metaloop])
        #print(f'len all_edges[{feature}]: ', len(all_edges[feature]))
        metaloop_features = None
        n_gc = gc.collect()
except Exception as e:
    print(e)
    IPython.embed()

# this will run only in the last one since for all of them is the same 
data = {}
for metaloop in edges_by_metaloop:
    this_loop = edges_by_metaloop[metaloop]
    if len(edges_by_metaloop[metaloop]) > 0 :
        data[metaloop] = int(''.join(str(int(n in this_loop)) for n in all_edges[feature]), 2)

        #print(f'data[{metaloop}]: ', data[metaloop])
metaloops = list(data.keys())

if metric == 'jaccard':
    distances = dict(((loop1, loop2),
                    popcount(data[loop1] & data[loop2]) / popcount(data[loop1] | data[loop2]))
                    for nloop, loop1 in enumerate(metaloops, 1)
                    for loop2 in metaloops[nloop:])
elif metric == 'overlap':
    distances = dict(((loop1, loop2),   
                    popcount(data[loop1] & data[loop2]) / min(popcount(data[loop1]), popcount(data[loop2])))
                    for nloop, loop1 in enumerate(metaloops, 1)
                    for loop2 in metaloops[nloop:])


distance_matrix = np.asarray([[1.0 - distances.get((n1, n2), distances.get((n2, n1), 1.0))
                                for n2 in metaloops] for n1 in metaloops])

## Apply clustering
D = np.asarray(distance_matrix)
# reverse and min/max normalize
D = (D - D.min()) / (D.max() - D.min())
#IPython.embed()
Y = sch.linkage(D, method='ward', metric='euclidean')

fig = plt.figure(figsize=(48, 56))
clusters, labels = plot_corr_mat(D, Y, labels=metaloops, max_clusters=100,
                                normed=True, cut_score='DB', fig=fig)


## convert linkage dendrogram to newick
def getNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick
        
tree = sch.to_tree(Y)


nw = getNewick(tree, "", tree.dist, labels)

tree_path = os.path.join(result_path, f'{data_type}_{metric}-distance-tree_new.nw')


nw_file = open(tree_path, 'w')
nw_file.writelines(nw)
nw_file.close()



tree_path = nw




def my_layout_clusters(node):
    node.img_style['size'] = 0

    if node.is_root():
        return

    if node.is_leaf():

        # each multi-loop represented as a sphere of a size proportional to
        # the number of edges (or nodes in case of bound-nodes).
        C = CircleFace(radius=20 * (popcount(data[node.name]) - min_node_size) / diff_node_size,
                    color='firebrick' if node.name in feature_loop['CTCF'] else 'royalblue',
                    style="circle")
        # Let's make the sphere transparent
        C.opacity = 0.9
        # And place as a float face over the tree
        faces.add_face_to_node(C, node, 0, position="float")
        return


    w = 10
    if node in clusters:
        node.img_style["hz_line_color"] = clusters[node]
        node.img_style["vt_line_color"] = clusters[node]
        node.img_style["vt_line_width"] = w
        node.img_style["hz_line_width"] = w



def my_layout_enrichments_CTCF(node):
    node.img_style['size'] = 0

    if node.is_root():
        return

    if node.is_leaf():
        if node.name in feature_loop['CTCF']:
            r = 0.85
            b = 0.3
        else:
            b = 0.85
            r = 0.3
        # each multi-loop represented as a sphere of a size proportional to
        # the number of edges (or nodes in case of bound-nodes).
        C = CircleFace(radius=20 * (popcount(data[node.name]) - min_node_size) / diff_node_size,
                    color=to_hex((r, 0.3, b)), style="circle")
        # Let's make the sphere transparent
        C.opacity = 0.9
        # And place as a float face over the tree
        faces.add_face_to_node(C, node, 0, position="float")
        return

    OR = node.OR
    r = 0.3
    g = 0.3
    b = 0.3
    w = 1
    if OR > 1 and node.pv < 0.05:
        r = max(0.35, min(0.95, OR - 1)**0.5)
        b = 0.3
        w = len(node)**0.5 / 2
    elif OR < 1 and node.pv < 0.05:
        b = max(0.35, min(0.95, (1 - OR)**0.5))
        r = 0.3
        w = len(node)**0.5 / 2

    else:
        w = max(3, w)
    node.img_style["hz_line_color"] = to_hex((r, 0.3, b))
    node.img_style["vt_line_color"] = to_hex((r, 0.3, b))
    node.img_style["vt_line_width"] = w
    node.img_style["hz_line_width"] = w


def my_layout_enrichments_GWAS(node):
    node.img_style['size'] = 0

    if node.is_root():
        return

    if node.is_leaf():
        if node.name in feature_loop['GWAS']:
            r = 0.85
            b = 0.3
        else:
            b = 0.85
            r = 0.3
        # each multi-loop represented as a sphere of a size proportional to
        # the number of edges (or nodes in case of bound-nodes).
        C = CircleFace(radius=20 * (popcount(data[node.name]) - min_node_size) / diff_node_size,
                    color=to_hex((r, 0.3, b)), style="circle")
        # Let's make the sphere transparent
        C.opacity = 0.9
        # And place as a float face over the tree
        faces.add_face_to_node(C, node, 0, position="float")
        return

    OR = node.OR
    r = 0.3
    g = 0.3
    b = 0.3
    w = 1
    if OR > 1 and node.pv < 0.05:
        r = max(0.35, min(0.95, OR - 1)**0.5)
        b = 0.3
        w = len(node)**0.5 / 2
    elif OR < 1 and node.pv < 0.05:
        b = max(0.35, min(0.95, (1 - OR)**0.5))
        r = 0.3
        w = len(node)**0.5 / 2

    else:
        w = max(3, w)
    node.img_style["hz_line_color"] = to_hex((r, 0.3, b))
    node.img_style["vt_line_color"] = to_hex((r, 0.3, b))
    node.img_style["vt_line_width"] = w
    node.img_style["hz_line_width"] = w

def my_layout_enrichments_feature(node):
    node.img_style['size'] = 0

    if node.is_root():
        return

    if node.is_leaf():
        if node.name in test_feature_loop:
            r = 0.85
            b = 0.3
        else:
            b = 0.85
            r = 0.3
        # each multi-loop represented as a sphere of a size proportional to
        # the number of edges (or nodes in case of bound-nodes).
        C = CircleFace(radius=20 * (popcount(data[node.name]) - min_node_size) / diff_node_size,
                    color=to_hex((r, 0.3, b)), style="circle")
        # Let's make the sphere transparent
        C.opacity = 0.9
        # And place as a float face over the tree
        faces.add_face_to_node(C, node, 0, position="float")
        return

    OR = node.OR
    r = 0.3
    g = 0.3
    b = 0.3
    w = 1
    if OR > 1 and node.pv < 0.05:
        r = max(0.35, min(0.95, OR - 1)**0.5)
        b = 0.3
        w = len(node)**0.5 / 2
    elif OR < 1 and node.pv < 0.05:
        b = max(0.35, min(0.95, (1 - OR)**0.5))
        r = 0.3
        w = len(node)**0.5 / 2

    else:
        w = max(3, w)
    node.img_style["hz_line_color"] = to_hex((r, 0.3, b))
    node.img_style["vt_line_color"] = to_hex((r, 0.3, b))
    node.img_style["vt_line_width"] = w
    node.img_style["hz_line_width"] = w



check = 0
#for enum_splits , nummber_splits in enumerate([4, 8, 16]):
for enum_splits , nummber_splits in enumerate([8]):
    tree = Tree(tree_path)
    #########################
   


    tree.add_feature('dist_root', 0)


    # get all possible cut lines to define different clusters
    # distances = set([tree.get_distance(n) for n in tree.iter_descendants()])
    distances = [2, 3, 4]

    # precompute distance to root for all nodes
    tree.add_feature('dist_root', 0)
    for n in tree.iter_descendants():
        if n.is_leaf():
            n.add_feature('dist_root', 1000)
            continue
        n.add_feature('dist_root', tree.get_distance(n))


    # get clusters
    splits = dict((2**d, [n for n in tree.iter_descendants()
                    if n.get_distance(tree, topology_only=True) >= d and
                n.up.get_distance(tree, topology_only=True) < d]) for d in sorted(distances))





    colors = ['firebrick', 'orange', 'royalblue', 'blue']

    if nummber_splits  > 4:
        colors = colors + ['olive', 'seagreen', 'mediumpurple', 'mediumorchid']
    if nummber_splits == 16:
        colors = colors + ['pink', 'green','gray', 'purple', 'brown', 'cyan', 'black', 'yellow']

    #metaloopcluster_path = f'{type_net}_loops_{type_enha}_tables_{type_filtered}_Custom_{ctcf_file}_{nummber_splits}_{data_type}_{metric}-metaloop_clusters.tsv'
    
    metaloopcluster_path = os.path.join(result_path, f'Custom_{ctcf_file}_{nummber_splits}_{data_type}_{metric}-metaloop_clusters.tsv')
    
    print(metaloopcluster_path)

    
    clusters = dict()
    ct = {}
    for feature in features:
        ct[feature] = sum(nn in feature_loop[feature] for nn in tree.iter_leaf_names())
    tt = len(tree)

    out = open(metaloopcluster_path, 'w')
    out_table = 'Gene\tcluster\t'
    
    
    annotation = {}
    for n, cluster in enumerate(splits[nummber_splits], 1):
        ORs = {}
        pvs = {}
        for feature in features:
            out_table += f'{feature}\t{feature}_OR\t{feature}_pv' + '\t'
            nfeature = sum(nn in feature_loop[feature] for nn in cluster.iter_leaf_names())
            lnode = len(cluster)
            
            lnode_no_feature = lnode - nfeature
            lnodes_out_cluster = tt - lnode
            lnodes_out_cluster_feature = ct[feature] - nfeature
            lnodes_out_cluster_no_feature = lnodes_out_cluster - lnodes_out_cluster_feature
            
            #if feature == "GWAS_Cat_prostate.carcinoma":
            #    IPython.embed()
            try:
                OR, pv = fisher_exact([[nfeature, lnode_no_feature], [lnodes_out_cluster_feature, lnodes_out_cluster_no_feature]])
            except:
                
                IPython.embed()

            #OR, pv = fisher_exact([[nfeature, lnode], [ct[feature], tt]])

            ORs[feature] = OR
            pvs[feature] = pv

            if ORs[feature] > 1 and pvs[feature] < 0.05:
                annotation[feature] = '+'
            elif ORs[feature] < 1 and pvs[feature] < 0.05:
                annotation[feature] = '-'
            else:
                annotation[feature] = '='
 
        out_table = out_table.strip() + "\n"
        clusters[cluster] = colors[n-1]

        out_line = '\n'.join(f'{leaf}\t{n}\t' + '\t'.join(f'{annotation[feature]}\t{ORs[feature]}\t{pvs[feature]}' 
        #out_line = '\n'.join(f'{leaf}\t{n}\t' + '\t'.join(f'{annotation[feature]}' 
            for feature in features) for leaf in cluster.iter_leaf_names())
        out_table += out_line.strip() + "\n"

        if check == 1:
            IPython.embed()
            check = 0
      
    out.write(out_table)
    out.close()

    min_node_size, max_node_size = np.percentile([popcount(data[node.name]) for node in tree], [1, 99])
    min_node_size = min_node_size * 0.8
    diff_node_size = max_node_size - min_node_size
    for n in tree.iter_descendants():
        # square root of distances for better visual
        n.dist = n.dist**0.5 * 10
        # shorten leave's branch-length for better visual too
        if n.is_leaf():
            n.dist=0.0



    ts = TreeStyle()
    ts.show_leaf_name = False
    # circular tree
    ts.mode = 'c'
    ts.root_opening_factor = 0.0
    ts.title.add_face(TextFace(f'Clustering of {data_type.replace("-", " ")} based on edge {metric}',
                            fsize=50), column=0)

    #metaloopcluster_tree_path = f'{type_net}_loops_{type_enha}_tables_{type_filtered}_Custom_clusters_{nummber_splits}_{data_type}_{metric}.png'
    metaloopcluster_tree_path = os.path.join(result_path, f'Custom_clusters_{nummber_splits}_{data_type}_{metric}.png')
    tree.render(metaloopcluster_tree_path, w=200, h=200, units="mm", tree_style=ts, layout=my_layout_clusters )

    
    for feature in features:
        if enum_splits != 0:
            break

        where = 'both' if feature == 'CTCF' else 'enhancer'
        where = 'both' if feature == 'CTCF' else 'promoter'

        ts = TreeStyle()
        ts.show_leaf_name = False
        # circular tree
        ts.mode = 'c'
        ts.root_opening_factor = 0.0
        ts.title.add_face(TextFace(f'{feature} on {where} clustering of {data_type.replace("-", " ")} based on edge {metric}',
                            fsize=50), column=0)
        for n in tree.iter_descendants():
            nfeature = sum(nn in feature_loop[feature] for nn in n.get_leaf_names())
            lnode = len(n)
            
            lnode_no_feature = lnode - nfeature
            lnodes_out_cluster = tt - lnode
            lnodes_out_cluster_feature = ct[feature] - nfeature
            lnodes_out_cluster_no_feature = lnodes_out_cluster - lnodes_out_cluster_feature
            try:
                OR, pv = fisher_exact([[nfeature, lnode_no_feature], [lnodes_out_cluster_feature, lnodes_out_cluster_no_feature]])
            except:
                print("Fishertest problem")
                IPython.embed()
            #OR, pv = fisher_exact([[nfeature, lnode], [ct[feature], tt]])

            try:
                n.del_feature('OR')
                n.del_feature('pv')
            except:
                pass

            n.add_feature('OR', OR)
            n.add_feature('pv', pv)

        my_layout_function = my_layout_enrichments_feature
        test_feature_loop = feature_loop[feature]

        
        #metaloopcluster_tree_path = f'{type_net}_loops_{type_enha}_tables_{type_filtered}_Custom_{feature}_{where}_enrichemnts_{data_type}_{metric}.png'
        metaloopcluster_tree_path = os.path.join(result_path, f'Custom_{feature}_{where}_enrichemnts_{data_type}_{metric}.png')

        tree.render(metaloopcluster_tree_path, w=200, h=200, units="mm", tree_style=ts, layout=my_layout_function )

