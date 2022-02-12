import networkx as nx
from numpy import *
from numpy import random

def get_angular(x,y):
    return sqrt(x**2+y**2), arctan2(y,x)

def get_cartesian(r,phi):
    return r*cos(phi), r*sin(phi)

def set_intervals(DG,n):
    if n == DG.graph['root']:
        DG.nodes[n]['lowerbound'] = 0
        DG.nodes[n]['upperbound'] = 2*pi
        #DG.nodes[n]['lowerbound'] = pi/12.
        #DG.nodes[n]['upperbound'] = pi-pi/12.
        DG.nodes[n]['position'] = 0

    if DG.nodes[n]['branchmass']>1:
        current_branchmass = 0.
        current_lowerbound = DG.nodes[n]['lowerbound']
        intervallength = DG.nodes[n]['upperbound'] - DG.nodes[n]['lowerbound']
        children = list(DG.successors(n))
        for c in children:
            current_branchmass += DG.nodes[c]['branchmass']
        for i in range(len(children)):
            DG.nodes[children[i]]['lowerbound'] = current_lowerbound
            DG.nodes[children[i]]['upperbound'] = current_lowerbound + DG.nodes[children[i]]['branchmass']/float(current_branchmass) * intervallength
            current_lowerbound  = DG.nodes[children[i]]['upperbound']
            set_intervals(DG,children[i]) 

def get_sorted_by_dist(nodes_,dists_):
    if len(nodes_)>0:
        dists, nodes = [ list(t) for t in zip(*sorted(zip(dists_,nodes_))) ]
        return nodes, dists
    else:
        return nodes_, dists_


def get_distributed_by_dist(nodes_, dists_):
    nodes, dists = get_sorted_by_dist(nodes_,dists_)
    n_l = len(nodes)
    if len(nodes)>1:
        nodes = nodes[0::2] + (nodes[1::2])[::-1]
        dists = dists[0::2] + (dists[1::2])[::-1]
    return nodes, dists

def get_inversely_distributed_by_dist(nodes_, dists_):
    nodes, dists = get_sorted_by_dist(nodes_, dists_)
    n_l = len(nodes)
    if len(nodes)>1:
        nodes = (nodes[0::2])[::-1] + (nodes[1::2])
        dists = (dists[0::2])[::-1] + (dists[1::2])
    return nodes, dists

def set_angular_positions(DG,n,dist_key):
    if n == DG.graph['root']:
        #DG.nodes[n]['lowerbound'] = pi/10.
        #DG.nodes[n]['upperbound'] = pi-pi/10.
        DG.nodes[n]['lowerbound'] = 0.
        DG.nodes[n]['upperbound'] = 2.*pi

    if DG.nodes[n]['branchmass']>1:
        current_branchmass = 0.
        current_lowerbound = DG.nodes[n]['lowerbound']
        intervallength = DG.nodes[n]['upperbound'] - DG.nodes[n]['lowerbound']

        children = DG.successors(n)
        for c in children:
            current_branchmass += DG.nodes[c]['branchmass']

        branches = [ c for c in children if DG.nodes[c]['branchmass']>1 ]
        leaves = list(set(children) - set(branches))
        branch_dists = [ 0.5*(min([ DG.nodes[c][dist_key] for c in DG.successors(b)])+DG.nodes[b][dist_key]) for b in branches ]
        leave_dists = [ DG.nodes[l][dist_key] for l in leaves ]

        if len(branches) == 0 and len(leaves) > 0:
            leaves, leave_dists = get_distributed_by_dist(leaves,leave_dists)
            for l in leaves:
                dp = 1./float(current_branchmass) * intervallength
                DG.nodes[l]['position'] = current_lowerbound + dp * 0.5
                current_lowerbound += dp
        elif len(leaves) == 0 and len(branches) > 0:
            branches, branch_dists = get_distributed_by_dist(branches,branch_dists)
            for b in branches:
                DG.nodes[b]['lowerbound'] = current_lowerbound
                DG.nodes[b]['upperbound'] = current_lowerbound + DG.nodes[b]['branchmass']/float(current_branchmass) * intervallength
                current_lowerbound  = DG.nodes[b]['upperbound']
                DG.nodes[b]['position'] = 0.5 * (DG.nodes[b]['upperbound'] + DG.nodes[b]['lowerbound'])
                set_angular_positions(DG,b,dist_key)
        else:
            #sort both lists for distance
            branches, branch_dists = get_sorted_by_dist(branches,branch_dists)
            leaves, leave_dists = get_sorted_by_dist(leaves,leave_dists)

            #fill the space underneath branches with leaves of lower distances
            #but maximally with a number of leaves that's lower than the number
            #of the branch's children
            leaves_of_branch = [[]]
            child_count_of_branch = [ len(DG.successors(b))-1 for b in branches ]
            bmark = 0
            lmark = 0
            while bmark<len(branches) and lmark<len(leaves):
                if leave_dists[lmark] < branch_dists[bmark] and \
                   child_count_of_branch[bmark] > len(leaves_of_branch[bmark]):
                    lmark += 1
                else:
                    bmark += 1
                    leaves_of_branch.append([])

            for bm in range(bmark+1,len(branches)):
                leaves_of_branch.append([])

            for bmark in range(len(branches)):
                child_count_of_branch[bmark] = min([ child_count_of_branch[bmark],
                                                     int(DG.nodes[branches[bmark]]['branchmass']/float(current_branchmass) * (lmark+1))])

            bmark = 0
            lmark = 0
            while bmark<len(branches) and lmark<len(leaves):
                if leave_dists[lmark] < branch_dists[bmark] and \
                   child_count_of_branch[bmark] > len(leaves_of_branch[bmark]):
                    leaves_of_branch[bmark].append(leaves[lmark])
                    lmark += 1
                else:
                    bmark += 1

            #now define angular fractions based on the branch masses for all branches and the remaining leaves
            remaining_leaves, remaining_dists = get_distributed_by_dist(leaves[lmark:], leave_dists[lmark:])
            #remaining_leaves, remaining_dists = get_inversely_distributed_by_dist(leaves[lmark:], leave_dists[lmark:])

            branch_masses = [ DG.nodes[b][dist_key] for b in branches ]
            branches_insert, branch_masses = get_inversely_distributed_by_dist(branches, branch_dists)
            n_l = len(remaining_leaves)
            current_branchmass -= (len(leaves) - len(remaining_leaves))
            current_nodes = remaining_leaves[:n_l//2] + branches_insert + remaining_leaves[n_l//2:]
            #current_dists = leave_distsbranch_dists + leave_dists[lmark:]
            #current_nodes, current_dists = get_inversely_distributed_by_dist(current_nodes,current_dists)
            for c in current_nodes:
                if c in branches_insert:
                    DG.nodes[c]['lowerbound'] = current_lowerbound
                    DG.nodes[c]['upperbound'] = current_lowerbound + DG.nodes[c]['branchmass']/float(current_branchmass) * intervallength
                    current_lowerbound  = DG.nodes[c]['upperbound']
                if c in leaves:
                    dp = DG.nodes[c]['branchmass']/float(current_branchmass) * intervallength
                    DG.nodes[c]['position'] = current_lowerbound + dp * 0.5
                    current_lowerbound += dp

            #arrange all branches and leaves and call self for the next branch
            for b in branches:

                ib = branches.index(b)
                dists = [ DG.nodes[l][dist_key] for l in leaves_of_branch[ib] ]
                leaves_of_branch[ib], dists = get_inversely_distributed_by_dist(leaves_of_branch[ib],dists)

                if len(leaves_of_branch[ib])==2:
                    print(leaves_of_branch[ib])
                    print(len(leaves_of_branch[ib]))
                    print("position: ", len(leaves_of_branch[ib])//2)
                leaves_of_branch[ib].insert(len(leaves_of_branch[ib])//2,b)
                if len(leaves_of_branch[ib])==3:
                    print(leaves_of_branch[ib])

                n_l = len(leaves_of_branch[ib])

                current_lowerbound = DG.nodes[b]['lowerbound']
                current_branchmass = float(len(leaves_of_branch[ib]))
                intervallength = DG.nodes[b]['upperbound'] - current_lowerbound
                for l in leaves_of_branch[ib]:
                    dp = 1./float(current_branchmass) * intervallength
                    DG.nodes[l]['position'] = current_lowerbound + dp * 0.5
                    current_lowerbound += dp

                set_angular_positions(DG,b,dist_key)

def prepare_tree(DG):
    nodes = DG.nodes()
    #nodes = list(set(nodes)-set([root]))
    #max_layer = 0
    for i in nodes:
        DG.nodes[i]['branchmass'] = len(nx.descendants(DG,i))+1
        #if i==root:
        #    DG.nodes[i]['parent'] = None
        #else:
        #    DG.nodes[i]['parent'] = DG.predecessors(i)[0]
        #DG.nodes[i]['layer']=len(nx.ancestors(DG,i))
        #max_layer = max([DG.nodes[i]['layer'],max_layer])

    set_intervals(DG,DG.graph['root'])

def get_initial_positions(DG,dist_key):
    nodes = DG.nodes()
    positions = {}
    for i in nodes:
        positions[i] = tuple(get_cartesian(DG.nodes[i][dist_key],
                             0.5*(DG.nodes[i]['upperbound'] + DG.nodes[i]['lowerbound']))
                        )
    return positions


def get_sophisticated_positions(DG,root,dist_key):

    #print "hello"
    set_angular_positions(DG,root,dist_key)
    positions = {}
    nodes = list(DG.nodes())
    print(nodes)
    print(DG.nodes['b'])
    for i in nodes:
        print(i)
        positions[i] = tuple(get_cartesian(DG.nodes[i][dist_key],
                                           DG.nodes[i]['position'])
                            )
        #positions[i] = tuple((DG.nodes[i][dist_key],
        #                     DG.nodes[i]['position']))

    #positions[root] = (0,pi)    
    return positions


def radial_distance_layout(tree,dist_key,mode='normal',save_data_to_tree=False):
    """
    Compute node positions of a tree where nodes are associated
    with a distance to the root node.

    Parameters
    ==========
    tree : networkx.DiGraph
        The tree for which to compute the layout
    dist_key : hashable type
        node property that carries the distance to the root
    mode : str
        can be 'normal' or 'soph'. Just go for 'normal',
        'soph' sucks.
    save_data_to_tree : bool, default = False
        whether or not to write the positions to
        the nodes of the tree

    Returns
    =======
    positions : dict
        keys are tree nodes, values are 2-tuples
        of x- and y-values
    """

    if save_data_to_tree:
        DG = tree
    else:
        DG = nx.DiGraph(tree)
    root = [ node for (node,ndict) in DG.nodes(data=True) if ndict[dist_key]==0 ][0]
    DG.graph['root'] = root
    prepare_tree(DG)

    if mode=='soph':
        return get_sophisticated_positions(DG,root,dist_key)
    else:
        return get_initial_positions(DG,dist_key)

if __name__=="__main__":

    import pylab as pl

    paths  = [ [ 'a','b','c'] ]
    paths += [ [ 'a','b','d'] ]
    paths += [ [ 'a','e','f','g'] ]
    paths += [ [ 'a','e','f','h'] ]
    paths += [ [ 'a','e','i'] ]
    paths += [ [ 'a','j','k'] ]
    paths += [ [ 'a','j','l'] ]

    dists = {'a': 0,
             'b':1.1, 'e': 1.2, 'j': 1.4,
             'c':2.1, 'd': 2.2, 'f': 2.1, 'i': 2.34, 'k':3.8, 'l':2.5,
             'g': 3.9, 'h': 3.8}
    T = nx.DiGraph()

    for p in paths:
        T.add_path(p)

    keystr = 'dist'

    nx.set_node_attributes(T,dists,keystr)

    fig,ax = pl.subplots(1,2,figsize=(15,8))

    #pos = radial_distance_layout(T,keystr,mode='soph')
    #nx.draw_networkx(T,pos,ax=ax[0])
    pos = radial_distance_layout(T,keystr,mode='normal')
    nx.draw_networkx(T,pos,ax=ax[1])
    pl.show()




