from collections import Counter
import numpy as np
from copy import copy
import math

#Copied from ToyFold 1D utils (Rhiju Das Matlab code originally)

def find_all(s, ch):
    if isinstance(ch, list):
        tmp_list=[]
        for c in ch:
            tmp_list.extend([i for i, ltr in enumerate(s) if ltr == c])
        return tmp_list
    elif isinstance(ch, str):
        return [i for i, ltr in enumerate(s) if ltr == ch]

def convert_structure_to_bps(secstruct):

    bps = []

    #Find other delimiters
    other_delimiters = [k for k in Counter(secstruct).keys() if k not in ['.','(',')','[',']','{','}',"<",">"]]

    for delim in other_delimiters:
        pos = find_all(secstruct, delim)
        assert(len(pos) % 2 == 0)

        N = int(len(pos)/2)
        i,j = pos[:N], pos[N:]

        if N > 1:
            assert(all(np.diff(i)==1))
            assert(all(np.diff(j)==1))

        for ind in range(N):
            bps.append([i[ind],j[-1-ind]])

    left_delimiters = ['(','{','[',"<"]
    right_delimiters = [')','}',']',">"]

    for (left_delim, right_delim) in list(zip(left_delimiters, right_delimiters)):

        left_list = []
        for i, char in enumerate(secstruct):
            if char == left_delim:
                left_list.append(i)

            elif char == right_delim:
                bps.append([left_list[-1],i])
                left_list = left_list[:-1]

        assert len(left_list)==0

    return bps

def parse_stems_from_bps(bps, debug=False):
    #Creates list of stems, where each stem is a list of base pairs.
    # i.e. '((.))' -> [[0,4],[1,3]]
    # '((.)).((.))' -> [[[0,4],[1,3]],[[6,10],[7,9]]]
    
    if debug: print(bps)

    if len(bps) == 0:
        stems = []
    else:
        nres = np.max(bps)
        stems = []

        while len(bps) > 0:
            bp = bps[0]
            bps = bps[1:]

            stem = [bp]

            if debug: print('stem init', stem)
            bp_next = copy(bp)
            if debug: print('bp_next', bp_next)

            # Check outward
            for i in list(reversed(range(bp[0]))):
                bp_next = [copy(bp_next)[0]-1,copy(bp_next)[1]+1]

                if debug: print('next_out', bp_next)
                if debug: print('bps here', bps)
                if len(bps) > 0:
                    gp = find_all([x[0] for x in bps], [bp_next[0]])
                    if debug: print('gp', gp)
                    if len(gp)>0:
                        if bps[gp[0]][1] == bp_next[1]: # found an extension
                            if debug: print('r')
                            stem.append(copy(bp_next))
                            del bps[gp[0]] # take out of bp list
                    else:
                        break

            bp_next = copy(bp)

            #Check inward
            for i in range(bp[0],nres+1):

                bp_next[0] = copy(bp_next)[0]+1
                bp_next[1] = copy(bp_next)[1]-1

                if debug: print('next_in', bp_next)
                if len(bps) > 0:
                    gp = find_all([x[0] for x in bps], [bp_next[0]])
                    if len(gp)>0:
                        if bps[gp[0]][1] == bp_next[1]: # found an extension
                            if debug: print('h')
                            stem = [copy(bp_next)]+copy(stem)
                            del bps[gp[0]] # take out of bp list
                        else:
                            break
            stems.append(stem)
            if debug: print('stem', stem)
    if debug: print('stems', stems)
    return stems

def parse_out_chainbreak(secstruct):
    secstruct_new = []
    is_chainbreak = []

    chainbreak_ctr=0

    for i, char in enumerate(secstruct):
        if char in [',','+',' ','&']:
            is_chainbreak.append(chainbreak_ctr)
        else:
            secstruct_new.append(char)
            chainbreak_ctr+=1
    return is_chainbreak, ''.join(secstruct_new)


def get_stem_assignment(secstruct):
    '''Returns vector length N_beads, 0 if not in a stem, otherwise assigned stem 1 to max(N_stems)
    Note basically not switched to zero-indexing for python version, unlike partner syntax'''

    is_chainbreak, secstruct = parse_out_chainbreak(secstruct)

    if not secstruct[0].isdigit():
        bps = convert_structure_to_bps(secstruct)

    else: # assume partner vector was given
        partner = secstruct
        bps = []
        for i in range(len(partner)):
            if partner[i] > i:
                bps.append([i, partner[i]])

    stems = parse_stems_from_bps(bps)

    stems = [x for x in sorted(stems)]

    stem_assignment = np.zeros([len(secstruct)])

    for i, stem in enumerate(stems):
        for bp in stem:
            stem_assignment[bp[0]] = i+1
            stem_assignment[bp[1]] = i+1

    return stem_assignment

def get_pairmap(secstruct):
    """
    (lifted from draw_rna)
    generates a list containing pair mappings

    args:
    secstruct contains secondary structure string

    returns:
    list with pair mappings
    """
    pair_stack = []
    end_stack = []
    pairs_array = []
    i_range = list(range(0,len(secstruct)))

    # initialize all values to -1, meaning no pair
    for ii in i_range:
        pairs_array.append(-1)

    left_delimiters = ['(','{','[',"<"]
    right_delimiters = [')','}',']',">"]

    # assign pairs based on secstruct
    for ii in i_range:
        if(secstruct[ii] in left_delimiters):
            pair_stack.append(ii)
        elif(secstruct[ii] in right_delimiters):
            if not pair_stack:
                end_stack.append(ii)
            else:
                index = pair_stack.pop()
                pairs_array[index] = ii
                pairs_array[ii] = index
    if len(pair_stack) == len(end_stack):
        n = len(pair_stack)
        for ii in range(n):
            pairs_array[pair_stack[ii]] = end_stack[-ii]
            pairs_array[end_stack[-ii]] = pair_stack[ii]
    else:
         print("ERROR: pairing incorrect %s" % secstruct)

    return pairs_array

def compose_structs(list_of_rg_objs, label_list=None):
    
    if label_list is not None:
        assert len(list_of_rg_objs)==len(label_list)
    rg1 = list_of_rg_objs[0]
    plot_nodes1 = [n for n in list(rg1.G.nodes) if isinstance(n,str)]
    subgraph1 = rg1.G.subgraph(plot_nodes1)
    new = nx.relabel.relabel_nodes(subgraph1, {k: '1%s' % k for k in subgraph1.nodes})

    for ind,rg2 in enumerate(list_of_rg_objs[1:]):
        plot_nodes2 = [n for n in list(rg2.G.nodes) if isinstance(n,str)]
        subgraph2 = rg2.G.subgraph(plot_nodes2)
        new_sub2 = nx.relabel.relabel_nodes(subgraph2, {k: '%d%s' % (ind+2, k) for k in subgraph2.nodes})

        new = nx.compose(new, new_sub2)

    G = new.to_undirected()
    pos =graphviz_layout(G,prog='neato')
    colors = [new[u][v]['color'] for u,v in new.edges()]
    colors = [new[u][v]['color'] for u,v in new.edges()]
    
    ax = plt.gca()
    for u,v in G.edges():
        
        x = [pos[u][0],pos[v][0]]
        y = [pos[u][1],pos[v][1]]
        l = Line2D(x,y, linewidth=2, solid_capstyle='round', color=G[u][v]['color'])
        ax.add_line(l)
    ax.autoscale()
    
    if label_list is not None:
    
        for k in range(len(list_of_rg_objs)):
            for x in pos.keys():
                if x.startswith("%d" % (k+1)):
                    ax.text(pos[x][0],pos[x][1], label_list[k])
                    break
    #colors = [new[u][v]['color'] for u,v in subgraph.edges()]
    #nx.draw(new, pos, width=3, node_size=0, edge_color = colors, arrows=False, solid_capstyle='round')

def _flip_helix(x,y,left,right):
    # get center line between 2 groups
    # left - list of indices 
    # right - list of indices 
    class1 = []
    class2 = []
    labels1 = []
    labels2 =[]
    for i in left:
        class1.append(np.array([x[i],y[i]]))
        labels1.append(-1)
    for i in right:
        class2.append(np.array([x[i],y[i]]))
        labels1.append(1)
    class1 = np.array(class1)
    class2 = np.array(class2)
    X = np.vstack((class1, class2))
    m = len(X)
    X = np.array([np.ones(m), X[:, 0], X[:, 1]]).T
    Y = np.concatenate((labels1, labels2)).T
    beta = np.linalg.inv(X.T @ X) @ (X.T @ Y)
    # 0 = b0 +b1x +b2y
    
    # reflect all points over center line
    for i in left+right:
        temp = -2*(beta[0]+beta[1]*x[i]+beta[2]*y[i])/(beta[1]**2+beta[2]**2)
        x[i] = temp*beta[1]+x[i]
        y[i] = temp*beta[2]+y[i]
    return x,y


def _move_group(x,y,offset,group):
    # group - list of indices 
    # offset - [x,y] for movement
    # move all items in the group
    for i in group:
        x[i] += offset[0]
        y[i] += offset[1]
    return x,y


def _rotate_group(x,y,angle,group):
    # group - list of indices 
    # angle in degrees
    rotation =  math.radians(angle)
    sum_x = 0
    sum_y = 0
    num_points = len(group)
    for i in group:
        sum_x += x[i]
        sum_y += y[i]
    centroid = [sum_x/num_points, sum_y/num_points]
    for i in group:
        x_orig = x[i]
        x[i] = centroid[0] + math.cos(rotation) * (x[i] - centroid[0]) - math.sin(rotation) * (y[i] - centroid[1])
        y[i] = centroid[1] + math.sin(rotation) * (x_orig - centroid[0]) + math.cos(rotation) * (y[i] - centroid[1])
    return x,y
