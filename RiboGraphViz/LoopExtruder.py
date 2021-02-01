from RiboGraphViz import RGV
import re
import numpy as np

def StackExtruder(sequence, structure, data=None, stack_size=0):
    assert len(sequence) == len(structure)
    mdl = RGV(structure, sequence=sequence)

    struct = '('*stack_size+' '+')'*stack_size

    full_list_motifs = []
    full_list_data = []

    for stem in mdl.stems:
        for stack_ind in range(len(stem)-stack_size+1):

            stack = stem[stack_ind:stack_ind+stack_size]
            
            side_1, side_2 = np.array(stack).T
            side_1 = list(reversed(side_1))
            side_2 = list(side_2)

            if data is not None:
                dat_vec = [data[x] for x in side_1]+[-1]+[data[x] for x in side_2]

            seq_vec = ''.join([mdl.sequence[x] for x in side_1]+[' ']+[mdl.sequence[x] for x in side_2])
            
            full_list_motifs.append("%s,%s" % (seq_vec,struct))

            if data is not None:
                full_list_data.append(np.array(dat_vec))
            #print(dat_vec)

    if data is not None:

        returning_dct = {}

        for i, k in enumerate(full_list_motifs):
            if k in returning_dct.keys():
                returning_dct[k].append(full_list_data[i])
            else:
                returning_dct[k] = [full_list_data[i]]

        return returning_dct
    else:
        return full_list_motifs

def LoopExtruder(sequence, structure, data=None, neighbor_bps=0):
    '''Given a sequence and structure, extract loops from the structure.
    If a data vector is provided, also extracts and returns the data values.

    Returns (if no data):
    List of motifs in format 'UUA CU,... ..' = 3x2 internal loop

    Returns (if data is provided):
    returning_dct(key, value): keys are motifs, values are list of associated per-nucleotide data.
    '''
    
    assert len(sequence) == len(structure)

    try:
        rg = RGV(structure, sequence=sequence)
    except:
        print('invalid RGV')
        raise RuntimeError('invalid RGV')

    nodes = [n for n in list(rg.G.nodes) if not isinstance(n,str)]
    string_assignment = ''.join(["%d" % int(x) for x in rg.stem_assignment])
    full_list_motifs = []
    full_list_data = []

    for i in range(1, max(nodes)+1):
        loop_motif=[]
        loop_seq=[]
        loop_data = []

        children = [x for x in rg.G[i]]
        #print(i,"has children", children)
        lst = [i]+children+[i]
        #print('working on node %d' % i)
        #print(lst)

        for j in range(len(lst)-1):
            start_stretch, end_stretch = lst[j], lst[j+1]
            if start_stretch == end_stretch:
                obj = re.search(r"%d00*%d" % (start_stretch, end_stretch), string_assignment) #in a hairpin
            else:
                obj = re.search(r"%d0*%d" % (start_stretch, end_stretch), string_assignment)

            start_ind, end_ind = obj.start(), obj.end()

            loop_motif.append(rg.secstruct[(start_ind-neighbor_bps+1):(end_ind+neighbor_bps-1)])
            loop_seq.append(rg.sequence[(start_ind-neighbor_bps+1):(end_ind+neighbor_bps-1)])
            
            if data is not None:
                loop_data.extend(data[(start_ind-neighbor_bps+1):(end_ind+neighbor_bps-1)])
                loop_data.extend([-1])

        motif = ' '.join(loop_seq) + ',' + ' '.join(loop_motif)
        if data is not None:
            del loop_data[-1]

            assert(len(' '.join(loop_seq)) == len(loop_data))
        assert(len(' '.join(loop_seq)) == len(' '.join(loop_motif)))

        full_list_motifs.append(motif)

        if data is not None:
            full_list_data.append(np.array(loop_data))

    if data is not None:
        returning_dct = {}

        for i, k in enumerate(full_list_motifs):
            if k in returning_dct.keys():
                returning_dct[k].append(full_list_data[i])
            else:
                returning_dct[k] = [full_list_data[i]]

        return returning_dct

    else:
        return full_list_motifs