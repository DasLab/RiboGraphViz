import networkx as nx
import RiboGraphViz.RG_utils as utils
from copy import copy
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.colors as mcolors
from matplotlib import cm
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Circle
import sys
import numpy as np

class RiboGraphViz(object):

    def __init__(self, secstruct, sequence=None):
        '''Create a NetworkX graph representing an RNA secondary structure, where nodes represent
        hairpin loops, internal loops, bulges, anything not a helix, and stems represent helices.

        Edge weights / lengths are equal to the number of base pairs present in the helix.
        Nodes are sized according to the number of unpaired bases in the loop.

        Input: secondary structure in dot-parentheses notation. Currently cannot handle pseudoknots.

        Attributes:

        G: NetworkX Directed Graph object. Amenable to much more analysis.

        n_hairpins, n_internal_loops, n_multiloops : number of loops
        n_helices : number of helices
        max_ladder_distance : maximum end_to_end distance of helices
         present in structure, not counting lengths of loops.

        loop_sizes (dictionary) : keys correspond to loop numbering (node numbering),
         values hold the number of unpaired bases present in each loop.
        '''

        self.chainbreak, self.secstruct = utils.parse_out_chainbreak(secstruct)

        if len(self.chainbreak) > 0:
            self.is_multi = True
        else:
            self.is_multi = False

        stems = utils.parse_stems_from_bps(utils.convert_structure_to_bps(self.secstruct))
        self.stems = [x for x in sorted(stems)]
        self.stem_assignment = utils.get_stem_assignment(self.secstruct)
        self.pairmap = utils.get_pairmap(self.secstruct)
        self.G = nx.DiGraph()
        self.loop_sizes = { 0:0 }
        self.helix_multiplier = 1
        self.helix_width = 1
        self.N = len(self.secstruct)
        self.sequence=sequence

        self.fiveprime_x, self.fiveprime_y, self.threeprime_x, self.threeprime_y = None, None, None, None

        self.setup_graph()
        self.calc_MLD()

        self.n_hairpins, self.n_internal_loops, self.n_3WJs, self.n_4WJs, self.n_5WJs_up = self.count_loops()

    def get_53_distance(self):

        plot_nodes = [n for n in list(self.G.nodes) if isinstance(n,str)]
        subgraph = self.G.subgraph(plot_nodes)
        subgraph = subgraph.to_undirected()

        if not self.is_multi:
            #locate 5' nucleotide
            if 'n0' in list(subgraph.nodes()):
                five_loc = 'n0'
            else:
                five_loc = 'h1b'

            if 'n%d' % (self.N-1) in list(subgraph.nodes()):
                three_loc = 'n%d' % (self.N-1)

    def calc_MLD(self):
        nodes = [n for n in list(self.G.nodes)] # if not isinstance(n,str)]
        subgraph = self.G.subgraph(nodes).to_undirected()

        first_helix_node = None
        for n in nodes:
            if isinstance(n, str):
                if n.startswith('h'):
                    first_helix_node = n
                    break


        if first_helix_node is None:
            self.MLD = 0
        else:
            node1 = [x for x in nx.traversal.bfs_edges(subgraph, first_helix_node)][-1][-1]
            node2 = [x for x in nx.traversal.bfs_edges(subgraph, node1)][-1][-1]

            self.MLD = nx.shortest_path_length(subgraph, node1, node2, weight='MLD_weight')

    def setup_graph(self):
        '''Create graph by reading pairmap array and recursively creating edges.'''
        
        self.G.add_node(0)
        jj = 0

        while (jj < len(self.pairmap)):

            if self.pairmap[jj] > jj:

                # we're in the structure here
                self.add_edges_r_(jj, self.pairmap[jj],0,0)
                jj = self.pairmap[jj]+1

            else:
                #we're in external loop here
                jj += 1

                self.loop_sizes[0] += 1

        stem_assignment_left = np.concatenate([np.array([-1]), self.stem_assignment[:-1]])
        stem_assignment_right = np.concatenate([self.stem_assignment[1:], np.array([-1])])

        for nod in range(len(self.stem_assignment)):
            if (nod,nod) in self.G.edges:

                self.G.remove_edge(nod,nod)

        for i, pair in enumerate(self.pairmap):
            if pair==-1:
                self.G.add_node('n%d' % i)

                if stem_assignment_left[i] > 0:
                    if self.secstruct[i-1]=='(':
                        letter='a'
                    elif self.secstruct[i-1]==')':
                        letter='b'
                    if i not in self.chainbreak:
                        self.G.add_edge('n%d' % i,'h%d%s' % (self.stem_assignment[i-1], letter), len=1.25, MLD_weight=0)

                    # add helix_a node here

                if stem_assignment_right[i] > 0:
                    #print(i)

                    if self.secstruct[i+1]==')':
                        letter='a'
                    elif self.secstruct[i+1]=='(':
                        letter='b'
                    if i+1 not in self.chainbreak:
                        self.G.add_edge('n%d' % i,'h%d%s' % (self.stem_assignment[i+1], letter), len=1.25, MLD_weight=0)

                    # add helix_b node here

        nuc_nodes = [n for n in list(self.G.nodes) if isinstance(n,str) and n.startswith('n')] # hacky
        for nuc in nuc_nodes:
            ind=int(nuc.replace('n',''))
            if 'n%d' % (ind-1) in nuc_nodes:
                # check for chainbreak

                if ind not in self.chainbreak:
                    self.G.add_edge('n%d' % (ind-1), 'n%d' % ind, len=1, MLD_weight=0)

                #print('n%d' % (ind-1), 'n%d' % ind,)


        for stem_ind in range(1,len(self.stems)+1):


            stem_length = len(self.stems[stem_ind-1])
            color_ind = np.where(self.stem_assignment==stem_ind)[0][0]

            self.G.add_edge('h%da' % (stem_ind),'h%db' % (stem_ind), len=self.helix_multiplier*(stem_length-0.99), MLD_weight=stem_length)

        for i in range(self.N-1):

            if self.stem_assignment[i+1] != self.stem_assignment[i]:

                if self.stem_assignment[i+1] != 0.0 and self.stem_assignment[i] != 0.0:
                    stem_ind_1 = self.stem_assignment[i]
                    stem_ind_2 = self.stem_assignment[i+1]
                    #color_ind = np.where(self.stem_assignment==stem_ind_1)[0][0]

                    #print("DEBUG", stem_ind_1, stem_ind_2)
                    #print("DEBUG", self.secstruct[i:i+2])

                    if self.secstruct[i:i+2] == '((':
                        self.G.add_edge('h%da' % stem_ind_1,'h%db' % stem_ind_2, len=1, weight=1, MLD_weight=0)
                        #self.G.add_edge(int(stem_ind_1), int(stem_ind_2)) # hkws hack dec 14
                    elif self.secstruct[i:i+2] == ')(':
                        if i+1 not in self.chainbreak:
                            self.G.add_edge('h%db' % stem_ind_1,'h%db' % stem_ind_2, len=1, weight=1, MLD_weight=0)
                            #self.G.add_edge(stem_ind_1, stem_ind_2)

                    elif self.secstruct[i:i+2] == '))':
                        self.G.add_edge('h%db' % stem_ind_1,'h%da' % stem_ind_2, len=1, weight=1, MLD_weight=0)
                        #self.G.add_edge(stem_ind_1, stem_ind_2)

    def add_edges_r_(self, start_index, end_index, last_helix, last_loop):
        '''Recursive method to add edges to graph.'''

        if (start_index > end_index):
            print('Error, start_index > end_index')
            sys.exit(0)
            
        if (self.pairmap[start_index] == end_index):
            #print("we're in helix %d! keep going" % self.stem_assignment[start_index])

            self.add_edges_r_(start_index+1, end_index-1, self.stem_assignment[start_index], last_loop)

        else:
            jj = start_index
            while jj <=end_index:
                if self.pairmap[jj] > jj:
                    #print('where are we?')
                    self.G.add_edge(int(last_loop), int(last_helix)) #int(self.stem_assignment[jj])) # needs more attributes?
                    last_loop = copy(last_helix)

                    self.add_edges_r_(jj, self.pairmap[jj],self.stem_assignment[jj], last_loop)
                    jj = self.pairmap[jj] + 1 #leaving helix
                else:
                    #print("DEBUG", jj,'in a loop!')
                    #print('DEBUG, Last helix was %d, last loop was %d' % (last_helix, last_loop))
                    #print('DEBUG', self.G.nodes)
                    if last_helix not in list(self.G.nodes):
                        #print(last_helix, len(stems[int(last_helix-1)]))
                        self.G.add_edge(int(last_loop), int(last_helix),
                                   len = len(self.stems[int(last_helix-1)]),
                                  weight = len(self.stems[int(last_helix-1)]), MLD_weight=0)

                        last_loop = copy(last_helix)

                    if int(last_helix) not in self.loop_sizes.keys():
                        self.loop_sizes[int(last_helix)] = 0
                        #self.G.add_edge('n%d'%jj, int(last_helix))

                    self.loop_sizes[int(last_helix)] += 1

                    jj+=1

    def get_coordinates(self, align=False, align_mode='COM'):
        '''
        align (bool): set first nucleotide at [0,0] and rotates structure according to align_mode.
        align_mode ("COM","end"): if 'COM', aligns center of mass to x axis. if "end", aligns 3' end to x axis.
        Return: array of x_coords, array of y_coords 
        '''

        if self.is_multi and align:
            raise RuntimeError('Alignment and multiple structures not yet implemented.')

        plot_nodes = [n for n in list(self.G.nodes) if isinstance(n,str)]
        subgraph = self.G.subgraph(plot_nodes)
        subgraph = subgraph.to_undirected()

        if not self.is_multi:
            #locate 5' nucleotide
            if 'n0' in list(subgraph.nodes()):
                subgraph.add_edge('n0',"5'",len=2)
            else:
                subgraph.add_edge('h1b',"5'",len=2)

            #locate 3' nucleotide. Only located if in a distinct string, otherwise can't distinguish end of outer stem.

            if 'n%d' % (self.N-1) in list(subgraph.nodes()):
                subgraph.add_edge('n%d' % (self.N-1),"3'",len=2)

        pos = graphviz_layout(subgraph, prog='neato')

        if not self.is_multi:
            self.fiveprime_x, self.fiveprime_y = pos["5'"]

            if "3'" in list(subgraph.nodes()):
                self.threeprime_x, self.threeprime_y = pos["3'"]


        for (u,v) in list(subgraph.edges()):
            if u.startswith('n') and v.startswith('n'):
                x1, x2, y1, y2 = pos[u][0], pos[v][0],pos[u][1], pos[v][1]
                break

        bond_width = np.sqrt((x2-x1)**2+(y2-y1)**2)
        node_positions_x, node_positions_y = {},{}

        for node in list(subgraph.nodes()):

            if node.startswith('n'):
                seq_ind = int(node[1:])
                node_positions_x[seq_ind] = pos[node][0]
                node_positions_y[seq_ind] = pos[node][1]

        for i, stem in enumerate(self.stems):

            start_x, start_y = pos['h%da' % (i+1)]
            fin_x, fin_y = pos['h%db' % (i+1)]

            x_dist = fin_x - start_x
            y_dist = fin_y - start_y

            path_length = np.sqrt((start_x - fin_x)**2 + (start_y - fin_y)**2)
            stem_angle = np.arctan2(y_dist, x_dist)

            x_diff = np.cos(stem_angle+np.pi/2)*0.5*bond_width*self.helix_width
            y_diff = np.sin(stem_angle+np.pi/2)*0.5*bond_width*self.helix_width
            # equidistant divide stem length, scatter points along it

            for j in range(len(stem)):
                x_midpoint = start_x+j* x_dist/(len(stem)-0.99)
                y_midpoint = start_y+j* y_dist/(len(stem)-0.99)

                if stem_angle < 0:
                    node_positions_x[stem[j][1]] = x_midpoint+x_diff
                    node_positions_x[stem[j][0]] = x_midpoint-x_diff
                    node_positions_y[stem[j][1]] = y_midpoint+y_diff
                    node_positions_y[stem[j][0]] = y_midpoint-y_diff

                else:
                    node_positions_x[stem[j][0]] = x_midpoint+x_diff
                    node_positions_x[stem[j][1]] = x_midpoint-x_diff
                    node_positions_y[stem[j][0]] = y_midpoint+y_diff
                    node_positions_y[stem[j][1]] = y_midpoint-y_diff

        if self.threeprime_x is not None:
            self.threeprime_x = self.threeprime_x - node_positions_x[0]
            self.threeprime_y = self.threeprime_y - node_positions_y[0]

        if self.fiveprime_x is not None:
            self.fiveprime_x = self.fiveprime_x - node_positions_x[0]
            self.fiveprime_y = self.fiveprime_y - node_positions_y[0]

        node_pos_list_x = [node_positions_x[i]-node_positions_x[0] for i in range(self.N)]
        node_pos_list_y = [node_positions_y[i]-node_positions_y[0] for i in range(self.N)]

        if align:
            if align_mode=='end':
                vec_01_x = node_pos_list_x[self.N-1]
                vec_01_y = node_pos_list_y[self.N-1]
            elif align_mode=='COM':
                vec_01_x = np.mean(node_pos_list_x)
                vec_01_y = np.mean(node_pos_list_y)
            elif isinstance(align_mode, int):
                vec_01_x = node_pos_list_x[align_mode]
                vec_01_y = node_pos_list_y[align_mode]

            curr_angle = np.arctan2(vec_01_y, vec_01_x)

            new_node_pos_list_x = [np.cos(-1*curr_angle)*node_pos_list_x[i] - np.sin(-1*curr_angle) *node_pos_list_y[i] for i in range(self.N)]
            new_node_pos_list_y = [np.sin(-1*curr_angle)*node_pos_list_x[i] + np.cos(-1*curr_angle) *node_pos_list_y[i] for i in range(self.N)]

            if self.fiveprime_x is not None:
                oldx, oldy = copy(self.fiveprime_x), copy(self.fiveprime_y)
                self.fiveprime_x = np.cos(-1*curr_angle)*oldx - np.sin(-1*curr_angle) *oldy
                self.fiveprime_y = np.sin(-1*curr_angle)*oldx + np.cos(-1*curr_angle) *oldy

            if self.threeprime_x is not None:
                oldx, oldy = copy(self.threeprime_x), copy(self.threeprime_y)
                self.threeprime_x = np.cos(-1*curr_angle)*oldx - np.sin(-1*curr_angle) *oldy
                self.threeprime_y = np.sin(-1*curr_angle)*oldx + np.cos(-1*curr_angle) *oldy

            if np.mean(new_node_pos_list_y) < 0:
                new_node_pos_list_y = [-1*y for y in new_node_pos_list_y]
                if self.threeprime_y is not None: self.threeprime_y *= -1
                if self.fiveprime_y is not None: self.fiveprime_y *= -1

            if new_node_pos_list_y[1] < new_node_pos_list_y[0]:
                new_node_pos_list_y = [-1*y for y in new_node_pos_list_y]
                if self.threeprime_y is not None: self.threeprime_y *= -1
                if self.fiveprime_y is not None: self.fiveprime_y *= -1

            node_pos_list_y = new_node_pos_list_y
            node_pos_list_x = new_node_pos_list_x
                

        if self.threeprime_x is not None:
            self.threeprime_x /= bond_width
            self.threeprime_y /= bond_width

        if self.fiveprime_x is not None:
            self.fiveprime_x /= bond_width
            self.fiveprime_y /= bond_width

        # Normalize to bond width
        node_pos_list_x /= bond_width
        node_pos_list_y /= bond_width


        return node_pos_list_x, node_pos_list_y



    def draw(self, label=None, struct_label=None, line=False,
        align=False, align_mode='COM', x0=0, y0=0,
         c = None, cmap='plasma', alpha=None, ax=None,
         vmin=None, vmax=None, fontsize=12):

        '''
        Draw structure using GraphViz layout.

        Input: RGV object or list of objects.

        label (str): list of labels for nucleotides.
        c: color. Can be single matplotlib color letter, string of matplotlib color letters,
                        or list-like object of values (to be used with cmap).
        alpha: Transparency value. can be single value or vector.
        
        ax: axis object to use for drawing.
        figsize: figure size.
        fontsize (int): fontsize for labels.

        Returns: ax object.
        '''
        if struct_label is not None:
            if self.is_multi:
                assert len(struct_label) == len(self.chainbreak)+1
            else:
                raise RuntimeError('struct_label present and not multiple structures.')


        x_coords, y_coords = self.get_coordinates(align=align, align_mode=align_mode)

        #fig = plt.gcf()
        ax = ax or plt.gca()
        ax.set_aspect('equal')

        if isinstance(c,str):
            if len(c) == 1:
                colors=[c]*self.N
            else:
                assert len(c) == self.N
                colors=c

        elif c is None:
            colors=['k']*self.N

        else:
            assert len(c) == self.N
            colormap = plt.get_cmap(cmap) 
            if vmin is None:
                vmin = np.min(c)
            if vmax is None:
                vmax = np.max(c)
            cNorm  = mcolors.Normalize(vmin=vmin, vmax=vmax) #3 for reac
            scalarMap = cm.ScalarMappable(norm=cNorm, cmap=colormap)
            colors = [scalarMap.to_rgba(val) for val in c]

        ax.set_xlim([np.min(x_coords)-5+x0, np.max(x_coords)+5+x0])
        ax.set_ylim([np.min(y_coords)-5+y0, np.max(y_coords)+5+y0])

        # Line mode: (coming)
        # for u,v in G.edges():
            
        #     x = [pos[u][0],pos[v][0]]
        #     y = [pos[u][1],pos[v][1]]
        #     l = Line2D(x,y, linewidth=2, solid_capstyle='round', color=G[u][v]['color'])
        #     ax.add_line(l)

        if alpha is None:
            alpha=[1]*self.N
        elif isinstance(alpha,float):
            alpha = [alpha]*self.N
        elif isinstance(alpha, list):
            assert len(alpha)==self.N

        for i in range(self.N):
            circ = plt.Circle((x_coords[i]+x0, y_coords[i]+y0), radius=1/2, color=colors[i], alpha=alpha[i],linewidth=0)
            ax.add_artist(circ)
            if label:
                plt.text(x_coords[i]+x0, y_coords[i]+y0, label[i],
                 horizontalalignment='center', verticalalignment='center', fontsize=fontsize)
        #plt.scatter(new_node_pos_list_x, new_node_pos_list_y, s=1, c = c, cmap=cmap)

        if self.fiveprime_x is not None:
            plt.text(self.fiveprime_x+x0, self.fiveprime_y+y0, "5'")

        if self.threeprime_x is not None:
            plt.text(self.threeprime_x+x0, self.threeprime_y+y0, "3'")

        if struct_label and self.is_multi:
            label_pos = [0]+self.chainbreak
            for i, text in enumerate(struct_label):
                plt.text(x_coords[label_pos[i]]+x0, y_coords[label_pos[i]]+y0, text)

        ax.axis('off')

        # nx.draw(subgraph, pos, node_size=node_sizes, labels = node_labels, node_color = sns.color_palette('deep',1), 
        #     width=1, arrows=False, solid_capstyle='round',font_size=20, **kwargs)
        return ax

    def count_loops(self):
        n_1, n_2, n_3, n_4, n_5 = 0,0,0,0,0
        nodes = [n for n in list(self.G.nodes) if not isinstance(n,str)]
        subgraph = self.G.subgraph(nodes)
        subgraph = subgraph.to_undirected()

        for x in list(subgraph.degree):
            if x[1]==1:
                n_1 +=1
            elif x[1] == 2:
                n_2 += 1
            elif x[1] == 3:
                n_3 += 1
            elif x[1] == 4:
                n_4 += 1
            elif x[1] > 4:
                n_5 += 1
        return n_1 - 1, n_2, n_3, n_4, n_5 #subtract off 1 to not count exterior loop as hairpin

    def get_info(self):
        print("Max ladder distance: %d" % self.MLD)
        print("n_hairpins: %d" % self.n_hairpins)
        print("n_internal_loops: %d" % self.n_internal_loops)
        print("n_3WJs: %d" % self.n_3WJs)
        print("n_4WJs: %d" % self.n_4WJs)
        print("n_5WJs_up: %d" % self.n_5WJs_up)

