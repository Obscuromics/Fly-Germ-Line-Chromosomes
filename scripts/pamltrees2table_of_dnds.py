#!/usr/bin/env python3

import os
from collections import defaultdict
# import numpy as np
from Bio import Phylo # type: ignore
import sys

# species in the analysis
tree_path = 'data/grc_dnds/newicks/'
outgroup = 'dmel'
sciarid_grcs = ['bcop_grc', 'bimp_grc', 'ling_grc']
sciarid_core = ['bcop_core', 'bimp_core', 'ling_core', 'phyg']
cecidomyiidae = ['aaphi', 'contarinia', 'orobi']
sciaridae = set(['bcop', 'ling', 'phyg','bimp'])

def is_just_grcs(clade):
    return(all(['grc' in tip.name for tip in clade.get_terminals()]))

def is_monophyletic_core_sciaridae(clade):
    return(all([any([tip.name.startswith(sci) for sci in sciarid_core]) for tip in clade.get_terminals()]))
    # return(all([tip.name.split('_')[0] in sciaridae for tip in clade]))

def is_monophyletic_cecidomyiidae(clade):
    return(all([tip.name.split('.')[0] in cecidomyiidae or 'grc' in tip.name for tip in clade.get_terminals()]))

def tip2sp_name(tip):
    return(tip.name.split('.')[0])

def tabulate_branches(og, outtable):
    dnfile = tree_path + og + '.dn.nwk'
    dsfile = tree_path + og + '.ds.nwk'
    dndsfile = tree_path + og + '.dnds.nwk'

    dntree = Phylo.read(dnfile, "newick")
    dstree = Phylo.read(dsfile, "newick")
    dndstree = Phylo.read(dndsfile, "newick")

    # root_by_dmel(dstree)
    # root_by_dmel(dntree)
    # root_by_dmel(dndstree)

    dsnodes = dstree.find_clades()
    dnnodes = list(dntree.find_clades())
    dndsnodes = list(dndstree.find_clades())

    # they ahve all the same index, so all I need to do is test individual nodes and print the table
    for idx, clade in enumerate(dsnodes):
        if clade.name:
            clade.name = tip2sp_name(clade)
        else:
            clade.name = str(idx)
        asn = 'other'
        if is_monophyletic_cecidomyiidae(clade):
            asn = 'ceci' # this is GRC + cecidomyiidae

        if is_just_grcs(clade):
            asn = 'GRC' # if it's just GRCs, override ceci assignment
        elif is_monophyletic_core_sciaridae(clade):
            asn = 'sci_core' # this is core sciaridae only

        type = 'branch'
        if clade.is_terminal():
            type = 'tip'
        
        if dndsnodes[idx].branch_length != None and dnnodes[idx].branch_length != None and clade.branch_length != None:
            outtable.write("\t".join([og, asn, type, str(dnnodes[idx].branch_length), str(clade.branch_length), str(dndsnodes[idx].branch_length), clade.name]) + '\n')
    return 0

# this was just for sanity checking if the dnds tree branch lengths are indeed ds/dn
# import numpy as np
# def get_banchlengths(tree):
#     bl = []
#     for idx, clade in enumerate(tree.find_clades()):
#         if clade.name:
#             clade.name = "%d_%s" % (idx, clade.name)
#         else:
#             clade.name = str(idx)
#         bl.append(clade.branch_length)
#     return bl
# 
# og = ogs_to_process[1000]
# og = "OG0008831"
# 
# dnfile = tree_path + og + '.dn.nwk'
# dsfile = tree_path + og + '.ds.nwk'
# dndsfile = tree_path + og + '.dnds.nwk'
# 
# dntree = Phylo.read(dnfile, "newick")
# dstree = Phylo.read(dsfile, "newick")
# dndstree = Phylo.read(dndsfile, "newick")
# 
# Phylo.draw(dntree)
# Phylo.draw(dstree)
# Phylo.draw(dndstree)
# 
# np.array(get_banchlengths(dntree), dtype=float)/ np.array(get_banchlengths(dstree), dtype=float) 
# np.array(get_banchlengths(dndstree), dtype=float)

# --------- this is in case we want to root the treees
# def root_by_dmel(tree):
#     dmel_tips = [tip.name for tip in tree.get_terminals() if "dmel" in tip.name]
#     return(tree.root_with_outgroup(dmel_tips))

# with open('tables/L-busco-phylogenies-summary.tsv', 'w') as tab:

all_files = os.listdir(tree_path)

ogs_to_process = list(set([f.split('.')[0] for f in all_files]))

with open('tables/DnDs_per_branch_summary.tsv', 'w') as tab:
    tab.write('orthogroup\tasn\ttype\tdn\tds\tdnds\tspecies\n')
    for og in ogs_to_process:
        #sys.stdout.write(file)
        #sys.stdout.write("\n")
        # file = 'OG0003003_tree.txt'
        # print(file)
        
        try:
            tabulate_branches(og, tab)
            sys.stderr.write(og + ": Done\n")
            # sys.stdout.write(gene + '\t' + tree2assigments(input_newick) + '\n')
        except IndexError:
            sys.stderr.write(og + ": failed to be processed")


