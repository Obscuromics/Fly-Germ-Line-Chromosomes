#!/usr/bin/env python3

import os
from collections import defaultdict
# import numpy as np
from Bio import Phylo
import sys

# species in the analysis
outgroup = 'dmel'
# sciarid_grcs = ['bcop_grc', 'bimp_grc', 'ling_grc']
# sciarid_core = ['bcop_core', 'bimp_core', 'ling_core', 'phyg']
cecidomyiidae = ['aaphi', 'contarinia', 'orobi']
sciaridae = set(['bcop', 'ling', 'phyg','bimp'])

def is_just_grcs(clade):
    return(all(['grc' in tip.name for tip in clade.get_terminals()]))

def is_monophyletic_sciaridae(clade):
    return(all([tip.name.split('_')[0] in sciaridae for tip in clade]))

def is_monophyletic_cecidomyiidae(clade):
    return(all([tip.name.split('_')[0] in cecidomyiidae or 'grc' in tip.name for tip in clade]))

def tip2sp_name(tip):
    return(tip.name.split('_')[0])

def tip2gene(tip):
    return("_".join(tip.name.split('_')[-1:]).rstrip("'"))

def print_assignment_group(tree_name, subtree, assignment):
    for tip in subtree:
        if 'grc' in tip.name:
            sp = tip2sp_name(tip)
            gene = tip2gene(tip)
            branch_length = str(tip.branch_length)
            sys.stdout.write("\t".join([tree_name, sp, gene, assignment, branch_length]) + '\n')

def print_no_assignment_orthogroup(tree_name):
    sys.stdout.write("\t".join([tree_name, 'NA', 'NA', 'other', 'NA']) + '\n')            

def find_family_splitting_clade(node_to_test, tree_name):
    if is_just_grcs(node_to_test.clades[0]): # if all in one branch are just GRC genes
        print_assignment_group(tree_name, node_to_test.clades[0], 'other') # those will be assinged as 'other'
        return(find_family_splitting_clade(node_to_test.clades[1], tree_name)) # and we continue the search for the node that splits families
    if is_just_grcs(node_to_test.clades[1]): # if all in the other branch are just GRC genes ...
        print_assignment_group(tree_name, node_to_test.clades[1], 'other')
        return(find_family_splitting_clade(node_to_test.clades[0], tree_name))
    return(node_to_test) # if not, then the split is separating families (or at least family members)

def tree2assigments(input_newick):
    tree_name = input_newick.split('/')[-1:][0].split('_')[0]
    tree = Phylo.read(input_newick, "newick")

    sp2node_name = defaultdict(list) # it's not really sp2tip anymore... instead clade2tip
    for tip in tree.get_terminals():
        sp = tip2sp_name(tip)
        is_grc = 'grc' in tip.name
        if sp in sciaridae and not is_grc:
            sp2node_name['sciaridae'].append(tip) # this will be a list of all non-grc sciarids
        elif sp in cecidomyiidae and not is_grc:
            sp2node_name['cecidomyiidae'].append(tip) # this will be a list of all cecidomyiidae species 
        else:
            sp2node_name[sp].append(tip) # this is just Dmel and grcs, I don't use GRCs anywhere at the moment (I could)

    if not sp2node_name[outgroup]:
        sys.stderr.write(tree_name + ': Outgroup (' + outgroup + ') absent\n')
        sys.stdout.write("\t".join([tree_name, 'NA', 'NA', 'other', 'NA']) + '\n')
        return(0)

    outgroup_node = sp2node_name[outgroup][0]
    tree.root_with_outgroup(outgroup_node)
    
    if tree.root.clades[0] != outgroup_node:
        gnat_subtree = tree.root.clades[0]
    else:
        gnat_subtree = tree.root.clades[1]

    gnat_subtree = find_family_splitting_clade(gnat_subtree, tree_name)

    gnats1 = gnat_subtree.clades[0].get_terminals()
    gnats2 = gnat_subtree.clades[1].get_terminals()
        
    if is_monophyletic_sciaridae(gnats1):
        sciaridae_subtree = gnats1
    elif is_monophyletic_cecidomyiidae(gnats1):
        cecidomyiidae_subtree = gnats1
    else:
        sys.stderr.write(tree_name + ': families are not monophyletic\n')
        print_no_assignment_orthogroup(tree_name)
        return(0)

    if is_monophyletic_sciaridae(gnats2):
        sciaridae_subtree = gnats2
    elif is_monophyletic_cecidomyiidae(gnats2):
        cecidomyiidae_subtree = gnats2
    else:
        sys.stderr.write(tree_name + ': families are not monophyletic\n')
        print_no_assignment_orthogroup(tree_name)
        return(0)

    print_assignment_group(tree_name, sciaridae_subtree, 'sciaridae')
    print_assignment_group(tree_name, cecidomyiidae_subtree, 'cecidomyiidae')

    return(0)


# input_dir = sys.argv[1]
input_dir = 'data/testing_trees'
# tree_files = [i for i in os.listdir(input_dir) if i.endswith('treefile')]
tree_files = os.listdir('data/testing_trees/')

# with open('tables/L-busco-phylogenies-summary.tsv', 'w') as tab:
sys.stdout.write('orthogroup\tspecies\tgene\tclassification\tbranch_lengths\n')

for file in tree_files:
    # file = 'OG0003003_tree.txt'
    # print(file)
    input_newick = input_dir + '/' + file
    tree2assigments(input_newick)
    # sys.stdout.write(gene + '\t' + tree2assigments(input_newick) + '\n')


