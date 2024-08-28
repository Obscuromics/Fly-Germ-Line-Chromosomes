import os
from collections import defaultdict
from Bio import Phylo
import sys

######### global constants
# these will be picked up from a file
outgroup = 'Sylvicola cinctus'
sciaridae = set(['Bradysia odoriphaga', 'Pseudolycoriella hygida', 'Phytosciara flavipes', 'Trichosia splendens', 'Bradysia coprophila', 'Bradysia impatiens', 'Lycoriella ingenua'])
cecidomyiidae = set(['Sitodiplosis mosellana', 'Mayetiola destructor', 'Obolodiplosis robiniae', 'Catotricha subobsoleta', 'Lestremia cinerea', 'Aphidoletes aphidimyza', 'Contarinia nasturtii', 'Resseliella maxima'])

############ FUNCTIONS

def is_just_grcs(clade):
    return(all([tip.name.startswith('L-') for tip in clade.get_terminals()]))

def is_monophyletic_sciaridae(clade):
    return(all([tip2species_name(tip) in sciaridae for tip in clade.get_terminals()]))

def is_monophyletic_cecidomyiidae(clade):
    return(all([tip2species_name(tip) in cecidomyiidae or tip.name.startswith('L-') for tip in clade.get_terminals()]))

def count_grcs(clade):
    return(sum([tip.name.startswith('L-') for tip in clade.get_terminals()]))

def tip2species_name(tip):
    sp_string = ' '.join(tip.name.split('_')[0:2])
    if '-' in sp_string:
        return(sp_string.split('-')[1])
    else:
        return(sp_string)

def tip2short_name(tip):
    name_in_sections = tip2species_name(tip).split(' ')
    return(name_in_sections[0][0] + name_in_sections[1][0:3])

######### file manipulation bit
# input_dir = sys.argv[1]
input_dir = 'data/testing_trees_busco'
tree_files = [i for i in os.listdir(input_dir) if i.endswith('treefile')]

# tree_summary = 'summary_all_trees.tsv'
# GRC_genes_summary = 'phylogenetic_placement_of_all_grc_genes.tsv'

# tree_file = tree_files[0]
# trees = [Phylo.read(input_dir + '/' + tree_file, "newick") for tree_file in tree_files]
# tree = trees[0]
# Phylo.draw(tree)

row_to_print = 'BUSCO_id\ttotal_genes\tmonophyletic_sci\tmonophyletic_ceci\tGRC_species\tGRCs_total\tGRCs_sci\tGRCs_ceci\n)'
sys.stdout.write(row_to_print)

for tree_file in tree_files:
    BUSCO_id = tree_file.rstrip('.treefile') ## 1
    input_newick = input_dir + '/' + tree_file
    tree = Phylo.read(input_newick, "newick")

    ### preprocessing
    # sp2tips = defaultdict(list) # I don't think I need this in the end
    present_sciaridae = set()
    present_cecidomyiidae = set()
    present_GRCs = set()
    total_tips = str(len(tree.get_terminals())) ## 2
    for tip in tree.get_terminals():
        sp_name = tip2species_name(tip)
        # sp2tips[sp_name].append(tip)
        if tip.name.startswith('L-'):
            present_GRCs.add(tip)
        elif sp_name in sciaridae:
            present_sciaridae.add(tip)
        if sp_name in cecidomyiidae:
            present_cecidomyiidae.add(tip)
        
    ### testing
    monophyletic_sci = is_monophyletic_sciaridae(tree.common_ancestor(present_sciaridae)) ## 3
    monophyletic_ceci = is_monophyletic_cecidomyiidae(tree.common_ancestor(present_cecidomyiidae)) ## 4

    GRCs_total = str(len(present_GRCs))
    GRC_species = ','.join(set([tip2short_name(tip) for tip in present_GRCs]))
    GRCs_sci = str(count_grcs(tree.common_ancestor(present_sciaridae))) if monophyletic_sci else 'NA'
    GRCs_ceci = str(count_grcs(tree.common_ancestor(present_cecidomyiidae))) if monophyletic_ceci else 'NA'

    row_to_print = '\t'.join([BUSCO_id, total_tips, str(monophyletic_sci), str(monophyletic_ceci), GRC_species, GRCs_total, GRCs_sci, GRCs_ceci]) + '\n'
    sys.stdout.write(row_to_print)


# cecido_branch = tree.common_ancestor(sp2tips['Mayetiola destructor'] + sp2tips['Resseliella maxima'])
# is_just_grcs(cecido_branch)
# is_monophyletic_sciaridae(cecido_branch)
# is_monophyletic_cecidomyiidae(cecido_branch)
# Ok, all the code seems to behave in the way I would expect it to