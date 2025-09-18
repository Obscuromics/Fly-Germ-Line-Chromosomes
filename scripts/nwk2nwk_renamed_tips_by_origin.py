from collections import defaultdict
from os import listdir
from os import path
from os import makedirs
from Bio import Phylo
# from sys import argv
from sys import stderr

input_dir = argv[1]
# input_dir = 'data/testing_trees_busco/'
ouput_dir = argv[2]
# ouput_dir = 'data/testing_trees_busco_renamed/'
busco_table_filename = argv[3]
# busco_table_filename = 'tables/busco_grc_classification_odb12_diptera_v2.tsv' #'tables/busco_grc_classification_odb12_diptera_v2.tsv' 


# PROBLEM-NA-Bradysia_impatiens:0.03838

if not path.exists(ouput_dir):
    stderr.write("Creating:" + ouput_dir + '\n')
    makedirs(ouput_dir)

origin_dict = defaultdict(lambda: 'NA')
which_GRC_dict = defaultdict(lambda: 'PROBLEM')

with open(busco_table_filename) as f:
    for line in f:
        BUSCO_line = line.strip().split()
        if BUSCO_line[0] == 'BUSCO_id':
            continue
        key = BUSCO_line[2] + BUSCO_line[0] + BUSCO_line[4]
        # sp_busco e.g. Ling_40951at7147
        which_GRC_dict[key] = BUSCO_line[3]
        if BUSCO_line[5] == "Cecidomyiidae":
            origin_dict[key] = 'c'
        else:
            origin_dict[key] = 's'

stderr.write("Extracted " + str(len(origin_dict.keys())) + ' busco + species + location -> origin records\n')

tree_files = [i for i in listdir(input_dir) if i.endswith('treefile')]
target_species = ['Bradysia_coprophila', 'Bradysia_impatiens', 'Lycoriella_ingenua']

stderr.write("Found: " + str(len(tree_files)) + ' trees in ' + input_dir + '\n')
stderr.write('Output trees will be written in ' + ouput_dir + '\n')

for tree_file in tree_files:
    busco = tree_file.rstrip('.treefile')
    input_newick = input_dir + '/' + tree_file
    output_newick = ouput_dir + '/' + tree_file
    tree = Phylo.read(input_newick, "newick")
    for i in range(len(tree.get_terminals())):
        terminal_orig_name = tree.get_terminals()[i].name
        species = "_".join(terminal_orig_name.split('_')[:2])
        if '-' in species:
            species = species.split('-')[1]
        if species in target_species:
            if terminal_orig_name[0] == 'A':
                new_tip_name = 'A-' + species
            else:
                sp_short = species.split("_")[0][0] + species.split("_")[1][0:3]
                location = terminal_orig_name.split('|')[1].split('_')[2]
                dict_key = sp_short + busco + location
                new_tip_name = which_GRC_dict[dict_key] + '-' + origin_dict[dict_key] + '-' + species
        else:
            new_tip_name = species
        tree.get_terminals()[i].name = new_tip_name
    Phylo.write(tree, output_newick, "newick")

stderr.write('All done\n')