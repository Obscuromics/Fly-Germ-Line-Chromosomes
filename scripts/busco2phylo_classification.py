from os import listdir
from Bio import Phylo
import sys

# this should be probably argparse object (package for named arguements), for now making it quick and dirty using sys
input_dir = sys.argv[1]
# input_dir = 'data/testing_trees_busco'

output_pattern = sys.argv[2] # will generate <output_pattern>_per_tree_summary.tsv and <output_pattern>_per_gene_summary.tsv
# output_pattern = 'data/busco'

meta_information_table_filename = sys.argv[3]

######### global constants
outgroup = set()
sciaridae = set()
cecidomyiidae = set()

with open(meta_information_table_filename, 'r') as meta_information_table:
    for line in meta_information_table:
        sp, family, is_outgroup = line.rstrip('\n').split('\t')
        if is_outgroup == 'yes':
            outgroup.add(sp)
        if family == 'Cecidomyiidae':
            cecidomyiidae.add(sp)
        if family == 'Sciaridae':
            sciaridae.add(sp)

sys.stderr.write("Outgroup: " + "\t".join(outgroup) + '\n')
sys.stderr.write("Sciaridae: " + "\t".join(sciaridae) + '\n')
sys.stderr.write("Cecidomyiidae: " + "\t".join(cecidomyiidae) + '\n')

# outgroup = 'Sylvicola cinctus'
# sciaridae = set(['Bradysia odoriphaga', 'Pseudolycoriella hygida', 'Phytosciara flavipes', 'Trichosia splendens', 'Bradysia coprophila', 'Bradysia impatiens', 'Lycoriella ingenua'])
# cecidomyiidae = set(['Sitodiplosis mosellana', 'Mayetiola destructor', 'Obolodiplosis robiniae', 'Catotricha subobsoleta', 'Lestremia cinerea', 'Aphidoletes aphidimyza', 'Contarinia nasturtii', 'Resseliella maxima'])

############ FUNCTIONS

def is_grc(tip):
    return(tip.name.startswith('L-'))

def is_just_grcs(clade):
    return(all([is_grc(tip) for tip in clade.get_terminals()]))

def is_monophyletic_sciaridae(clade):
    return(all([tip2species_name(tip) in sciaridae for tip in clade.get_terminals()]))

def is_monophyletic_cecidomyiidae(clade):
    return(all([tip2species_name(tip) in cecidomyiidae or is_grc(tip) for tip in clade.get_terminals()]))

def count_grcs(clade):
    return(sum([is_grc(tip) for tip in clade.get_terminals()]))

def tip2species_name(tip):
    sp_string = ' '.join(tip.name.split('_')[0:2])
    if '-' in sp_string:
        return(sp_string.split('-')[1])
    else:
        return(sp_string)

def tip2short_name(tip):
    name_in_sections = tip2species_name(tip).split(' ')
    return(name_in_sections[0][0] + name_in_sections[1][0:3])

def grc_tip2sp_scf_loc(tip):
    sp = tip2short_name(tip)
    if not "|" in tip.name:
        tipname_split = tip.name.split('_')
        chrom = tipname_split[2]
        loc = tipname_split[3]
    else:
        tipname_split = tip.name.split('|')[1].split('_')
        chrom = '_'.join(tipname_split[0:2])
        loc = tipname_split[2]
    return(sp, chrom, loc)

######### file manipulation bit

tree_files = [i for i in listdir(input_dir) if i.endswith('treefile')]

# tree_summary = 'summary_all_trees.tsv'
# GRC_genes_summary = 'phylogenetic_placement_of_all_grc_genes.tsv'

# tree_file = tree_files[1]
# trees = [Phylo.read(input_dir + '/' + tree_file, "newick") for tree_file in tree_files]
# tree = trees[1]
# Phylo.draw(tree)

######### Per TREE summary
tree_summary_filename = output_pattern + '_per_tree_summary.tsv'
gene_summary_filename = output_pattern + '_per_gene_summary.tsv'

with open(tree_summary_filename, 'w') as tree_summary, open(gene_summary_filename, 'w') as gene_summary:
    sys.stderr.write('Running the classification analysis of ...\n')

    ### headers
    row_to_print = 'BUSCO_id\ttotal_genes\tmonophyletic_sci\tmonophyletic_ceci\tGRC_species\tGRCs_total\tGRCs_sci\tGRCs_ceci\n)'
    tree_summary.write(row_to_print)
    row_to_print = 'BUSCO_id\tsp\tchromosome\tlocation\tGRC_closest_relative\tGRC_bootstrap\tnon-GRC_closest_relative\tnon-GRC_bootstrap\tbranch_length\n)'
    gene_summary.write(row_to_print)

    for tree_file in tree_files:
        BUSCO_id = tree_file.rstrip('.treefile') ## 1
        sys.stderr.write('\t' + BUSCO_id + '\n')
        input_newick = input_dir + '/' + tree_file
        tree = Phylo.read(input_newick, "newick")

        ### rooting of the tree
        outgroup_tips = [tip for tip in tree.get_terminals() if tip2species_name(tip) in outgroup] # subselects all the terminal branches that match the specified outgroup in the config file
        if outgroup_tips: # if there is the outgroup
            tree.root_with_outgroup(outgroup_tips) # root the tree by that outgroup
        
        ### preprocessing
        present_sciaridae = set()
        present_cecidomyiidae = set()
        present_GRCs = set()
        total_tips = str(len(tree.get_terminals())) ## 2
        for tip in tree.get_terminals():
            sp_name = tip2species_name(tip)
            if is_grc(tip):
                present_GRCs.add(tip)
            elif sp_name in sciaridae:
                present_sciaridae.add(tip)
            if sp_name in cecidomyiidae:
                present_cecidomyiidae.add(tip)
                
        ######### Per Gene summary
        sciaridae_ancestor = tree.common_ancestor(present_sciaridae)
        cecidomyiidae_ancestor = tree.common_ancestor(present_cecidomyiidae)

        monophyletic_sci = is_monophyletic_sciaridae(sciaridae_ancestor) ## 3
        monophyletic_ceci = is_monophyletic_cecidomyiidae(cecidomyiidae_ancestor) ## 4

        GRCs_total = str(len(present_GRCs))
        GRC_species = ','.join(set([tip2short_name(tip) for tip in present_GRCs]))
        GRCs_sci = str(count_grcs(sciaridae_ancestor)) if monophyletic_sci else 'NA'
        GRCs_ceci = str(count_grcs(cecidomyiidae_ancestor)) if monophyletic_ceci else 'NA'

        row_to_print = '\t'.join([BUSCO_id, total_tips, str(monophyletic_sci), str(monophyletic_ceci), GRC_species, GRCs_total, GRCs_sci, GRCs_ceci]) + '\n'
        tree_summary.write(row_to_print)

        ######### Per gene summary
        for grc in present_GRCs:
            sys.stderr.write('\t\t' + grc.name + '\n')
            sp, chrom, loc = grc_tip2sp_scf_loc(grc) ### 2,3,4 of the reported values: species, chromsome, location
            grc2root = list(reversed(tree.get_path(grc)[:-1])) # I am storing the whole path from root to the investigated termial branch (excluding that branch, thjat's [:-1] which removes the last value)
            non_grc_bootstraps = []
            grc_bootstraps = []

            for target_node in grc2root: # target_node is is the last common ancestor (node) of the GRC gene and any non-GRC gene
                if target_node.name: 
                    non_grc_bootstraps.append(target_node.name) # this will record the whole line of bootstraps
                if not is_just_grcs(target_node): # breaks on the fist node with non-GRC ancestors OR at the root
                    break
            non_grc_ancestor = target_node 
            closest_non_grc_relatives = ';'.join([tip2short_name(tip) for tip in non_grc_ancestor.get_terminals() if not is_grc(tip)]) ### 7
            if non_grc_bootstraps:
                closest_non_grc_relatives_bs = ';'.join(non_grc_bootstraps) ### 8
            else:
                closest_non_grc_relatives_bs = 'NA'

            for target_node in grc2root: # target_node is is the last common ancestor (node) of the GRC gene and any other GRC gene
                if target_node.name:
                    grc_bootstraps.append(target_node.name) # this will record the whole line of bootstraps, necesarily will overlap with the others
                if [tip for tip in target_node.get_terminals() if is_grc(tip) and tip != grc]: # breaks on the fist node with GRC ancestors OR at the root
                    break
            grc_ancestor = target_node 
            closest_grc_relatives = ';'.join([tip2short_name(tip) for tip in grc_ancestor.get_terminals() if is_grc(tip) and not tip == grc]) ### 5
            if non_grc_bootstraps:
                closest_grc_relatives_bs = ';'.join(grc_bootstraps) ### 6
            else:
                closest_grc_relatives_bs = 'NA'
            
            row_to_print = '\t'.join([BUSCO_id, sp, chrom, loc, closest_grc_relatives, closest_grc_relatives_bs, closest_non_grc_relatives, closest_non_grc_relatives_bs, str(grc.branch_length)]) + '\n'
            gene_summary.write(row_to_print)

            # Phylo.write(non_grc_ancestor, sys.stdout, "newick") # it get's messy with Newick

