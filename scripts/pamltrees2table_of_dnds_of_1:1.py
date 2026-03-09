#!/usr/bin/env python3

import os
# from collections import defaultdict
# import numpy as np
from Bio import Phylo # type: ignore
import matplotlib.pyplot as plt
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
    

def is_monophyletic_sciaridae(clade):
    return(all([tip.name.split('_')[0] in sciaridae for tip in clade.get_terminals()]))


def is_monophyletic_cecidomyiidae(clade):
    return(all([tip.name.split('.')[0] in cecidomyiidae or 'grc' in tip.name for tip in clade.get_terminals()]))

def tip2sp_name(tip):
    return(tip.name.split('.')[0])

def tip2gene_name(tip):
    return(tip.name.split('.')[1])

def plot_trees(tree1, tree2, tree3, filename=None, titles=None):

    trees = [tree1, tree2, tree3]

    if titles is None:
        titles = ["Tree 1", "Tree 2", "Tree 3"]

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    for ax, tree, title in zip(axes, trees, titles):

        def branch_label(clade):
            if clade.branch_length is not None:
                return f"{clade.branch_length:.3f}"
            return None

        Phylo.draw(
            tree,
            axes=ax,
            do_show=False,
            branch_labels=branch_label
        )

        ax.set_title(title)

    # plt.tight_layout()

    if filename is not None:
        plt.savefig(filename, bbox_inches="tight", dpi=300)

    # plt.show()

def tabulate_branches(og, outtable):
    dnfile = tree_path + og + '.dn.nwk'
    dsfile = tree_path + og + '.ds.nwk'
    dndsfile = tree_path + og + '.dnds.nwk'

    dntree = Phylo.read(dnfile, "newick")
    dstree = Phylo.read(dsfile, "newick")
    dndstree = Phylo.read(dndsfile, "newick")

    dmel_tips = [tip.name for tip in dntree.get_terminals() if "dmel" in tip.name]
    if len(dmel_tips) == 0:
        sys.stderr.write("No dmel tip found in " + og + '\n')
        return 0
    if len(dmel_tips) == 1:
        dntree.root_with_outgroup(dmel_tips[0])
        dmel_tips = [tip.name for tip in dstree.get_terminals() if "dmel" in tip.name]
        dstree.root_with_outgroup(dmel_tips[0])
        dmel_tips = [tip.name for tip in dndstree.get_terminals() if "dmel" in tip.name]
        dndstree.root_with_outgroup(dmel_tips[0])
    else: 
        mrca = dntree.common_ancestor(dmel_tips)
        dntree.root_with_outgroup(mrca)
        dmel_tips = [tip.name for tip in dstree.get_terminals() if "dmel" in tip.name]
        mrca = dstree.common_ancestor(dmel_tips)
        dstree.root_with_outgroup(mrca)
        dmel_tips = [tip.name for tip in dndstree.get_terminals() if "dmel" in tip.name]
        mrca = dndstree.common_ancestor(dmel_tips)
        dndstree.root_with_outgroup(mrca)

    # plot_trees(dntree, dstree, dndstree, filename='figures/DnDs_trees/python/' + og + '_trees.png', titles=['dN tree', 'dS tree', 'dN/dS tree'])

    # dsnodes = dstree.find_clades()
    # dnnodes = list(dntree.find_clades())
    # dndsnodes = list(dndstree.find_clades())

    all_tips = dndstree.get_terminals()
    tip_sp = [tip.name.split('.')[0] for tip in all_tips]
    
    # all_sciaridae_tips = [tip for tip in all_tips if tip.name.split('.')[0][0:4] in sciaridae]
    all_sciaridae_core_tips =  [tip for tip in all_tips if any([tip.name.startswith(sci) for sci in sciarid_core])]

    sciaridae_ancestor = dndstree.common_ancestor(all_sciaridae_core_tips) #all_sciaridae_tips

    # tip_to_process = []
    for sp in ['bcop', 'bimp', 'ling']:
        # sys.stderr.write('\t processing ... ' + sp + '\n')
        core_tips = [tip_name.startswith(sp + "_core") for tip_name in tip_sp]
        grc_tips = [tip_name.startswith(sp + "_grc") for tip_name in tip_sp]
        if sum(core_tips) == 1 and sum(grc_tips) == 1: # if there is one core and one grc tip, then we can process this tree
            core_tip = [tip for idx, tip in enumerate(all_tips) if core_tips[idx]][0]
            grc_tip = [tip for idx, tip in enumerate(all_tips) if grc_tips[idx]][0]
        else:
            continue # if not, we skip species in this tree, because we can't be sure which tips to process

        # now I have both GRC and core tips to process        
        core_gene = tip2gene_name(core_tip)
        dn_branch_core = 0
        ds_branch_core = 0
        dnds_branch_core = 0
        dn_tip = dntree.find_any(core_tip.name)
        if dn_tip == None:
            sys.stderr.write("Malformated dn tree when processing " + core_tip.name + " in " + og + '\n')
            break
        ds_tip = dstree.find_any(core_tip.name)
        if ds_tip == None:
            sys.stderr.write("Malformated ds tree when processing " + core_tip.name + " in " + og + '\n')
            break
        dn_path = list(reversed(dntree.get_path(dn_tip)))
        ds_path = list(reversed(dstree.get_path(ds_tip)))
        for idx, clade in enumerate(reversed(dndstree.get_path(core_tip))):
            dn_node = dn_path[idx]
            ds_node = ds_path[idx]
            if clade == sciaridae_ancestor:
                sys.stderr.write("Cutting off ancestral Sciaridae node in " + og + '\n')
                break
            if is_monophyletic_sciaridae(clade):
                dnds_branch_core += clade.branch_length if clade.branch_length else 0
                dn_branch_core += dn_node.branch_length if dn_node.branch_length else 0
                ds_branch_core += ds_node.branch_length if ds_node.branch_length else 0
            else:
                break

        grc_gene = tip2gene_name(grc_tip)
        dn_branch_grc = 0
        ds_branch_grc = 0
        dnds_branch_grc = 0
        dn_path = list(reversed(dntree.get_path(dntree.find_any(grc_tip.name))))
        ds_path = list(reversed(dstree.get_path(dstree.find_any(grc_tip.name))))
        for idx, clade in enumerate(reversed(dndstree.get_path(grc_tip))):
            dn_node = dn_path[idx]
            ds_node = ds_path[idx]
            if is_just_grcs(clade):
                dnds_branch_grc += clade.branch_length if clade.branch_length else 0
                dn_branch_grc += dn_node.branch_length if dn_node.branch_length else 0
                ds_branch_grc += ds_node.branch_length if ds_node.branch_length else 0
            else:
                break
        
        if dn_branch_core != 0 and ds_branch_core != 0 and dn_branch_grc != 0 and ds_branch_grc != 0:
            outtable.write("\t".join([og, sp, core_gene, str(dn_branch_core), str(ds_branch_core), str(dnds_branch_core), grc_gene, str(dn_branch_grc), str(ds_branch_grc), str(dnds_branch_grc)]) + '\n')

    return 0

all_files = os.listdir(tree_path)
# all_plots = os.listdir('figures/DnDs_trees/python/')

ogs_to_process = list(set([f.split('.')[0] for f in all_files]))
# ogs_to_skip = set([f.split('.')[0] for f in all_plots])

with open('tables/DnDs_1to1_grc_core_orthologs.tsv', 'w') as tab:
    tab.write('orthogroup\tspecies\tcore_gene\tcore_dn\tcore_ds\tcore_dnds\tgrc_gene\tgrc_dn\tgrc_ds\tgrc_dnds\n')
    for og in ogs_to_process:
        # sys.stdout.write(og +'\n')
        # if og in ogs_to_skip:
        #     sys.stderr.write(og + ": already processed, skipping\n")
        #     continue
        try:
            tabulate_branches(og, tab)
            sys.stderr.write(og + ": Done\n")
            # sys.stdout.write(gene + '\t' + tree2assigments(input_newick) + '\n')
        except IndexError:
            sys.stderr.write(og + ": failed to be processed")


