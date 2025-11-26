

dnds_table <- read.table('tables/DnDs_per_branch_summary.tsv', header=TRUE, sep='\t')
dnds_table$log10dnds <- log10(dnds_table$dnds)

pal = c('dmel' = '#f4f26d92', 'sci_core' = '#12ebaa57', 'sci_core_branch' = '#8cd1bc97',  
        'ceci' = '#f4263a67', 'ceci_branch' = '#e28a9367', 'GRC' = '#ffbb0096',  'GRC_branch' = '#e9ca77a3')

# plot(log10(dnds_table[, 'ds']), log10(dnds_table[, 'dn']), xlab = 'log10(ds)', ylab = 'log10(dn)', main = 'dn vs ds per branch')
# plot(log10(dnds_table[, 'ds']), log10(dnds_table[, 'dnds']), xlab = 'log10(ds)', ylab = 'log10(dn/ds)', main = 'dN/dS vs dS per branch')

png('figures/paml_dnds_vs_ds_per_tip_cores.png', width=800, height=600)
    plot(log10(dnds_table[, 'ds']), log10(dnds_table[, 'dnds']), xlab = 'log10(ds)', ylab = 'log10(dnds)', main = 'dnds vs ds per tip')
    points(log10(dnds_table[dnds_table$species == 'dmel', 'ds']), log10(dnds_table[dnds_table$species == 'dmel', 'dnds']), col=pal['dmel'], pch=19)
    points(log10(dnds_table[dnds_table$asn == 'sci_core' & dnds_table$type == 'tip', 'ds']), log10(dnds_table[dnds_table$asn == 'sci_core' & dnds_table$type == 'tip', 'dnds']), col=pal['sci_core'], pch=19)
    points(log10(dnds_table[dnds_table$asn == 'ceci' & dnds_table$type == 'tip', 'ds']), log10(dnds_table[dnds_table$asn == 'ceci' & dnds_table$type == 'tip', 'dnds']), col=pal['ceci'], pch=19)
    legend('topright', legend = c('all tips', 'D. melanogaster tips', 'Sciaridae core tips', 'Cecidomyiidae tips'), 
        col = c('black', pal['dmel'], pal['sci_core'], pal['ceci']), pch=19, bty = 'n')
dev.off()

png('figures/paml_dnds_vs_ds_per_branch_Sci_GRC.png', width=800, height=600)
    plot(log10(dnds_table[dnds_table$type == 'branch', 'ds']), log10(dnds_table[dnds_table$type == 'branch', 'dnds']), xlab = 'log10(ds)', ylab = 'log10(dnds)', main = 'dnds vs ds per branch')
    points(log10(dnds_table[dnds_table$asn == 'sci_core' & dnds_table$type == 'branch', 'ds']), log10(dnds_table[dnds_table$asn == 'sci_core' & dnds_table$type == 'branch', 'dnds']), col=pal['sci_core_branch'], pch=19)
    points(log10(dnds_table[dnds_table$asn == 'GRC' & dnds_table$type == 'branch', 'ds']), log10(dnds_table[dnds_table$asn == 'GRC' & dnds_table$type == 'branch', 'dnds']), col=pal['GRC_branch'], pch=19)
    legend('topright', legend = c('all branches', 'Sciaridae core branches', 'GRC branches'), 
        col = c('black', pal['sci_core_branch'], pal['GRC_branch']), pch=19, bty = 'n')
dev.off()

png('figures/paml_dnds_vs_ds_per_tip_Sci_GRC.png', width=800, height=600)
    plot(log10(dnds_table[dnds_table$type == 'tip', 'ds']), log10(dnds_table[dnds_table$type == 'tip', 'dnds']), xlab = 'log10(ds)', ylab = 'log10(dnds)', main = 'dnds vs ds per tip')
    points(log10(dnds_table[dnds_table$asn == 'sci_core' & dnds_table$type == 'tip', 'ds']), log10(dnds_table[dnds_table$asn == 'sci_core' & dnds_table$type == 'tip', 'dnds']), col=pal['sci_core_branch'], pch=19)
    points(log10(dnds_table[dnds_table$asn == 'GRC' & dnds_table$type == 'tip', 'ds']), log10(dnds_table[dnds_table$asn == 'GRC' & dnds_table$type == 'tip', 'dnds']), col=pal['GRC_branch'], pch=19)
    legend('topright', legend = c('all tips', 'Sciaridae core tips', 'GRC tips'), 
        col = c('black', pal['sci_core_branch'], pal['GRC_branch']), pch=19, bty = 'n')
dev.off()

png('figures/paml_ds_histogram_all_data.png', width=800, height=600)
        hist(log10(dnds_table[, 'ds'])) # < 0.01
dev.off()

hist(dnds_table[, 'ds'], breaks = 100) # > 10

hist(log10(dnds_table[, 'dn'])) # < 0.001

dnds_table_ds_filtered <- dnds_table[dnds_table$ds > 0.01 & dnds_table$ds < 10, ] 
hist(log10(dnds_table_ds_filtered[, 'dn'])) 

dnds_table_filtered <- dnds_table[dnds_table$ds > 0.01 & dnds_table$dn > 0.001 & dnds_table$ds < 10, ] 
hist(log10(dnds_table_filtered[, 'dnds']))


png('figures/paml_dnds_vs_ds_per_tip_ds_filtered.png', width=800, height=600)
    plot(log10(dnds_table_ds_filtered[dnds_table_ds_filtered$type == 'tip', 'ds']), log10(dnds_table_ds_filtered[dnds_table_ds_filtered$type == 'tip', 'dnds']), xlab = 'log10(ds)', ylab = 'log10(dnds)', main = 'dnds vs ds per tip (ds filtered)')
    points(log10(dnds_table_ds_filtered[dnds_table_ds_filtered$asn == 'sci_core' & dnds_table_ds_filtered$type == 'tip', 'ds']), log10(dnds_table_ds_filtered[dnds_table_ds_filtered$asn == 'sci_core' & dnds_table_ds_filtered$type == 'tip', 'dnds']), col=pal['sci_core'], pch=19)
    points(log10(dnds_table_ds_filtered[dnds_table_ds_filtered$asn == 'ceci' & dnds_table_ds_filtered$type == 'tip', 'ds']), log10(dnds_table_ds_filtered[dnds_table_ds_filtered$asn == 'ceci' & dnds_table_ds_filtered$type == 'tip', 'dnds']), col=pal['ceci'], pch=19)
    legend('topright', legend = c('all tips', 'Sciaridae core tips', 'Cecidomyiidae tips'), 
        col = c('black', pal['sci_core'], pal['ceci']), pch=19, bty = 'n')
dev.off()

png('figures/paml_dnds_vs_ds_per_branch_ds_filtered.png', width=800, height=600)
    plot(log10(dnds_table_ds_filtered[dnds_table_ds_filtered$type == 'branch', 'ds']), log10(dnds_table_ds_filtered[dnds_table_ds_filtered$type == 'branch', 'dnds']), xlab = 'log10(ds)', ylab = 'log10(dnds)', main = 'dnds vs ds per branch (ds filtered)')
    points(log10(dnds_table_ds_filtered[dnds_table_ds_filtered$asn == 'sci_core' & dnds_table_ds_filtered$type == 'branch', 'ds']), log10(dnds_table_ds_filtered[dnds_table_ds_filtered$asn == 'sci_core' & dnds_table_ds_filtered$type == 'branch', 'dnds']), col=pal['sci_core_branch'], pch=19)
    points(log10(dnds_table_ds_filtered[dnds_table_ds_filtered$asn == 'GRC' & dnds_table_ds_filtered$type == 'branch', 'ds']), log10(dnds_table_ds_filtered[dnds_table_ds_filtered$asn == 'GRC' & dnds_table_ds_filtered$type == 'branch', 'dnds']), col=pal['GRC_branch'], pch=19)
    legend('topright', legend = c('all branches', 'Sciaridae core branches', 'GRC branches'), 
        col = c('black', pal['sci_core_branch'], pal['GRC_branch']), pch=19, bty = 'n')
dev.off()

png('figures/paml_dnds_boxplots_ds_filtering.png', width=800, height=600)
boxplot(log10dnds ~ type + asn, data = dnds_table_ds_filtered, 
        main = 'log10(dN/dS) across branch assignments and types',
        xlab = 'Branch Assignment and Type',
        ylab = 'log10(dN/dS)',
        col = pal[c('GRC_branch', 'GRC', 'ceci', 'ceci_branch',  'grey',  'grey', 'sci_core_branch', 'sci_core')]
        )
dev.off()

png('figures/paml_dnds_boxplots_ds_filtering_by_clade.png', width=800, height=600)
boxplot(log10dnds ~ asn, data = dnds_table_ds_filtered, 
        main = 'log10(dN/dS) across branch assignments and types',
        xlab = 'Branch Assignment',
        ylab = 'log10(dN/dS)',
        col = pal[c('GRC', 'ceci', 'grey', 'sci_core')]
        )
dev.off()

png('figures/paml_dnds_boxplots_dn_and_ds_filtering.png', width=800, height=600)
boxplot(log10dnds ~ asn + type, data = dnds_table_filtered, 
        main = 'log10(dN/dS) across branch assignments and types',
        xlab = 'Branch Assignment and Type',
        ylab = 'log10(dN/dS)',
        col = pal[c('GRC_branch', 'ceci_branch', 'grey', 'sci_core_branch', 'GRC', 'ceci', 'grey','sci_core')]
        )
dev.off()

wilcox.test( dnds_table_filtered[dnds_table_filtered$asn == 'sci_core' & dnds_table_filtered$type == 'tip', 'dnds'],
             dnds_table_filtered[dnds_table_filtered$asn == 'GRC' & dnds_table_filtered$type == 'tip', 'dnds'],
             alternative = 'two.sided'
           )

median( dnds_table_filtered[dnds_table_filtered$asn == 'sci_core' & dnds_table_filtered$type == 'tip', 'dnds'])
median( dnds_table_filtered[dnds_table_filtered$asn == 'GRC' & dnds_table_filtered$type == 'tip', 'dnds'])


median( dnds_table_filtered[dnds_table_filtered$asn == 'sci_core' & dnds_table_filtered$type == 'branch', 'dnds'])
median( dnds_table_filtered[dnds_table_filtered$asn == 'GRC' & dnds_table_filtered$type == 'branch', 'dnds'])

median( dnds_table_filtered[dnds_table_filtered$asn == 'sci_core', 'dnds'])
median( dnds_table_filtered[dnds_table_filtered$asn == 'GRC', 'dnds'])

wilcox.test( dnds_table_filtered[dnds_table_filtered$asn == 'sci_core', 'dnds'],
             dnds_table_filtered[dnds_table_filtered$asn == 'GRC', 'dnds'],
             alternative = 'two.sided'
           )
