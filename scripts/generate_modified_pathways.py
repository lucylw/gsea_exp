import os
import copy
import tqdm
import random


def write_gmt_file(output_file, gs_output):
    """
    Write a GSEA gmt file
    :param output_file: name of output file
    :param gene_sets: list of gene sets; each entry is (gene_set_name, gene_set_origin, list of gene symbols)
    :return:
    """
    assert not os.path.exists(output_file)

    with open(output_file, 'w') as f:
        for gs_name, gs_entry in gs_output.items():
            f.write('{}\t{}\t{}\n'.format(
                gs_name,
                gs_entry['origin'],
                '\t'.join(gs_entry['genes'])
            ))


# load original data (MSigDB Reactome pathways)
base_dir = os.path.expanduser('~/git/synthetic_gsea/data/')
orig_gmt_file = os.path.join(base_dir, 'c2.cp.reactome.v6.2.symbols.gmt')

gene_sets = dict()

with open(orig_gmt_file, 'r') as f:
    contents = f.readlines()
    for line in contents:
        parts = line.strip().split('\t')
        gs_name = parts[0]
        entry = {
            'origin': parts[1],
            'genes': parts[2:]
        }
        gene_sets[gs_name] = entry

uniq_genes = [entry['genes'] for _, entry in gene_sets.items()]
uniq_genes = [g for gs in uniq_genes for g in gs]
uniq_genes = list(set(uniq_genes))
print('Unique genes: {}'.format(len(uniq_genes)))

# redundant percentages
percs = [0.0, 0.04, 0.08, 0.12, 0.16, 0.2]

for perc_redundant in tqdm.tqdm(percs):

    modified_gene_sets = copy.copy(gene_sets)

    redundant_genes = random.sample(uniq_genes, int(perc_redundant * len(uniq_genes)))

    for gene in redundant_genes:
        including_gsets = [
            gs_name for gs_name, gs_entry in modified_gene_sets.items()
            if gene in gs_entry['genes']
        ]
        new_gene_name = gene + '_REDUNDANT'
        mod_gsets = random.sample(including_gsets, int(0.5 * len(including_gsets)))

        for gs in mod_gsets:
            orig_genes = modified_gene_sets[gs]['genes']
            modified_gene_sets[gs]['genes'] = [
                new_gene_name if g == gene else g for g in orig_genes
            ]

    outfile = os.path.join(base_dir, 'reactome_gene_sets_{0:.2f}.gmt'.format(perc_redundant))
    write_gmt_file(outfile, modified_gene_sets)



