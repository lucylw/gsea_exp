import os
import shutil
import copy
import tqdm
import random
import pickle

import gseapy
from collections import defaultdict

import pandas as pd


def process_all_leading_genes(f_path):
    """
    Read all leading edge genes from output file
    :param f:
    :return:
    """
    with open(f_path, 'r') as f:
        contents = f.read()
    parts = contents.strip().split('\t')
    genes = parts[2:]
    return genes


def process_leading_genes(f_path):
    """
    Process genes from leading genes file
    :param f_path:
    :return:
    """
    with open(f_path, 'r') as f:
        lines = f.readlines()

    gene_list = lines[2].strip().split('\t')[2:]
    gene_sets = []

    by_gene = defaultdict(list)

    for l in lines[3:]:
        parts = l.strip().split('\t')
        gene_sets.append(parts[0])
        for p, gene in zip(parts[2:], gene_list):
            by_gene[gene].append(int(p))

    totals = [(k, sum(v)) for k, v in by_gene.items()]
    totals.sort(key=lambda x: x[1], reverse=True)

    return totals


def process_results_file(f_path):
    """
    Process GSEA summary results file
    :param f_path:
    :return:
    """
    results = pd.read_csv(f_path, sep='\t', header=0)
    keep_cols = {'GS', 'SIZE', 'ES', 'NES', 'p-val'}
    results = results[:20].filter(keep_cols)
    return results


def write_cls_file(output_file, classes, samples):
    """
    Generate a GSEA cls file
    :param output_file: name of output file
    :param classes: list of class names
    :param samples: number of samples per class in order given in classes
    :return:
    """
    assert (len(samples) == len(classes))
    assert (len(set(classes)) == len(classes))

    total_samples = sum(samples)
    total_classes = len(classes)

    with open(output_file, 'w') as f:
        f.write('{} {} 1\n'.format(total_samples, total_classes))
        f.write('# {}\n'.format('\t'.join(classes)))
        for s_num, cl_name in zip(samples, classes):
            for i in range(s_num):
                f.write('{} '.format(cl_name))
        f.write('\n')


def write_gct_file(output_file, class_names, class_counts, expression_matrix):
    """
    Generate a GSEA gct file
    :param output_file: name of output file
    :param expression_matrix: dataframe containing expression levels
    :return:
    """
    total_genes = len(expression_matrix)
    first_key = list(expression_matrix.keys())[0]
    total_samples = len(expression_matrix[first_key])

    headers = ['NAME', 'DESCRIPTION']

    for c_name, c_count in zip(class_names, class_counts):
        for i in range(c_count):
            headers.append('{}_{}'.format(c_name, i + 1))

    with open(output_file, 'w') as f:
        f.write('#1.2\n')
        f.write('{} {}\n'.format(total_genes, total_samples))
        f.write('\t'.join(headers))
        f.write('\n')

        for g_name, values in expression_matrix.items():
            f.write(g_name)
            f.write('\tna\t')
            f.write('\t'.join(
                ['{0:.2f}'.format(v) for v in values]
            ))
            f.write('\n')


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


# set directories
base_dir = os.path.expanduser('~/git/synthetic_gsea/')
gene_exp_file = os.path.join(base_dir, 'data', 'tcga_luad_gene_exp.tsv')
pt_info_file = os.path.join(base_dir, 'data', 'tcga_luad_pt_info.tsv')
orig_gmt_file = os.path.join(base_dir, 'data', 'c2.cp.reactome.v6.2.symbols.gmt')

# # read clinical and exposure data from file
# pt_data = pd.read_csv(pt_info_file, sep='\t')
# controls = pt_data[pt_data.tissue_type == 'control'].file_name.tolist()
# tumors = pt_data[pt_data.tissue_type == 'tumor'].file_name.tolist()
#
# # read gene expression from file
# gene_exp = pd.read_csv(gene_exp_file, sep='\t', index_col=0)
#
# # write class file
# cls_file = os.path.join(base_dir, 'output', 'gsea_exp.cls')
# cls_names = ['CONTROL', 'TUMOR']
# cls_counts = [len(controls), len(tumors)]
# write_cls_file(cls_file, cls_names, cls_counts)
#
# # write gene expression file
# gct_file = os.path.join(base_dir, 'output', 'gsea_exp.gct')
# control_data = gene_exp.loc[controls]
# tumor_data = gene_exp.loc[tumors]
# final_data = pd.concat([control_data, tumor_data])
#
# gene_names = final_data.columns.tolist()
# exp_data = final_data.values.T
# exp_matrix = {g_name: g_exp for g_name, g_exp in zip(gene_names, exp_data)}
# write_gct_file(gct_file, cls_names, cls_counts, exp_matrix)


# read gene sets from original file
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
iterations = 10

for perc_redundant in tqdm.tqdm(percs):
    print('Perc redundant: {}'.format(perc_redundant))

    for i in range(iterations):
        print('\ti = {}'.format(i))

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

        # write modified gene sets to disk
        gmt_file = os.path.join(
            base_dir,
            'output',
            'gsea_{}'.format(perc_redundant),
            'reactome_gene_sets_{0:.2f}.gmt'.format(perc_redundant)
        )

        write_gmt_file(gmt_file, modified_gene_sets)

        # run GSEA
        cls_file = os.path.join(base_dir, 'output', 'gsea_exp.cls')
        gct_file = os.path.join(base_dir, 'output', 'gsea_exp.gct')

        gsea_dir = os.path.join(base_dir, 'output', 'gsea_{}'.format(perc_redundant), 'gsea_output')
        shutil.rmtree(gsea_dir)
        os.mkdir(gsea_dir)
        gseapy.gsea(
            data=gct_file, gene_sets=gmt_file, cls=cls_file, outdir=gsea_dir
        )

        # gsea output files to process
        tumor_all_leading_genes_file = os.path.join(
            gsea_dir,
            'syngsea.all.leading.genes.TUMOR.gmt'
        )

        tumor_leading_genes_file = os.path.join(
            gsea_dir,
            'syngsea.leading.genes.TUMOR.gct'
        )

        tumor_summary_results_file = os.path.join(
            gsea_dir,
            'syngsea.SUMMARY.RESULTS.REPORT.TUMOR.txt'
        )

        tumor_leading_genes = process_all_leading_genes(tumor_all_leading_genes_file)
        tumor_leading_gene_occurrences = process_leading_genes(tumor_leading_genes_file)
        tumor_summary_results = process_results_file(tumor_summary_results_file)

        gsea_output_dict = {
            'leading_genes': tumor_leading_genes,
            'leading_genes_by_occurrence': tumor_leading_gene_occurrences,
            'summary': tumor_summary_results,
            'gene_sets': modified_gene_sets
        }

        # save to pickle
        gsea_pickle_file = os.path.join(
            base_dir,
            'output',
            'gsea_{}'.format(perc_redundant),
            'trial_{}.pkl'.format(i)
        )

        pickle.dump(gsea_output_dict, open(gsea_pickle_file, 'wb'))