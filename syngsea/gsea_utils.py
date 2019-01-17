import os
from typing import List

import pandas as pd


def generate_cls_file(output_file: str, classes: List, samples: List):
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
    return output_file


def generate_gct_file(output_file: str, expression_matrix: pd.DataFrame):
    """
    Generate a GSEA gct file
    :param output_file: name of output file
    :param expression_matrix: dataframe containing expression levels
    :return:
    """
    total_genes = len(expression_matrix)
    total_samples = len(expression_matrix.columns) - 2

    with open(output_file, 'w') as f:
        f.write('#1.2\n')
        f.write('{} {}\n'.format(total_genes, total_samples))
        f.write('\t'.join(expression_matrix.columns))
        f.write('\n')
        for index, row in expression_matrix.iterrows():
            values = row.tolist()
            f.write('\t'.join(values[0:2]))
            f.write('\t')
            values = ["{0:.2f}".format(v) for v in values[2:]]
            f.write('\t'.join(values))
            f.write('\n')

    return output_file


def generate_gmt_file(output_file: str, gene_sets: List):
    """
    Generate a GSEA gmt file
    :param output_file: name of output file
    :param gene_sets: list of gene sets; each entry is (gene_set_name, gene_set_origin, list of gene symbols)
    :return:
    """
    with open(output_file, 'w') as f:
        for gs_name, gs_origin, symbols in gene_sets:
            f.write('{}\t{}\t{}\n'.format(
                gs_name,
                gs_origin,
                '\t'.join(symbols)
            ))
    return output_file
