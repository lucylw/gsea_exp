import os
import glob

import pandas as pd
import gseapy


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


# set directories
base_dir = os.path.expanduser('~/git/synthetic_gsea/')
gene_exp_file = os.path.join(base_dir, 'data', 'tcga_luad_gene_exp.tsv')
pt_info_file = os.path.join(base_dir, 'data', 'tcga_luad_pt_info.tsv')

# read clinical and exposure data from file
pt_data = pd.read_csv(pt_info_file, sep='\t')
controls = pt_data[pt_data.tissue_type == 'control'].file_name.tolist()
tumors = pt_data[pt_data.tissue_type == 'tumor'].file_name.tolist()

# read gene expression from file
gene_exp = pd.read_csv(gene_exp_file, sep='\t', index_col=0)

# write class file
cls_file = os.path.join(base_dir, 'output', 'gsea_exp.cls')
cls_names = ['CONTROL', 'TUMOR']
cls_counts = [len(controls), len(tumors)]
write_cls_file(cls_file, cls_names, cls_counts)

# write gene expression file
gct_file = os.path.join(base_dir, 'output', 'gsea_exp.gct')
control_data = gene_exp.loc[controls]
tumor_data = gene_exp.loc[tumors]
final_data = pd.concat([control_data, tumor_data])

gene_names = final_data.columns.tolist()
exp_data = final_data.values.T
exp_matrix = {g_name: g_exp for g_name, g_exp in zip(gene_names, exp_data)}
write_gct_file(gct_file, cls_names, cls_counts, exp_matrix)

# # run GSEA
# cls_file = os.path.join(base_dir, 'output', 'gsea_exp.cls')
# gct_file = os.path.join(base_dir, 'output', 'gsea_exp.gct')
# gene_set_files = glob.glob(os.path.join(base_dir, 'data', 'reactome_gene_sets_*.gmt'))
#
# for gmt_file in gene_set_files:
#     gmt_name = os.path.basename(gmt_file)
#     increment = gmt_name[:-4].split('_')[-1]
#     output_dir = os.path.join(base_dir, 'output', 'gsea_{}'.format(increment))
#     os.mkdir(output_dir)
#     gseapy.gsea(
#         data=gct_file, gene_sets=gmt_file, cls=cls_file, outdir=output_dir
#     )
