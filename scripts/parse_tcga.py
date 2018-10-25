import os
import csv
import json
import glob
import gzip
import tqdm

import pandas as pd

from collections import OrderedDict, defaultdict
from pybiomart import Server

base_path = os.path.expanduser('~/git/TCGA_data/')

luad_data_path = os.path.join(base_path, 'LUAD')
clinical_file = os.path.join(base_path, 'clinical.tsv')
exposure_file = os.path.join(base_path, 'exposure.tsv')
metadata_file = os.path.join(base_path, 'metadata.cart.2018-10-22.json')

# read clinical and exposure data from file
clinical_data = pd.read_csv(clinical_file, sep='\t', index_col='submitter_id')
exposure_data = pd.read_csv(exposure_file, sep='\t', index_col='submitter_id')

clinical_data = pd.concat((clinical_data, exposure_data), axis=1, join='outer')

# read and parse metadata from file
metadata_content = json.load(open(metadata_file, 'r'))

metadata_data = []

for entry in metadata_content:
    f_id = entry['file_id']
    f_name = entry['file_name']
    assoc_entity = entry['associated_entities'][0]
    entity_id = assoc_entity['entity_id']
    entity_submitter_id = assoc_entity['entity_submitter_id']

    id_parts = entity_submitter_id.split('-')
    participant_id = '-'.join(id_parts[:3])
    sample_id = int(id_parts[3][:2])
    if sample_id < 10:
        tissue_type = 'tumor'
    else:
        tissue_type = 'control'

    metadata_data.append(OrderedDict({
        'submitter_id': participant_id,
        'tissue_type': tissue_type,
        'file_name': f_name,
        'entity_id': entity_id,
        'entity_submitter_id': entity_submitter_id
    }))

metadata = pd.DataFrame(metadata_data)

# combine metadata and clinical data
pt_data = pd.merge(metadata, clinical_data, on='submitter_id')

pt_outfile = '/Users/lwang/git/synthetic_gsea/data/tcga_luad_pt_info.tsv'
pt_data.to_csv(pt_outfile, sep='\t', index=False)

# extract gene expression data
gene_exp = dict()

for file_path in tqdm.tqdm(
        glob.glob(os.path.join(luad_data_path, '*', '*.FPKM.txt.gz'))
):
    parts = file_path.split('/')
    f_name = parts[-1]
    exp_dict = dict()
    with gzip.open(file_path, 'rt', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        for gene_id, expression in reader:
            exp_dict[gene_id] = float(expression)
    gene_exp[f_name] = exp_dict

print('Extract and sort Ensembl gene ids...')
all_genes = [list(exp.keys()) for k, exp in gene_exp.items()]
all_genes = [item for sublist in all_genes for item in sublist]
all_genes = list(set(all_genes))
all_genes.sort()

# map ensemble identifiers to gene names
server = Server(host='http://www.ensembl.org')

ensembl_dataset = (
    server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl']
)

results = ensembl_dataset.query(
    attributes=['ensembl_gene_id', 'external_gene_name']
)

ens_dict = dict(zip(results['Gene stable ID'], results['Gene name']))

# map ensemble genes to gene names
all_genes = [ens_id for ens_id in all_genes if ens_id.split('.')[0] in ens_dict]

# create sorted gene expression dictionary
sorted_gene_exp = dict()

for f_name, exp_dict in tqdm.tqdm(gene_exp.items()):
    get_gene_exp = lambda k: exp_dict[k] if k in exp_dict else None
    new_dict = OrderedDict({
        ens_dict[gene_id.split('.')[0]]: get_gene_exp(gene_id)
        for gene_id in all_genes
    })
    sorted_gene_exp[f_name] = new_dict

converted_gene_exp = pd.DataFrame.from_dict(sorted_gene_exp, orient='index')

gexp_outfile = '/Users/lwang/git/synthetic_gsea/data/tcga_luad_gene_exp.tsv'
converted_gene_exp.to_csv(gexp_outfile, sep='\t')
