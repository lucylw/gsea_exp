import os
import tqdm
import pandas as pd
import requests
import pickle

from bioservices import BioDBNet
from pybiomart import Server

import syngsea.gsea_utils as gsea_utils
from syngsea.paths import GSEAFilePath


paths = GSEAFilePath()

study_dict = {
    'ACT': 'ACT_study_data',
    'AMP': 'AMP_AD',
    'LUAD': 'TCGA_LUAD',
    'HNSCC': 'TCGA_HNSCC'
}

study_paths = {
    k: os.path.join(paths.data_dir, v)
    for k, v in study_dict.items()
}


def get_ensembl_dict():
    """
    Create BioMart server and get Ensembl id mapping dict
    :return:
    """
    # map ensemble identifiers to gene names
    server = Server(host='http://www.ensembl.org')

    ensembl_dataset = (
        server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl']
    )

    results = ensembl_dataset.query(
        attributes=['ensembl_gene_id', 'external_gene_name']
    )

    ens_dict = dict(zip(results['Gene stable ID'], results['Gene name']))

    return ens_dict


def map_ensembl_ids(ensembl_ids, lookup_dict, tmp_file=None):
    """
    Map Ensembl identifiers to Gene symbol using BioDBNet
    :param ensembl_ids:
    :param lookup_dict:
    :param tmp_file:
    :return:
    """
    input_db = 'Ensembl Gene ID'
    output_db = ['Gene Symbol']
    species = 9606
    service = BioDBNet()

    gene_symbols = []

    buffer = []
    missing = 0
    num_run = 200

    for count, ens_id in enumerate(tqdm.tqdm(ensembl_ids)):
        if ens_id in lookup_dict:
            gene_symbols.append(lookup_dict[ens_id])
        else:
            buffer.append(ens_id)
            if len(buffer) == num_run:
                results = service.db2db(input_db, output_db, buffer, species)

                for e_id, g_symb in zip(results.index, results['Gene Symbol']):
                    if g_symb != '-':
                        gene_symbols.append(g_symb)
                        lookup_dict[e_id] = g_symb
                    else:
                        gene_symbols.append(e_id)
                        lookup_dict[e_id] = e_id
                        missing += 1

                buffer = []

                if tmp_file:
                    pickle.dump(lookup_dict, open(tmp_file, 'wb'))

    # process remainder
    if buffer:
        results = service.db2db(input_db, output_db, buffer, species)

        for e_id, g_symb in zip(results.index, results['Gene Symbol']):
            if g_symb != '-':
                gene_symbols.append(g_symb)
                lookup_dict[e_id] = g_symb
            else:
                gene_symbols.append(e_id)
                lookup_dict[e_id] = e_id
                missing += 1

    print('Total missed: {}'.format(missing))

    if tmp_file:
        pickle.dump(lookup_dict, open(tmp_file, 'wb'))

    return gene_symbols, lookup_dict


def map_ensembl_ids_to_hgnc_name(ensembl_ids, lookup_dict, tmp_file=None):
    """
    Map Ensembl IDs to HGNC Gene names
    :param ensembl_ids:
    :param lookup_dict:
    :param temporary_file:
    :return:
    """
    ensembl_server = "https://rest.ensembl.org"
    hgnc_names = []

    print('Looking up {} Ensembl Ids'.format(len(ensembl_ids)))

    try:
        for count, ens_id in enumerate(tqdm.tqdm(ensembl_ids)):
            if ens_id in lookup_dict:
                hgnc_names.append(lookup_dict[ens_id])
            else:
                ext = "/xrefs/id/{}?".format(ens_id)

                r = requests.get(ensembl_server + ext, headers={
                    "external_db": "HGNC",
                    "Content-Type": "application/json",
                    "all_levels": "1"
                })

                add_id = ens_id

                if r.ok:
                    decoded = r.json()
                    for entry in decoded:
                        if entry['dbname'] == 'HGNC':
                            add_id = entry['display_id']

                hgnc_names.append(add_id)
                lookup_dict[ens_id] = add_id

            if count % 1000 == 0:
                if tmp_file:
                    pickle.dump(lookup_dict, open(tmp_file, 'wb'))

    except Exception:
        if tmp_file:
            pickle.dump(lookup_dict, open(tmp_file, 'wb'))
        raise Exception
    except KeyboardInterrupt:
        if tmp_file:
            pickle.dump(lookup_dict, open(tmp_file, 'wb'))
        raise KeyboardInterrupt

    if tmp_file:
        pickle.dump(lookup_dict, open(tmp_file, 'wb'))

    return hgnc_names, lookup_dict


def load_act():
    """
    Parse ACT data to GSEA format
    :return:
    """
    print('Processing ACT data...')

    normalized_data_file = os.path.join(
        study_paths['ACT'], 'fpkm_table_normalized.csv'
    )
    gene_id_file = os.path.join(
        study_paths['ACT'], 'rows-genes.csv'
    )
    patient_id_file = os.path.join(
        study_paths['ACT'], 'columns-samples.csv'
    )
    metadata_file = os.path.join(
        study_paths['ACT'], 'DonorInformation.csv'
    )

    # load gene identifier data
    gene_data = pd.read_csv(gene_id_file)
    keep_gene_cols = ['gene_id', 'gene_symbol']
    gene_data = gene_data.filter(items=keep_gene_cols)

    # load patient identifier data
    patient_data = pd.read_csv(patient_id_file)
    keep_pt_cols = ['rnaseq_profile_id', 'donor_id', 'specimen_id', 'structure_name']
    patient_data = patient_data.filter(items=keep_pt_cols)

    metadata = pd.read_csv(metadata_file)
    keep_md_cols = ['donor_id', 'apo_e4_allele', 'ever_tbi_w_loc', 'act_demented']
    metadata = metadata.filter(items=keep_md_cols)

    pt = patient_data.merge(metadata, on='donor_id')

    # load gene expression data
    exp_data = pd.read_csv(normalized_data_file)
    exp_data.rename(index=str, columns={"gene_id \ rnaseq_profile_id": "gene_id"}, inplace=True)

    df = gene_data.merge(exp_data, on='gene_id', how='inner')
    df = df.drop(columns=['gene_id'])
    df = df.rename(index=str, columns={"gene_symbol": "rnaseq_profile_id"})
    df.set_index('rnaseq_profile_id', inplace=True)
    df = df.T
    df.index = df.index.astype('int64')
    df = pt.merge(df, left_on='rnaseq_profile_id', right_index=True)

    structures = [
        'hippocampus (hippocampal formation)',
        'parietal neocortex',
        'temporal neocortex',
        'white matter of forebrain'
    ]

    structure_abbrev = {
        'hippocampus (hippocampal formation)': 'hippo',
        'parietal neocortex': 'p_neo',
        'temporal neocortex': 't_neo',
        'white matter of forebrain': 'fore'
    }

    partitions = {
        structure_abbrev[k]: df.loc[df['structure_name'] == k]
        for k in structures
    }

    partitions = {
        k: v.drop(columns=['specimen_id', 'structure_name', 'apo_e4_allele', 'ever_tbi_w_loc'])
            .sort_values(['act_demented'], ascending=[True])
        for k, v in partitions.items()
    }

    for k, v in partitions.items():
        print('\tProcessing {}...'.format(k))
        cls_file = os.path.join(paths.processed_data_dir, '{}_{}.cls'.format(
            'ACT', k
        ))
        gct_file = os.path.join(paths.processed_data_dir, '{}_{}.gct'.format(
            'ACT', k
        ))

        dementia = v.loc[v['act_demented'] == 'Dementia']
        controls = v.loc[v['act_demented'] == 'No Dementia']

        gene_exp = v.drop(columns=['rnaseq_profile_id', 'act_demented'])
        gene_exp.set_index('donor_id', inplace=True)
        gene_exp = gene_exp.T
        gene_exp.index.name = 'NAME'
        gene_exp.insert(loc=0, column='DESCRIPTION', value=['na']*len(gene_exp))
        gene_exp.reset_index(level=0, inplace=True)

        column_headers = ['AD_{}'.format(i) for i in range(1, len(dementia) + 1)] \
                         + ['CONTROL_{}'.format(i) for i in range(1, len(controls) + 1)]
        gene_exp.columns = ['NAME', 'DESCRIPTION'] + column_headers

        gsea_utils.generate_cls_file(cls_file, ['AD', 'CONTROL'], [len(dementia), len(controls)])
        gsea_utils.generate_gct_file(gct_file, gene_exp)


def load_amp(ensembl_lookup_dict):
    """
    Parse AMP data into GSEA format
    :param ensembl_lookup_dict:
    :return:
    """
    brodmann_files = {
        'BM10': 'AMP-AD_MSBB_MSSM_BM_10.normalized.sex_race_age_RIN_PMI_exonicRate_rRnaRate_batch_adj.tsv',
        'BM22': 'AMP-AD_MSBB_MSSM_BM_22.normalized.sex_race_age_RIN_PMI_exonicRate_rRnaRate_batch_adj.tsv',
        'BM36': 'AMP-AD_MSBB_MSSM_BM_36.normalized.sex_race_age_RIN_PMI_exonicRate_rRnaRate_batch_adj.tsv',
        'BM44': 'AMP-AD_MSBB_MSSM_BM_44.normalized.sex_race_age_RIN_PMI_exonicRate_rRnaRate_batch_adj.tsv'
    }

    brodmann_map = {
        'BM10': 'anterior_prefrontal_cortex',
        'BM22': 'left_posterior_superior_temporal_gyrus',
        'BM36': 'temporal/fusiform_gyrus',
        'BM44': 'pars opercularis/brocas'
    }

    clinical_file = os.path.join(study_paths['AMP'], 'metainfo', 'MSBB_clinical.csv')
    rnaseq_id_file = os.path.join(study_paths['AMP'], 'metainfo', 'MSBB_RNAseq_covariates_November2018Update.csv')

    clinical_data = pd.read_csv(clinical_file, sep=',', quotechar='"')
    keep_cl_cols = ['individualIdentifier', 'CDR', 'NP.1']
    clinical_data = clinical_data.filter(items=keep_cl_cols)
    clinical_data['has_AD'] = clinical_data['NP.1'] > 1

    rnaseq_ids = pd.read_csv(rnaseq_id_file, sep=',')
    keep_rna_cols = ['sampleIdentifier', 'BrodmannArea', 'individualIdentifier']
    rnaseq_ids = rnaseq_ids.filter(items=keep_rna_cols)
    rnaseq_ids = rnaseq_ids.drop_duplicates()

    pt = rnaseq_ids.merge(clinical_data, on='individualIdentifier')
    pt = pt.sort_values(['has_AD', 'BrodmannArea'], ascending=[False, True])

    temp_file = os.path.join(paths.base_dir, 'temp', 'gene_lookup_dict.pickle')
    if os.path.exists(temp_file):
        gene_dict = pickle.load(open(temp_file, 'rb'))
    else:
        gene_dict = dict()

    for bm_area, bm_file in brodmann_files.items():
        print(bm_area)
        cls_file = os.path.join(paths.processed_data_dir, '{}_{}.cls'.format(
            'AMP', bm_area
        ))
        gct_file = os.path.join(paths.processed_data_dir, '{}_{}.gct'.format(
            'AMP', bm_area
        ))

        pt_bm = pt.loc[pt['BrodmannArea'] == bm_area]
        keep_cols = ['sampleIdentifier', 'has_AD']
        pt_bm = pt_bm.filter(items=keep_cols)

        file_path = os.path.join(study_paths['AMP'], bm_file)
        df = pd.read_csv(file_path, sep='\t')

        df = df.T
        df['sampleIdentifier'] = [index.split('.')[1] for index in df.index]
        df = pt_bm.merge(df, on='sampleIdentifier')

        dementia = df.loc[df['has_AD'] == True]
        controls = df.loc[df['has_AD'] == False]

        gene_exp = df.drop(columns=['sampleIdentifier', 'has_AD'])
        gene_exp = gene_exp.T

        # map gene expression index to
        ensembl_ids = gene_exp.index

        # map ensemble genes to gene names
        all_genes = [ens_id for ens_id in ensembl_ids if ens_id.split('.')[0] in ensembl_lookup_dict]
        print('Keep {} genes out of {}'.format(len(all_genes), len(ensembl_ids)))

        keep_rows = [ind in all_genes for ind in gene_exp.index]
        gene_exp = gene_exp[keep_rows]
        gene_names = [ensembl_lookup_dict[ens_id] for ens_id in gene_exp.index]

        gene_exp.index.name = 'NAME'
        gene_exp.insert(loc=0, column='NAME', value=gene_names)
        gene_exp.insert(loc=1, column='DESCRIPTION', value=['na']*len(gene_exp))
        gene_exp.reset_index(level=0, inplace=True, drop=True)

        column_headers = ['AD_{}'.format(i) for i in range(1, len(dementia) + 1)] \
                         + ['CONTROL_{}'.format(i) for i in range(1, len(controls) + 1)]
        gene_exp.columns = ['NAME', 'DESCRIPTION'] + column_headers

        gsea_utils.generate_cls_file(cls_file, ['AD', 'CONTROL'], [len(dementia), len(controls)])
        gsea_utils.generate_gct_file(gct_file, gene_exp)


def load_tcga(header):
    """
    Parse TCGA data to GSEA format
    :param header:
    :return:
    """
    print('TCGA_{}'.format(header))

    cls_file = os.path.join(paths.processed_data_dir, '{}.cls'.format(header))
    gct_file = os.path.join(paths.processed_data_dir, '{}.gct'.format(header))

    pt_info_file = os.path.join(study_paths[header], 'tcga_{}_pt_info.tsv'.format(header.lower()))
    gene_exp_file = os.path.join(study_paths[header], 'tcga_{}_gene_exp.tsv'.format(header.lower()))

    pt_info = pd.read_csv(pt_info_file, sep='\t')
    pt_info = pt_info[['tissue_type', 'file_name']]

    gene_exp = pd.read_csv(gene_exp_file, sep='\t')
    gene_exp = gene_exp.rename(index=str, columns={'Unnamed: 0': 'file_name'})

    df = pt_info.merge(gene_exp, on='file_name')
    df.sort_values(['tissue_type'], ascending=[False])

    tumor = df.loc[pt_info['tissue_type'] == 'tumor']
    control = df.loc[pt_info['tissue_type'] == 'control']

    gene_exp = df.drop(columns=['tissue_type', 'file_name'])
    gene_exp = gene_exp.T

    gene_exp.index.name = 'NAME'
    gene_exp.insert(loc=0, column='NAME', value=gene_exp.index)
    gene_exp.insert(loc=1, column='DESCRIPTION', value=['na'] * len(gene_exp))
    gene_exp.reset_index(level=0, inplace=True, drop=True)

    column_headers = ['TUMOR_{}'.format(i) for i in range(1, len(tumor) + 1)] \
                     + ['CONTROL_{}'.format(i) for i in range(1, len(control) + 1)]
    gene_exp.columns = ['NAME', 'DESCRIPTION'] + column_headers

    gsea_utils.generate_cls_file(cls_file, ['TUMOR', 'CONTROL'], [len(tumor), len(control)])
    gsea_utils.generate_gct_file(gct_file, gene_exp)



ensembl_dict = get_ensembl_dict()

load_act()
load_amp(ensembl_dict)
load_tcga('LUAD')
load_tcga('HNSCC')
