import os
import csv
import json
import glob
import gzip
import tqdm

import pandas as pd

from collections import OrderedDict
from pybiomart import Server


class TCGAParser:
    def __init__(
            self,
            data_path,
            clinical_file,
            exposure_file,
            metadata_file,
            pt_outfile,
            gexp_outfile
    ):
        self.data_path = data_path
        self.clinical_file = clinical_file
        self.exposure_file = exposure_file
        self.metadata_file = metadata_file
        self.pt_outfile = pt_outfile
        self.gexp_outfile = gexp_outfile

    def _parse_metadata(self):
        """
        Parse clinical data and metadata
        :return:
        """
        # read clinical and exposure data from file
        clinical_data = pd.read_csv(self.clinical_file, sep='\t', index_col='submitter_id')
        exposure_data = pd.read_csv(self.exposure_file, sep='\t', index_col='submitter_id')

        clinical_data = pd.concat((clinical_data, exposure_data), axis=1, join='outer')

        # read and parse metadata from file
        metadata_content = json.load(open(self.metadata_file, 'r'))

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

        return pt_data

    def _parse_gexp_data(self):
        """
        Parse gene expression data
        :return:
        """
        # extract gene expression data
        gene_exp = dict()

        for file_path in tqdm.tqdm(
                glob.glob(os.path.join(self.data_path, '*', '*.FPKM.txt.gz'))
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

        return converted_gene_exp

    def parse(self):
        """
        Parse TCGA gene expression data
        :return:
        """
        print('Parsing metadata...')
        pt_data = self._parse_metadata()
        pt_data.to_csv(self.pt_outfile, sep='\t', index=False)

        print('Parsing gene expression data')
        gexp_data = self._parse_gexp_data()
        gexp_data.to_csv(self.gexp_outfile, sep='\t')


if __name__ == '__main__':
    base_path = os.path.expanduser('~/git/')

    # # process LUAD data
    # print('---LUAD---')
    # luad_tcga_path = os.path.join(base_path, 'TCGA_data', 'TCGA-LUAD')
    # luad_out_path = os.path.join(base_path, 'gsea_exp', 'data', 'TCGA_LUAD')
    # luad_header = 'luad'
    #
    # luad_data_path = os.path.join(luad_tcga_path, 'gene_exp')
    # luad_clinical_file = os.path.join(luad_tcga_path, 'clinical.tsv')
    # luad_exposure_file = os.path.join(luad_tcga_path, 'exposure.tsv')
    # luad_metadata_file = os.path.join(luad_tcga_path, 'metadata.cart.2018-10-22.json')
    #
    # luad_pt_outfile = os.path.join(luad_out_path, 'tcga_{}_pt_info.tsv'.format(luad_header))
    # luad_gexp_outfile = os.path.join(luad_out_path, 'tcga_{}_gene_exp.tsv'.format(luad_header))
    #
    # parser = TCGAParser(
    #     data_path=luad_data_path,
    #     clinical_file=luad_clinical_file,
    #     exposure_file=luad_exposure_file,
    #     metadata_file=luad_metadata_file,
    #     pt_outfile=luad_pt_outfile,
    #     gexp_outfile=luad_gexp_outfile
    # )
    # parser.parse()

    # process HNSCC data
    print('---HNSCC---')
    hnscc_tcga_path = os.path.join(base_path, 'TCGA_data', 'TCGA-HNSC')
    hnscc_out_path = os.path.join(base_path, 'gsea_exp', 'data', 'TCGA_HNSCC')
    hnscc_header = 'hnscc'

    hnscc_data_path = os.path.join(hnscc_tcga_path, 'gene_exp')
    hnscc_clinical_file = os.path.join(hnscc_tcga_path, 'clinical.tsv')
    hnscc_exposure_file = os.path.join(hnscc_tcga_path, 'exposure.tsv')
    hnscc_metadata_file = os.path.join(hnscc_tcga_path, 'metadata.cart.2019-01-23.json')

    hnscc_pt_outfile = os.path.join(hnscc_out_path, 'tcga_{}_pt_info.tsv'.format(hnscc_header))
    hnscc_gexp_outfile = os.path.join(hnscc_out_path, 'tcga_{}_gene_exp.tsv'.format(hnscc_header))

    parser = TCGAParser(
        data_path=hnscc_data_path,
        clinical_file=hnscc_clinical_file,
        exposure_file=hnscc_exposure_file,
        metadata_file=hnscc_metadata_file,
        pt_outfile=hnscc_pt_outfile,
        gexp_outfile=hnscc_gexp_outfile
    )
    parser.parse()
