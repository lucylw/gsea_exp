import os
import shutil
import copy
import tqdm
import random
import pickle

import rpy2.robjects as robjects
from collections import defaultdict

import pandas as pd
from multiprocessing import Pool


class GSEAExperiment:
    def __init__(self):
        # set directories
        self.base_dir = os.path.expanduser('~/git/synthetic_gsea/')
        self.gsea_r_location = os.path.join(self.base_dir, 'GSEA', 'GSEA.1.0.R')
        self.gene_exp_file = os.path.join(self.base_dir, 'data', 'tcga_luad_gene_exp.tsv')
        self.pt_info_file = os.path.join(self.base_dir, 'data', 'tcga_luad_pt_info.tsv')
        self.orig_gmt_file = os.path.join(self.base_dir, 'data', 'c2.cp.reactome.v6.2.symbols.gmt')
        # self._read_pt_data()
        self.gene_sets, self.uniq_genes = self._read_gset_file()
        self.iterations = 10

    @staticmethod
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

    @staticmethod
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

    @staticmethod
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

    @staticmethod
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

    @staticmethod
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

    @staticmethod
    def write_gmt_file(output_file, gs_output):
        """
        Write a GSEA gmt file
        :param output_file: name of output file
        :param gene_sets: list of gene sets; each entry is (gene_set_name, gene_set_origin, list of gene symbols)
        :return:
        """
        with open(output_file, 'w') as f:
            for gs_name, gs_entry in gs_output.items():
                f.write('{}\t{}\t{}\n'.format(
                    gs_name,
                    gs_entry['origin'],
                    '\t'.join(gs_entry['genes'])
                ))

    def _read_pt_data(self):
        # read clinical and exposure data from file
        pt_data = pd.read_csv(self.pt_info_file, sep='\t')
        controls = pt_data[pt_data.tissue_type == 'control'].file_name.tolist()
        tumors = pt_data[pt_data.tissue_type == 'tumor'].file_name.tolist()

        # read gene expression from file
        gene_exp = pd.read_csv(self.gene_exp_file, sep='\t', index_col=0)

        # write class file
        cls_file = os.path.join(self.base_dir, 'output', 'gsea_exp.cls')
        cls_names = ['CONTROL', 'TUMOR']
        cls_counts = [len(controls), len(tumors)]
        self.write_cls_file(cls_file, cls_names, cls_counts)

        # write gene expression file
        gct_file = os.path.join(self.base_dir, 'output', 'gsea_exp.gct')
        control_data = gene_exp.loc[controls]
        tumor_data = gene_exp.loc[tumors]
        final_data = pd.concat([control_data, tumor_data])

        gene_names = final_data.columns.tolist()
        exp_data = final_data.values.T
        exp_matrix = {g_name: g_exp for g_name, g_exp in zip(gene_names, exp_data)}
        self.write_gct_file(gct_file, cls_names, cls_counts, exp_matrix)
        return

    def _read_gset_file(self):
        # read gene sets from original file
        gene_sets = dict()

        with open(self.orig_gmt_file, 'r') as f:
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
        return gene_sets, uniq_genes

    def _run_gsea(self, gct_file, gmt_file, cls_file, gsea_dir):
        """
        Run GSEA using rpy2
        :param gct_file:
        :param gmt_file:
        :param cls_file:
        :param gsea_dir:
        :return:
        """
        r = robjects.r
        r.source(self.gsea_r_location)
        r("""GSEA(                             # Input/Output Files :-------------------------------------------
                 input.ds =  "{}",               # Input gene expression Affy dataset file in RES or GCT format
                 input.cls = "{}",               # Input class vector (phenotype) file in CLS format
                 gs.db =     "{}",    # Gene set database in GMT format
                 output.directory      = "{}/",            # Directory where to store output and results (default: "")
                #  Program parameters :----------------------------------------------------------------------------------------------------------------------------
                 doc.string            = "syngsea",     # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
                 non.interactive.run   = F,               # Run in interactive (i.e. R GUI) or batch (R command line) mode (default: F)
                 reshuffling.type      = "sample.labels", # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels" 
                 nperm                 = 1000,            # Number of random permutations (default: 1000)
                 weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
                 nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
                 fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
                 fdr.q.val.threshold   = 0.25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
                 topgs                 = 20,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
                 adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
                 gs.size.threshold.min = 10,              # Minimum size (in genes) for database gene sets to be considered (default: 25)
                 gs.size.threshold.max = 500,             # Maximum size (in genes) for database gene sets to be considered (default: 500)
                 reverse.sign          = F,               # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
                 preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
                 random.seed           = 111,             # Random number generator seed. (default: 123456)
                 perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
                 fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
                 replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)
                 save.intermediate.results = F,           # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
                 OLD.GSEA              = F,               # Use original (old) version of GSEA (default: F)
                 use.fast.enrichment.routine = T          # Use faster routine to compute enrichment for random permutations (default: T)
            )""".format(gct_file, cls_file, gmt_file, gsea_dir))

        r("""GSEA.Analyze.Sets(
               directory = "{}/",        # Directory where to store output and results (default: "")
               topgs = 20,            # number of top scoring gene sets used for analysis
               height = 16,
               width = 16
            )""".format(gsea_dir))

    def run_gsea_experiments(self, perc_redundant):
        """
        Run GSEA experiments for given redundancy percent and interation
        :param perc_redundant:
        :param iterations: number of experiments to run
        :param gsets: original gene sets from file
        :param uniqs: unique genes
        :return:
        """
        print('Perc redundant: {}'.format(perc_redundant))

        for i in range(self.iterations):
            print('\ti = {}'.format(i))

            modified_gene_sets = copy.copy(self.gene_sets)

            redundant_genes = random.sample(self.uniq_genes, int(perc_redundant * len(self.uniq_genes)))

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
                self.base_dir,
                'output',
                'gsea_{0:.2f}'.format(perc_redundant),
                'reactome_gene_sets_{0:.2f}.gmt'.format(perc_redundant)
            )

            self.write_gmt_file(gmt_file, modified_gene_sets)

            # run GSEA
            cls_file = os.path.join(self.base_dir, 'output', 'gsea_exp.cls')
            gct_file = os.path.join(self.base_dir, 'output', 'gsea_exp.gct')

            gsea_dir = os.path.join(self.base_dir, 'output', 'gsea_{0:.2f}'.format(perc_redundant), 'gsea_output')
            shutil.rmtree(gsea_dir)
            os.mkdir(gsea_dir)

            self._run_gsea(gct_file, gmt_file, cls_file, gsea_dir)

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

            tumor_leading_genes = self.process_all_leading_genes(tumor_all_leading_genes_file)
            tumor_leading_gene_occurrences = self.process_leading_genes(tumor_leading_genes_file)
            tumor_summary_results = self.process_results_file(tumor_summary_results_file)

            gsea_output_dict = {
                'leading_genes': tumor_leading_genes,
                'leading_genes_by_occurrence': tumor_leading_gene_occurrences,
                'summary': tumor_summary_results,
                'gene_sets': modified_gene_sets
            }

            # save to pickle
            gsea_pickle_file = os.path.join(
                self.base_dir,
                'output',
                'gsea_{0:.2f}'.format(perc_redundant),
                'trial_{}.pkl'.format(i)
            )

            pickle.dump(gsea_output_dict, open(gsea_pickle_file, 'wb'))


if __name__ == '__main__':

    gsea = GSEAExperiment()

    percs = [0.00, 0.04, 0.08, 0.12, 0.16, 0.20]

    # start processes
    p = Pool(6)
    p.map(gsea.run_gsea_experiments, percs)