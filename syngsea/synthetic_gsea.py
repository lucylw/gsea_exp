import os
import random
from datetime import datetime
import shutil
import pandas as pd
from collections import defaultdict

from numpy.random import normal, poisson
from syngsea.paths import SynGSEAFilePath

# generator class wrapper for numpy.random.normal
class NormalNoiseGenerator:
    def __init__(self, mu=0., sigma=1.):
        self.mu = mu
        self.sigma = sigma

    def __next__(self):
        return normal(self.mu, self.sigma)

    def __iter__(self):
        return self


# generator class wrapper for numpy.random.poisson
class PoissonNoiseGenerator:
    def __init__(self, lam=1.):
        self.lam = lam

    def __next__(self):
        return poisson(self.lam)

    def __iter__(self):
        return self


# class for conducting synthetic GSEA experiments
class SyntheticGSEA:
    def __init__(
            self,
            perc_redundant,
            num_patients=100,
            total_genes=2000,
            expressed_genes=500,
            perc_diff_expressed=0.1,
            perc_positive=0.5,
            num_gene_sets=100,
            min_gs_size=10,
            max_gs_size=100
    ):
        paths = SynGSEAFilePath()
        self.path = paths.base_dir
        self.output_path = paths.output_dir

        self.perc_redundant = perc_redundant
        self.num_patients = num_patients
        self.total_genes = total_genes
        self.expressed_genes = expressed_genes
        self.perc_diff_expressed = perc_diff_expressed
        self.perc_positive = perc_positive
        self.num_gene_sets = num_gene_sets

        self.min_gs_size = min_gs_size
        self.max_gs_size = max_gs_size

        assert self.perc_redundant >= 0.
        assert self.perc_redundant <= 1.
        assert self.num_patients > 0
        assert self.total_genes > 0
        assert self.expressed_genes <= self.total_genes
        assert self.perc_diff_expressed >= 0.
        assert self.perc_diff_expressed <= 1.
        assert self.perc_positive >= 0.1
        assert self.perc_positive <= 0.9
        assert self.num_gene_sets > 0
        assert self.min_gs_size >= 5
        assert self.max_gs_size <= 200

        self.synthetic_genes = []
        self.identifier_to_gene_dict = dict()
        self.gene_to_identifier_dict = defaultdict(list)
        self.expression_profiles_dict = dict()

        self.sys_noise_generator = NormalNoiseGenerator(0., 25.)
        self.exp_noise_generator = NormalNoiseGenerator(0., 50.)
        self.fold_diff_generator = PoissonNoiseGenerator(5.)

        self.name_pos_class = 'POS'
        self.name_neg_class = 'NEG'
        self.num_pos_class = int(self.perc_positive * self.num_patients)
        self.num_neg_class = self.num_patients - self.num_pos_class

    def _generate_cls_file(self, output_file, classes, samples):
        """
        Generate a GSEA cls file
        :param output_file: name of output file
        :param classes: list of class names
        :param samples: number of samples per class in order given in classes
        :return:
        """
        assert not os.path.exists(output_file)
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

    def _generate_gct_file(self, output_file, expression_matrix):
        """
        Generate a GSEA gct file
        :param output_file: name of output file
        :param expression_matrix: dataframe containing expression levels
        :return:
        """
        assert not os.path.exists(output_file)

        total_genes = len(expression_matrix)
        total_samples = len(expression_matrix.columns) - 2

        with open(output_file, 'w') as f:
            f.write('#1.2\n')
            f.write('{} {}\n'.format(total_genes, total_samples))
            f.write('\t'.join(list(expression_matrix.keys())))
            f.write('\n')
            for index, row in expression_matrix.iterrows():
                values = row.tolist()
                f.write('\t'.join(values[0:2]))
                f.write('\t')
                values = ["{0:.2f}".format(v) for v in values[2:]]
                f.write('\t'.join(values))
                f.write('\n')

        return output_file

    def _generate_gmt_file(self, output_file, gene_sets):
        """
        Generate a GSEA gmt file
        :param output_file: name of output file
        :param gene_sets: list of gene sets; each entry is (gene_set_name, gene_set_origin, list of gene symbols)
        :return:
        """
        assert not os.path.exists(output_file)

        with open(output_file, 'w') as f:
            for gs_name, gs_origin, symbols in gene_sets:
                f.write('{}\t{}\t{}\n'.format(
                    gs_name,
                    gs_origin,
                    '\t'.join(symbols)
                ))

        return output_file

    def _generate_synthetic_genes(self):
        """
        Generate synthetic genes and mapping dicts
        :return:
        """
        # generate all genes and add initial identifier to id dict
        for i in range(self.total_genes):
            gene_id = 'GENE' + str(i + 1)
            uid = 'UID' + str(i + 1)
            self.synthetic_genes.append(gene_id)
            self.identifier_to_gene_dict[uid] = gene_id
            self.gene_to_identifier_dict[gene_id].append(uid)

        # sample random percent for redundant identifiers
        num_redundant = int(self.perc_redundant * self.total_genes)
        redundant_sample = random.sample(self.synthetic_genes, num_redundant)

        # add redundant identifiers
        max_id = self.total_genes + 1
        for r in redundant_sample:
            uid = 'UID' + str(max_id)
            self.identifier_to_gene_dict[uid] = r
            self.gene_to_identifier_dict[r].append(uid)
            max_id += 1

        return

    def _generate_expression_profile(self):
        """
        Generate synthetic expression profile for normal and positive phenotypes
        :param norm_list: list of not differentially expressed genes
        :param diff_list: list of differentially expressed genes
        :return:
        """
        # initialize
        detected_subset = random.sample(self.synthetic_genes, self.expressed_genes)

        # convert to gene identifiers
        detected_gene_ids = [random.sample(
            self.gene_to_identifier_dict[gene], 1
        )[0] for gene in detected_subset]

        # initialize expression dict
        neg_exp = dict()
        for g_id in detected_gene_ids:
            neg_exp[g_id] = next(self.exp_noise_generator)

        # generate differential expression profile
        pos_exp = dict(neg_exp)
        diff_exp_subset = random.sample(detected_gene_ids, int(self.perc_diff_expressed * len(detected_subset)))
        for g_id in diff_exp_subset:
            pos_exp[g_id] = next(self.fold_diff_generator) * neg_exp[g_id]

        print('Differentially expressed: ')
        print(diff_exp_subset)

        return pos_exp, neg_exp

    def _generate_experimental_data(self, pos_exp, neg_exp):
        """
        Generate synthetic expression data for patients
        :param pos_exp: positive expression levels
        :param neg_exp: negative expression levels
        :return:
        """
        expression_dict = dict()

        gene_names = neg_exp.keys()

        expression_dict['NAME'] = list(gene_names)
        expression_dict['DESCRIPTION'] = ['na'] * len(gene_names)

        name_tag = [self.name_pos_class] * self.num_pos_class \
            + [self.name_neg_class] * self.num_neg_class

        for i, tag in enumerate(name_tag):
            name = '{}_{}'.format(tag, i+1)
            if tag == self.name_pos_class:
                use_exp = pos_exp
            else:
                use_exp = neg_exp

            values = []
            for g_id in gene_names:
                values.append(use_exp[g_id] + next(self.sys_noise_generator))
            expression_dict[name] = values

        expression = pd.DataFrame(data=expression_dict)

        return expression

    def _generate_gene_sets(self):
        """
        Generate synthetic gene sets from synthetic genes
        :return:
        """
        g_sets = []

        for i in range(self.num_gene_sets):
            gs_name = 'GSET_{}'.format(i)
            gs_size = random.randint(self.min_gs_size, self.max_gs_size)
            genes_in_set = random.sample(self.synthetic_genes, gs_size)
            gene_ids = [random.sample(
                self.gene_to_identifier_dict[gene], 1
            )[0] for gene in genes_in_set]
            g_sets.append((gs_name, 'SYNGSEA', tuple(gene_ids)))

        return g_sets

    def _write_files(self, data, gsets):
        """
        Write cls, gct, and gmt files to output folder
        :param data: gene expression data
        :param gsets: gene set data
        :return:
        """
        output_experiment_dir = os.path.join(
            self.output_path,
            '{}-{}'.format('syndata',
                           datetime.now().strftime('%Y-%m-%d'))
        )

        try:
            # make output directory
            if not(os.path.exists(output_experiment_dir)):
                os.makedirs(output_experiment_dir)

            os.makedirs(os.path.join(output_experiment_dir, 'gsea_output'))

            # write class information to cls file
            print('Generating GSEA cls file...')
            cls_file = os.path.join(output_experiment_dir, 'gsea_exp.cls')
            self._generate_cls_file(
                cls_file,
                (self.name_pos_class, self.name_neg_class),
                (self.num_pos_class, self.num_neg_class)
            )

            # write gene expression data to gct file
            print('Generating GSEA gct file...')
            gct_file = os.path.join(output_experiment_dir, 'gsea_exp.gct')
            self._generate_gct_file(gct_file, data)

            # write gene sets to gmt file
            print('Generating GSEA gmt file...')
            gmt_file = os.path.join(output_experiment_dir, 'gsea_exp.gmt')
            self._generate_gmt_file(gmt_file, gsets)

            print('done.')

        except BaseException:
            print('Problem generating data, deleting incomplete data...')
            shutil.rmtree(output_experiment_dir)

        return

    def generate_all_files(self):
        """
        Generate all synthetic data files from input specifications
        :return:
        """
        # generate synthetic genes
        self._generate_synthetic_genes()

        # generate expression profile from gene list
        positive_expression, negative_expression = self._generate_expression_profile()

        # simulate positive/negative phenotype expression data
        expression_data = self._generate_experimental_data(
            positive_expression, negative_expression
        )

        # generate synthetic gene sets
        gene_sets = self._generate_gene_sets()

        # write all files to output dir
        self._write_files(expression_data, gene_sets)

        return

    def run_gsea(self):
        """
        Run GSEA on synthetic data
        :return:
        """
        return