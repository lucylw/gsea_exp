import os
import csv
import json

from collections import defaultdict
import pandas as pd

from syngsea.libs.rbo import rbo
from syngsea.paths import GSEAFilePath
import syngsea.constants as constants


class GSEAProcessor:
    def __init__(self):
        paths = GSEAFilePath()

        self.gsea_baseline_path = os.path.join(paths.output_dir, 'gsea_baseline')
        self.gsea_normalized_path = os.path.join(paths.output_dir, 'gsea_normalized')
        self.comp_output_path = os.path.join(paths.output_dir, 'comparative_eval')
        self.viz_output_path = os.path.join(paths.output_dir, 'visualization')

        self.tissues = [
            'ACT_fore',
            'ACT_hippo',
            'ACT_p_neo',
            'ACT_t_neo',
            'AMP_BM10',
            'AMP_BM22',
            'AMP_BM36',
            'AMP_BM44',
            'LUAD',
            'HNSCC'
        ]

        self.pw_root = ['PW_0000002', 'PW_0000013', 'PW_0000754', 'PW_0000004', 'PW_0000003']

        pw_path = os.path.join(
            paths.base_dir,
            '..', 'pathhier',
            'data',
            'pathway_ontology',
            'pw.json'
        )
        self.pw_tree = self.construct_pw_tree(pw_path)

    def construct_pw_tree(self, pw_path):
        """
        Construct a tree from PW classes
        :return:
        """
        pw_ont = json.load(open(pw_path, 'r'))
        pw_ont = {k.split('/')[-1]: v for k, v in pw_ont.items()}

        for k, v in pw_ont.items():
            if v['subClassOf']:
                v['subClassOf'] = [uid.split('/')[-1] for uid in v['subClassOf']]
                v['subClassOf'] = [uid for uid in v['subClassOf'] if uid.startswith('PW')]
            if v['part_of']:
                v['part_of'] = [uid.split('/')[-1] for uid in v['part_of']]
                v['part_of'] = [uid for uid in v['part_of'] if uid.startswith('PW')]

        return pw_ont

    def add_all_parents(self, ns, ls, start_id, ont):
        """
        Add all parent nodes and edges to graph
        :param ns: nodes
        :param ls: links
        :param start_id:
        :param ont:
        :return:
        """
        parents = ont[start_id]['subClassOf']
        for par_id in parents:
            if par_id not in ns and par_id in ont:
                par_name = ont[par_id]['name'].replace(' ', '_').upper()
                ns[par_id] = (par_name, -1.0)
            if (par_id, start_id) not in ls:
                ls.add((par_id, start_id))
            ns, ls = self.add_all_parents(ns, ls, par_id, ont)

        return ns, ls

    def construct_child_tree(self, current_nodes, ns, ls):
        """
        Construct a json tree
        :param current_node:
        :param ns:
        :param ls:
        :return:
        """
        node_list = []

        for node in current_nodes:
            if node in ns:
                parents = [par for par, cld in ls if cld == node]
                if len(parents) == 0:
                    parent_id = "PW_0000001"
                else:
                    parent_id = parents[0]

                children = [cld for par, cld in ls if par == node]

                node_list.append({
                    "uid": node,
                    "name": ns[node][0],
                    "score": ns[node][1],
                    "parent": parent_id,
                    "children": self.construct_child_tree(children, ns, ls)
                })

        return node_list

    def generate_json_file(self, pw_paths, pw, outfile):
        """
        Extract PW class names and generate json file for visualization
        :param pw_paths:
        :param pw:
        :param outfile:
        :return:
        """
        nodes = dict()
        links = set()

        for es, path_name in zip(pw_paths.ES_normalized, pw_paths.GS_normalized):
            path_parts = path_name.split('_')
            pw_id = '_'.join(path_parts[:2])
            if pw_id in pw:
                pw_name = pw[pw_id]['name'].replace(' ', '_').upper()
                nodes[pw_id] = (pw_name, es)
                nodes, links = self.add_all_parents(nodes, links, pw_id, pw)

        root_dict = [
            {
                "name": "PATHWAY",
                "score": -1.0,
                "parent": "null",
                "children": self.construct_child_tree(self.pw_root, nodes, links)
            }
        ]

        # dump json to file
        json.dump(
            root_dict,
            open(outfile, 'w'),
            sort_keys=True,
            indent=4
        )
        print('Wrote graph to json.')

    def process_results_file(self, f_path):
        """
        Process GSEA summary results file
        :param f_path:
        :return:
        """
        keep_cols = ['Term', 'es', 'nes', 'pval', 'ledge_genes']
        new_col_names = ['GS', 'ES', 'NES', 'pval', 'ledge']

        df = pd.read_csv(f_path, sep=',', header=0)
        df = df[keep_cols]
        df.columns = new_col_names
        df = df.sort_values(by=['NES'], ascending=False)
        df = df[:constants.KEEP_TOP_N_GSEA_RESULTS].filter(new_col_names)
        df = df.reset_index(drop=True)

        return df

    def process_r_results(self, f_path):
        """
        Process GSEA summary results from R output
        :param f_path:
        :return:
        """
        keep_cols = ['GS', 'ES', 'NES', 'p-val']
        new_col_names = ['GS', 'ES', 'NES', 'pval']

        df = pd.read_csv(f_path, sep='\t', header=0)
        df = df[keep_cols]
        df.columns = new_col_names
        df = df.sort_values(by=['ES'], ascending=False)
        df = df[:constants.KEEP_TOP_N_GSEA_RESULTS].filter(keep_cols)
        df = df.reset_index(drop=True)

        return df

    def _count_leading_edge(self, ledge_genes):
        """
        Count leading edge genes
        :param ledge_genes:
        :return:
        """
        gene_dict = defaultdict(int)

        for gene_list in ledge_genes:
            genes = gene_list.split(';')
            for gene in genes:
                gene_dict[gene] += 1

        return gene_dict

    def compute_ledge_rbo(self, dat1, dat2):
        """
        Compute rank biased overlap between two datasets
        :param dat1:
        :param dat2:
        :return:
        """
        g1_genes = self._count_leading_edge(dat1.ledge)
        g2_genes = self._count_leading_edge(dat2.ledge)

        g1_ordered = [(k, v) for k, v in g1_genes.items()]
        g1_ordered.sort(key=lambda x: x[1], reverse=True)
        g1_ordered = [entry[0] for entry in g1_ordered]

        g2_ordered = [(k, v) for k, v in g2_genes.items()]
        g2_ordered.sort(key=lambda x: x[1], reverse=True)
        g2_ordered = [entry[0] for entry in g2_ordered]

        rbo_val = rbo(g1_ordered, g2_ordered, p=0.9)
        jaccard = len(set(g1_ordered).intersection(set(g2_ordered))) / len(set(g1_ordered).union(set(g2_ordered)))

        print('Baseline: {}...'.format(', '.join(g1_ordered[:10])))
        print('Normalized: {}...'.format(', '.join(g2_ordered[:10])))
        print('RBO: {:.3f}'.format(rbo_val['ext']))
        print('Jaccard: {:.3f}'.format(jaccard))

    def read_r_ledge_file(self, ledge_file):
        """
        Read data from R leading edge file
        :param ledge_file:
        :return:
        """
        with open(ledge_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader)
            next(reader)
            header = next(reader)
            genes = header[2:]
            results = defaultdict(int)

            for l in reader:
                for gene, in_le in zip(genes, l[2:]):
                    if in_le == '1':
                        results[gene] += 1

        results_ordered = [(k, v) for k, v in results.items()]
        results_ordered.sort(key=lambda x: x[1], reverse=True)
        results_ordered = [entry[0] for entry in results_ordered]

        return results_ordered

    def compute_r_ledge_rbo(self, ledge_file1, ledge_file2):
        """
        Compute leading edge stats for files generated in R
        :param ledge_file1:
        :param ledge_file2:
        :return:
        """
        ledge1 = self.read_r_ledge_file(ledge_file1)
        ledge2 = self.read_r_ledge_file(ledge_file2)

        rbo_val = rbo(ledge1, ledge2, p=0.9)
        jaccard = len(set(ledge1).intersection(set(ledge2))) / len(set(ledge1).union(set(ledge2)))

        print('Baseline: {}...'.format(', '.join(ledge1[:10])))
        print('Normalized: {}...'.format(', '.join(ledge2[:10])))
        print('RBO: {:.3f}'.format(rbo_val['ext']))
        print('Jaccard: {:.3f}'.format(jaccard))

    def process_gsea_output(self):
        """
        Compare GSEA output for each tissue
        :return:
        """
        for t_name in self.tissues:
            print(t_name)

            # AMP files were processed with R
            if t_name.startswith('AMP'):
                baseline_report = os.path.join(self.gsea_baseline_path, t_name, 'GSEA.analysis.SUMMARY.RESULTS.REPORT.AD.txt')
                baseline_data = self.process_r_results(baseline_report)

                normalized_report = os.path.join(self.gsea_normalized_path, t_name, 'GSEA.analysis.SUMMARY.RESULTS.REPORT.AD.txt')
                normalized_data = self.process_r_results(normalized_report)

                baseline_ledge_file = os.path.join(self.gsea_baseline_path, t_name, 'GSEAanalysis.leading.genes.AD.gct')
                normalized_ledge_file = os.path.join(self.gsea_baseline_path, t_name, 'GSEAanalysis.leading.genes.AD.gct')
                self.compute_r_ledge_rbo(baseline_ledge_file, normalized_ledge_file)

            # Remainder were processed with gseapy
            else:
                baseline_report = os.path.join(self.gsea_baseline_path, t_name, 'gseapy.gsea.gene_set.report.csv')
                baseline_data = self.process_results_file(baseline_report)

                normalized_report = os.path.join(self.gsea_normalized_path, t_name, 'gseapy.gsea.gene_set.report.csv')
                normalized_data = self.process_results_file(normalized_report)

                self.compute_ledge_rbo(baseline_data, normalized_data)

            baseline_data.columns = [str(col) + '_baseline' for col in baseline_data.columns]
            normalized_data.columns = [str(col) + '_normalized' for col in normalized_data.columns]

            combined = pd.merge(baseline_data, normalized_data, left_index=True, right_index=True)
            keep_cols = ['NES_baseline', 'GS_baseline', 'NES_normalized', 'GS_normalized']
            combined = combined.filter(keep_cols)

            output_file = os.path.join(self.comp_output_path, 'comparative_{}.tsv'.format(t_name))
            header_row = ['Baseline_NES', 'Baseline_GS', 'Normalized_NES', 'Normalized_GS']
            combined.index += 1

            combined.to_csv(
                open(output_file, 'w'),
                sep='\t',
                float_format='%.3f',
                header=header_row,
                index=True
            )

            output_json_file = os.path.join(self.viz_output_path, '{}.json'.format(t_name))
            self.generate_json_file(normalized_data, self.pw_tree, output_json_file)


if __name__ == '__main__':
    data_processor = GSEAProcessor()
    data_processor.process_gsea_output()


