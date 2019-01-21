import os
import glob
import gseapy

from syngsea.paths import GSEAFilePath


class GSEAExperiment:
    def __init__(
            self,
            data_dir: str,
            gene_set_file: str,
            output_dir: str
    ):
        self.data_dir = data_dir
        self.gmt_use = gene_set_file
        self.output_dir = output_dir
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)

    @staticmethod
    def run_gsea(
            cls_file: str,
            gct_file: str,
            gmt_file: str,
            save_dir: str
    ):
        """
        Run GSEA
        :param cls_file:
        :param gct_file:
        :param gmt_file:
        :param save_dir:
        :return:
        """
        assert os.path.exists(cls_file)
        assert os.path.exists(gct_file)
        assert os.path.exists(gmt_file)
        assert os.path.exists(save_dir)

        gseapy.gsea(
            data=gct_file,
            gene_sets=gmt_file,
            cls=cls_file,
            outdir=save_dir,
            verbose=True
        )

    def iterate_data(self):
        """
        Iterate through data and run GSEA
        :return:
        """
        cls_files = glob.glob(
            os.path.join(self.data_dir, '*.cls')
        )
        cls_files.sort()

        for cls_file in cls_files:
            fname_w_ext = os.path.basename(cls_file)
            file_name, file_ext = os.path.splitext(fname_w_ext)

            print('Running GSEA on {}...'.format(file_name))

            gct_file = os.path.join(self.data_dir, '{}.{}'.format(file_name, 'gct'))
            out_dir = os.path.join(self.output_dir, file_name)

            if not os.path.exists(out_dir):
                os.mkdir(out_dir)

            self.run_gsea(
                cls_file=cls_file,
                gct_file=gct_file,
                gmt_file=self.gmt_use,
                save_dir=out_dir
            )


if __name__ == '__main__':
    paths = GSEAFilePath()
    gset_file = os.path.join(paths.data_dir, 'gene_sets', 'c2.cp.v6.2.symbols.gmt')
    output_dir = os.path.join(paths.output_dir, 'gsea_baseline')
    exp = GSEAExperiment(
        data_dir=paths.processed_data_dir,
        gene_set_file=gset_file,
        output_dir=output_dir
    )
    exp.iterate_data()