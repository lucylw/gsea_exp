import os

INIT_DIR = '/Users/lwang/git/gsea_exp'


class GSEAFilePath(object):
    source_folder = 'syngsea'
    output_folder = 'output'
    data_folder = 'data'
    gsea_folder = 'GSEA'

    def __init__(self, base_dir=INIT_DIR):
        """Set self.base_dir.
        """
        self.base_dir = base_dir

    @property
    def data_dir(self):
        return os.path.join(
            self.base_dir, self.data_folder
        )

    @property
    def processed_data_dir(self):
        return os.path.join(
            self.data_dir, 'processed_data'
        )

    @property
    def output_dir(self):
        return os.path.join(
            self.base_dir, self.output_folder
        )

    @property
    def gsea_dir(self):
        return os.path.join(
            self.base_dir, self.gsea_folder
        )