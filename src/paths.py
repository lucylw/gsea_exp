import os

INIT_DIR = '/Users/lwang/git/synthetic_gsea'

class SynGSEAFilePath(object):
    syngsea_source_folder = 'src'
    syngsea_output_folder = 'output'

    def __init__(self, base_dir=INIT_DIR):
        """Set self.base_dir.
        """
        self.base_dir = base_dir

    @property
    def output_dir(self):
        return os.path.join(
            self.base_dir, self.syngsea_output_folder
        )