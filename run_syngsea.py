#!/usr/bin/python
from syngsea.synthetic_gsea import SyntheticGSEA
from subprocess import call

gsea = SyntheticGSEA(0.2)
gsea.generate_all_files()

call(["Rscript", "GSEA/run_gsea.R"])