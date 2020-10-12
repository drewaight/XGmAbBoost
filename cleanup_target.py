import sys
import numpy as np
import pandas as pd
import re
import subprocess
from pathlib import Path
import fileinput
from biopandas.pdb import PandasPdb

def main(pdbfile, outfile):

    ppdb = PandasPdb()
    ppdb.read_pdb(pdbfile)

    chains_A = ppdb.df['ATOM']['chain_id'].unique()
    chain_num_A = len(chains_A)
    print("Number of Antigen chains:", chain_num_A)
    if chain_num_A > 1:
        sys.exit("Ya ok this only works with one Antigen chain")

    ppdb.df['ATOM'].chain_id = 'A'
    ppdb.df['ATOM'].b_factor = 1.0
    ppdb.df['ATOM'] = ppdb.df['ATOM'][ppdb.df['ATOM']['element_symbol'] != 'H']
    del ppdb.df['ANISOU']
    del ppdb.df['HETATM']
    ppdb.to_pdb(outfile, records=None)

if __name__ == "__main__":
    main(pdbfile="/media/hdd1/XGmAbBoost/5vkj_CD22_glycan/5vkj_g_M.pdb",
         outfile='/media/hdd1/XGmAbBoost/5vkj_CD22_glycan/5vjk_test.pdb')