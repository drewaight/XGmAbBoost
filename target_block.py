import sys
import os
import numpy as np
import pandas as pd
import re
import subprocess
from pathlib import Path
import fileinput
from biopandas.pdb import PandasPdb

def main(path, pdb, reslist, outfile):

    blist = []
    with open(os.path.join(path, reslist), 'r') as r:
        for line in r:
            line = line.strip('\n')
            resnum = line.split("-")
            resnum = [int(i) for i in resnum]
            blist.append(resnum)
    ppdb = PandasPdb()
    ppdb.read_pdb(os.path.join(path, pdb))
    for block in blist:
        if len(block) < 2:
            ppdb.df['ATOM'].loc[ppdb.df['ATOM'].residue_number == block[0], 'residue_name'] = 'BLK'
        else:
            for i in range(block[0], block[1]+1):
                ppdb.df['ATOM'].loc[ppdb.df['ATOM'].residue_number == i, 'residue_name'] = 'BLK'
    
    ppdb.to_pdb(os.path.join(path, outfile), records=None)

if __name__ == "__main__":
 main(path="/media/hdd1/XGmAbBoost/5vkj_CD22_glycan",
      pdb="5vkj_gR.pdb", reslist='5vkj_block', outfile='5vkj_gR_blocking.pdb')
