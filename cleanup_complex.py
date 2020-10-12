import sys
import numpy as np
import pandas as pd
import re
import subprocess
from pathlib import Path
import fileinput
from biopandas.pdb import PandasPdb

def main(pdbfile, epfile, outfile):
    
    ppdb = PandasPdb()
    ppdb.read_pdb(pdbfile)
    ppdb2= PandasPdb()
    ppdb2.read_pdb(epfile)
   
    # ppdb.df['OTHERS'] = ppdb.df['OTHERS'].copy()
    # hetdf = ppdb2.df['OTHERS'][ppdb2.df['OTHERS'].record_name == "HETNAM"].copy()
    # hetdf.line_idx = np.arange(1, len(hetdf)+1)
    # ppdb.df['OTHERS'].line_idx = ppdb.df['OTHERS'].line_idx + len(hetdf)
    # ppdb.df['HETATM'].line_idx = ppdb.df['HETATM'].line_idx + len(hetdf)
    # ppdb.df['ATOM'].line_idx = ppdb.df['ATOM'].line_idx + len(hetdf)
    # ppdb.df['OTHERS'] = hetdf.append(ppdb.df['OTHERS'])

    df =ppdb.df['ATOM'].copy()
    hdf = df[(df['chain_id'] == 'H')].reset_index(drop=True)
    ldf = df[(df['chain_id'] == 'L')].reset_index(drop=True)
    adf = df[(df['chain_id'] == 'A')].reset_index(drop=True)

    idx = int(len(ppdb.df['ATOM'])+2)

    remap = pd.DataFrame(np.arange(2, idx), columns=['atom_number'])
    ppdb.df['ATOM'].loc[(ppdb.df['ATOM'].chain_id == 'A'), 'atom_number'] = remap.atom_number

    ppdb.to_pdb(outfile, records=None)     

if __name__ == "__main__":
    main(pdbfile="//media/hdd1/XGmAbBoost/5vkj_CD22_glycan/inotuzumab_5vkj_gR/megadock/complex_0.pdb",
         epfile="/media/hdd1/XGmAbBoost/5vkj_CD22_glycan/inotuzumab_5vkj_gR/megadock/5vkj_gR_1.pdb" )


