import sys
import numpy as np
import pandas as pd
import re
import subprocess
from pathlib import Path
import fileinput
from biopandas.pdb import PandasPdb
import os.path

def main(pdbfile, outfile):

    ppdb = PandasPdb()
    ppdb.read_pdb(pdbfile)

    ppdb.df['ATOM']['chain_id'] = 'A'
    idx = int(len(ppdb.df['ATOM'])+1)
    remap = pd.DataFrame(np.arange(1, idx), columns=['atom_number'])
    ppdb.df['ATOM'].loc[(ppdb.df['ATOM'].chain_id == 'A'), 'atom_number'] = remap.atom_number

    ppdb.to_pdb(outfile, records=['ATOM'])
    if not os.path.exists('complex.fasta'):
        seq = ppdb.amino3to1()
        seq = str(''.join(seq.loc[seq['chain_id'] == 'A', 'residue_name']))
        with open('complex.fasta', 'w') as f:
            f.write('>complex'+'\n')
            f.write(seq)

if __name__ == "__main__":
    main(pdbfile=sys.argv[1],
         outfile=sys.argv[2])