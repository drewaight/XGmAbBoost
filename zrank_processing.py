import numpy as np
import pandas as pd
import fileinput
import sys
from pathlib import Path
import subprocess
import re

def main(megadockfile, ligfile, recfile, path, t_block):
    ligh = ligfile+'.h'
    rech = recfile+'.h'
    with open(ligh, 'w') as redlig:
        subprocess.run(['reduce', ligfile], stdout=redlig)
    with open(rech, 'w') as redrec:
        subprocess.run(['reduce', recfile], stdout=redrec)

    md_recname = ''.join([Path(recfile).stem, '_blocking.pdb'])
    if t_block:
        md_ligname = ''.join([Path(ligfile).stem, '_blocking.pdb'])
    else:
        md_ligname = ligfile
    
    for line in fileinput.input(megadockfile, inplace=1):
        if md_ligname in line:
            line = line.replace(md_ligname, ligh)
        if md_recname in line:
            line = line.replace(md_recname, rech)
        sys.stdout.write(line)

    subprocess.run(['mpiexec','-n','16', 'python', path+'zrank_mpi.py', megadockfile, path, '12000'])
    # subprocess.run(['mpiexec','-n','16', 'python', path+'zrank_mpi.py', megadockfile, path, '2000'])

    zrankfile = Path(megadockfile).stem+'_unsorted.zr.out'
    df1 = pd.read_csv(megadockfile, header= 3, sep='\t', names=['rot1','rot2','rot3','vox1','vox2','vox3','score'])
    df2 = pd.read_csv(zrankfile, sep='\t', names=['id','zrank'])
    df1 = df1.join(df2)
    df1['id'] = df1['id'].values.astype(np.int64)
    df1['score'] = df1['score'].round(2)
    df1.sort_values(by=['zrank'], inplace=True)
    df1 = df1.reset_index(drop=True)

    return(df1)

if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

