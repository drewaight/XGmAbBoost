import re
import os
import sys
import subprocess
import numpy as np
import pandas as pd
from pathlib import Path
import fileinput
import shutil

def main(dfout, epitope, num_rank, path):
     decoygen_path = path
     num_rank = int(num_rank)
     dfout = dfout.copy()
     for i in dfout.index:
          dock_list = []
          if i <= num_rank:
               ligout_file = ''.join(['lig_', dfout.loc[i]['model_name'], '.', str(i), '.pdb'])
               dfout.loc[i, 'ligout'] = ligout_file
               subprocess.run([os.path.join(decoygen_path, 'decoygen'), ligout_file, epitope,
                              ''.join([dfout.loc[i]['model_name'], '_complex_zranked.out']),
                                   str(dfout.loc[i]['index'])], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
          else:
               continue
     outlist = []
     for i in dfout.index:
          if i <= num_rank:
               outfilename = ''.join(['complex_', str(i), '.pdb'])
               with open(outfilename, 'wb') as outFile:
                    with open(''.join([dfout.loc[i]['model_name'], '.pdb']), 'rb') as rec, \
                         open(dfout.loc[i]['ligout'], 'rb') as lig:
                              shutil.copyfileobj(rec, outFile)
                              shutil.copyfileobj(lig, outFile)
               outlist.append(outfilename)
          else:
               return(outlist)

if __name__ == "__main__":
     main(sys.argv[1],sys.argv[2],sys.argv[3])
   