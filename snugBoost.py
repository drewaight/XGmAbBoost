#########################################################################################################
#   XGmAbBoost 0.3 - XGBoost pyRosetta ProQ3D docking program 06/01/2020 Andrew Waight
#   Takes Fab/Fv and Epitope pdb files and returns most accurately docked poses, if there are any.  
#########################################################################################################
#   ENVIRONMENT VARIABLES

def setvars():

    exec_dir = '/media/hdd1/XGmAbBoost/drewdock_exec/'
    ProQ3_dir = '/home/drewaight/proq3'
    decoygen_dir = '/home/drewaight/drewdock/'

    return(exec_dir, ProQ3_dir, decoygen_dir)
#
#########################################################################################################

import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-a", "--antibody", dest = "antibody", help="Antibody structure file in PDB format")
parser.add_argument("-t", "--target", dest = "target", help="Target antigen file in PDB format")
parser.add_argument("-l", "--list", dest = "str_list", help="Antibody list of PDB files")
parser.add_argument("-n", "--number", dest = "num", default = 96, help="Number of decoys to generate")
parser.add_argument("-c", "--cutoff", dest = "cutoff", default = 0.9, help="Probability cutoff")
parser.add_argument("-d", "--directory", dest = "str_dir", default = './', help="Structure directory")
parser.add_argument("-b", "--target_block", dest = "t_block", default = None, help="Target Residue Blocking File")
parser.add_argument("-r", "--repeat", dest = "repeat", default="1", help="Number of Relax Repeats")

args = parser.parse_args()

import os
import shutil
import re
import sys
import fileinput
import subprocess
from pathlib import Path
import cleanup_Fab
import cleanup_target
import fab_blocking_extract
import pandas_anarci_numberer
import target_block
import bfactor_to_blocking_biopandas
import decoygen
import structure_predictor
import zrank_processing
import numpy as np
import pandas as pd
import chain_check
import collector
import natural_sort
import timeit
import drew_relax
import progressbar
import proQ_prepare

def main(str_dir, exec_dir, ProQ3_dir, decoynum, cutoff, t_block):
    
    # os.chdir(str_dir)
    cwd = os.getcwd()
    print(cwd)
    jsonfile= "relaxed_list.json"
    with open('relax_list.txt') as f:
        complex_list = [line.rstrip() for line in f]
  
    # The relax_mpi.py program uses mpi to run the drew_relax.py program over 16 cores. The output is a json file of the
    # relaxed and interface parameters of each complex
    # print(''.join(['Running Relax: ', cwd]))

    # bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength).start()
    # p = subprocess.Popen(['mpiexec', '-n', '16', 'python', os.path.join(exec_dir, 'relax_mpi.py'),
    #                       'relax_list.txt', cwd, '1'], 
    #                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # result = []
    # while p.stdout is not None:
    #     bar.update()
    #     line = p.stdout.readline()
    #     result.append(line.decode('UTF-8').rstrip('\r'))
    #     if not line:
    #         print("\n")
    #         p.stdout.flush()
    #         break  
    # with open("relax_mpi.log", 'w') as f:
    #     f.write(''.join(result))
   
    # print('Align:')
    # subprocess.run(['mpiexec', '-n', '16', 'python', exec_dir+'align_mpi.py', cwd, jsonfile], 
    #                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # Read back in the output json file from the above mpi programs
    ros_df = pd.read_json(os.path.join(cwd, jsonfile))
    ros_df = ros_df.drop(['dG_cross', 'dG_cross/dSASAx100', 'packstat'], axis=1)

    # We want ProQ3D to use our relaxed structures write the list for it to process
    # with open('proQ3_list.txt', 'w') as p:
    #     for cmplx in complex_list:
    #         proQ_prepare.main(''.join([Path(cmplx).stem, '_relaxed.pdb']), 
    #                           ''.join([Path(cmplx).stem, '_relaxed_reres.pdb']))
    #         p.write(Path(cmplx).stem + '_relaxed_reres.pdb'+'\n')

    # ProQ3D is a big complicated program. It is called using GNU parallel through the proq3_all.sh script, if there is 
    # not yet a profile one is created in the top level working directory. All ProQ3D output is put into the run_all 
    # directory
    # if os.path.exists('run_all'):
    #     shutil.rmtree('run_all')
    # os.makedirs('run_all')
    # try:
    #     print('ProQ3 processing')
    #     subprocess.run([os.path.join(ProQ3_dir,'proq3_all.sh'), './complex.fasta', 'proQ3_list.txt', 
    #                     os.path.join(str_dir, ''.join([str_dir, '_profile'])), 'run_all', '16'], 
    #                     stdout=subprocess.PIPE)
    # except:
    #     print('ProQ3 parallel error')
    
    # os.chdir('run_all')
    # ProQ_list = []
    # for file in os.listdir():
    #     if file.endswith(".pdb"):
    #         ProQ_list.append(file)
    # ProQ_list_sorted = natural_sort.main(ProQ_list)
 
    # The collector.py program scrapes the scores from the proq3.global files and returns a dataframe
    # total_sc = collector.main()
    # Join the rosetta scores
    # total_sc = total_sc.join(ros_df)
    # total_sc.to_json('total_sc.json')
    total_sc = pd.read_json(os.path.join('run_all','total_sc.json'))
    # Merge the zrank data into the total scores
    # total_sc = total_sc.merge(dfout, how='left', on='complex') 
    # for i in abdf.index:
    #     abdf.loc[i,'model_name'] = Path(abdf.loc[i, 'model_name']).stem
    # Merge the antibody information the single statistics for the antibody including CDR classifications are populated
    # to all the complexes
    # total_sc = total_sc.merge(abdf, how='right', on='model_name')
    # Shuffle important columns to the front
    # col = total_sc.pop("zrank")
    # total_sc.insert(0, col.name, col)
    # col = total_sc.pop("model_name")
    # total_sc.insert(0, col.name, col)
    # stop_time = timeit.default_timer()
    os.chdir(str_dir)
    outfp = 'snugdock_TOTAL.csv'
    total_sc = total_sc.sort_values(by=['complex'])
    total_sc.to_csv(outfp, sep='\t', mode='a')
    total_sc.to_json('snugdock_TOTAL.json')

    # xgb_clf_file=os.path.join(exec_dir, 'final_xgb_noscale.dat')
    # xgb_reg_file=os.path.join(exec_dir, 'final_reg_xgb_unscaled.dat')

    # cwd = os.getcwd()
    # result = structure_predictor.main(cwd, total_sc, cutoff, xgb_clf_file, xgb_reg_file)
    # if result.empty:
    #     print("No suitable docking poses found with %3.2f cutoff. Sorry" % (cutoff))
    # else:
    #     os.makedirs('docked')
    #     for i in result.index:
    #         shutil.copy(os.path.join('megadock', ''.join([result.iloc[i]['name'],'_relaxed.pdb'])), 
    #                     os.path.join(cwd, ''.join([result.iloc[i]['name'], '_docked.pdb'])))  

    #     result.to_csv(''.join([directory, '_SCORE.csv']), sep='\t', mode='a')
    #     result.to_json(''.join([directory, '_SCORE.json']))


if __name__ == '__main__':
    exec_dir, ProQ3_dir, decoygen_dir = setvars()
    main(exec_dir = exec_dir, str_dir = args.str_dir, ProQ3_dir = ProQ3_dir,
    decoynum = args.num, cutoff= args.cutoff, t_block= args.t_block)
    