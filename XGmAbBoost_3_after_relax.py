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

def main(mab, epitope, str_list, str_dir, exec_dir, ProQ3_dir, decoynum, cutoff, t_block):
    # The .theano directory gets very large RAM usage over time using ProQ3D
    # subprocess.run(['rm', '-rf', '/home/drewaight/.theano'])
    # start_time = timeit.default_timer()
    
    #Parse the input
    ab_list = []
    if str_list:
        name = str_list
        with open(str_list) as l:
            ab_list = l.read().splitlines()
    else:
        name = mab
        ab_list.append(mab)
    #Working directory where you want all the results stored
    os.chdir(str_dir)
    topdir = os.getcwd()
    # scheme = 'aho'

    # # # Make a folder and copy the target structures into it
    directory = ''.join([Path(name).stem, "_", Path(epitope).stem])
    # if os.path.exists(directory):
    #     shutil.rmtree(directory)
    # os.makedirs(directory)
    # os.chdir(str_dir)
    # cwd = os.getcwd()
    # print(cwd)
    # for ab in ab_list:
    #     shutil.copy(os.path.join(str_dir, ab), directory)
    # shutil.copy(os.path.join(str_dir, epitope), directory)
    os.chdir(os.path.join(str_dir,directory))
  
    # # The biopandas_cleanup module selects the best vL vH Epitope complex from a structure. It uses Bfactors as criteria
    # # the best complex is relabelled 'L' 'H' 'A' and the antibodies are renumbered by Anarci with the input scheme
    # # and returns dataframes with the heavy and light stats. See biopandas_cleanup.py for more information

    # epitope_clean = ''.join([Path(epitope).stem, '_1.pdb'])
    # try:
    #     cleanup_target.main(epitope, epitope_clean)
    # except:
    #     sys.exit("One of your pdb files somehow does not fit the criteria")
    
    # ab_list_clean = []
    # for ab in ab_list:
    #     struct_clean = ''.join([Path(ab).stem, '_1.pdb'])
    #     ab_list_clean.append(struct_clean)
    #     try:
    #         h_stats, l_stats = cleanup_Fab.main(ab, scheme, struct_clean)
    #     except:
    #         sys.exit("One of your pdb files somehow does not fit the criteria")

    # # RaBdaB takes IgG Kappa definition pulled here from the anarci stats
    # if re.search(r"IGK", str(l_stats['v_gene'])):
    #      light_chain = 'kappa'
    # elif re.search(r"IGL", str(l_stats['v_gene'])):
    #     light_chain = 'lambda'
    # else:
    #     light_chain = 'undefined'

    cwd = os.getcwd()
    # The program antibody_info_and_dock.prep.py calls structure_extract_relax.py which uses pyrosetta to run packing 
    # and minimization on the input structure and extracts the releveant antibody CDR statistics and Dunbrack clusters 
    # into a returned dataframe. After this the Aho numbered CDRs + 3 stem residues are set to 1.0 in the Bfactor column
    # and all other residues are set to 0.0 the antibody_info_and_dock.prep then separates the LH and A molecules 
    # randomizes the orientation and slides them back together in preparation for docking. 
    # abdf = pd.DataFrame()
    # for ab in ab_list_clean:
    #     try:
    #         df = fab_blocking_extract.main(cwd, ab, light_chain)
    #         df['model_name'] = ab
    #         abdf = abdf.append(df)
    #     except:
    #         sys.exit("Rosetta doesn't like your structure(s)")

    # abdf = abdf.reset_index(drop=True)
    abdf = pd.read_json('abdf.json')
    # abdf.to_json('abdf.json')    
    # h_stats.to_json('h_stats.json')
    # l_stats.to_json('l_stats.json')

    # At this point we want to fork the two prepared structures into the broad CDR definition above, and a targeted 
    # defintion using Parapred. The randomized LH structure is copied with the ppd extension to denote Parapred scoring
    # dock_list = []
    # for ab in ab_list_clean:
    #     dock_list.append([ab, epitope_clean])
    #     print([ab, epitope_clean])
    
    # # The bfactor_to_blocking_biopandas.py program simply converts all residues with bfactor == 0 to the residue BLK 
    # # which avoids rigid body docking on these residues (ZDOCK and Megadock) standard
    # for ab in dock_list:
    #     bfactor_to_blocking_biopandas.main(ab[0])
    #     block_cdr = Path(ab[0]).stem + "_blocking.pdb"
    #     ab.append(block_cdr)

    # # Writing the jobs table for the megadock program"

    # if t_block == None:
    #     for ab in dock_list:
    #         ab.append(ab[1])
    # else:
    #     blocked_target = ''.join([Path(epitope_clean).stem, '_blocking.pdb'])
    #     target_block.main(cwd, epitope_clean, os.path.join('..', t_block), blocked_target)
    #     for ab in dock_list:
    #         ab.append(blocked_target)
    
    # with open("SAMPLE.table", "w") as f:
    #     f.write("TITLE= sample jobs\n")
    #     f.write("PARAM= -R $1 -L $2\n")
    #     for ab in dock_list:
    #         f.write(ab[2]+"\t"+ab[3]+"\n")

    # Execute megadock GPU docking on all inputs
    # subprocess.run(['/usr/bin/mpiexec', '-n', '16', '--use-hwthread-cpus', 'megadock-gpu-dp',
    #                 '-v', '1.0', '-D', '-t', '3', '-N', '54000', '-tb', 'SAMPLE.table'])
    # subprocess.run(['/usr/bin/mpiexec', '-n', '16', '--use-hwthread-cpus', 'megadock-gpu-dp',
    #                '-tb', 'SAMPLE.table'])
    
    # Everything gets moved into the megadock folder
    # if os.path.exists('megadock'):
    #     shutil.rmtree('megadock')
    # os.makedirs('megadock')
    
    # shutil.move(epitope_clean, 'megadock')
    # if t_block != None:
    #     shutil.move(''.join([Path(epitope_clean).stem, '_blocking.pdb']), 'megadock')

    # for ab in dock_list:
    #     shutil.move(ab[0], 'megadock')
    #     shutil.move(ab[2], 'megadock')
        
    # dock_out_list = []
    # for ab in dock_list:
    #     outfile = ''.join([Path(ab[0]).stem, '_complex.out'])
    #     shutil.move(''.join([Path(ab[2]).stem, '-', Path(ab[3]).stem, '.out']), os.path.join('megadock', outfile))
    #     ab.append(outfile)
                  
    os.chdir('megadock')

    # for row in dock_list:
    #     print(row)

    # The zrank_processing.py module rescores the megadock output files using the zrank program, it calls the mpi 
    # enabled module zrank_mpi.py and returns a dataframe of the results sorted by zscore

    # df_zsc = pd.DataFrame()
    # for ab in dock_list:
    #     print(ab)
    # for ab in dock_list:   
    #     infile = ab[4]
    #     df = zrank_processing.main(infile, ab[1], ab[0], exec_dir, t_block)
    #     df = df.sort_values(by=['zrank'])
    #     df = df.reset_index(drop=True)
    #     df['index'] = df.index
    #     df.drop(df.index[0], inplace=True)
    #     df['model_name'] = Path(ab[0]).stem
    #     df_zsc = df_zsc.append(df)

    # df1 = df_zsc.sort_values(by=['zrank'])
    # df1 = df1.reset_index(drop=True)

    # df1['rot1'] = df1['rot1'].map(lambda x: '%.6f' % x)
    # df1['rot2'] = df1['rot2'].map(lambda x: '%.6f' % x)
    # df1['rot3'] = df1['rot3'].map(lambda x: '%.6f' % x)
    # df1['score'] = df1['score'].map(lambda x: '%.2f' % x)
    # df1['zrank'] = df1['zrank'].map(lambda x: '%.5f' % x)

    # Reformat the merged zranked dataframes back into megadock format
    # for ab in dock_list:
    #     name = Path(ab[0]).stem
    #     outfile = ''.join([Path(ab[0]).stem, '_complex_zranked.out'])
    #     outfp = open(outfile, 'w')
    #     with open(ab[4], 'r') as fp:
    #         for i, line in enumerate(fp):
    #             if i < 4:
    #                 outfp.write(line)
    #     df = df1[df1.model_name == name].copy()
    #     df2 = df.drop(['model_name','id','zrank', 'index'],axis=1)
    #     df2.to_csv(outfp, sep='\t', mode='a', header=False, index=False)

    # We want to keep track of the zranked stats, this data is split off as dfout, the rest gets re-written to the 
    # outfile for decoy generation
    # df3 = df1.drop(['model_name','id','zrank'],axis=1)
    # df3.to_csv(outfp, sep='\t', mode='a', header=False, index=False)
    # dfout = df1.drop(['rot1','rot2','rot3','vox1','vox2','vox3','score'],axis=1)
    # dfout['complex'] = dfout.index
    # print(dfout.head())
    # dfout.to_csv('zrank.csv', sep='\t')
    # dfout.to_json('zrank.json')
    dfout = pd.read_json('zrank.json')

    # The decoygen.py program converts the zranked list of rotational coordinates to back into pdb complexes. The number
    # of complexes returned is set by the decoynum variable
    # zinfile = outfile
    # complexes = decoygen.main(dfout, epitope_clean, decoynum, decoygen_dir)
    # complex_list = natural_sort.main(complexes)
    # with open('relax_list.txt', 'w') as r:
    #     for cmplx in complex_list:
    #         r.write(cmplx+'\n')

    cwd = os.getcwd()
    jsonfile= "relaxed_list.json"
    with open('relax_list.txt') as f:
        complex_list = [line.rstrip() for line in f]
    for cmplx in complex_list:
        print(cmplx)
    # The relax_mpi.py program uses mpi to run the drew_relax.py program over 16 cores. The output is a json file of the
    # relaxed and interface parameters of each complex
    # print(''.join(['Running Relax: ', directory]))

    # bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength).start()
    # p = subprocess.Popen(['mpiexec', '-n', '16', 'python', os.path.join(exec_dir, 'relax_mpi.py'),
    #                       'relax_list.txt', cwd], 
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
   
    # The rmsd_align_mpi.py program uses mpi to run rmsd.py which aligns all of the relaxed complexes by their Antigen 
    # chain
    print('Align:')
    subprocess.run(['mpiexec', '-n', '16', 'python', exec_dir+'align_mpi.py', cwd, jsonfile], 
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # Read back in the output json file from the above mpi programs
    ros_df = pd.read_json(os.path.join(cwd, jsonfile))
    ros_df = ros_df.drop(['dG_cross', 'dG_cross/dSASAx100', 'packstat'], axis=1)

    # We want ProQ3D to use our relaxed structures write the list for it to process
    with open('proQ3_list.txt', 'w') as p:
        for cmplx in complex_list:
            proQ_prepare.main(''.join([Path(cmplx).stem, '_relaxed.pdb']), 
                              ''.join([Path(cmplx).stem, '_relaxed_reres.pdb']))
            p.write(Path(cmplx).stem + '_relaxed_reres.pdb'+'\n')

    # ProQ3D is a big complicated program. It is called using GNU parallel through the proq3_all.sh script, if there is 
    # not yet a profile one is created in the top level working directory. All ProQ3D output is put into the run_all 
    # directory
    if os.path.exists('run_all'):
        shutil.rmtree('run_all')
    os.makedirs('run_all')
    try:
        print('ProQ3 processing')
        subprocess.run([os.path.join(ProQ3_dir,'proq3_all.sh'), './complex.fasta', 'proQ3_list.txt', 
                        os.path.join(topdir, ''.join([directory, '_profile'])), 'run_all', '16'], 
                        stdout=subprocess.PIPE)
    except:
        print('ProQ3 parallel error')
    
    os.chdir('run_all')
    ProQ_list = []
    for file in os.listdir():
        if file.endswith(".pdb"):
            ProQ_list.append(file)
    ProQ_list_sorted = natural_sort.main(ProQ_list)
 
    # The collector.py program scrapes the scores from the proq3.global files and returns a dataframe
    total_sc = collector.main()
    # Join the rosetta scores
    total_sc = total_sc.join(ros_df)
    total_sc.to_json('total_sc.json')
    # Merge the zrank data into the total scores
    total_sc = total_sc.merge(dfout, how='left', on='complex') 
    for i in abdf.index:
        abdf.loc[i,'model_name'] = Path(abdf.loc[i, 'model_name']).stem
    # Merge the antibody information the single statistics for the antibody including CDR classifications are populated
    # to all the complexes
    total_sc = total_sc.merge(abdf, how='right', on='model_name')
    # Shuffle important columns to the front
    col = total_sc.pop("zrank")
    total_sc.insert(0, col.name, col)
    col = total_sc.pop("model_name")
    total_sc.insert(0, col.name, col)
    stop_time = timeit.default_timer()
    os.chdir(os.path.join(topdir,directory))
    outfp = ''.join([directory,'_TOTAL.csv'])
    total_sc = total_sc.sort_values(by=['complex'])
    total_sc.to_csv(outfp, sep='\t', mode='a')
    total_sc.to_json(''.join([directory,'_TOTAL.json']))

    xgb_clf_file=os.path.join(exec_dir, 'final_xgb_noscale.dat')
    xgb_reg_file=os.path.join(exec_dir, 'final_reg_xgb_unscaled.dat')

    cwd = os.getcwd()
    result = structure_predictor.main(cwd, total_sc, cutoff, xgb_clf_file, xgb_reg_file)
    if result.empty:
        print("No suitable docking poses found with %3.2f cutoff. Sorry" % (cutoff))
    else:
        os.makedirs('docked')
        for i in result.index:
            shutil.copy(os.path.join('megadock', ''.join([result.iloc[i]['name'],'_relaxed.pdb'])), 
                        os.path.join(cwd, ''.join([result.iloc[i]['name'], '_docked.pdb'])))  

        result.to_csv(''.join([directory, '_SCORE.csv']), sep='\t', mode='a')
        result.to_json(''.join([directory, '_SCORE.json']))

    stop_time = timeit.default_timer()
    time_tot = str(stop_time-start_time)
    with open ('time.txt', 'w') as t:
        t.write(time_tot)
    os.chdir(str_dir)

if __name__ == '__main__':
    exec_dir, ProQ3_dir, decoygen_dir = setvars()
    main(mab= args.antibody, epitope=args.target, str_list=args.str_list,
    exec_dir = exec_dir, str_dir = args.str_dir, ProQ3_dir = ProQ3_dir,
    decoynum = args.num, cutoff= args.cutoff, t_block= args.t_block)
    