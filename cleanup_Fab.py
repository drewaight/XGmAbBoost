import sys
import numpy as np
import pandas as pd
import re
import subprocess
from pathlib import Path
import fileinput
from biopandas.pdb import PandasPdb

# ----------------------------------------------------------Anarci--------------------------------------------------

def anarci(seq,scheme):
    out =subprocess.run(['anarci', '--sequence', seq, '--scheme', scheme, '--assign_germline'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = out.stdout.decode("utf-8").splitlines()
    stats = {}
    hd = []
    anarci_dict = {'chain':[], 'new_residue_number':[], 'new_insertion':[]}
    for line_i, line in enumerate(output):
        l = line.strip()
        li = line.split()
        if line.startswith('#'):
            hspl= line.split('|')
            hd.append(hspl)
        else:
            if len(li) == 3:
                if li[2] == '-':
                    continue
                else:
                    anarci_dict['chain'].append(li[0])
                    anarci_dict['new_residue_number'].append(int(li[1]))
                    anarci_dict['new_insertion'].append('')
                    residue = li[2]
                    
            elif len(li) == 4:
                anarci_dict['chain'].append(li[0])
                anarci_dict['new_residue_number'].append(int(li[1]))
                anarci_dict['new_insertion'].append(str(li[2]))
                residue = li[2]   

    stats = {'chain' : hd[5][2], 'start': hd[5][5], 'stop': hd[5][6], 'species' : hd[8][1], 
            'v_gene' : hd[8][2], 'v_id' : hd[8][3], 'j_gene' : hd[8][4], 'j_id' : hd[8][5]} 
    df_renum = pd.DataFrame.from_dict(anarci_dict)
    return(stats, df_renum)

# ------------------------------------------------Functions--------------------------------------------------

def seq_order(df):
    from collections import OrderedDict
    df['residue_insertion'] = df['residue_number'].astype(str)+df['insertion'].astype(str)
    ordered_seq = list(OrderedDict.fromkeys(df['residue_insertion']))
    seq_dict = {ordered_seq[i]: i+1 for i in range(0, len(ordered_seq))}
    df['residue_insertion'] = df['residue_insertion'].map(seq_dict)
    df['residue_number'] = df['residue_insertion']
    df.drop(['residue_insertion'], axis=1, inplace = True)
    df['insertion'] = ''
    return(df)

def atom_renum(df):
    from collections import OrderedDict
    atom_list = list(OrderedDict.fromkeys(df['atom_number']))
    atom_dict = {atom_list[i]: i+1 for i in range(0, len(atom_list))}
    df['atom_number'] = df['atom_number'].map(atom_dict)
    return(df)

# ------------------------------------------------BioPandas Cleanup--------------------------------------------------

def main(pdbfile, scheme, outfile):
    
    ppdb = PandasPdb()
    ppdb.read_pdb(pdbfile)

    chains = ppdb.df['ATOM']['chain_id'].unique()
    chain_num = len(chains)
    print("Number of Fv chains:", chain_num)
    if chain_num > 2:
        sys.exit("Ya ok this only works with one Fab or ScFv molecule")

    del ppdb.df['HETATM']
    del ppdb.df['ANISOU']
    del ppdb.df['OTHERS']

    def df_to_seq(chain):
        seq = ppdb.amino3to1()
        seq = str(''.join(seq.loc[seq['chain_id'] == chain, 'residue_name']))
        return(seq)

    chain_info = []
    df =ppdb.df['ATOM'].copy()
    ppdb.df['ATOM'].drop(ppdb.df['ATOM'].index, inplace=True)
    for chain in df['chain_id'].unique():
        cdf =  chain+'df'
        cdf = df[(df['chain_id'] == chain)].reset_index(drop=True)
        cdf = seq_order(cdf)
        ppdb.df['ATOM'] = ppdb.df['ATOM'].append(cdf).reset_index(drop=True)
        seq = df_to_seq(cdf['chain_id'].unique()[0])
        chain_info.append([chain, chain+'df', seq])

    df =ppdb.df['ATOM'].copy()

    for chain in chain_info:
        stats, num = anarci(chain[2], scheme)
        if stats['chain'] == 'H':
            print('Chain %s is a Heavy chain' % chain[0])
            df.loc[df.chain_id == chain[0], 'chain_id'] = 'H'
            hdf_stats = stats
            h_seq = chain[2]
            hdf_renum = num
            hdf = df[(df['chain_id'] == 'H')].reset_index(drop=True)
            hdf_renum.drop(['chain'], inplace = True, axis=1)
            hdf_renum['residue_number'] = np.arange(len(hdf_renum)) + 1
            hdf.drop(hdf[hdf['residue_number'] > int(hdf_stats['stop'])].index, inplace=True)
            hdf = hdf.merge(hdf_renum, how='left') \
            .drop(columns=['residue_number', 'insertion']) \
            .rename(columns={'new_residue_number': 'residue_number','new_insertion': 'insertion'})
        elif stats['chain'] == 'K':
            print('Chain %s is a Kappa Light chain' % chain[0])
            df.loc[df.chain_id == chain[0], 'chain_id'] = 'L'
            ldf_stats = stats
            l_seq = chain[2]
            ldf_renum = num
            ldf = df[(df['chain_id'] == 'L')].reset_index(drop=True)
            ldf_renum.drop(['chain'], inplace = True, axis=1)
            ldf_renum['residue_number'] = np.arange(len(ldf_renum)) + 1
            ldf.drop(ldf[ldf['residue_number'] > int(ldf_stats['stop'])].index, inplace=True)
            ldf = ldf.merge(ldf_renum, how='left') \
            .drop(columns=['residue_number', 'insertion']) \
            .rename(columns={'new_residue_number': 'residue_number','new_insertion': 'insertion'})
        elif stats['chain'] == 'L':
            print('Chain %s is a Lambda Light chain' % c_id)
            df.loc[df.chain_id == chain[0], 'chain_id'] = 'L'
            ldf_stats = stats
            l_seq = chain[2]
            ldf_renum = num
            ldf = df[(df['chain_id'] == 'L')].reset_index(drop=True)
            ldf_renum.drop(['chain'], inplace = True, axis=1)
            ldf_renum['residue_number'] = np.arange(len(ldf_renum)) + 1
            ldf.drop(ldf[ldf['residue_number'] > int(ldf_stats['stop'])].index, inplace=True)
            ldf = ldf.merge(ldf_renum, how='left') \
            .drop(columns=['residue_number', 'insertion']) \
            .rename(columns={'new_residue_number': 'residue_number','new_insertion': 'insertion'})
        else:
            sys.exit("Chain %s is not an antibody chain, sorry" % chain[0])  

    ppdb.df['ATOM'].drop(ppdb.df['ATOM'].index, inplace=True)
    hdf = hdf.append(ldf, sort=False)
    hdf = atom_renum(hdf)
    hdf.line_idx = np.arange(1, len(hdf)+1)
    ppdb.df['ATOM'] = ppdb.df['ATOM'].append(hdf, sort=False).reset_index(drop=True)

    ppdb.to_pdb(outfile, records=None)  

    h_stats = pd.DataFrame.from_records([hdf_stats])
    l_stats = pd.DataFrame.from_records([ldf_stats])
    return(h_stats, l_stats)

if __name__ == "__main__":
    main(pdbfile=sys.argv[1],
         scheme=sys.argv[2], 
         outfile=sys.argv[3])


