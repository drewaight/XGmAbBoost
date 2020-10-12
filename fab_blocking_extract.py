from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.protocols.rosetta_scripts import *
from pyrosetta.rosetta.protocols.antibody import *
from pyrosetta.rosetta.protocols.antibody.design import *
from pyrosetta.rosetta.core.scoring import *
import pandas_anarci_numberer
from pyrosetta.rosetta.utility import *
import os
import sys
import timeit
import numpy as np
import pandas as pd
import re
from pathlib import Path

def main(path, structure, light_chain):
    init('-no_fconfig @/media/hdd1/XGmAbBoost/drewdock_exec/common')
    pose = pose_from_pdb(os.path.join(path, structure))
    scorefxn = get_score_function()
    original_pose = pose.clone()

    abinfo = AntibodyInfo(original_pose, AHO_Scheme, North)
    if light_chain == "kappa":
        abinfo.set_light_chain_type(kappa)
    elif light_chain == "lambda":
        abinfo.set_light_chain_type(LightChainTypeEnum_start)
    print(abinfo)    

    pdbinfo = pose.pdb_info()
    chains = ("H","L")
    for res in range(1, pose.total_residue()+1):
        for atom in range(1, pdbinfo.natoms(res)+1):
            pdbinfo.temperature(res,atom, 0.0)
                      
    for num in range(1,3):        
        for res in range(pose.conformation().chain_begin(num), pose.conformation().chain_begin(num)+4):
            for atom in range(1, pdbinfo.natoms(res)+1):
                pdbinfo.temperature(res,atom, 1.0)
        
    for chain in chains:
        for res in range(pdbinfo.pdb2pose(chain, 73), pdbinfo.pdb2pose(chain, 93)):
            for atom in range(1, pdbinfo.natoms(res)+1):
                pdbinfo.temperature(res,atom, 1.0)
    
    for i in range(1, 7):
        CDR_range = abinfo.get_CDR_start(CDRNameEnum(i), pose), abinfo.get_CDR_end(CDRNameEnum(i), pose)
        for r in range(CDR_range[0]-3, CDR_range[1]+4):
            for atom in range(1, pdbinfo.natoms(r)+1):
                pdbinfo.temperature(r,atom, 1.0)

    pose.dump_pdb(os.path.join(path, structure))
    
    list = []
    for i in range(1, 7):
        list.append(abinfo.get_cluster_name(abinfo.get_CDR_cluster(CDRNameEnum(i)).cluster()))
        list.append(abinfo.get_CDR_cluster(CDRNameEnum(i)).normalized_distance_in_degrees())
        list.append(abinfo.get_CDR_sequence_with_stem(CDRNameEnum(i), pose, 3, 3))
    abdf_T = pd.DataFrame(list)
    abdf = pd.DataFrame(abdf_T.T)
    abdf.columns=['H1_cluster', 'H1_distance', 'H1_sequence', 'H2_cluster', 'H2_distance', 'H2_sequence',
                'H3_cluster', 'H3_distance', 'H3_sequence', 'L1_cluster', 'L1_distance', 'L1_sequence',
                'L2_cluster', 'L2_distance', 'L2_sequence', 'L3_cluster', 'L3_distance', 'L3_sequence']
        
    print(abdf)
    return(abdf)

if __name__ == '__main__':
    main(path='/media/hdd1/XGmAbBoost/daratumumab7/daramodels', structure='model_341_886.pdb',light_chain='kappa')

