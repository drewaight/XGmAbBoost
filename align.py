from pyrosetta import *
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.protocols.docking import *
from rosetta.protocols.rosetta_scripts import *
from pyrosetta.rosetta.protocols.simple_moves import *
import numpy as np
import pandas as pd
import os


def main(path, cmplx):
    scores = []
    init()
    pose = pose_from_pdb(os.path.join(path, cmplx))
    reference = pose_from_pdb(os.path.join(path, 'complex_0_relaxed.pdb'))
    acm = AlignChainMover()
    acm.source_chain(3)
    acm.target_chain(3)
    acm.pose(reference)
    acm.apply(pose)
    pose.dump_pdb(os.path.join(path, cmplx))

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
