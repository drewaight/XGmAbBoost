import os
import pandas_anarci_numberer_LH
from pathlib import Path
os.chdir('/media/hdd1/XGmAbBoost/daratumumab4/snugdock/snugdock_models')
files = os.listdir()

with open('relax_list.txt') as f:
    complex_list = [line.rstrip() for line in f]
for cmplx in complex_list:
    pandas_anarci_numberer_LH.main(cmplx, "aho", cmplx)