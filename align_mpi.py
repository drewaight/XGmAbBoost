import sys
import fileinput
import concurrent.futures
import timeit
import subprocess
from random import random
from mpi4py import MPI
import numpy as np
import pandas as pd
from pathlib import Path
from communication import *
import align
import json

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def align(path, cmplx):
  start_time = timeit.default_timer()
  with concurrent.futures.ThreadPoolExecutor() as executor:
    future = executor.submit(align.main, path, cmplx)
    rms = future.result()
    executor.shutdown(wait=True)
  print('Processor %d finished'%rank)
  return(rms)

def main(path, jsonfile):
  ros_df = pd.read_json(path + '/' + jsonfile)
  print(ros_df)
  z = None
  nproc = size
  ncmplx = np.int64(len(ros_df.index))

  if rank == 0:
    x = np.arange(0, ncmplx, dtype=np.int32)
    print('Test started with %d processors.'%size)
    z = np.zeros([ncmplx, 2], dtype=np.float64)
    for x in x:
      zr = z[x]
      zr[0] = x

  z_scatt = scatter_array(z)

  for z in z_scatt:
    zi = int(z[0])
    zin = str(zi)
    loc_cmplx = ('complex_'+zin+'_relaxed.pdb')
    rms = align(path, loc_cmplx)
    z[1] = rms

  gath = np.zeros([ncmplx, 2], dtype=np.float64)
  gather_array(gath, z_scatt)
 
  if rank == 0:
    print("All Aligned")

  # MPI.Finalize()

if __name__ == '__main__':
  # main(path = '/media/hdd1/roproQ3drew/2xwt_nopack/megadock', jsonfile="relaxed_list.json", reference='complex_0_relaxed.pdb')
  main(sys.argv[1],sys.argv[2],sys.argv[3])