import os, sys, argparse
from os import listdir, path
from os.path import isfile, join, isdir
import numpy as np
import requests
import Bio
from Bio.PDB import PDBList
import pandas as pd
import urllib.request

def read_file(file):
  in_file = open(file, 'r')

  header_title = in_file.readline().strip()

  new_lines = []

  for line in in_file:
    line = line.strip().upper()
    new_lines.append(line)

  return header_title, new_lines

def download_pdb(pdb_lst, args, cd):
  missing_pdbs = []

  for i in range(len(pdb_lst)):
    print(i, pdb_lst[i])
    
    try:
      url = 'https://files.rcsb.org/download/' + pdb_lst[i].lower() + '.pdb'
      pdb = pdb_lst[i].lower() + '.pdb'
      urllib.request.urlretrieve(url, pdb)
    except:
      print("Could not download PDB: " + pdb_lst[i])
      if pdb_lst[i] not in missing_pdbs:
        missing_pdbs.append(pdb_lst[i])
      continue

  # Dataframe for PDBs that need to be manually downloaded as cif or other file formats
  missing_df = pd.DataFrame()
  missing_df['pdb'] = missing_pdbs

  os.chdir(cd)
  missing_df.to_csv('unfinished_downloaded_pdbs.txt', index = False, sep = ' ')

# python3 download_pdbs.py -idir 'pdb_files' -i 'pdbs_to_download.txt' -pdbc 'pdbid'
def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-idir', '--input_directory', dest='idir', nargs='?', default="pdb_files", type=str, help='default pdb_files, options: pdb_ent_files')
  parser.add_argument('-i', '--input', dest='i', nargs='?', default="pdb_files.txt", type=str, help='default pdb_files.txt')
  parser.add_argument('-pdbc', '--pdb_column', dest='pdbc', nargs='?', default="PDB", type=str, help='default PDB; options pdbid')
  parser.add_argument('-wo', '--write_outfile', dest='wo', nargs='?', default="True", type=str, help='default True')
  args = parser.parse_args()

  # Move into data folder 
  os.chdir('data')
  cd = os.getcwd()

  # Read Input CSV with column of individual PDB ids ['6MHC', '39P0']
  pdb_df = pd.read_csv(args.i, sep = '\t')
  pdbs = list(pdb_df[args.pdbc].unique())
  print("Total PDBs found in input file: " + str(len(pdbs)))

  # Check PDB directory exists
  path = cd + '/' + args.idir
  os.makedirs(path, exist_ok = True)
 
  # Check which PDB files have already been downloaded
  os.chdir(path) 
  files = [f for f in listdir(os.getcwd()) if isfile(join(os.getcwd(), f))]

  # Iterate through input PDBs to only download those which have not already been downloaded
  unfinished = []
  for i in range(len(pdbs)):
    filename = pdbs[i].lower() + '.pdb'
    if (filename not in files) & (filename not in unfinished):
      unfinished.append(pdbs[i])

  print("Total PDBs from input file that need to be downloaded: " + str(len(unfinished)))
  
  if args.wo == "True":
    os.chdir(path) 
    # Download PDBs via https://files.rcsb.org/download/4hhb.pdb
    download_pdb(unfinished, args, cd)

if __name__ == "__main__":
  main()
