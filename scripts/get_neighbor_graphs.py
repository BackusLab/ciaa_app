import Bio
from Bio.PDB import *
import numpy, os, sys, subprocess
from os.path import isfile, join
import sys, os, argparse, subprocess, time
import numpy as np
import pandas as pd

def read_reference_csv_file(file):
  in_file = open(file, 'r')
  pdb_list = []

  for line in in_file:
    line = line.strip().split(',')
    pdb = line[0]

    # reformat pdb ids
    if '+' in pdb:
      new_pdb = str(pdb[0]) + 'E' + str(pdb[-2:])
    elif '-' in pdb:
      new_pdb = str(pdb[0]) + pdb[-3:].upper()
    else:
      new_pdb = pdb

    # update pdb dictionary
    if new_pdb not in pdb_list:
      pdb_list.append(new_pdb)

  in_file.close()
  return pdb_list

def read_reference_file(file):
  in_file = open(file, 'r')
  pdb_list = []

  for line in in_file:
    line = line.strip()
    pdb = line.upper()

    if pdb not in pdb_list:
      pdb_list.append(pdb)

  in_file.close()
  return pdb_list

def get_pdb_details(pdb_files, rd, nd, aa, sa, at, dist, wo):
  count = 0

  for i in range(len(pdb_files)):

    pdb = pdb_files[i]
    new_pdb = pdb.replace('reduced_', '').replace('.pdb', '').upper()

    parser = Bio.PDB.PDBParser(QUIET=True) 
    try: 
      structures = parser.get_structure(pdb_files[i][:-4], pdb)
    except:
      print('Error', pdb)
      continue

    if len(structures) > 1:
      print("multiple structures")
      print(pdb, len(structures))

    try: 
      model = structures[0]
    except KeyError:
      print('KeyError', pdb)
      continue

    chain_list = Selection.unfold_entities(model, "C")
    
    master_list = []
    count = 0
    for chain in chain_list:
      chainid = chain.get_id()
      res_dict = {}
      aa_list = []
      residues = chain.get_residues()

      for residue in residues:
        segid = list(residue.id)[1]
        res_dict[segid] = residue.get_resname()
        if (residue.get_resname() == aa):
          aa_list.append(segid)
          temp_list = get_neighbors(new_pdb, model, chainid, segid, sa, at, dist)
          master_list = master_list + temp_list

    header = 'identifier,pdb,resname,segid,atomid,vector,distance'

    os.makedirs(nd, exist_ok = True)
    os.chdir(nd)
    if wo == "True":
      write_file(pdb_files[i][:-4] + '_neighbors.csv', header, master_list)
    os.chdir(rd)

def angle_between_vectors(v1, v2):
    dot_product = np.dot(v1, v2)
    
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    
    cosine_angle = dot_product / (norm_v1 * norm_v2)
    
    cosine_angle = np.clip(cosine_angle, -1, 1)
    
    angle = np.arccos(cosine_angle)
    angle_degrees = np.degrees(angle)
    
    return angle_degrees

def get_neighbors(pdb, model, chainid, segid, sa, atomid, dist):

  try: 
    target_atom = model[chainid][segid][atomid]
  except KeyError:
    return []
  
  atoms  = Bio.PDB.Selection.unfold_entities(model, "A")
  ns = Bio.PDB.NeighborSearch(atoms)

  main_identifier = str(pdb) + '_' + str(chainid) + '_' + sa +  str(segid)

  close_atoms_details = []
  target_vector =  str(target_atom.get_vector()).replace(',', ' ').replace("Vector ", "").strip()
  close_atoms_details.append([main_identifier, pdb, 'CYS', segid, 'SG', target_vector[1:-1], 0, 0])
  
  close_atoms = ns.search(target_atom.coord, float(dist))
  close_atoms.remove(target_atom)

  for atom in close_atoms:
    dist = numpy.linalg.norm(target_atom.get_vector() - atom.get_vector())
    angle = angle_between_vectors(target_atom.get_vector(), atom.get_vector())
    
    # reformat entry
    identifier = str(pdb) + '_' + str(chainid) + '_' + sa +  str(segid)
    combo_vector = str(atom.get_vector()).replace(',', ' ').replace("Vector ", "").strip()
    combo_parent = str(atom.get_parent())[9:-9]
    residue = combo_parent[:3].strip()
    resseq = combo_parent.split("=")[-1]

    details = [identifier, pdb, residue, resseq, atom.get_name(), combo_vector[1:-1], dist, angle]
    if details not in close_atoms_details:
      close_atoms_details.append(details)
  return close_atoms_details

def write_file(file, header, close_atoms):
  out_file = open(file, 'w')

  for i in range(len(close_atoms)):
    current = close_atoms[i]
    st = ""
    for j in range(len(current)):
      st += str(current[j]) + ','
    out_file.write(st[:-1] + '\n')

  out_file.close()

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--input_pdb_file', dest='i', nargs='?', default="pdb_files.txt", type=str, help='default pdb_files.txt')
  parser.add_argument('-pc', '--protein_column', dest='pc', nargs='?', default="PDB", type=str, help='default 0; options = pdbid')
  parser.add_argument('-aa', '--amino_acid', dest='aa', nargs='?', default="CYS", type=str, help='default CYS')
  parser.add_argument('-sa', '--single_amino_acid', dest='sa', nargs='?', default="C", type=str, help='default C')
  parser.add_argument('-at', '--atom', dest='at', nargs='?', default="SG", type=str, help='default SG')
  parser.add_argument('-d', '--distance', dest='d', nargs='?', default="10", type=str, help='default 10; options = 5')
  parser.add_argument('-o', '--output', dest='o', nargs='?', default="output", type=str, help='default output')
  parser.add_argument('-wo', '--write_outfile', dest='wo', nargs='?', default="True", type=str, help='default True')
  args = parser.parse_args()

  cd = os.getcwd()

  pdb_df = pd.read_csv('data/' + args.i, sep = '\t')
  pdb_list = list(pdb_df[args.pc].unique())
  print(len(pdb_list))
  
  rd = cd + '/data/reduced_files'
  try:
      os.chdir(rd)
  except Exception as e:
      raise RuntimeError("Failed to change directory: {}".format(e))
      sys.exit()

  os.chdir(rd)
  reduced_files = [f for f in os.listdir('.') if f.startswith("reduced_")]

  os.chdir(cd)
  nd = cd + '/data/neighbor_files_' + args.aa.lower()

  try:
    os.mkdir(nd)
    print("Directory " , nd ,  " Created ") 
  except FileExistsError:
    print("Directory " , nd ,  " already exists")

  os.chdir(nd)
  completed_files = [f for f in os.listdir('.') if f.endswith(".csv")]
  print(len(completed_files))
  
  os.chdir(rd)
  target_pdbs = []

  for i in range(len(pdb_list)):
    reduced_pdb = 'reduced_' + pdb_list[i].lower() + '.pdb'
    if (reduced_pdb in reduced_files) & (reduced_pdb.replace('.pdb', '_neighbors.csv') not in completed_files):
      target_pdbs.append(reduced_pdb)

  print("Number of target pdbs: " + str(len(target_pdbs)))

  get_pdb_details(target_pdbs, rd, nd, args.aa, args.sa, args.at, args.d, args.wo)

if __name__ == "__main__":
  main()