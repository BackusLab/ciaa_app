# Date Created: 230321
# Date Modified: 240220
# Author: Lisa Boatner
# Description: Compute descriptors from PDBs for ML.

import Bio
from Bio.PDB import *
import numpy, os, sys, subprocess, argparse
from os.path import isfile, join
from Bio.PDB.DSSP import DSSP
from pathlib import Path

class Neighbor:
  def __init__ (self, resname, segid, atomid, vector, distance):
    self.resname = resname
    self.segid = segid
    self.atomid = atomid
    self.vector_st = vector.split('  ')
    self.vector_lst = [float(item) for item in self.vector_st]
    self.vector = numpy.array(self.vector_lst)
    self.distance = distance
    self.identifier = self.resname + '_' + self.segid

    self.aa_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    self.charge_dict = {'CYS': 0, 'ASP': -1, 'SER': 0, 'GLN': 0, 'LYS': 1,
     'ILE': 0, 'PRO': 0, 'THR': 0, 'PHE': 0, 'ASN': 0, 
     'GLY': 0, 'HIS': 1, 'LEU': 0, 'ARG': 1, 'TRP': 0, 
     'ALA': 0, 'VAL': 0, 'GLU': -1, 'TYR': 0, 'MET': 0}

    # van der waals radii
    # https://pubs.acs.org/doi/pdf/10.1021/j100785a001
    self.vdwr_dict = {'H': 1.2, 'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8}

    # amino acid categories
    # https://pubmed.ncbi.nlm.nih.gov/25030274/
    self.acidity_dict = {'CYS': 'P', 'ASP': 'A', 'SER': 'P', 'GLN': 'P', 'LYS': 'B',
     'ILE': 'ALI', 'PRO': 'P', 'THR': 'P', 'PHE': 'ARO', 'ASN': 'P', 
     'GLY': 'ALI', 'HIS': 'B', 'LEU': 'ALI', 'ARG': 'B', 'TRP': 'ARO', 
     'ALA': 'ALI', 'VAL':'ALI', 'GLU': 'A', 'TYR': 'ARO', 'MET': 'ALI'}

    # kyte-doolittle scale
    # http://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/650/Hydrophobicity_scales.html
    self.hydrophobicity_dict = {'CYS': 2.5, 'ASP': -3.5, 'SER': -0.8, 'GLN': -3.5, 'LYS': -3.9,
     'ILE': 4.5, 'PRO': -1.6, 'THR': -0.7, 'PHE': 2.8, 'ASN': -3.5, 
     'GLY': -0.4, 'HIS': -3.2, 'LEU': 3.8, 'ARG': -4.5, 'TRP': -0.9, 
     'ALA': 1.8, 'VAL':4.2, 'GLU': -3.5, 'TYR': -1.3, 'MET': 1.9}

    # GRaSP parameters for atomic interactions
    # https://pubmed.ncbi.nlm.nih.gov/33381849/
    self.atomic_dict = {'ALA_N': 'DON', 'ALA_CA': '-', 'ALA_C': '-', 'ALA_O': 'ACP', 'ALA_CB': 'HPB',
      'ARG_N': 'DON', 'ARG_CA': '-', 'ARG_C': '-', 'ARG_O': 'ACP', 'ARG_CB': 'HPB', 'ARG_CG': 'HPB', 'ARG_CD': '-', 'ARG_NE': 'POS,DON', 'ARG_CZ': 'POS', 'ARG_NH1': 'POS,DON', 'ARG_NH2': 'POS,DON',
      'ASN_N': 'DON', 'ASN_CA': '-', 'ASN_C': '-', 'ASN_O': 'ACP', 'ASN_CB': 'HPB', 'ASN_CG': '-', 'ASN_OD1': '-', 'ASN_ND2': 'DON',
      'ASP_N': 'DON', 'ASP_CA': '-', 'ASP_C': '-', 'ASP_O': 'ACP', 'ASP_CB': 'HPB', 'ASP_CG': 'HPB', 'ASP_OD1': 'NEG,ACP', 'ASP_OD2': 'NEG,ACP',
      'CYS_N': 'DON', 'CYS_CA': '-', 'CYS_C': '-', 'CYS_O': 'ACP', 'CYS_CB': 'HPB', 'CYS_SG': 'DON,ACP,SSB',
      'GLN_N': 'DON', 'GLN_CA': '-', 'GLN_C': '-', 'GLN_O': 'ACP', 'GLN_CB': 'HPB', 'GLN_CG': 'HPB', 'GLN_CD': '-', 'GLN_OE1': 'ACP', 'GLN_NE2': 'DON',
      'GLU_N': 'DON', 'GLU_CA': '-', 'GLU_C': '-', 'GLU_O': 'ACP', 'GLU_CB': 'HPB', 'GLU_CG': 'HPB', 'GLU_CD': '-', 'GLU_OE1': 'NEG,ACP', 'GLU_OE2': 'NEG,ACP',
      'GLY_N': 'DON', 'GLY_CA': '-', 'GLY_C': '-', 'GLY_O': 'ACP',
      'HIS_N': 'DON', 'HIS_CA': '-', 'HIS_C': '-', 'HIS_O': 'ACP', 'HIS_CB': 'HPB', 'HIS_CG': 'ARM', 'HIS_ND1': 'ARM,POS,DON,ACP', 'HIS_CD2': 'ARM', 'HIS_CE1': 'ARM', 'HIS_NE2': 'ARM,POS,DON,ACP',  
      'ILE_N': 'DON', 'ILE_CA': '-', 'ILE_C': '-', 'ILE_O': 'ACP', 'ILE_CB': 'HPB', 'ILE_CG1': 'HPB', 'ILE_CG2': 'HPB', 'ILE_CD1': 'HPB',
      'LEU_N': 'DON', 'LEU_CA': '-', 'LEU_C': '-', 'LEU_O': 'ACP', 'LEU_CB': 'HPB', 'LEU_CG': 'HPB', 'LEU_CD1': 'HPB', 'LEU_CD2': 'HPB',
      'LYS_N': 'DON', 'LYS_CA': '-', 'LYS_C': '-', 'LYS_O': 'ACP', 'LYS_CB': 'HPB', 'LYS_CG': 'HPB', 'LYS_CD': 'HPB', 'LYS_CE': '-', 'LYS_NZ': '-',
      'MET_N': 'DON', 'MET_CA': '-', 'MET_C': '-', 'MET_O': 'ACP', 'MET_CB': 'HPB', 'MET_CG': 'HPB', 'MET_SD': 'ACP', 'MET_CE': 'HPB',
      'PHE_N': 'DON', 'PHE_CA': '-', 'PHE_C': '-', 'PHE_O': 'ACP', 'PHE_CB': 'HPB', 'PHE_CG': 'HPB', 'PHE_CD1': 'HPB,ARM', 'PHE_CD2': 'HPB,ARM', 'PHE_CE1': 'HPB,ARM', 'PHE_CE2': 'HPB,ARM', 'PHE_CZ': 'HPB,ARM',
      'PRO_N': 'DON', 'PRO_CA': '-', 'PRO_C': '-', 'PRO_O': 'ACP', 'PRO_CB': 'HPB', 'PRO_CG': 'HPB', 'PRO_CD': '-',
      'SER_N': 'DON', 'SER_CA': '-', 'SER_C': '-', 'SER_O': 'ACP', 'SER_CB': '-', 'SER_OG': 'DON,ACP',
      'THR_N': 'DON', 'THR_CA': '-', 'THR_C': '-', 'THR_O': 'ACP', 'THR_CB': '-', 'THR_OG1': 'DON,ACP', 'THR_CG2': 'HPB',
      'TRP_N': 'DON', 'TRP_CA': '-', 'TRP_C': '-', 'TRP_O': 'ACP', 'TRP_CB': 'HPB', 'TRP_CG': 'HPB,ARM', 'TRP_CD1': 'ARM', 'TRP_CD2': 'HPB,ARM', 'TRP_NE1': 'ARM,DON', 'TRP_CE2': 'ARM', 'TRP_CE3': 'HPB,ARM', 'TRP_CZ2': 'HPB,ARM', 'TRP_CZ3': 'HPB,ARM', 'TRP_CH2': 'HPB,ARM',       
      'TYR_N': 'DON', 'TYR_CA': '-', 'TYR_C': '-', 'TYR_O': 'ACP', 'TYR_CB': 'HPB', 'TYR_CG': 'HPB,ARM', 'TYR_CD1': 'HPB,ARM', 'TYR_CD2': 'HPB,ARM', 'TYR_CE1': 'HPB,ARM', 'TYR_CE2': 'HPB,ARM', 'TYR_CZ': 'ARM', 'TYR_OH': 'DON,ACP',
      'VAL_N': 'DON', 'VAL_CA': '-', 'VAL_C': '-', 'VAL_O': 'ACP', 'VAL_CB': 'HPB', 'VAL_CG1': 'HPB', 'VAL_CG2': 'HPB'
      }

    self.aa = self.aa_dict[self.resname]
    self.charge = self.charge_dict[self.resname]
    self.acidity = self.acidity_dict[self.resname]
    self.vdwr = self.vdwr_dict[self.atomid[0]]
    self.hydrophobicity = self.hydrophobicity_dict[self.resname]

    # reformat to match atomic_dict keys
    if self.atomid == 'OXT':
      self.atomic_identifier = self.resname + '_O'
    else:
      self.atomic_identifier = self.resname + '_' + self.atomid

    if self.atomic_identifier in self.atomic_dict.keys():
      self.atomic_details = self.atomic_dict[self.atomic_identifier]
    else:
      self.atomic_details = ''
    
def read_neighbors_file(file, dist_cutoff):
  in_file = open(file, 'r')
  header = ['identifier', 'pdb', 'resname', 'segid', 'atomid', 'vector', 'distance']
  # header = in_file.readline().strip().split(',')
 
  neighbors_dict = {}

  # skip nucleic acids
  bp = ['A', 'C', 'G', 'U']
  
  # skip modified neighbors
  aa = ['CYS', 'ASP', 'SER', 'GLN', 'LYS', 'ILE', 'PRO', 'THR', 'PHE', 'ASN', 
  'GLY', 'HIS', 'LEU', 'ARG', 'TRP', 'ALA', 'VAL', 'GLU', 'TYR', 'MET']

  for line in in_file:
    line = line.strip().split(',')
  
    identifier = line[header.index('identifier')]

    # skip bp/nucleic acids
    if line[header.index('resname')] in bp:
      continue 

    # skip modified neighbors
    if line[header.index('resname')] not in aa:
      continue 
    
    # skip atoms farther than constraint
    if float(line[header.index('distance')]) > dist_cutoff:
      continue

    # skip atoms that are too close
    if float(line[header.index('distance')]) <= 0:
      continue

    # skip atoms in current cysteine residue
    if (line[header.index('resname')] == 'CYS') & (line[header.index('segid')] == identifier.split('_C')[-1]):
      continue

    if identifier not in neighbors_dict.keys():
      neighbor = Neighbor(line[header.index('resname')], line[header.index('segid')], line[header.index('atomid')], line[header.index('vector')], line[header.index('distance')])
      neighbors_dict[identifier] = [neighbor]
    else:
      current = neighbors_dict[identifier]
      neighbor = Neighbor(line[header.index('resname')], line[header.index('segid')], line[header.index('atomid')], line[header.index('vector')], line[header.index('distance')])
      neighbors_dict[identifier] = current + [neighbor]   

  in_file.close()

  return neighbors_dict

class Cysteine:
  def __init__ (self, model, dssp, pdbid, chainid, residue, segid, seq, identifier, neighbors):

    self.model = model
    self.dssp = dssp

    self.pdbid = pdbid
    self.chainid = chainid
    self.residue = residue
    self.segid = segid
    self.identifier = identifier

    self.sequence = seq
    self.target_atom = self.model[self.chainid][self.segid]['SG']
    self.vector = self.target_atom.get_vector()
    self.b_factor = self.target_atom.get_bfactor()
    self.disorder = self.residue.is_disordered()
    if dssp != '':
      self.secondary_structure = self.dssp[2]
      self.saas = self.dssp[3]
    else:
      self.secondary_structure = ''
      self.saas = ''
    self.vdwr = 1.8

    # proximity to other cysteine descriptors
    self.c_idxs = []

    # 2D
    self.up_2D_dist = 0
    self.down_2D_dist = 0
    self.up_2D_3D_dist = 0
    self.down_2D_3D_dist = 0

    # 3D
    self.sorted_3D_neighbors = []
    self.closest_3D_2D_dist = 0
    self.closest_3D_dist = 0

    self.get_c_idxs()
    if len(self.c_idxs) >= 2:
      self.get_c_2D_neighbors()
      self.get_c_3D_neighbors()

    # neighborhood descriptors
    self.neighborslist = neighbors
    self.neighborsresnames = []
    self.neighborcount = 0
    self.atomcount = 0
    self.netcharge = 0
    self.hydrophobic = 0
    self.hydrophobicity_kd = 0
    self.steric_list = []
    self.steric_count = 0

    self.aa_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    self.aa_count_dict = {'CYS': 0, 'ASP': 0, 'SER': 0, 'GLN': 0, 'LYS': 0,
     'ILE': 0, 'PRO': 0, 'THR': 0, 'PHE': 0, 'ASN': 0, 
     'GLY': 0, 'HIS': 0, 'LEU': 0, 'ARG': 0, 'TRP': 0, 
     'ALA': 0, 'VAL': 0, 'GLU': 0, 'TYR': 0, 'MET': 0}

    self.charge_dict = {'CYS': 0, 'ASP': -1, 'SER': 0, 'GLN': 0, 'LYS': 1,
     'ILE': 0, 'PRO': 0, 'THR': 0, 'PHE': 0, 'ASN': 0, 
     'GLY': 0, 'HIS': 1, 'LEU': 0, 'ARG': 1, 'TRP': 0, 
     'ALA': 0, 'VAL': 0, 'GLU': -1, 'TYR': 0, 'MET': 0}

    self.acidity_dict = {'CYS': 'P', 'ASP': 'A', 'SER': 'P', 'GLN': 'P', 'LYS': 'B',
     'ILE': 'ALI', 'PRO': 'P', 'THR': 'P', 'PHE': 'ARO', 'ASN': 'P', 
     'GLY': 'ALI', 'HIS': 'B', 'LEU': 'ALI', 'ARG': 'B', 'TRP': 'ARO', 
     'ALA': 'ALI', 'VAL':'ALI', 'GLU': 'A', 'TYR': 'ARO', 'MET': 'ALI'}

    self.hydrophobicity_dict = {'CYS': 2.5, 'ASP': -3.5, 'SER': -0.8, 'GLN': -3.5, 'LYS': -3.9,
     'ILE': 4.5, 'PRO': -1.6, 'THR': -0.7, 'PHE': 2.8, 'ASN': -3.5, 
     'GLY': -0.4, 'HIS': -3.2, 'LEU': 3.8, 'ARG': -4.5, 'TRP': -0.9, 
     'ALA': 1.8, 'VAL':4.2, 'GLU': -3.5, 'TYR': -1.3, 'MET': 1.9}

    self.acidity_count_dict = {'P': 0, 'ARO': 0, 'ALI': 0, 'A': 0, 'B':0}

    self.grasp_count_dict = {'HPB': 0, 'DON': 0, 'ACP': 0, 'ARM': 0, 'POS': 0, 'NEG': 0, 'SSB': 0}
    self.grasp_list_dict = {'HPB':{}, 'DON':{}, 'ACP':{}, 'ARM':{}, 'POS':{}, 'NEG':{}, 'SSB':{}}
    self.grasp_intxn_count_dict = {'ARO_STACK': 0, 'HB': 0, 'HYDROPHOBIC': 0, 'REPULSIVE': 0, 'SALT_BRIDGES': 0}
    self.grasp_intxn_min = {'ARO_STACK': 1.5, 'HB': 2.0, 'HYDROPHOBIC': 2.0, 'REPULSIVE': 2.0, 'SALT_BRIDGES': 2.0}
    self.grasp_intxn_max = {'ARO_STACK': 3.5, 'HB': 3.0, 'HYDROPHOBIC': 3.8, 'REPULSIVE': 6.0, 'SALT_BRIDGES': 6.0}

    self.update()

  def update(self):
    # iterate through list of neighbor atoms
    for i in range(len(self.neighborslist)):
      neighbor = self.neighborslist[i]
      identifier = neighbor.resname + '_' + neighbor.segid

      # calculate residue level counts and interactions
      # skip atoms in residues already accounted in calculations
      if identifier not in self.neighborsresnames:
        self.netcharge += neighbor.charge
        self.hydrophobic += neighbor.hydrophobicity
        self.aa_count_dict[neighbor.resname] +=1 
        self.acidity_count_dict[self.acidity_dict[neighbor.resname]] += 1
        self.neighborsresnames.append(identifier)

      # calucate atomic level sterics
      self.get_sterics(neighbor)
      self.update_neighbors(neighbor)

    if len(self.neighborslist) > 0:
      # average hydrophobicity of residues within cutofoff distance
      self.hydrophobicity_kd = self.hydrophobic / len(self.neighborsresnames)
      # Interactions between atoms based on pre-determined distances
      self.update_interactions()

    # update overall counts
    self.steric_count = len(self.steric_list)
    self.neighborcount = len(self.neighborsresnames)
    self.atomcount = len(self.neighborslist)

  def get_c_idxs(self):
    for i in range(len(self.sequence)):
      if self.sequence[i] == 'C':
        self.c_idxs.append(i+1)

  def get_c_2D_neighbors(self):
    try: 
      current = self.c_idxs.index(self.segid)

      if (current-1) != 0:
        upstream = self.c_idxs[current-1]
      else:
        upstream = self.c_idxs[0]

      if (current+1) != len(self.c_idxs):
        downstream = self.c_idxs[current+1]
      else:
        downstream = self.c_idxs[-1]

      self.up_2D_dist = self.segid - upstream
      self.down_2D_dist = downstream - self.segid
    except:
      self.up_2D_dist = ''
      self.down_2D_dist = ''

    try:
      self.upstream_target_atom = self.model[self.chainid][upstream]['SG']
      self.up_2D_3D_dist = numpy.linalg.norm(self.target_atom.get_vector() - self.upstream_target_atom.get_vector())
    except:
      self.upstream_target_atom = ''

    try:
      self.downstream_target_atom = self.model[self.chainid][downstream]['SG']
      self.down_2D_3D_dist = numpy.linalg.norm(self.target_atom.get_vector() - self.downstream_target_atom.get_vector())
    except:
      self.downstream_target_atom = ''

  def get_c_3D_neighbors(self):
    try:
      current = self.c_idxs.index(self.segid)
      c_3D_dict = {}
      for i in range(len(self.c_idxs)):
        if self.c_idxs[i] != self.segid:
          try:
            new_target_atom = self.model[self.chainid][self.c_idxs[i]]['SG']
          except:
            continue
          dist = numpy.linalg.norm(self.vector - new_target_atom.get_vector())
          c_3D_dict[self.c_idxs[i]] = dist

      if len(c_3D_dict.keys()) > 0:
        self.sorted_3D_neighbors = sorted(c_3D_dict.items(), key=lambda kv: kv[1])
        self.closest_3D_2D_dist = self.sorted_3D_neighbors[0][0]
        self.closest_3D_dist = self.sorted_3D_neighbors[0][1]

    except:
      self.sorted_3D_neighbors = []
      self.closest_3D_2D_dist = ''
      self.closest_3D_dist = ''

  # def update_neighbors(self, neighbor, central_sg):
  def update_neighbors(self, neighbor):
    details = neighbor.atomic_details
    if (details != '') & (details != '-'):
      if ',' in details:
        details = details.split(',')
        for i in range(len(details)):
          self.grasp_count_dict[details[i]] += 1
          current = self.grasp_list_dict[details[i]]
          current[neighbor.atomic_identifier] = neighbor
          self.grasp_list_dict[details[i]] = current

      elif len(details) == 3:
        self.grasp_count_dict[details] += 1
        current = self.grasp_list_dict[details]
        current[neighbor.atomic_identifier] = neighbor
        self.grasp_list_dict[details] = current

  def update_interactions(self):
    self.get_pairs(True, 'ARM', 'ARM', 'ARO_STACK')
    self.get_pairs(True, 'POS', 'POS', 'REPULSIVE')
    self.get_pairs(True, 'NEG', 'NEG', 'REPULSIVE')
    self.get_pairs(True, 'HPB', 'HPB', 'HYDROPHOBIC')
    self.get_pairs(False, 'DON', 'ACP', 'HB')
    self.get_pairs(False, 'POS', 'NEG', 'SALT_BRIDGES')
    
  def get_pairs(self, single, k1, k2, k3):

    if single == True:
      current = list(self.grasp_list_dict[k1].keys())
      
      # output is a list of tuples
      list_of_pairs = [(current[p1], current[p2]) for p1 in range(len(current)) for p2 in range(p1+1,len(current))]
    else:
      lst1 = self.grasp_list_dict[k1].keys()
      lst2 = self.grasp_list_dict[k2].keys()

      # output is a list of lists
      list_of_pairs = [[a, b] for a in lst1 for b in lst2 if a != b] 
    
    count = self.check_interactions(list_of_pairs, self.grasp_list_dict[k1], self.grasp_list_dict[k2], self.grasp_intxn_min[k3], self.grasp_intxn_max[k3])
    self.grasp_intxn_count_dict[k3] = self.grasp_intxn_count_dict[k3] + count

  def check_interactions(self, pairs, kv_dict1, kv_dict2, min_dist, max_dist):
    count = 0
    for i in range(len(pairs)):
      current = pairs[i]
      if kv_dict1[current[0]].segid != kv_dict2[current[1]].segid:
        dist = numpy.linalg.norm(kv_dict1[current[0]].vector - kv_dict2[current[1]].vector) 
        if (dist >= min_dist) & (dist <= max_dist):
          count += 1

    return count

  # find steric interactions based on van der waals radii
  # http://www.ambrsoft.com/TrigoCalc/Sphere/TwoSpheres/Intersection.htm
  def get_sterics(self, neighbor):
    if self.segid == neighbor.segid:
      return
    else:
      dist = numpy.linalg.norm(self.vector - neighbor.vector)
      if (dist < (self.vdwr + neighbor.vdwr)) & (neighbor.identifier not in self.steric_list):
        self.steric_list.append(neighbor.identifier)

  def output(self):
    return self.identifier + ',' + str(self.b_factor) + ',' + str(self.disorder) + ',' + self.secondary_structure + ',' + str(self.down_2D_dist) + ',' + str(self.down_2D_3D_dist) + ',' +  str(self.up_2D_dist) + ',' + str(self.up_2D_3D_dist) + ',' + str(self.closest_3D_2D_dist) + ',' + str(self.closest_3D_dist) + ',' + str(self.saas) + ',' + str(self.atomcount) + ',' + str(self.neighborcount) + ',' + str(self.netcharge) + ',' + str(self.hydrophobicity_kd) + ',' + str(self.steric_count) + ',' + self.list_to_string(self.steric_list) + ',' + self.dict_to_string(self.aa_count_dict) + ',' + self.dict_to_string(self.acidity_count_dict) + ',' + self.dict_to_string(self.grasp_count_dict) + ',' + self.dict_to_string(self.grasp_intxn_count_dict)

  def dict_to_string(self, dct):
    st = ''
    for k in dct:
      st += str(dct[k]) + ','
    return st[:-1]

  def list_to_string(self, lst):
    return (';'.join([str(elem) for elem in lst]))

# read pdb file
def get_pdb_details(pdb_files, descriptor_files, dist_cutoff, cd, rd):
  
  # master_list = []
  count = 0

  for i in range(len(pdb_files)):
    master_list = []

    # format 'reduced_pdbid.pdb'
    old_pdb = pdb_files[i]
    # new_pdb = old_pdb.replace('reduced_pdb', '').replace('.pdb', '').upper()
    new_pdb = old_pdb.replace('reduced_', '').replace('.pdb', '').upper()

    # skip completed files
    # if (old_pdb.replace('.pdb', '_descriptors.csv') in descriptor_files):
    if (old_pdb.replace('.pdb', '_' + str(dist_cutoff) + '_descriptors.csv') in descriptor_files):
      continue
    else:
      print(pdb_files[i])

    parser = Bio.PDB.PDBParser(QUIET=True) 

    try:
      structures = parser.get_structure(pdb_files[i][:-4], old_pdb)
    except UnicodeDecodeError:
      print('UnicodeDecodeError ' + old_pdb)
      continue
    except:
      print('Error for strucutre in ' + old_pdb)
      continue

    if len(structures) > 1:
      print("multiple structures")
      print(old_pdb, len(structures))

    try:
      model = structures[0] 
    except:
      print('KeyError for model in ' + old_pdb)
      continue
    chain_list = Selection.unfold_entities(model, "C")

    # access dssp for secondary structure and saas
    try:
      dssp = DSSP(model, old_pdb, dssp='/Users/lisamarieboatner/anaconda3/bin/mkdssp')

    except:
      print('Error for DSSP in ' + old_pdb)
      dssp = ''
    
    count = 0

    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
      'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
      'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
      'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    for chain in chain_list:
      chainid = chain.get_id()
      res_dict = {}
      cys_list = []
      residues = chain.get_residues()
      res_dict = {}
      chain_dict = {}

      for residue in residues:
        resname = residue.get_resname()
        segid = list(residue.id)[1]

        # skip amino acids with negative indices
        if (residue.get_resname() == "CYS") & (segid >= 0):
          cys_list.append(residue)

        if resname in d.keys():
          res_dict[segid] = d[resname]
        else:
          res_dict[segid] = '-'
          # print("not found", resname)

      chain_dict[chainid] = res_dict

      seq = get_sequence(res_dict)

      cysteine_list = get_cysteines(old_pdb, new_pdb, model, dssp, chainid, seq, cys_list, dist_cutoff, cd)
      master_list = master_list + cysteine_list

    dd = cd + '/descriptor_files_' + str(dist_cutoff)

    os.chdir(dd)

    if len(master_list) != 0:
      write_file(old_pdb.replace('.pdb', '_' + str(dist_cutoff) + '_descriptors.csv'), master_list, True)

    os.chdir(rd)
    
  return master_list

def get_sequence(res_dict):
  seq_len = max(res_dict.keys())
  seq = ['-'] * (seq_len + 1)

  for k in res_dict:
    # skip negative amino acids
    if k < 0:
      continue
    if k >= 0 | k <= seq_len:
      seq[k] = res_dict[k]

  str_seq = list_to_string(seq[1:], '')
  return str_seq

def get_cysteines(old_pdb, new_pdb, model, dssp, chainid, seq, cys_list, dist_cutoff, cd):
  cysteine_list = []
  for i in range(len(cys_list)):
    residue = cys_list[i]
    segid = list(residue.id)[1]
    identifier = new_pdb + '_' + chainid + '_C' + str(segid)

    # nd = cd + '/neighbor_files'
    nd = cd + '/neighbor_files_cys'

    os.chdir(nd)

    nf = old_pdb.replace('.pdb', '') + '_neighbors.csv'
    nf_path = Path(nd + '/' + nf)

    if nf_path.is_file():
      neighbors_dict = read_neighbors_file(old_pdb.replace('.pdb', '') + '_neighbors.csv', dist_cutoff)
    else:
      neighbors_dict = {}
      print('neighbors file not found')

    os.chdir(cd + '/reduced_files')

    if identifier in neighbors_dict.keys():
      current_neighbors = neighbors_dict[identifier]
    else:
      current_neighbors = []
      print('identifier not found in neighbors', identifier)

    try:
      target_atom = model[chainid][segid]['SG']
    except:
      print('identifier not found in pdb', identifier)
      continue

    try:
      dssp_info = dssp[chainid,(' ', segid, ' ')]
    except:
      dssp_info = ''
    cysteine = Cysteine(model, dssp_info, new_pdb, chainid, residue, segid, seq, identifier, current_neighbors)
    cysteine_list.append(cysteine)

  return cysteine_list

def list_to_string(lst, symbol):
  return (symbol.join([str(elem) for elem in lst]))

# def write_file(file, header, misc, special):
def write_file(file, misc, special):
  out_file = open(file, 'w', encoding='utf-8')

  print('writing output')

  aa_count_keys = 'CYS,ASP,SER,GLN,LYS,ILE,PRO,THR,PHE,ASN,GLY,HIS,LEU,ARG,TRP,ALA,VAL,GLU,TYR,MET'
  acidity_count_keys = 'P,ARO,ALI,A,B'
  grasp_count_keys = 'HPB,DON,ACP,ARM,POS,NEG,SSB'
  grasp_intxn_keys = 'ARO_STACK,HB,HYDROPHOBIC,REPULSIVE,SALT_BRIDGES'
  header = 'identifier,b_factor,disorder,secondary_structure,2D_down,2D_down_dist,2D_up,2D_up_dist,3D_near,3D_dist,saas,atom_count,neighbor_count,net_charge,hydrophobicity_kd,steric_count,steric_aas'
  header = header + ',' + aa_count_keys + ',' + acidity_count_keys + ',' + grasp_count_keys + ',' + grasp_intxn_keys
  
  out_file.write(header + '\n')

  if isinstance(misc, dict):
    write_dict_file(out_file, misc, special)

  elif isinstance(misc, list):
    write_list_file(out_file, misc, special)

  else:
    print("Cannot write file, unrecognized format. " + str(misc))

  out_file.close()

def write_dict_file(out_file, misc):
  for k in misc:
    out_file.write(k + ',' + str(misc[k]) + '\n')

def write_list_file(out_file, misc, special):
  for i in range(len(misc)):
    current = misc[i]
    if special == True:
      st = current.output()
      out_file.write(st + '\n')
    else:
      st = ''
      for j in range(len(current)):
        st +=  str(current[j]) + ','
      out_file.write(st[:-1] + '\n')

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

    reduced_pdb = 'reduced_pdb' + new_pdb.lower() + '.pdb'

    # update pdb dictionary
    if reduced_pdb not in pdb_list:
      pdb_list.append(reduced_pdb)

  in_file.close()
  return pdb_list

def main():
  cd = os.getcwd()
  rd = cd + '/reduced_files'

  parser = argparse.ArgumentParser()
  parser.add_argument('-d', '--distance', dest='d', nargs='?', default=10, type=str, help='default 10')

  args = parser.parse_args()

  dist_cutoff = int(args.d) 

  print('processing residues')
  os.chdir(rd)
  reduced_files = [f for f in os.listdir('.') if f.startswith("reduced_")]

  dd = cd + '/descriptor_files_' + str(args.d)

  try:
    # Create target Directory
    os.mkdir(dd)
    print("Directory " , dd ,  " Created ") 
  except FileExistsError:
    print("Directory " , dd ,  " already exists")
  
  os.chdir(dd)
  descriptor_files = [f for f in os.listdir('.') if f.endswith(".csv")]
  os.chdir(rd)

  master_list = get_pdb_details(reduced_files, descriptor_files, dist_cutoff, cd, rd)

if __name__ == "__main__":
  main()