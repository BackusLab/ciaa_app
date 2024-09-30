# Date Created: 230424
# Date Modified: 240325
# Author: Lisa Boatner
# Description: OHE descriptors for ML.

import glob
import re
import sys, os, argparse

import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
from statistics import mean
import matplotlib.pyplot as plt
from scipy import stats

from numpy import array
from numpy import argmax
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder

def encode_and_bind(original_dataframe, feature_to_encode):
    dummies = pd.get_dummies(original_dataframe[[feature_to_encode]])
    res = pd.concat([original_dataframe, dummies], axis=1)
    return(res)

def calculate_sterics(df, cutoff):
    vals = []
    
    for index, row in df.iterrows():
        if (row['steric_count_' + str(cutoff)] != 0) & (row['neighbor_count_' + str(cutoff)] != 0):
            try:
              vals.append((float(row['steric_count_' + str(cutoff)])/float(row['neighbor_count_' + str(cutoff)]))* 100)
            except:
              print("Encountered 0 by division error.")
              print(row['steric_count_' + str(cutoff)], row['neighbor_count_' + str(cutoff)])
              vals.append(0)
        else: 
            vals.append(0)
            
    df['steric_P_' + str(cutoff)] = vals

    return df

def get_ss_type(df, helix_codes, beta_codes, loop_codes):
    tmp = df.copy()
    tmp['ss_helix'] = np.where(tmp['secondary_structure'].isin(helix_codes), 1, 0)
    tmp['ss_beta'] = np.where(tmp['secondary_structure'].isin(beta_codes), 1, 0)
    tmp['ss_loop'] = np.where(tmp['secondary_structure'].isin(loop_codes), 1, 0)
    
    return tmp

def get_aa_percentage(df, aa, cutoff):
    tmp = df.copy()
    
    tmp["CYS_P" + "_" + str(cutoff)] = (tmp['CYS'].astype(int)/tmp['neighbor_count_' + str(cutoff)].astype(int))*100
    tmp["ASP_P" + "_" + str(cutoff)] = (tmp['ASP'].astype(int)/tmp['neighbor_count_' + str(cutoff)].astype(int))*100
    tmp["SER_P" + "_" + str(cutoff)] = (tmp['SER'].astype(int)/tmp['neighbor_count_' + str(cutoff)].astype(int))*100
    tmp["GLN_P" + "_" + str(cutoff)] = (tmp['GLN'].astype(int)/tmp['neighbor_count_' + str(cutoff)].astype(int))*100
    tmp["LYS_P" + "_" + str(cutoff)] = (tmp['LYS'].astype(int)/tmp['neighbor_count_' + str(cutoff)].astype(int))*100
    tmp["ILE_P" + "_" + str(cutoff)] = (tmp['ILE'].astype(int)/tmp['neighbor_count_' + str(cutoff)].astype(int))*100
    tmp["PRO_P" + "_" + str(cutoff)] = (tmp['PRO'].astype(int)/tmp['neighbor_count_' + str(cutoff)].astype(int))*100
    tmp["THR_P" + "_" + str(cutoff)] = (tmp['THR'].astype(int)/tmp['neighbor_count_' + str(cutoff)].astype(int))*100
    tmp["PHE_P" + "_" + str(cutoff)] = (tmp['PHE'].astype(int)/tmp['neighbor_count_' + str(cutoff)].astype(int))*100
    tmp["ASN_P" + "_" + str(cutoff)] = (tmp['ASN'].astype(int)/tmp['neighbor_count_' + str(cutoff)].astype(int))*100
    tmp["GLY_P" + "_" + str(cutoff)] = (tmp['GLY'].astype(int)/tmp['neighbor_count_' + str(cutoff)].astype(int))*100
    tmp["HIS_P" + "_" + str(cutoff)] = (tmp['HIS'].astype(int)/tmp['neighbor_count_' + str(cutoff)].astype(int))*100
    tmp["LEU_P" + "_" + str(cutoff)] = (tmp['LEU'].astype(int)/tmp['neighbor_count_' + str(cutoff)].astype(int))*100
    tmp["ARG_P" + "_" + str(cutoff)] = (tmp['ARG'].astype(int)/tmp['neighbor_count_' + str(cutoff)].astype(int))*100
    tmp["TRP_P" + "_" + str(cutoff)] = (tmp['TRP'].astype(int)/tmp['neighbor_count_' + str(cutoff)].astype(int))*100
    tmp["ALA_P" + "_" + str(cutoff)] = (tmp['ALA'].astype(int)/tmp['neighbor_count_' + str(cutoff)].astype(int))*100
    tmp["VAL_P" + "_" + str(cutoff)] = (tmp['VAL'].astype(int)/tmp['neighbor_count_' + str(cutoff)].astype(int))*100
    tmp["GLU_P" + "_" + str(cutoff)] = (tmp['GLU'].astype(int)/tmp['neighbor_count_' + str(cutoff)].astype(int))*100
    tmp["TYR_P" + "_" + str(cutoff)] = (tmp['TYR'].astype(int)/tmp['neighbor_count_' + str(cutoff)].astype(int))*100
    tmp["MET_P" + "_" + str(cutoff)] = (tmp['MET'].astype(int)/tmp['neighbor_count_' + str(cutoff)].astype(int))*100

    return tmp

def get_relative_percentages(df, cols, name, cutoff):
    tmp = df.copy()
    
    tmp[name + '_sum'] = tmp[cols].astype(int).sum(axis=1)
    
    for i in range(len(cols)):
        tmp[cols[i] + '_P_' + str(cutoff)] = (tmp[cols[i]].astype(int) / tmp[name + '_sum'])*100
        
    tmp = tmp.drop(columns = [name + '_sum'])
    return tmp

def get_environment_percentages(df, cols, sum_col_name, cutoff):
    tmp = df.copy()
        
    for i in range(len(cols)):
        tmp[cols[i] + '_P_' + str(cutoff)] = (tmp[cols[i]].astype(int) / tmp[sum_col_name].astype(int))*100
        
    return tmp

def get_1D_min_value(df):
    tmp = df.copy()
    
    min_dists = []
    min_values = []
    left_sides = []
    right_sides = []
    
    for index, row in tmp.iterrows():
        up_dist = row['2D_up_dist']
        down_dist = row['2D_down_dist']
        up = int(row['2D_up'])
        down = int(row['2D_down'])
        
        min_dist = 0
        min_value = 0
        
        # checked
        if ((up_dist == 0) & (down_dist == 0)):
            min_dist = 0
            min_value = 0
            left_sides.append(0)
            right_sides.append(0)
            
        elif ((up_dist == down_dist) & (up_dist != 0)):
            min_dist = up_dist
            min_value = abs(up)
            left_sides.append(1)
            right_sides.append(1)
            
        elif((up_dist == 0) & (down_dist != 0)):
            min_dist = down_dist
            min_value = abs(down)
            left_sides.append(0)
            right_sides.append(1)

        elif((up_dist != 0) & (down_dist == 0)):
            min_dist = up_dist
            min_value = abs(up)
            left_sides.append(1)
            right_sides.append(0)
            
        elif up_dist < down_dist:
            min_dist = up_dist
            min_value = abs(up)
            left_sides.append(1)
            right_sides.append(0)
            
        elif up_dist > down_dist:
            min_dist = down_dist
            min_value = abs(down)
            left_sides.append(0)
            right_sides.append(1) 
            
        min_dists.append(min_dist)
        min_values.append(min_value)
        
    return min_dists, min_values, left_sides, right_sides
            
def get_3D_min_value(df):
    tmp = df.copy()
    min_dists = []
    min_values = []
    left_sides = []
    right_sides = []
    
    for index, row in tmp.iterrows():
        if 'PDB' in row['identifier']:
            identifier = row['identifier'].replace('PDB', '').split('_')
        else:
            identifier = row['identifier'].split('_')

        cys = int(identifier[2].replace('C', ''))
        
        min_dist = 0
        min_value = 0

        # case: no cysteines nearby
        if row['3D_dist'] == 0:
            left_sides.append(0)
            right_sides.append(0)
        # case: cysteines nearby
        else:
            min_dist = row['3D_dist']
            
            difference = int(row['3D_near']) - cys
            min_value = abs(difference)
            
            if difference < 0:
                left_sides.append(1)
                right_sides.append(0)
            elif difference > 0:
                left_sides.append(0)
                right_sides.append(1)
        
        min_dists.append(min_dist)
        min_values.append(min_value)
        
    return min_dists, min_values, left_sides, right_sides

def get_new_columns(df, cutoff, rename_columns):
    
    old_names = df.columns.to_list()
    new_names = []
    
    for i in range(len(old_names)):
        if old_names[i] not in rename_columns:
            new_names.append(old_names[i])
        else:
            new_names.append(old_names[i] + '_' + str(cutoff))
            
    return new_names

def main():

  parser = argparse.ArgumentParser()
  parser.add_argument('-idir', '--input_directory', dest='idir', nargs='?', default="/Users/lisamarieboatner/dropbox/cysteine_reactivity_data'", type=str, help='default current_directory, options: /Users/lisamarieboatner/dropbox/cysteine_reactivity_data')
  parser.add_argument('-i', '--input', dest='i', nargs='?', default="input", type=str, help='default input')
  parser.add_argument('-o', '--output', dest='o', nargs='?', default="output", type=str, help='default output_modified_descriptors_10.csv')
  parser.add_argument('-c', '--cutoff', dest='c', nargs='?', default="10", type=str, help='default 10; options 5')
  parser.add_argument('-wo', '--write_outfile', dest='wo', nargs='?', default="True", type=str, help='default True')
  args = parser.parse_args()

  # os.chdir(args.idir)
  cutoff = int(args.c)

  # os.chdir('data')
  descr_10_df = pd.read_csv(args.i)
  descr_10_df = descr_10_df[descr_10_df['identifier'] != 'identifier']

  print("Read input.")

  descr_10_df['pdb'] = descr_10_df['identifier'].map(lambda x: str(x).split('_')[0].replace('PDB', ''))  
  descr_10_df['pdbid'] = descr_10_df['identifier'].map(lambda x: str(x).split('_')[0].replace('PDB', '') + '_' + str(x).split('_')[-1])
  descr_10_df['chainid'] = descr_10_df['identifier'].map(lambda x: str(x).split('_')[1])
  descr_10_df['structureid'] = descr_10_df['identifier'].map(lambda x: str(x).split('_')[0] + '_' + str(x).split('_')[1] + '_' + str(x).split('_')[-1].replace('C', ''))  
  print("Identifiers input.")

  descr_10_df = descr_10_df.rename(columns = {'steric_count': 'steric_count_' + str(cutoff), 'atom_count': 'atom_count_' + str(cutoff), 'neighbor_count': 'neighbor_count_' + str(cutoff)})
  descr_10_df['saas_P_' + str(cutoff)] = descr_10_df['saas'].map(lambda x: float(x) * 100)
  # descr_10_df['saas_P_' + str(cutoff)] = 0
  descr_10_df['secondary_structure'] = descr_10_df['secondary_structure'].fillna('-')
  descr_10_df = descr_10_df.drop(columns = ['steric_aas', 'saas'])
  print("SASA finished.")

  # descr_10_df['steric_P_' + str(cutoff)] = (descr_10_df['steric_count'].astype(float) / descr_10_df['neighbor_count'].astype(float)) * 100
  descr_10_df = calculate_sterics(descr_10_df, cutoff)
  descriptors_10_df = descr_10_df.copy()
  print("Sterics finished.")

  # convert ss to categories: helix, beta sheet, loop
  helix_codes = ['H', 'G', 'I']
  beta_codes = ['B', 'E', 'b']
  loop_codes = ['T', 'S', 'C']

  # ohe ss
  wc_df = encode_and_bind(descriptors_10_df, 'secondary_structure')
  wc_df = get_ss_type(wc_df, helix_codes, beta_codes, loop_codes)

  print("SS finished.")

  # convert raw counts to percentages
  aa = ['CYS', 'ASP', 'SER', 'GLN', 'LYS', 'ILE', 'PRO', 'THR', 'PHE', 'ASN', 
  'GLY', 'HIS', 'LEU', 'ARG', 'TRP', 'ALA', 'VAL', 'GLU', 'TYR', 'MET']
  wc_df = get_aa_percentage(wc_df, aa, cutoff)

  print('AA percentages finished.')

  # convert counts of interactions, acidity, and residues to percentages
  # atom descriptors
  atoms = [
    'HPB',
    'DON',
    'ACP',
    'ARM',
    'POS',
    'NEG',
    'SSB',
  ]

  # atom interactions
  interactions = [
    'ARO_STACK',
    'HB',
    'HYDROPHOBIC',
    'REPULSIVE',
    'SALT_BRIDGES'
  ]

  # residue descriptors
  residues = [
    'P',
    'ARO',
    'ALI',
    'A',
    'B'
  ]

  wc_df = get_environment_percentages(wc_df, atoms, 'atom_count_' + str(cutoff), cutoff)
  wc_df = get_environment_percentages(wc_df, interactions, 'atom_count_' + str(cutoff), cutoff)
  wc_df = get_environment_percentages(wc_df, residues, 'neighbor_count_' + str(cutoff), cutoff)
  wc_df = wc_df.fillna(0)

  print("Interactions finished.")

  # consolidate 2D_up, 2D_down, 2D_down_dist, 2D_up_dist
  min_dists, min_values, left_sides, right_sides = get_1D_min_value(wc_df)
  prox_df = wc_df.copy()
  prox_df['1D_prox_3D'] = min_dists
  prox_df['1D_prox_1D'] = min_values
  prox_df['1D_prox_left'] = left_sides
  prox_df['1D_prox_right'] = right_sides

  # convert 3D_near and 3D_dist descriptors
  min_3D_dist, min_3D_values, left_3D_sides, right_3D_sides = get_3D_min_value(wc_df)

  prox_3D_df = prox_df.copy()
  prox_3D_df['3D_prox_3D'] = min_3D_dist
  prox_3D_df['3D_prox_1D'] = min_3D_values
  prox_3D_df['3D_prox_left'] = left_3D_sides
  prox_3D_df['3D_prox_right'] = right_3D_sides

  print("Proximity finished.")

  desc_10_df = prox_3D_df.copy()

  dropped_columns_10 = [
    # 'saas',
    # 'saas_P',
    '2D_down',
    '2D_down_dist',
    '2D_up',
    '2D_up_dist',
    '3D_near',
    '3D_dist',
    'secondary_structure',
    'secondary_structure_-',
  ]

  desc_10_df = desc_10_df.drop(columns = dropped_columns_10)

  rename_columns = [
    'atom_count',
	# 'neighbor_count',
	'net_charge',
	'hydrophobicity_kd',
	# 'steric_count',
	'P',
	'ARO',
	'ALI',
	'A',
	'B',
	'HPB',
	'DON',
	'ACP',
	'ARM',
	'POS',
	'NEG',
	'SSB',
	'ARO_STACK',
	'HB',
	'HYDROPHOBIC',
	'REPULSIVE',
	'SALT_BRIDGES',
	'1D_prox_3D',
	'1D_prox_1D',
	'1D_prox_left',
	'1D_prox_right',
	'3D_prox_3D',
	'3D_prox_1D',
	'3D_prox_left',
	'3D_prox_right'
  ]

  new_columns_10 = get_new_columns(desc_10_df, args.c, rename_columns)
  d10_df = desc_10_df.copy()
  d10_df.columns = new_columns_10

  empty_df = d10_df[d10_df.isna().any(axis=1)]
  print(empty_df.shape)

  d10_df.to_csv(args.o + '_modified_descriptors_' + str(args.c) + '.csv', index = False)

if __name__ == "__main__":
  main()