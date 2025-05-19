"""
Read experimental isoTOP-ABPP excel and check corresponding PDBs to make the final dataset.
"""

import re
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot
from prody import parsePDB, fetchPDB, confProDy, buildBiomolecules, writePDB
from MDAnalysis import Universe
from MDAnalysis.lib.util import convert_aa_code

ADDITIVES = ('SO4', 'GOL', 'EDO', 'SDS', 'TRS', 'LDA', 'DMS', 'PO4', 'ACT',
             'EPE', 'NH4', 'PEG', 'EOH', 'BME', 'FMT', 'PYR', 'MES', 'BCT',
             'MPD', 'TAR', 'NO3', 'IOD', 'FLC', 'LMR', 'SIN', 'BTB', 'CAC',
             'MLI', 'POG', 'EOH', 'IPA', 'DIO', 'IMD', 'HEZ', 'JEF', 'TBU')

def sort_isotop_excel_file(df):
  cleaned = []
  dirty = []

  protein_groups = df.groupby(['protein'])

  for protein_name, protein_group in protein_groups:
    peptide_groups = protein_group.groupby(['peptides'])

    for peptide_name, peptide_group in peptide_groups:        
      for index, peptide in peptide_group.iterrows():
        clean_flag = True

        if ';' in peptide['peptides']:
          sequences = peptide['peptides'].split(';')
        else:
          sequences = [peptide['peptides']]

        uniprot_resid = peptide['identifier'].split('_')[1]
        if uniprot_resid.startswith('C'):
          uniprot_resid = uniprot_resid[1:]

        cleaned_sequences = []

        for sequence in sequences:
          positions = np.array([m.start() for m in re.finditer('*', sequence)])
          positions -= np.arange(1, len(positions) + 1)
          sequence = sequence.replace('*', '')

          cleaned_sequences.append((sequence, positions))

          if positions.shape[0] > 1:
            clean_flag = False

        if clean_flag:
          peptide['uniprot_resid'] = uniprot_resid
          peptide['peptides'] = cleaned_sequences[0][0]
          peptide['position'] = cleaned_sequences[0][1][0]
          cleaned.append(peptide)
        else:
          dirty.append(peptide)

  df_cleaned = pd.DataFrame(cleaned).sort_values(by=['protein', 'identifier'])
  df_dirty = pd.DataFrame(dirty).sort_values(by=['protein', 'identifier'])

  columns = ['identifier', 'protein', 'uniprot_resid', 'peptides', 'position', 
              '2022_count', '2019_count', '2010_count', 'experiment_count',
              'isotop-3_median',
              'isotop-11_median',
              'isotop-6_median',
              'isotop-5_median',
              'isotop-12_median',
              'isotop-9_median',
              'isotop-7_median',
              'isotop-10_median',
              'isotop-2_median',
              'isotop-1_median',
              'isotop-13_median',
              'isotop-8_median',
              'isotop-4_median',
              'isotop-14_median',
              'isotop-15_median',
              'isotop-16_median',
              'isotop-17_median',
              'isotop-18_median',
              'isotop-19_median',
              'isotop-20_median',
              'isotop-21_median',
              'isotop-22_median',
              'isotop-23_median',
              'isotop-24_median',
              'mean', 'sd', 'pdb']

  df_cleaned = df_cleaned[columns]

  df_cleaned.to_excel('220505_isotop_pdb_cleaned.xlsx', engine='openpyxl', index=False)
  df_dirty.to_excel('220505_isotop_pdb_dirty.xlsx', engine='openpyxl', index=False)
  
  return df_cleaned, df_dirty


def get_peptides_pdbs(df):
  all_peptides = []
  all_pdbs = []

  error_count = 0

  for protein_name, protein_group in df_cleaned.groupby(by=['protein']):
      
    try:
        pdbs = sorted(set([p.strip() for p in ';'.join(protein_group['pdb'].values).split(';')]))
    except TypeError:
        error_count += 1
        continue
    
    # PDB
    data = [(protein_name, pdb) for pdb in pdbs]
    all_pdbs.extend(data)
    tmp = pd.DataFrame(data=data, columns=['protein', 'pdb'])
    
    ratio_count = protein_group['experiment_count']
    protein_group['ratio_count'] = ratio_count
    protein_group.rename(columns={'peptides': 'peptide', 'mean': 'ratio_mean', 'sd': 'ratio_sd'}, inplace=True)
    
    columns = ['protein', 'uniprot_resid', 'peptide', 'position', 'ratio_mean', 'ratio_sd', 'ratio_count']
    all_peptides.append(protein_group[columns].convert_dtypes())

  df_pdbs = pd.DataFrame(data=all_pdbs, columns=['protein', 'pdb'])

  # We don't need the same (protein, pdb) combination
  df_pdbs.drop_duplicates(inplace=True)
  df_pdbs.to_csv('list_all_pdbs.csv', index=False)

  df_peptides = pd.concat(all_peptides)
  df_peptides.to_csv('list_all_peptides.csv', index=False)


def get_peptides():
  df_pdbs = pd.read_csv('list_all_pdbs.csv')
  df_peptides = pd.read_csv('list_all_peptides.csv')

  data = []
  no_data = []

  for protein in df_peptides['protein'].unique():
      print("Protein           : %s" % protein)
      protein_pdbs = df_pdbs[df_pdbs['protein'] == protein]
      protein_peptides = df_peptides[df_peptides['protein'] == protein]
      
      for pdb in protein_pdbs['pdb']:
          print("PDB               : %s" % pdb)
          
          try:
            p = parsePDB('isotop_pdbs/pdb%s.ent' % pdb.lower())
          except OSError:
            no_data.append(pdb)
            pdb_file = fetchPDB(pdb.lower(), compressed=False)
            if pdb_file is not None:
              os.system('mv ' + '%s.pdb' % pdb.lower() + ' isotop_pdbs/pdb%s.ent' % pdb.lower())
              try:
                p = parsePDB(pdb_dir + '/pdb%s.ent' % pdb.lower())
              except:
                print("Could not download nor parse pdb.")
                continue
            else:
                continue
          
          i = 0
          j = 0

          for row in protein_peptides.itertuples():
          # We take only a fragment of the peptide i - 3:i:i + 3 because longer is 
          # the peptide, higher are our chance to get have a missing residue in the PDB
          start = row.position - 3 if row.position - 3 > 0 else 0
          end = row.position + 3
          selected_peptide = row.peptide[start:end + 1]
          position_in_selected_peptide = row.position - start
          
          sel = p.select('ca sequence "%s"' % selected_peptide)
          
          if sel is not None:
            # The sequence can match in multiple chainids
            current_chainids = list(dict.fromkeys(sel.getChids()))
            i += 1
            
            for chainid in current_chainids:
              tmp_sel = sel.select('chid %s' % chainid)
              
              if tmp_sel.getResnames()[position_in_selected_peptide] == 'CYS':
                  data.append((protein, row.uniprot_resid, pdb, row.peptide, row.position,
                               tmp_sel.getResnums()[position_in_selected_peptide], 
                               tmp_sel.getChids()[position_in_selected_peptide],
                               row.ratio_mean, row.ratio_sd, row.ratio_count))
      
                  j += 1
              else:
                  print('ERROR             : position is not Cys for %s' % pdb)
      
          print('Peptide identified: %d (in %d chain(s))' % (i, j))
    
    print("")

  df.sort_values(by=['protein', 'uniprot_resid', 'pdb', 'resid', 'chainid', 'ratio_mean'], inplace=True)
  df = df[~df[['protein', 'uniprot_resid', 'pdb', 'resid', 'chainid']].apply(frozenset, axis=1).duplicated()]
  df.to_csv('list_found_peptides.csv', index=False)
  
  nonparsed_pdbs_df = pd.DataFrame(nonparsed_pdbs, columns = ['pdb'])
  nonparsed_pdbs_df.to_csv('nonparsed_pdbs_df.csv', index=False)


def ss_dis():
  # Source: https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/primary-sequences-and-the-pdb-format

  data_ss_dis = []
  data = ""

  with open('ss_dis.txt') as f:
    lines = f.readlines()
  
    for line in lines:
      if line.startswith('>'):
        if data:
          data_ss_dis.append((pdb[1:], chainid, info_type, data))
          data = ""
          
        pdb, chainid, info_type = line.strip().split(':')
      else:
        data += line.strip()

  columns = ('pdb', 'chainid', 'info_type', 'data')
  df_ss_dis = pd.DataFrame(data=data_ss_dis, columns=columns)
  df_ss_dis.to_csv('ss_dis.csv', index=False)


def get_complete(current_dir, pdb_dir):
  df_ss_dis = pd.read_csv('ss_dis.csv')
  df_found_peptides = pd.read_csv('list_found_peptides.csv')

  data_pdb = []
  not_complete = []
  written_pdbs = []
  columns = ('protein', 'pdb', 'n_biomolecule', 'total_biomolecule', 'n_segid', 'total_segid', 
             'n_terminus_disorder', 'c_terminus_disorder', 'total_disorder',
             'experiment', 'resolution', 'mutation', 'engineered', 'modified', 'complete_seq')

  for protein in df_found_peptides['protein'].unique():
    print("Protein           : %s" % protein)
    
    protein_peptides = df_found_peptides[df_found_peptides['protein'] == protein]
    pdbs = protein_peptides['pdb'].unique()
    
    for pdb in pdbs:
      print("PDB               : %s" % pdb)
      
      # Get sequence info (missing residues, secondary structures, seqres)
      ss_dis = df_ss_dis[df_ss_dis['pdb'] == pdb]
      
      try:
        p, header = parsePDB('isotop_pdbs/pdb%s.ent' % pdb.lower(), header=True, biomol=True)
      except:
        # Prody fails on very big complex, so we don't construct the biological unit
        p, header = parsePDB('isotop_pdbs/pdb%s.ent' % pdb.lower(), header=True)

      polymers = header['polymers']
      chemicals = header['chemicals']
      
      if not isinstance(p, list):
        p = [p]
      
      modified = any([x.modified for x in polymers])
      mutation = any([x.mutation for x in polymers])
      engineered = any([x.engineered for x in polymers])
      experiment = header['experiment']
      try:
        resolution = header['resolution']
      except KeyError:
        # There is no resolution header for NMR structures
        resolution = 'nan'
      
      """
      print('MODIFIED  ', modified)
      print('MUTATION  ', mutation)
      print('CHID      ', polymers[0].chid)
      print('EC        ', polymers[0].ec)
      print('FRAGMENT  ', polymers[0].fragment)
      print('ENGINEERED', engineered)
      print('DBREFS    ', polymers[0].dbrefs)
      print('COMMENTS  ', polymers[0].comments)
      """
      
      for i, b in enumerate(p):
        n_terminus_disorder = 0
        c_terminus_disorder = 0
        total_disorder = 0
        complete = False
        chainids_corrected = None
        n_duplicate_chainids = 0
        ca = b.select('ca')
        
        chids = b.getChids()
        chindices = b.getChindices()
        
        try:
          # Get the chainids in the current biological unit
          current_chainids = list(dict.fromkeys(ca.getChids()))
        except AttributeError:
          # This happens for 3VYY, because the biological unit contains only the DNA
          # So we skip it
          continue
        
        # Reconstruct the full missing residues sequence (in the same chainid order)
        full_disorder = []
        for chainid in current_chainids:
          tmp = ss_dis[ss_dis['segid'] == chainid]
            
          try:
            full_disorder.append(tmp[tmp['info_type'] == 'disorder']['data'].values[0])
          except IndexError:
            # It means that the additionnal chainid is something very small like a peptide 
            # and was not considered as a protein by the PDB
            # So we ignore it, it is not part of the full sequence
            pass
        
        # Check for disordered residues in N and C-ter and the rest
        # If residues are only missing in the N and C-ter, we won't consider them as incomplete
        for d in full_disorder:
          total_disorder += d.count('X')
          
          # Look for X on the Nter
          for r in d:
            if r != 'X':
              break
            n_terminus_disorder += 1
        
          # Look for X on the Cter
          for r in d[::-1]:
            if r != 'X':
              break
            c_terminus_disorder += 1
        
        # Correct potential duplicated chainids
        # It means that there are duplicates chainids
        if np.unique(chids).size < np.unique(chindices).size:
          try:
            new_chainids = [string.ascii_uppercase[chindice] for chindice in chindices]
            b.setChids(new_chainids)
            chainids_corrected = True
          except IndexError:
            chainids_corrected = False
          
          n_duplicate_chainids = int(np.unique(chindices).size / np.unique(chids).size)
              
        # Remove segnames because MDAnalysis prefers to use segments than chains
        b.setSegnames('')
        
        if (n_terminus_disorder + c_terminus_disorder) == total_disorder:
          print('COMPLETE!     (%d, %d, %d)' % (n_terminus_disorder, c_terminus_disorder, total_disorder))
          complete = True
        else:
          print('NOT COMPLETE! (%d, %d, %d)' % (n_terminus_disorder, c_terminus_disorder, total_disorder))
        
        # Write PDB biological unit
            
        # Separate water and crystallographic additives from the rest
        additives_selection_str = ' or '.join(['resname %s' % a for a in ADDITIVES])
        not_water_additives = b.select('not (resname WAT or resname HOH or %s)' % additives_selection_str)
        
        Write PDBs without water and additives
        writePDB('%s_%d.pdb' % ( pdb.lower(), i + 1), not_water_additives)
        
        Write water molecules
        water = b.select('resname WAT or resname HOH')
        try:
          writePDB('%s_%d_water.pdb' % ( pdb.lower(), i + 1), water)
        
        except TypeError:
          pass
            
        Write additives
        additives = b.select(additives_selection_str)
        try:
          writePDB('%s_%d_additives.pdb' % ( pdb.lower(), i + 1), additives)
            
        except TypeError:
          pass
    
        data_pdb.append((protein, pdb,
                         i + 1, len(p),
                         len(current_chainids), len(polymers), n_duplicate_chainids, chainids_corrected,
                         n_terminus_disorder, c_terminus_disorder, total_disorder, complete,
                         experiment, resolution, mutation, engineered))
      
    print('')
    
  df_data_pdb = pd.DataFrame(data=data_pdb, columns=columns)
  df_data_pdb.to_csv('list_all_pdb.csv', index=False)


def get_environment(current_dir, pdb_dir):
  df_data_pdb = pd.read_csv('list_all_complete_pdb.csv')
  df_found_peptides = pd.read_csv('list_found_peptides.csv')

  df_complete_pdb = df_data_pdb[df_data_pdb['complete_seq'] == True]
  df_not_complete_pdb = df_data_pdb[df_data_pdb['complete_seq'] == False]

  df_complete_pdb.to_csv('list_all_complete_pdb.csv', index=False)
  df_not_complete_pdb.to_csv('list_all_not_complete_pdb.csv', index=False)

  data = []
  columns = ['protein', 'pdb', 'n_biomolecule', 'peptide', 'resid', 'chainid', 'reactive', 'ratio_mean', 'ratio_sd', 'ratio_count',
             'complete_res', 'non_std_aa', 'hetero', 'ion', 'zinc', 'altloc', 'nucleic']

  # Selection of the complete PDB files
  complete_pdbs = df_data_pdb[(df_data_pdb['mutation'] == False) & (df_data_pdb['chainids_corrected'] != False)]

  non_complete_pdbs = []

  AA = {'MET': 8, 'THR': 7, 'GLN': 9, 'ILE': 8, 'ASP': 8, 
        'LEU': 8, 'GLY': 4, 'PRO': 7, 'ASN': 8, 'TRP': 14, 
        'PHE': 11, 'ARG': 11, 'VAL': 7, 'LYS': 9, 'GLU': 9, 
        'SER': 6, 'ALA': 5, 'CYS': 6, 'HIS': 10, 'GLN': 9, 
        'SEC': 6, 'TYR': 12}


  for pdb in complete_pdbs.itertuples():
    print("Protein (PDB - Biomol) : %s (%s - %d)" % (pdb.protein, pdb.pdb, pdb.n_biomolecule))
    
    os.chdir(biomol_dir)
    
    pdb_filename = '%s_%d.pdb' % (pdb.pdb.lower(), pdb.n_biomolecule)
    
    try:
      p = parsePDB(pdb_filename, secondary=False)
    except:
      non_complete_pdbs.append(pdb.pdb + '_' + str(pdb.n_biomolecule))
      continue

    # We want only the PDBs that correspond to this protein
    # Sometimes a PDB can contain multiple protein (protein != chain/segment)
    df_reactive_cysteines = df_found_peptides[(df_found_peptides['protein'] == pdb.protein) & (df_found_peptides['pdb'] == pdb.pdb)]
    
    try:
      # Get all the cysteines that are in the same chain(s)
      chains_str = ' or '.join('chain %s' % c for c in df_reactive_cysteines['chainid'].unique())
      all_cysteines = p.select('name SG and resname CYS and (%s)' % chains_str)
    except:
      continue
    
    # This happens when there was a CYS in the original PDB, but it was composed of multiple biological units
    # in which one does not contains a CYS (example: DNA)
    if all_cysteines is not None:
      print('CYSTEINES (R / NR)     : %d / %d' % (len(df_reactive_cysteines['resid']), len(all_cysteines.getResnums())))
      
      for cys in all_cysteines:
        n_hetero = 0
        n_ion = 0
        n_zinc = 0
        n_altloc = 0
        peptide = ''
        close_to_nucleic = False
        uniprot_resid = ''
        ratio_mean = ''
        ratio_sd = ''
        ratio_count = ''
        reactive = False
        complete = False
        duplicate_chainid = False
        non_std_aa = False
        resid = cys.getResnum()
        chid = cys.getChid()
        
        df_selected_cysteines = df_reactive_cysteines[(df_reactive_cysteines['resid'] == resid) & (df_reactive_cysteines['chainid'] == chid)]

        if not df_selected_cysteines.empty:
          reactive = True
          uniprot_resid = df_reactive_cysteines['uniprot_resid']
          ratio_mean = df_selected_cysteines['ratio_mean'].values[0]
          ratio_sd = df_selected_cysteines['ratio_sd'].values[0]
          ratio_count = df_selected_cysteines['ratio_count'].values[0]
          peptide = df_selected_cysteines['peptide'].values[0]
        else:
          resids = range(resid - 4, resid + 4 + 1)
          selection_str = ' or '.join(['(resnum %s and chain %s and ca)' % (x, chid) for x in resids])
          
          try:
            tmp = p.select(selection_str)
            peptide = ''.join([convert_aa_code(x) for x in tmp.getResnames()])
            
            # It means that we have duplicate chainids
            if len(peptide) > 9:
                peptide = peptide[:int(len(peptide) / 2)]
          except:
            print(selection_str)
            print('ERROR: Something is wrong with the selection.')
            continue

        cysteine_selection_str = '(name SG and resname CYS and resnum %d and chain %s)' % (resid, chid)
        hetero = p.select('(not water and not ion) hetero within 6 of %s' % cysteine_selection_str)
        ion = p.select('ion within 6 of %s' % cysteine_selection_str)
        altloc = p.select('(protein or nucleic) and altloc A within 6 of %s' % cysteine_selection_str)
        nucleic = p.select('nucleic within 6 of %s' % cysteine_selection_str)
        residues = p.select('protein within 6 of %s' % cysteine_selection_str)

        if hetero is not None:
          n_hetero = len(hetero.getResnums())
        if ion is not None:
          n_ion = len(ion.getResnums())
          try:
              n_zinc = Counter(ion.getResnames())['ZN']
          except KeyError:
              pass
        if altloc is not None:
          n_altloc = len(set([(x.getResname(), x.getResnum(), x.getChid()) for x in altloc]))
        if nucleic is not None:
          close_to_nucleic = True
        if residues is not None:
          # We need to do a new selection in order to get all the atoms in each residue
          residues_list = set([(x.getResname(), x.getResnum(), x.getChid()) for x in residues])
          selection_str = ' or '.join(['(resnum %s and chain %s and not hydrogen)' % (x[1], x[2]) for x in residues_list])

          try:
            residues = p.select(selection_str)
          except:
            print(selection_str)
            print('ERROR: Something is wrong with the selection.')
            continue

          current_n_atom = len(residues.getResnums())
          expected_n_atom = 0
          
          for x in residues_list:
            try:
              expected_n_atom += AA[x[0]]
            except KeyError:
              print("ERROR: unknown amino acid %s. This PDB should have been tagged as modified!" % x[0])
              continue
          
          # We don't penalize if more atoms than expected
          if current_n_atom >= expected_n_atom:
            complete = True
          
          if [res.getResname() for res in residues if not res.getResname() in AA]:
            non_std_aa = True

        if complete:
          info_str = 'COMPLETE     - '
        else:
          info_str = 'NOT COMPLETE - '

        info_str += 'HET (%d) ION (%s / %s) ALT (%s) NUC (%s) NSTD (%s)' % (n_hetero, n_ion, n_zinc, n_altloc, close_to_nucleic, non_std_aa)
        if reactive:
          print('CYS %4d - %s  REACTIVE   : %s' % (resid, chid, info_str))
        else:
          print('CYS %4d - %s             : %s' % (resid, chid, info_str))

        data.append((pdb.protein, pdb.pdb, pdb.n_biomolecule, peptide, resid, chid, reactive, ratio_mean, ratio_sd, ratio_count,
                     complete, non_std_aa, n_hetero, n_ion, n_zinc, n_altloc, close_to_nucleic))

      print("")
    
  df_cys_env_pdb = pd.DataFrame(data=data, columns=columns)
  df_cys_env_pdb.to_csv('list_cysteine_env_pdb.csv', index=False)


def select():
  df_data_pdb = pd.read_csv('list_all_pdb.csv', na_filter=False)
  df_cys_env_pdb = pd.read_csv('list_cysteine_env_pdb.csv', na_filter=False)

  tmp = df_data_pdb.merge(df_cys_env_pdb, on=['protein', 'pdb', 'n_biomolecule'])
  tmp['index'] = tmp.index

  df_final = tmp[(tmp['mutation'] == False) & (tmp['chainids_corrected'] != False)]
  df_complete_res = tmp[(tmp['mutation'] == False) & (tmp['chainids_corrected'] != False) & (tmp['complete_res'] == True)]

  # Without Zinc
  df_selection_extended_no_zinc = tmp[(tmp['mutation'] == False) & (tmp['chainids_corrected'] != False)
                     & (tmp['hetero'] == 0) & (tmp['ion'] == 0) & (tmp['nucleic'] == False) & (tmp['mutation'] == False)]

  # df_selection_extended_no_zinc.to_csv('final_selection_extended_no_zinc.csv', index=False)
  df_selection_selected_no_zinc = tmp[(tmp['complete_seq'] == True) & (tmp['mutation'] == False) & (tmp['chainids_corrected'] != False)
                      & (tmp['hetero'] == 0) & (tmp['ion'] == 0) & (tmp['nucleic'] == False) 
                      & (tmp['experiment'] == 'X-RAY DIFFRACTION') 
                      & (tmp['n_terminus_disorder'] <= 10) & (tmp['c_terminus_disorder'] <= 10)]

  # df_selection_selected_no_zinc.to_csv('final_selection_selected_no_zinc.csv', index=False)
  df_selection_complete_no_zinc = tmp[(tmp['complete_seq'] == True) & (tmp['mutation'] == False) & (tmp['chainids_corrected'] != False)
                      & (tmp['hetero'] == 0) & (tmp['ion'] == 0) & (tmp['nucleic'] == False) 
                      & (tmp['experiment'] == 'X-RAY DIFFRACTION') 
                      & (tmp['n_terminus_disorder'] == 0) & (tmp['c_terminus_disorder'] == 0)]

  # df_selection_selected_no_zinc.to_csv('final_selection_complete_no_zinc.csv', index=False)
  df_selection_extended_zinc = tmp[(tmp['mutation'] == False) & (tmp['chainids_corrected'] != False)
                     & (tmp['hetero'] == 0) & (tmp['nucleic'] == False) 
                     & (tmp['zinc'] > 0)
                     & (tmp['experiment'] == 'X-RAY DIFFRACTION')]
  df_selection_extended_zinc = df_selection_extended_zinc[df_selection_extended_zinc['ion'] == df_selection_extended_zinc['zinc']]

  # With Zinc
  # df_selection_extended_zinc.to_csv('final_selection_extended_zinc.csv', index=True)
  df_selection_selected_zinc = tmp[(tmp['complete_seq'] == True) & (tmp['mutation'] == False) & (tmp['chainids_corrected'] != False)
                      & (tmp['hetero'] == 0) & (tmp['zinc'] > 0) & (tmp['nucleic'] == False) 
                      & (tmp['experiment'] == 'X-RAY DIFFRACTION') 
                      & (tmp['n_terminus_disorder'] <= 10) & (tmp['c_terminus_disorder'] <= 10)]
  df_selection_selected_zinc = df_selection_selected_zinc[df_selection_selected_zinc['ion'] == df_selection_selected_zinc['zinc']]

  # df_selection_selected_zinc.to_csv('final_selection_selected_zinc.csv', index=True)
  df_selection_complete_zinc = tmp[(tmp['complete_seq'] == True) & (tmp['mutation'] == False) & (tmp['chainids_corrected'] != False)
                      & (tmp['hetero'] == 0) & (tmp['zinc'] > 0) & (tmp['nucleic'] == False) 
                      & (tmp['experiment'] == 'X-RAY DIFFRACTION') 
                      & (tmp['n_terminus_disorder'] == 0) & (tmp['c_terminus_disorder'] == 0)]
  df_selection_complete_zinc = df_selection_complete_zinc[df_selection_complete_zinc['ion'] == df_selection_complete_zinc['zinc']]

  # df_selection_complete_zinc.to_csv('final_selection_complete_zinc.csv', index=True)

  # Merge
  df_selection_extended = pd.concat([df_selection_extended_no_zinc, df_selection_extended_zinc])
  df_selection_extended = df_selection_extended.drop_duplicates()
  df_selection_selected = pd.concat([df_selection_selected_no_zinc, df_selection_selected_zinc])
  df_selection_selected = df_selection_selected.drop_duplicates()
  df_selection_complete = pd.concat([df_selection_complete_no_zinc, df_selection_complete_zinc])
  df_selection_complete = df_selection_complete.drop_duplicates()
  
  extended_ids = list(df_selection_extended['index'].unique())
  selected_ids = list(df_selection_selected['index'].unique())
  complete_ids = list(df_selection_complete['index'].unique())
  complete_res_ids = list(df_complete_res['index'].unique())
  final_ids = list(df_final['index'].unique())
  
  tmp['extended'] = np.where(tmp['index'].isin(extended_ids), 1, 0)
  tmp['selected'] = np.where(tmp['index'].isin(selected_ids), 1, 0)
  tmp['complete'] = np.where(tmp['index'].isin(complete_ids), 1, 0)
  tmp['residue_selection'] = np.where(tmp['index'].isin(complete_res_ids), 1, 0)
  tmp['final_selection'] = np.where(tmp['index'].isin(final_ids), 1, 0)
  tmp = tmp.drop(columns = ['index'])
  tmp.to_csv('final_selection.csv', index = False)
  
  not_selected_permissive_tmp = tmp[tmp['extended'] == 0]
  not_selected_permissive_selected_tmp = not_selected_permissive_tmp[not_selected_permissive_tmp['final_selection'] == 1]
  not_selected_permissive_selected_tmp = not_selected_permissive_selected_tmp.drop_duplicates()
  not_selected_permissive_selected_tmp.to_csv('selection_unfinished.csv', index = False)


def main():
  # Read experimental isoTOP cysteine reactivity measurements
  df = pd.read_excel('isotop_pdb.xlsx', engine='openpyxl')

  # Extract data from excel file
  df_cleaned, df_dirty = sort_isotop_excel_file(df)

  # Get all peptides and PDBs
  get_peptides_pdbs(df_cleaned)

  # Search for peptides in PDB files
  get_peptides(current_dir, pdb_dir)

  # Prepare ss_dis file
  ss_dis()

  # Check if sequences are complete in PDB files
  get_complete(current_dir, pdb_dir)

  # Check cysteine environment
  get_environment(current_dir, pdb_dir)

  # Final selection
  select()

if __name__ == "__main__":
  main()
