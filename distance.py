import biotite.database.rcsb as rcsb
import datetime
#import pymol
from pymol import *
import collections
#from collections import defaultdict
from tqdm.auto import tqdm
import os
import matplotlib.pyplot as plt
import numpy as np
#import csv

import streamlit as st

import pandas as pd
import shutil

def get_spike_ids(uniprot_id="P0DTC2", min_weight=400, max_resolution=4.0):
    """
    get all pdbs with defined weight and resolution,
    input the uniprot_id (the default is spike), min_weight , and max_resolution
    """
    # uniprot_id = "P0DTC2" #spike in Sars-cov-2
    # max_resolution = 4.0
    # min_weight =400
    """
    in Da, structure min mass to get rid of rbd only structures,
    Spike mass is 429 Da.
    """
    query_by_uniprot_id = rcsb.FieldQuery(
        "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
        exact_match=uniprot_id,
    )
    today = datetime.datetime.now()

    query_by_resolution = rcsb.FieldQuery(
        "rcsb_entry_info.resolution_combined", less_or_equal=max_resolution
    )

    query_by_polymer_weight = rcsb.FieldQuery(
        "rcsb_entry_info.molecular_weight", greater=min_weight
    )

    query = rcsb.CompositeQuery(
        [
            query_by_uniprot_id,
            query_by_resolution,
            # query_by_polymer_count,
            query_by_polymer_weight,
        ],
        "and",
    )
    pdb_ids = rcsb.search(query)
    #print(f"Number of spike structures on  {today.year}-{today.month}-{today.day} with "
    #      f"resolution less than or equal to {max_resolution} with mass more than or equal to {min_weight}: {len(pdb_ids)}")
    #print("Selected PDB IDs:\n", *pdb_ids)
    st.write(f"Number of spike structures on  {today.year}-{today.month}-{today.day} with "
          f"resolution less than or equal to {max_resolution}A with mass more than or equal to {min_weight}: {len(pdb_ids)}")
    return (pdb_ids)

def pdb_load_check(pdb_ids):
    # pdb_ids = ['6xm0', '6xr8']
    if not os.path.exists('PDB'):
        os.mkdir('PDB')
    error_pdbs = []
    chains = ['A', 'B', 'C']

    for i in tqdm(pdb_ids):
        i = i.lower()
        cmd.delete("*")
        #cmd.fetch(i, path='./PDB/')
        if os.path.exists('./PDB/' + i + '.pdb'):
            cmd.load("./PDB/" + i + ".pdb")
        elif os.path.exists('./PDB/' + i + '.cif') :
            cmd.load("./PDB/" + i + ".cif")
        else:
            file = cmd.fetch(i, path='./PDB/', type='pdb')
            #st.write(file)
            if  file != i:
                cmd.fetch(i, path='./PDB/')
            #st.write("cif instead of pdb is fetched")

        #else:
        #    st.write(f'could not load {i}' )

        for item in chains:
            # sanity check that the numbering is correct, check if 1000 is ARG in chain A B C
            p1 = cmd.select("p1", f'chain {item} and i. 1000 and r. ARG and n. CA')
            if p1 != 1:
                error_pdbs.append(i.upper())
                break

    new_pdb_ids = [pdb for pdb in pdb_ids if pdb not in error_pdbs]

    st.header('**Incorrectly numbered pdbs**')
    st.write(f'There are {len(error_pdbs)} structures with 1000 ARG error in at least one chain:', str(error_pdbs))
    st.write(f'There are {len(new_pdb_ids)} structures without 1000 ARG error')
    f = open("IncorrectNumberingPDB.txt", "w")
    f.write(str(error_pdbs))
    f.close()

    return new_pdb_ids


def distance_dif(pdb_ids, resid_1,  resid_2, atom_1, atom_2):
        if not os.path.exists('error_residue'):
            os.mkdir('error_residue')
        #chains_ordered = ['A', 'B', 'C']

        # make a chain order

        dist_list = collections.defaultdict(list)  # empty dictionary for future rmsd
        dist_list_reverse = collections.defaultdict(list) # empty dictionary for future rmsd from reverse
        missing_residue = []
        missing_residue_reverse = []
        for i in tqdm(pdb_ids):
            cmd.delete("*")
            chains_ordered = ['A']

            i = i.lower()

            try:
                cmd.load("./PDB/" + i + ".pdb")
            except CmdException:
                cmd.load("./PDB/" + i + ".cif")
                #st.write('load cif from folder')

            dist1 = cmd.get_distance(atom1=f'chain A and i. 971 and n. CA',
                                     atom2=f'chain B and i. 752 and n. CA')

            dist2 = cmd.get_distance(atom1=f'chain A and i. 971 and n. CA',
                                     atom2=f'chain C and i. 752 and n. CA')

            if min(dist1, dist2) == dist1:
                chains_ordered.append('B')
                chains_ordered.append('C')
            elif min(dist1, dist2) == dist2:
                chains_ordered.append('C')
                chains_ordered.append('B')

            for j in range(0,3):
                try:
                    if j == 2 :
                        dist = cmd.get_distance(atom1=f'chain {chains_ordered[j]} and i. {resid_1} and n. {atom_1}',
                                                atom2=f'chain {chains_ordered[0]} and i. {resid_2} and n. {atom_2}')
                        dist_list[i].append(dist)
                    else:

                        dist = cmd.get_distance(atom1=f'chain {chains_ordered[j]} and i. {resid_1} and n. {atom_1}',
                                            atom2=f'chain {chains_ordered[j + 1]} and i. {resid_2} and n. {atom_2}')
                        dist_list[i].append(dist)
                except CmdException:
                    missing_residue.append(i.upper())
                    #break

            #test for the second plot in other direction
            for j in range(0, 3):
                try:
                    if j == 2:
                        dist_reverse = cmd.get_distance(atom1=f'chain {chains_ordered[0]} and i. {resid_1} and n. {atom_1}',
                                                atom2=f'chain {chains_ordered[j]} and i. {resid_2} and n. {atom_2}')
                        dist_list_reverse[i].append(dist_reverse)
                    else:

                        dist_reverse = cmd.get_distance(atom1=f'chain {chains_ordered[j + 1]} and i. {resid_1} and n. {atom_1}',
                                                atom2=f'chain {chains_ordered[j]} and i. {resid_2} and n. {atom_2}')
                        dist_list_reverse[i].append(dist_reverse)
                except CmdException:
                    missing_residue_reverse.append(i.upper())
                    # break

        st.header('**Incorrectly numbered pdbs**')
        st.write(f'There are {len(missing_residue)} chains with a problem in chosen residue:', str(missing_residue))
        error_file_name = './error_residue/errorsPDB_' + str(resid_1) + str(atom_1) + '_' + str(resid_2) + str(atom_2) + '.txt'
        f = open(error_file_name, "w")
        f.write(str(missing_residue))
        f.write(str(missing_residue_reverse))
        f.close()
        return dist_list, dist_list_reverse

def distance_same(pdb_ids, resid_1,  resid_2, atom_1, atom_2, flag):
    dist_list = collections.defaultdict(list)  # empty dictionary for future rmsd
    missing_residue = []
    if not os.path.exists('error_residue'):
        os.mkdir('error_residue')
    for i in tqdm(pdb_ids):
        cmd.delete("*")

        try:
            cmd.load("./PDB/" + i + ".pdb")
        except CmdException:
            cmd.load("./PDB/" + i + ".cif")
            #st.write('load cif from folder')

        for chain in ['A', 'B', 'C']:
            try:
                dist = cmd.get_distance(atom1=f'chain {chain} and i. {resid_1} and n. {atom_1}',
                                            atom2=f'chain {chain} and i. {resid_2} and n. {atom_2}')
                dist_list[i].append(dist)
            except CmdException:
                missing_residue.append(i)

    st.header('**Incorrectly numbered pdbs**')
    st.write(f'There are {len(missing_residue)} chains with a problem in chosen residue:', str(missing_residue))
    error_file_name = './error_residue/errorsPDB_' + str(resid_1) + str(atom_1) + '_' + str(resid_2) + str(atom_2) + '.txt'
    f = open(error_file_name, "w")
    f.write(str(missing_residue))
    f.close()
    return dist_list


def analysis(distancesDict, resid_1, atom_1, resid_2, atom_2, flag):
    """
    plot the histogram of distances
    """
    distances_only = list(distancesDict.values())

    st.write(f'Number of structures with at least one analyzed chain {len(distancesDict)}')

    if not os.path.exists('plots'):
        os.mkdir('plots')

    if not os.path.exists('distances'):
        os.mkdir('distances')

    fig = plt.figure(figsize=(15, 7.5))
    plt.hist(np.hstack(distances_only), bins=100, color="skyblue", edgecolor='white')
    plt.hist(np.hstack(distances_only), bins=100, color="skyblue", edgecolor='white')
    # sns.histplot(data=distances_only , binwidth=0.2)
    plt.xlabel('distance, A', fontsize=20)
    plt.ylabel("Count", fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    st.pyplot(fig)

    df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in distancesDict.items()])).transpose()
    name = str(resid_1) + str(atom_1) + '_' + str(resid_2) + str(atom_2)
    df_name = './distances/distance_' + name + '_' + flag +'.csv'
    df.to_csv(df_name)
    st.write(df)
    plt_name = './plots/distance_' + name + '_' + flag +'.png'
    plt.savefig(plt_name, bbox_inches='tight')
    #print(f'number of corrected numbered and analyzed structures {len(distancesDict)}')
    #return distances_only

def delete_PDB_folder():
        shutil.rmtree('./PDB')

#pdb_ids = get_spike_ids()
##dist = distance(pdb_ids, 318, 'PHE', 292,  'ALA')
#dist = distance(pdb_ids, 1000, 'ARG', 975,  'SER')
#analysis(dist)

