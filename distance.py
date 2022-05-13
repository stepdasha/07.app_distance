import biotite.database.rcsb as rcsb
import datetime
import pymol
from pymol import *
import collections
from collections import defaultdict
from tqdm.auto import tqdm
import os
import matplotlib.pyplot as plt
import numpy as np
import csv

import streamlit as st

import pandas as pd

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

def distance(pdb_ids, resid_1, resName_1, resid_2, resName_2):

    # pdb_ids = ['6xm0', '6xr8']
    if not os.path.exists('PDB'):
        os.mkdir('PDB')
    dist318_292 = collections.defaultdict(list)  # empty dictionary for future rmsd
    error_pdbs = []

    for i in tqdm(pdb_ids):
        cmd.delete("*")
        cmd.fetch(i, path='./PDB/')
        for chain in ['A', 'B', 'C']:
            try:
                dist = cmd.get_distance(atom1=f'chain {chain} and i. {resid_1} and n. CA and r. {resName_1}',
                                        atom2=f'chain {chain} and i. {resid_2} and n. CA and r. {resName_2}')
                dist318_292[i].append(dist)
            except CmdException:
                error_pdbs.append(i)
                print('error with pdb numbering', i)

    #print(f'There are {len(error_pdbs)} incorrectly numbered PDBs: ', str(error_pdbs) )
    st.header('**Incorrectly numbered pdbs**')
    st.write(f'There are {len(error_pdbs)} chains with measurement error:', str(error_pdbs))
    f = open("IncorrectNumberingPDB.txt", "w")
    f.write(str(error_pdbs))
    f.close()
    return dist318_292


def analysis(distancesDict):
    """
    plot the histogram of distances
    """
    distances_only = list(distancesDict.values())

    st.write(f'Number of structures with at least one analyzed chain {len(distancesDict)}')

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
    df.to_csv("distancesDF.csv")
    st.write(df)
    plt.savefig('distance.png', bbox_inches='tight')
    #print(f'number of corrected numbered and analyzed structures {len(distancesDict)}')
    #return distances_only


#pdb_ids = get_spike_ids()
##dist = distance(pdb_ids, 318, 'PHE', 292,  'ALA')
#dist = distance(pdb_ids, 1000, 'ARG', 975,  'SER')
#analysis(dist)

