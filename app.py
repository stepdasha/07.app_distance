import streamlit as st
import pandas as pd
from PIL import Image
import subprocess
import os
import base64
import pickle

import distance
from distance import *

# Logo image
#image = Image.open('logo.png')

#st.image(image, use_column_width=True)

# Page title
st.markdown("""
# Distance measurement app for SARS-CoV-2 Spike PDB structures.

This app measures a distance between two residues CA in all SARS-CoV-2 Spike structures deposited to PDB.

---
""")

# Sidebar
with st.sidebar.header('Enter residues between which you measure distance. Note: residues should be in the same chain'):
    chain_1 = st.sidebar.text_input("Input chain of residue 1 (A, B or C)")
    resid_1 = st.sidebar.text_input("Input residue 1 id")
    #resName_1 = st.sidebar.text_input("Input residue 1 name")
    chain_2 = st.sidebar.text_input("Input chain of residue 2 (A, B or C)")
    resid_2 = st.sidebar.text_input("Input residue 2 id")
    #resName_2 = st.sidebar.text_input("Input residue 2 name")
    st.sidebar.markdown("""
""")

if st.sidebar.button('Measure'):
    load_data = get_spike_ids()
    st.header('**Available pdb structures**')
    st.write(load_data)

    with st.spinner("Measuring distance"):
        dist = distance(load_data, chain_1, resid_1, chain_2, resid_2)


    # Read in calculated descriptors and display the dataframe
    st.header('**Histogram of distances**')
    analysis(dist)
   # st.write(desc)
    #st.write(desc.shape)

else:
    st.info('Enter residues in the sidebar to start!')
