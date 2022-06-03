

from distance import *

# Logo image
#image = Image.open('logo.png')

#st.image(image, use_column_width=True)

# Page title
st.markdown("""
# Distance measurement program for SARS-CoV-2 Spike PDB structures.

This program measures a distance between two residues CA atoms in all SARS-CoV-2 (Uniprot ID P0DTC2) Spike structures deposited to PDB.

---
""")

# Sidebar
with st.sidebar.header('Enter residues between which you measure distance.'):
    resid_1 = st.sidebar.text_input("Input residue 1 id")
    atom_1 = st.sidebar.text_input("Input residue 1 atom name")
    resid_2 = st.sidebar.text_input("Input residue 2 id")
    atom_2 = st.sidebar.text_input("Input residue 2 atom name")

    st.sidebar.markdown("""
""")

with st.sidebar.header(''):
    option = st.selectbox(
        'Are residues in the same chain?',
        ('same', 'different'))
    st.sidebar.markdown("""
""")

if st.sidebar.button('Measure'):
    pdb_ids = get_spike_ids()
    st.header('**Available pdb structures**')
    st.write(pdb_ids)

    with st.spinner("Loading and checking PDBs"):
        pdb_ids_updated = pdb_load_check(pdb_ids)

    with st.spinner("Measuring distance"):
        if option == 'different':
            dist, dist_reverse = distance_dif(pdb_ids_updated, resid_1, resid_2, atom_1, atom_2)
            st.header('**Histogram of distances**')
            analysis(dist, resid_1, atom_1, resid_2, atom_2, flag = 'diff1')
            analysis(dist_reverse, resid_1, atom_1, resid_2, atom_2, flag = 'diff2')

        elif option == 'same':
            dist = distance_same(pdb_ids_updated, resid_1, resid_2, atom_1, atom_2, flag = 'same')
            st.header('**Histogram of distances**')
            analysis(dist, resid_1, atom_1, resid_2, atom_2, flag= 'same')

    #uncomment if want to delete all PDB files
    #delete_PDB_folder()

    # Read in calculated descriptors and display the dataframe
    #st.header('**Histogram of distances**')
    #analysis(dist, resid_1, atom_1, resid_2, atom_2)
   # st.write(desc)
    #st.write(desc.shape)

else:
    st.info('Enter residues in the sidebar to start!')
