This porgram calculates distributions of distances of PDB structures of SARS-CoV-2 Spike available on PDB.
If you want to use it, pull it from github using the following commadn:

	git clone git@github.com:stepdasha/07.app_distance.git
Create the environment called 'structure' from the file environment.yml file with the following command:

	conda env create --file environment.yml
ACtivate the environment with the following command:

	conda activate structure
Once the environment is activated one can run this little programm with the follwong streamlit command:
	
	streamlit run app.py 

