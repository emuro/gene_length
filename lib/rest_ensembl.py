import requests, sys
import pandas as pd

SERVER ="https://rest.ensembl.org"

def get_REST_ensembl_species():
    '''Acces to the REST ensembl server to retrieve the species it annotates
    input: None
    output:
        version (int)
        number_of_species (int)
        pandas dataframe: the different species sorted by taxonId
            taxonId, species formal name, species common name
    '''

    ext = "/info/species?"

    # get the species from rest.ensembl.org
    #
    r = requests.get(SERVER+ext,headers={"Content-Type" : "application/json"})
    if not r.ok:
      r.raise_for_status()
      sys.exit()
    decoded = r.json() #print(repr(decoded))
    ##print(type(decoded))
    ##print(decoded.keys())
    #
    # get the ensembl release where the data is obtained from
    number_of_species = len(decoded["species"])
    ensembl_release=decoded["species"][0]["release"]
    ##print("The data corresponds to the release:", ensembl_release, "from ensembl")
    ##print("There are", number_of_species, "different species in this release")

    # get the species in a dataframe
    # (sorted by taxon_id)
    for i in range(len(decoded["species"])):
        ##print(decoded["species"][i], end="\n\n") # print the raw data (species)
        pass
    species_labels=["taxonId", "name", "commonName"]
    #print("columns of df_species:", species_labels)
    df_species = pd.DataFrame (columns=species_labels)
    for i in range(len(decoded["species"])):
        species_elem=[decoded["species"][i]["taxon_id"], decoded["species"][i]["name"],
              decoded["species"][i]["common_name"]]
        df_species_aux = pd.DataFrame([species_elem], columns=species_labels)
        df_species=pd.concat([df_species,df_species_aux], ignore_index=True)
    df_species.taxonId=df_species.taxonId.astype(int)
    df_species_sorted=df_species.sort_values(by="taxonId", ascending=True)

    ##print(df_species.describe())
    ##print(df_species_sorted)
    return ensembl_release, number_of_species, df_species_sorted



