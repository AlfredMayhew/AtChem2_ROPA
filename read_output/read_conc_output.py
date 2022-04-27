#imports
import pandas as pd

def conc_output(out_path="/home/alfie/AtChem2/model/output/speciesConcentrations.output",
                species="ALL", **kwargs):
    """Function to convert AtChem2 concentration output files into dictionary 
    indexed by timestep"""
    
    #Create dataframe from file in out_path
    data=pd.read_csv(out_path, delim_whitespace=True, **kwargs)
    
    #select only the species specified in the "species" variable
    #
    if species=="ALL":
        pass
    
    elif type(species)== list:
        time_species=["t"]+species
        
        data= data[time_species]
    
    elif type(species)== str:
        time_species=["t"]+species.split()
        
        data= data[time_species]
   
    else:
        raise TypeError("""Species must be a list or string of species names to 
                        select from the model output file (or left blank to 
                        return all species present)""")
        
    #convert dataframe to dictionary indexed by time
    dict_data=data.set_index("t").to_dict()
    
    return dict_data