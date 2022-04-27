def time_list(out_path):
    """Function to produce list of time steps from AtChem2 rate output file"""
    
    with open(out_path,"r") as file:
        lines=file.readlines()
        lines=lines[1:] #remove first header line
        
    #remove new line characters and split each line by whitespace
    for i,l in enumerate(lines):
        lines[i]=l.strip("\n").split()

    #initialise empty set to add times to
    timeset=set()
    for l in lines:
        timeset.add(float(l[0]))
    
    times=sorted(timeset) #turn unordered set into ordered list

    return times
    
def reaction_dict(out_path):
    """Function to produce a dictionary of reactions in an AtChem2 rate output
    file indexed by reactionNumber"""
    with open(out_path,"r") as file:
        lines=file.readlines()
        lines=lines[1:] #remove first header line
        
    #remove new line characters and split each line by whitespace
    for i,l in enumerate(lines):
        lines[i]=l.strip("\n").split()

    #make dictionary of reaction definitions for each reactionNumber
    reactions={}
    for l in lines:
        reactions[int(l[3])]=l[5]
    
    return reactions

def read_p_rates(out_path, species="ALL", drop_0=True, drop_net_0=True, 
                 drop_rev=False, error_for_non_species = True):
    """function to convert production rate output files from AtChem2
    into a nested dictionary containing lists of rates for each reaction at 
    each timestep"""
    #read lines of rate output file from atchem2
    with open(out_path,"r") as file:
        lines=file.readlines()
        lines=lines[1:] #remove first header line
    
    #remove new line characters and split each line by whitespace
    for i,l in enumerate(lines):
        lines[i]=l.strip("\n").split()
    
    #create dictionary of rates indexed by species and reactions number
    rates={}
    for l in lines: #add species index
        rates[l[2]]={}

    for l in lines: #add reactionNumber index and initialise with empty list
        rates[l[2]][int(l[3])]=[]

    for l in lines: #loop through list appending reaction rate at each timestep
        rates[l[2]][int(l[3])]=rates[l[2]][int(l[3])]+[float(l[4])]
    
  

    #create dictionary of reactions indexed by reaction number
    reactions={}
    for l in lines:
        reactions[int(l[3])]=l[5].split("=") #create list of two lists (reactants and products)
    
    for r in reactions: #split reactants and products by
        for i,n in enumerate(reactions[r]):
            reactions[r][i]=n.split("+")



    #filter to only contain species specified in "species" input including 
    #families of compounds e.g. NOx, where rates will be summed for each species
    #produced in the family      
    if species=="ALL":
        pass
    
    elif species=="NOx": 
        fam_spec=["NO","NO2","NO3","N2O5"] #list of species in family
        temp_rates={} #create temporary dict to store summed NOx rates
        temp_rates["NOx"]={}
        #create dict of reaction multipliers for each reaction
        r_mult={}
        for i in reactions:
            nox_count=0 #count number of nox species for each reactant
            for r in set(reactions[i][1]): #set because don't want rates to be 
                                           #multiplied by two e.g. if two NO2 
                                           #are present in products (this is already done by AtChem2)
                if r in fam_spec:
                    if r == "N2O5": #N2O5 contains 2 N atoms, so must be counted twice
                        nox_count=nox_count+2
                    else:
                        nox_count=nox_count+1
            r_mult[i]=nox_count
            
        #multiply rates of each reaction by the nox count for that reaction,
        #reactions corresponding to multiple species will be overwritten by the
        #subsequent loops, but shouldn't be a problem as the rates should be 
        #the same each time
        for s in rates:
            for r in rates[s]:
                temp_rates["NOx"][r]=[x*r_mult[r] for x in rates[s][r]]
        
        #set temp_rates to be the new rates dict
        rates=temp_rates
                    
    elif type(species)==list: #dictionary comprehension replaces rates with a 
                              #new dict containing only the required entries
        new_rates = {}
        for s in species:
            try:
                new_rates[s] = rates[s]
            except KeyError as error:
                if error_for_non_species:
                    raise error
                else: #don't raise the key error for a missing species if specified
                    pass
        rates = new_rates

    elif type(species)==str:
        rates={ s: rates[s] for s in species.split() }
        
    else:
        raise TypeError("""Species must be a list or string of species names to 
                        select from the model output file (or left blank to 
                        return all species present)""")



    #initialise list to store indexes to remove from dict due to drop_0 and 
    #drop_net_0
    remove=[] 



    #drop reactions where rate=0 throughout if specified in "drop_0" input   
    if drop_0==True:
        for s in rates:
            for r in rates[s]:
                if rates[s][r].count(0)==len(rates[s][r]): #if all values in list are 0
                    remove.append([s,r])                    
                else:
                    pass
            
    elif drop_0==False:
        pass
      
    else:
        raise TypeError("""drop_0 must be input as boolean. Defaults to False""")



    #drop any reactions where a species is both a product AND a reactant if 
    #specified in "drop_net_0" input
    if drop_net_0==True:
        #make dictionary of reaction defniitions for each reactionNumber
        for r in reactions:            
            for s in reactions[r][1]: #iterate over all species in reaction products
                #if number of that species in reactants= number in products
                if reactions[r][0].count(s)==reactions[r][1].count(s):
                    remove.append([s,r])
    
    elif drop_0==False:
        pass
    
    else:
        raise TypeError("""drop_net_0 must be input as boolean. Defaults to False""")



    #drop reactions where another reaction is present in the reactions list
    #that is the reverse, if specified in "drop_rev" input. e.g NO3+NO2=N2O5
    if drop_rev==True:
        for s in rates:
            for r in rates[s]:
                for i in reactions:
                    r_eq_p=set(reactions[r][0])==set(reactions[i][1])
                    p_eq_r=set(reactions[r][1])==set(reactions[i][0])
                    if r_eq_p and p_eq_r:
                            remove.append([s,r])
   
    elif drop_rev==False:
        pass
    
    else:
        raise TypeError("""drop_rev must be input as boolean. Defaults to False""")



    #remove indexes listed in "remove" from rates dict
    for s, r in remove: 
        try:
            del rates[s][r]
        
        except KeyError: #reactions may be added to "remove" list twice, and 
                         #raise a KeyError in second deletion
            pass

    
    #return dict of rates indexed by species and reaction number
    return rates
        
        
def read_l_rates(out_path, species="ALL", drop_0=True, drop_net_0=True, 
                 drop_rev=False, error_for_non_species = True):
    """function to convert loss rate output files from AtChem2
    into a nested dictionary containing lists of rates for each reaction at 
    each timestep"""
    #read lines of rate output file from atchem2
    with open(out_path,"r") as file:
        lines=file.readlines()
        lines=lines[1:] #remove first header line
    
    #remove new line characters and split each line by whitespace
    for i,l in enumerate(lines):
        lines[i]=l.strip("\n").split()
    
    #create dictionary of rates indexed by species and reactions number
    rates={}
    for l in lines: #add species index
        rates[l[2]]={}

    for l in lines: #add reactionNumber index and initialise with empty list
        rates[l[2]][int(l[3])]=[]

    for l in lines: #loop through list appending reaction rate at each timestep
        rates[l[2]][int(l[3])]=rates[l[2]][int(l[3])]+[float(l[4])]
    
  

    #create dictionary of reactions indexed by reaction number
    reactions={}
    for l in lines:
        reactions[int(l[3])]=l[5].split("=") #create list of two lists (reactants and products)
    
    for r in reactions: #split reactants and products by
        for i,n in enumerate(reactions[r]):
            reactions[r][i]=n.split("+")



    #filter to only contain species specified in "species" input including 
    #families of compounds e.g. NOx, where rates will be summed for each species
    #produced in the family      
    if species=="ALL":
        pass
    
    elif species=="NOx": 
        fam_spec=["NO","NO2","NO3","N2O5"] #list of species in family
        temp_rates={} #create temporary dict to store summed NOx rates
        temp_rates["NOx"]={}
        #create dict of reaction multipliers for each reaction
        r_mult={}
        for i in reactions:
            nox_count=0 #count number of nox species for each reactant
            for r in set(reactions[i][0]): #set because don't want rates to be 
                                           #multiplied by two e.g. if two NO2 
                                           #are present in products (this is already done by AtChem2)
                if r in fam_spec:
                    if r == "N2O5": #N2O5 contains 2 N atoms, so must be counted twice
                        nox_count=nox_count+2
                    else:
                        nox_count=nox_count+1
            r_mult[i]=nox_count
            
        #multiply rates of each reaction by the nox count for that reaction,
        #reactions corresponding to multiple species will be overwritten by the
        #subsequent loops, but shouldn't be a problem as the rates should be 
        #the same each time
        for s in rates:
            for r in rates[s]:
                temp_rates["NOx"][r]=[x*r_mult[r] for x in rates[s][r]]
        
        #set temp_rates to be the new rates dict
        rates=temp_rates
                    
    elif type(species)==list: #dictionary comprehension replaces rates with a 
                              #new dict containing only the required entries
        new_rates = {}
        for s in species:
            try:
                new_rates[s] = rates[s]
            except KeyError as error:
                if error_for_non_species:
                    raise error
                else: #don't raise the key error for a missing species if specified
                    pass
        rates = new_rates
        
    elif type(species)==str:
        rates={ s: rates[s] for s in species.split() }
        
    else:
        raise TypeError("""Species must be a list or string of species names to 
                        select from the model output file (or left blank to 
                        return all species present)""")



    #initialise list to store indexes to remove from dict due to drop_0 and 
    #drop_net_0
    remove=[] 



    #drop reactions where rate=0 throughout if specified in "drop_0" input   
    if drop_0==True:
        for s in rates:
            for r in rates[s]:
                if rates[s][r].count(0)==len(rates[s][r]): #if all values in list are 0
                    remove.append([s,r])                    
                else:
                    pass
            
    elif drop_0==False:
        pass
      
    else:
        raise TypeError("""drop_0 must be input as boolean. Defaults to True""")



    #drop any reactions where a species is both a product AND a reactant if 
    #specified in "drop_net_0" input
    if drop_net_0==True:
        #make dictionary of reaction defniitions for each reactionNumber
        for r in reactions:            
            for s in reactions[r][1]: #iterate over all species in reaction products
                #if number of that species in reactants= number in products
                if reactions[r][0].count(s)==reactions[r][1].count(s):
                    remove.append([s,r])
    
    elif drop_net_0==False:
        pass
    
    else:
        raise TypeError("""drop_net_0 must be input as boolean. Defaults to True""")



    #drop reactions where another reaction is present in the reactions list
    #that is the reverse, if specified in "drop_rev" input. e.g NO3+NO2=N2O5
    if drop_rev==True:
        for s in rates:
            for r in rates[s]:
                for i in reactions:
                    r_eq_p=set(reactions[r][0])==set(reactions[i][1])
                    p_eq_r=set(reactions[r][1])==set(reactions[i][0])
                    if r_eq_p and p_eq_r:
                            remove.append([s,r])
   
    elif drop_rev==False:
        pass
    
    else:
        raise TypeError("""drop_rev must be input as boolean. Defaults to False""")



    #remove indexes listed in "remove" from rates dict
    for s, r in remove: 
        try:
            del rates[s][r]
        
        except KeyError: #reactions may be added to "remove" list twice, and 
                         #raise a KeyError in second deletion
            pass
    
        
    #return dict of rates indexed by species and reaction number
    return rates       
        
        
        
        
        
        
        
        
        
    
    