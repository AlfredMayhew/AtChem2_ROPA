"""Script to plot stacked plot of production and loss rates of species from
AtChem2 output files"""
#imports
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
from read_output.read_rate_output import time_list, reaction_dict, read_l_rates, read_p_rates
import pandas as pd
from datetime import date
import sys
import ast

#define functions
def string_to_bool(string):
    """function to evaluate strings of 'True' and 'False' into bools"""
    if string.casefold() == "True".casefold():
        return True
    elif string.casefold() == "False".casefold():
        return False
    else:
        raise Exception("Provided bool does not match True or False")
        
def get_keys(dictionary,val): 
    """function to return keys for any value"""
    keys=[]
    for key, value in dictionary.items(): 
         if val == value: 
             keys.append(key)
    if keys:
        return keys
    else:
        #if value isn't found the raise exception
        raise Exception(f"Specified value ({val}) doesn't exist in the dictionary")

#read in arguments from command line
args=sys.argv
kwarg_dict={"drop_rev":True,"title_page_text":"","remove_l_reactions":[],
            "remove_p_reactions":[],"exclusive_l_reactions":[],
            "exclusive_p_reactions":[],"lump_l_reactions":{},"lump_p_reactions":{}}

if len(args) < 6:
    raise Exception("""Must provide at least 5 arguments (in this order):
                     - Model Output Path
                     - Species of Interest (Comma Separated List e.g. NO2,O3,NO3)
                     - Number of Reactions to List (int, e.g. 10 for Top 10 Reactions, with the rest being lumped into "Other")
                     - Start Time (In Model Time or "START")
                     - End Time (In Model Time or "END")
                     
                     Additional key word arguments (e.g. drop_rev=True) are:
                     - drop_rev (bool, Ignore reversible reactions? Default: True)
                     - title_page_text (string, Text to add to the title page of the pdf output.)
                     - remove_l_reactions (Comma Separated List, loss reactions to exclude from the plot)
                     - remove_p_reactions (Comma Separated List, production reactions to exclude from the plot)
                     - exclusive_l_reactions (Comma Separated List, the only reactions to be included in the loss plots)
                     - exclusive_p_reactions (Comma Separated List, the only reactions to be included in the production plots)
                     - lump_l_reactions (Dictionary for each species (e.g. "{'NO2':{'NOx':['NO2+O3=NO3','HO2+NO2=OH+NO2']},'HO2':{'HOx':'HO2+O3=OH'}}"), loss reactions to add together into one catagory)
                     - lump_p_reactions (Dictionary for each species (e.g. "{'NO2':{'NOx':['NO2+O3=NO3','HO2+NO2=OH+NO2']},'HO2':{'HOx':'HO2+O3=OH'}}"), production reactions to add together into one catagory)
                     """)
else: #if 5 or more arguments passed then get the initial 5
    print(args)    

    out_path=args[1]
    
    species=args[2].strip("[]").split(",") #strip [] in case it was entered as python syntax list
    #turn "species" into a string (as opposed to a list) if only one specified
    if len(species) == 1:
        species = species[0]
    
    top_n=int(args[3])
    
    if args[4].casefold() == "start".casefold():
        start = "START"
    else:
        try:
            start = float(args[4])
        except ValueError:
            raise Exception("Start time must be numeric or 'START'")
    
    if args[5].casefold() == "end".casefold():
        end = "END"
    else:
        try:
            end = float(args[5])
        except ValueError:
            raise Exception("Start time must be numeric or 'START'")    

if len(args) > 6: #if there are more than 5 arguments then process the kwargs
    for kwarg in args[6:]:
        kw=kwarg.split("=")[0]
        arg=kwarg.split("=",1)[1] #only split at the first =, may be subsequent = from reaction definitions
        
        if (kw == "drop_rev"): #bools
            kwarg_dict[kw] = string_to_bool(arg)
        elif kw == "title_page_text": #strings
            kwarg_dict[kw] = arg
        elif (kw == "remove_l_reactions") or (kw == "remove_p_reactions") or (kw == "exclusive_l_reactions") or (kw == "exclusive_p_reactions"):  #lists
            kwarg_dict[kw] = arg.strip("[]").split(",") #strip [] in case it was entered as python syntax list
        else: #dicts
            kwarg_dict[kw] = ast.literal_eval(arg)

#read rates output files and record production and loss rates
l_rates=read_l_rates(out_path+"/lossRates.output",species=species,
                     drop_rev=kwarg_dict["drop_rev"], error_for_non_species=False)
p_rates=read_p_rates(out_path+"/productionRates.output",species=species,
                     drop_rev=kwarg_dict["drop_rev"], error_for_non_species=False)


#get list of times and definitions of reactionNumbers
times=time_list(out_path+"/productionRates.output")
l_reactions=reaction_dict(out_path+"/lossRates.output")
p_reactions=reaction_dict(out_path+"/productionRates.output")


#filter rates to only the specified times
if start != "START":
    absolute_difference_function = lambda list_value : abs(list_value - start)
    closest_value = min(times, key=absolute_difference_function)
    
    startindex= times.index(closest_value)
else:
    startindex = 0

if end != "END":
    absolute_difference_function = lambda list_value : abs(list_value - end)
    closest_value = min(times, key=absolute_difference_function)
    
    endindex = times.index(closest_value)
else:
    endindex = -1

times = times[startindex:endindex]
for k1 in p_rates:
    for k2 in p_rates[k1]:
        p_rates[k1][k2] = p_rates[k1][k2][startindex:endindex]
for k1 in l_rates:
    for k2 in l_rates[k1]:
        l_rates[k1][k2] = l_rates[k1][k2][startindex:endindex]


#if individual reaction(s) specified in "exclusive_reactions" then select only 
#that reactions
if kwarg_dict["exclusive_l_reactions"] != []:
    exclloss=[]
    for r in kwarg_dict["exclusive_l_reactions"]:
        exclloss += [x for x in get_keys(l_reactions,r)]
    for s in l_rates:
        l_rates[s]={k:l_rates[s][k] for k in l_rates[s].keys() if k in exclloss}

if kwarg_dict["exclusive_p_reactions"] != []:
    exclprod=[]
    for r in kwarg_dict["exclusive_p_reactions"]:
        exclprod += [x for x in get_keys(p_reactions,r)]
    for s in p_rates:
        p_rates[s]={k:p_rates[s][k] for k in p_rates[s].keys() if k in exclprod}
            
    
#remove reactions that don't want plotted (specified in remove_x_reactions variables)
remloss=[]
for r in kwarg_dict["remove_l_reactions"]:
    remloss += [x for x in get_keys(l_reactions,r)]
for s in l_rates:
    try:
        l_rates[s]={k:l_rates[s][k] for k in l_rates[s].keys() if k not in remloss} #replace l_rates with dictionary containing only desired reactions
    except KeyError: #KeyError raised if one of the exclusive reactions doesn't include one of the species
        pass

remprod=[]
for r in kwarg_dict["remove_p_reactions"]:
    remprod += [x for x in get_keys(p_reactions,r)]
for s in p_rates:
    try:
        p_rates[s]={k:p_rates[s][k] for k in p_rates[s].keys() if k not in remprod} #replace l_rates with dictionary containing only desired reactions
    except KeyError: #KeyError raised if one of the exclusive reactions doesn't include one of the species
        pass
         
#
#use dataframes to sort dictionary and select top n
pdrop=[] #empty dictionaries raise errors in plotting, so remove them (i.e. if there are no production and/or loss reactions of a species)
ldrop=[]
for s in l_rates:
    try:
        l_rates_df=pd.DataFrame.from_dict(l_rates[s]).reindex(pd.DataFrame.from_dict(l_rates[s]).median().sort_values(ascending=False).index, axis=1)
        l_rates[s]=l_rates_df.iloc[:,0:top_n].to_dict("list")
        if l_rates[s]:
            l_rates[s]["Other"] = l_rates_df.iloc[:,top_n:].sum(axis=1).to_list()
    except KeyError as e:
        if e.args[0] == 0: #if the KeyError is due to an empty dict, append species to list of keys to be removed
            ldrop.append(s)
l_rates={k:l_rates[k] for k in l_rates.keys() if (l_rates[k]) and (k not in ldrop)}

for s in p_rates:
    try:
        p_rates_df=pd.DataFrame.from_dict(p_rates[s]).reindex(pd.DataFrame.from_dict(p_rates[s]).median().sort_values(ascending=False).index, axis=1)
        p_rates[s]=p_rates_df.iloc[:,0:top_n].to_dict("list")
        if p_rates[s]:
            p_rates[s]["Other"] = p_rates_df.iloc[:,top_n:].sum(axis=1).to_list()
    except KeyError as e:
        if e.args[0] == 0: #if the KeyError is due to an empty dict, append species to list of keys to be removed
            pdrop.append(s)
p_rates={k:p_rates[k] for k in p_rates.keys() if (p_rates[k]) and (k not in pdrop)}


#
#if lumping reactions by reaction type then sum relavent reactions
if kwarg_dict["lump_l_reactions"] != {}:
    lump_l_rxns = kwarg_dict["lump_l_reactions"]
    new_l_rates={} #initialise dictionary with same structure as l_rates dict
    for s in lump_l_rxns.keys():
        new_l_rates[s]={}
        for cat in lump_l_rxns[s]:
            reaction_keys = []
            for reaction in lump_l_rxns[s][cat]:
                reaction_keys += get_keys(l_reactions,reaction)
                
            new_l_rates[s][cat]={k:[] for k in reaction_keys} 
            new_l_rates[s][cat]=[sum(x) for x in zip(*[l_rates[s][k] for k in reaction_keys])] #for each species and lump of reactions, add all required reaction rates
    
    old_l_rates=l_rates
    l_rates=new_l_rates #overwrite l_rates dict
    
if kwarg_dict["lump_p_reactions"] != {}:
    lump_p_rxns = kwarg_dict["lump_p_reactions"]
    new_p_rates={} #initialise dictionary with same structure as l_rates dict
    for s in lump_p_rxns.keys():
        new_p_rates[s]={}
        for cat in lump_p_rxns[s]:
            reaction_keys = []
            for reaction in lump_p_rxns[s][cat]:
                reaction_keys += get_keys(p_reactions,reaction)
                
            new_p_rates[s][cat]={k:[] for k in reaction_keys} 
            new_p_rates[s][cat]=[sum(x) for x in zip(*[p_rates[s][k] for k in reaction_keys])] #for each species and lump of reactions, add all required reaction rates
    
    old_p_rates=p_rates
    p_rates=new_p_rates #overwrite p_rates dict



#calculate procution/loss rates as a % of loss/producion from the topn reactions
pcl_rates={}
tot_l_rate={k:[] for k in l_rates.keys()}
for s in tot_l_rate.keys():
    tot_l_rate[s]=[sum(x) for x in zip(*[l_rates[s][k] for k in l_rates[s].keys()])]
for s in l_rates.keys():
    pcl_rates[s]={}
    for r in l_rates[s].keys():
        pcl_rates[s][r]=[] #list of percentage loss rates
        for rate,tot in zip(l_rates[s][r],tot_l_rate[s]):
            if tot != 0: #avoid divide by zero error
                pcl_rates[s][r].append((rate/tot)*100)
            else:
                pcl_rates[s][r].append(0)
pcp_rates={}
tot_p_rate={k:[] for k in p_rates.keys()}
for s in tot_p_rate.keys():
    tot_p_rate[s]=[sum(x) for x in zip(*[p_rates[s][k] for k in p_rates[s].keys()])]
for s in p_rates.keys():
    pcp_rates[s]={}
    for r in p_rates[s].keys():
        pcp_rates[s][r]=[] #list of percentage loss rates
        for rate,tot in zip(p_rates[s][r],tot_p_rate[s]):
            if tot != 0: #avoid divide by zero error
                pcp_rates[s][r].append((rate/tot)*100)
            else:
                pcp_rates[s][r].append(0)


#
#Plotting
pp=PdfPages("temp_rates_plot.pdf")
cols=cm.get_cmap("tab20")
#count number of species in rates (for plotting)
n_species=len(l_rates)

#labels for stacked plots
l_labels={}
for s in l_rates:
    l_labels[s]=[]
    for k in l_rates[s].keys():
        if not kwarg_dict["lump_l_reactions"]: #if not using lumped reactions then label is reaction
            if k != "Other":
                l_labels[s].append(l_reactions[k])
            else:
                l_labels[s].append(k)
        else: #if using lumped reactions, then label is the name of the "lump" i.e. what type of reaction
            l_labels[s].append(k)

p_labels={}
for s in p_rates:
    p_labels[s]=[]
    for k in p_rates[s].keys():
        if not kwarg_dict["lump_p_reactions"]: #if not using lumped reactions then label is reaction
            if k != "Other":
                p_labels[s].append(p_reactions[k])
            else:
                p_labels[s].append(k)
        else: #if using lumped reactions, then label is the name of the "lump" i.e. what type of reaction
            p_labels[s].append(k)
            
#create title page for pdf output
confirstPage = plt.Figure()
plottxt=f"""Production and loss rates of: {l_rates.keys()}

Reversible Reactions Removed : {kwarg_dict['drop_rev']}

Exclusive Reactions : {kwarg_dict['exclusive_p_reactions']+kwarg_dict['exclusive_l_reactions']}
Removed Reactions : {kwarg_dict['remove_p_reactions']+kwarg_dict['remove_l_reactions']}

Lumped Prodution reactions : {kwarg_dict['lump_p_reactions']}
Lumped Loss reactions : {kwarg_dict['lump_l_reactions']}

{kwarg_dict['title_page_text']}

Plotted on {date.today()}"""
confirstPage.text(0.5,0.2,plottxt, size=10, ha="center", wrap=True)
pp.savefig(confirstPage)
        
#Plot loss rates as function of time
lfig=plt.figure(figsize=(10,5*n_species))

for i,s in enumerate(l_rates):
    #plot temporary plot without "other" to get axis limits
    temp_fig=plt.figure()
    temp_ax=temp_fig.add_subplot(111)
    temp_ax.stackplot(times,[v for v in l_rates[s].values()][:-1])
    ylims=temp_ax.get_ylim()
        
    #plot actual plot
    ax=lfig.add_subplot(n_species,1,i+1)
    ax.stackplot(times,[v for v in l_rates[s].values()],labels=l_labels[s],
                 colors=cols([i for i,x in enumerate(l_labels[s])]))
    
    ax.set_ylim(ylims)
    ax.set_ylabel(f"k[{s}] / molecule cm-3 s-1")
    ax.set_xlabel("time / s")
    ax.set_title(f"{s} Loss")
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

lfig.suptitle("Loss Reactions")
pp.savefig(lfig,bbox_inches="tight")


#Plot production rates as function of time
pfig=plt.figure(figsize=(10,5*n_species))

for i,s in enumerate(p_rates):
    #plot temporary plot without "other" to get axis limits
    temp_fig=plt.figure()
    temp_ax=temp_fig.add_subplot(111)
    temp_ax.stackplot(times,[v for v in p_rates[s].values()][:-1])
    ylims=temp_ax.get_ylim()

    #plot actual plot
    ax=pfig.add_subplot(n_species,1,i+1)
    ax.stackplot(times,[v for v in p_rates[s].values()],labels=p_labels[s],
                 colors=cols([i for i,x in enumerate(p_labels[s])]))
    
    ax.set_ylim(ylims)
    ax.set_ylabel(f"k[{s}] / molecule cm-3 s-1")
    ax.set_xlabel("time / s")
    ax.set_title(f"{s} Production")
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

pfig.suptitle("Production Reactions")
pp.savefig(pfig,bbox_inches="tight")
    
#plot %loss rates vs time
pclfig=plt.figure(figsize=(10,5*n_species))
for i,s in enumerate(pcl_rates):
    ax=pclfig.add_subplot(n_species,1,i+1)
    ax.stackplot(times,[v for v in pcl_rates[s].values()],labels=l_labels[s],
                 colors=cols([i for i,x in enumerate(l_labels[s])]))
    ax.set_ylim(0,100)
    
    ax.set_ylabel("Proportion of loss / %")
    ax.set_xlabel("time / s")
    ax.set_title(f"{s} Loss")
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

pfig.suptitle("% Loss Reactions")
pp.savefig(pclfig,bbox_inches="tight")


#plot %production rates vs time
pcpfig=plt.figure(figsize=(10,5*n_species))
for i,s in enumerate(pcp_rates):
    ax=pcpfig.add_subplot(n_species,1,i+1)
    ax.stackplot(times,[v for v in pcp_rates[s].values()],labels=p_labels[s],
                 colors=cols([i for i,x in enumerate(p_labels[s])]))
    ax.set_ylim(0,100)
    
    ax.set_ylabel("Proportion of production / %")
    ax.set_xlabel("time / s")
    ax.set_title(f"{s} production")
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

pfig.suptitle("% production Reactions")
pp.savefig(pcpfig,bbox_inches="tight")


pp.close()

