from read_output.read_rate_output import reaction_dict, read_p_rates, read_l_rates, time_list
from statistics import mean
from itertools import islice
import sys

#function to select first n objects in iter (dict here)
def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))

#put command-line args into list
args=sys.argv

if len(args)<6: #if not enough args provided
        raise Exception("""Requred arguments are (in this order): 
                        - path to model output file
                        - species of interest (list separated by commas or "ALL") 
                        - top 'n' reactions to list for each species (e.g. '10' for top 10 most important reactions for each species)
                        - start point (in model time) for averaging over, or "START" 
                        - end point (in model time) for averaging over, or "END" """)

else: #note args[0] will be the name of this script
    out_path=args[1] #give path to model output directory
    species=args[2] #define species to be picked out
    top_n=int(args[3]) #define number of top reactions to list in output
    start_t=args[4]
    end_t=args[5]
    
#convert lists of species to list type
if "," in species:
    species=species.split(",")

#convert "start" and "end" into indexes (based on the timesteps of the model)
#to reference for averaging later
timelist=time_list(out_path+"/productionRates.output") #create list of the model timesteps

try:
    start_i=timelist.index(int(start_t))
except ValueError:
    if start_t=="START":
        start_i=0
    else:
        raise Exception("start value must be integer or 'START'")
#convert "end" into integer
try:
    end_i=timelist.index(int(end_t))
except ValueError:
    if end_t=="END":
        end_i=-1
    else:
        raise Exception("end value must be integer or 'END'")



#read production rates output file and record NO3 production rates
p_rates=read_p_rates(out_path+"/productionRates.output", species=species)
l_rates=read_l_rates(out_path+"/lossRates.output",species=species)

#get definitions of reactionNumbers
p_reactions=reaction_dict(out_path+"/productionRates.output")
l_reactions=reaction_dict(out_path+"/lossRates.output")



#calculate net production/loss for reversible reactions.

#Produce reaction dictionaries split into lists of products and reactants
split_p_reactions={}
for rp in p_reactions: #iterate through production reactions
    split_p_reactions[rp]=p_reactions[rp].split("=") #split string into reactants and products
    for i,s in enumerate(split_p_reactions[rp]):
        split_p_reactions[rp][i]=s.split("+") #split reactants and products into individual species
split_l_reactions={}
for rl in l_reactions: #iterate through loss reactions
    split_l_reactions[rl]=l_reactions[rl].split("=") 
    for i,s in enumerate(split_l_reactions[rl]):
        split_l_reactions[rl][i]=s.split("+")

#go through production reactions and look for reverse loss reaction
rev_rxns=[]
for rp in split_p_reactions:
    for rl in split_l_reactions:
        p_to_l=set(split_p_reactions[rp][0])==set(split_l_reactions[rl][1]) #do reactants of production reaction= products of loss?
        l_to_p=set(split_p_reactions[rp][1])==set(split_l_reactions[rl][0]) #do reactants of products production reaction= reactants of loss?
        if p_to_l and l_to_p:
            rev_rxns.append([rp,rl])

r_rates={} #create new dict to record reversible reactions
for s in p_rates: #initialise with keys for all species of interest
    r_rates[s]={}
for s in l_rates:
    r_rates[s]={}

pl_remove=[] #list to record indexes to remove from p_ and l_rates (any reactions added to r_rates)
for v in rev_rxns:
    for s in r_rates:
        try:
            prod=p_rates[s][v[0]]
            loss=l_rates[s][v[1]]
            diffs=[a-b for a, b in zip(prod, loss)]
            r_rates[s][v[0]]=diffs #record net production/loss (with the production reaction as key)
            
            pl_remove.append([s,v[0],v[1]]) #add species and reaction numbers to list for later removal
            
        except KeyError: #Reaction often won't exist for the species "s"
            pass 
#remove reactions in r_rates from p_rates and l_rates
for s,r1,r2 in pl_remove: 
    try:
        del p_rates[s][r1]
        del l_rates[s][r2]
    except KeyError:
        pass


#remove empty keys from r_rates
r_remove=[] #list to record empty indexes to remove from r_rates
for s in r_rates:
    if bool(r_rates[s].keys())==False:#see if dict has no keys (i.e. is empty)
        r_remove.append(s) #if so then remove it from the r_rates dict
    elif bool(r_rates[s].keys()):
        pass

for s in r_remove: #remove empty dict keys from r_rates
    del r_rates[s]
    
        
#create dictionary with same structure as production/loss/reversible rates 
#dict, but with average rates created from the list provided
avg_p_rates={}
for s in p_rates:
    avg_p_rates[s]={}
    for r in p_rates[s]:
        avg_p_rates[s][r]=mean(p_rates[s][r][start_i:end_i])

avg_l_rates={}
for s in l_rates:
    avg_l_rates[s]={}
    for r in l_rates[s]:
        avg_l_rates[s][r]=mean(l_rates[s][r][start_i:end_i])

avg_r_rates={}
for s in r_rates:
    avg_r_rates[s]={}
    for r in r_rates[s]:
        avg_r_rates[s][r]=mean(r_rates[s][r][start_i:end_i])


#sort the average rates by value
for s in avg_p_rates: #dictionary comprehension to sort dict by value
        avg_p_rates[s]={k: v for k, v in sorted(avg_p_rates[s].items(), 
                                                key=lambda item: item[1],
                                                reverse=True)}     
for s in avg_l_rates:
        avg_l_rates[s]={k: v for k, v in sorted(avg_l_rates[s].items(), 
                                                key=lambda item: item[1],
                                                reverse=True)}
for s in avg_r_rates: 
        avg_r_rates[s]={k: v for k, v in sorted(avg_r_rates[s].items(), 
                                                key=lambda item: abs(item[1]),
                                                reverse=True)}


#print the average rates for each species
print("-"*10)
print(f"Showing top {top_n} rate(s)")
print("-"*10)
print("PRODUCTION")
print()
for s in avg_p_rates:
    print()
    print(f"{s} Production")
    top=take(top_n, avg_p_rates[s].items()) #select top "n" reactions 
    total_sum = sum(avg_p_rates[s].values())
    for i,r in top:
        print(f"{p_reactions[i]}: average rate= {r:.3e} molecule cm-3 s-1 ({(r/total_sum)*100:.3g}%)")

print("-"*10)
print("LOSS")
print()
for s in avg_l_rates:
    print()
    print(f"{s} Loss")
    top=take(top_n, avg_l_rates[s].items()) #select top "n" reactions
    total_sum = sum(avg_l_rates[s].values())
    for i,r in top:
        print(f"{l_reactions[i]}: average rate= {r:.3e} molecule cm-3 s-1 ({(r/total_sum)*100:.3g}%)")

print("-"*10)
print("REVERSIBLE (+ve= net production, -ve= net loss)")
print()
for s in avg_r_rates:
    print()
    print(f"{s}")
    top=take(top_n, avg_r_rates[s].items()) #select top "n" reactions 
    for i,r in top: #uses p_reactions
        print(f"{p_reactions[i]}: average rate= {r:.3e} molecule cm-3 s-1")

