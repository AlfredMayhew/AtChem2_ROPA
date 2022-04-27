"""script to produce the total rate of reaction for a species (in s-1) from 
AtChem2 output files"""
#imports
from read_output.read_rate_output import read_p_rates, read_l_rates, time_list
from read_output.read_conc_output import conc_output
import sys
import pandas as pd

#put command-line args into list
args=sys.argv

if len(args)<3: #if not enough args provided
        raise Exception("""Requred arguments are (in this order): 
                        - path to model output file
                        - species of interest (individual species only)""")

else: #note args[0] will be the name of this script
    out_path=args[1] #give path to model output directory
    species=args[2] #define species to be picked out


#make list of model time steps
times=time_list(out_path+"/productionRates.output")

#read production rates output file and record NO3 production rates
p_rates=read_p_rates(out_path+"/productionRates.output", species=species)
l_rates=read_l_rates(out_path+"/lossRates.output",species=species)

#get concentration for species at each timestep (used to convert rate constants 
#from molecule_cm-3_s-1 to s-1)
concs=conc_output(out_path+"/speciesConcentrations.output",species)

#sum all production rates for the species of interest (for each timestep)
total_prod=pd.DataFrame(index=times,columns=[species])
for i,t in enumerate(times):
    temp_sum=0
    for r in p_rates[species]:
        temp_sum=temp_sum+p_rates[species][r][i]
    total_prod.loc[t,species]=temp_sum/concs[species][t] #divide by conc to convert to s-1
#sum all loss rates for the species of interest (for each timestep)
total_loss=pd.DataFrame(index=times)
for i,t in enumerate(times):
    temp_sum=0
    for r in l_rates[species]:
        temp_sum=temp_sum+l_rates[species][r][i]
    total_loss.loc[t,species]=temp_sum/concs[species][t] #divide by conc to convert to s-1

#take difference between production and loss to get overall rate (+ve=production, -ve=loss)
rate_diff=total_prod-total_loss

#output to file
rate_diff.to_csv("temp_total_rates", sep=" ")
