#!/usr/bin/env python
# coding: utf-8


import numpy as np
import operator as op
import random as rd
import matplotlib.pyplot as plt 
import pandas as pd
from collections import Counter


#     1. GENETIC DRIFT

# a) Trace the allele frequencies over time


pA=0.4
pB=1-pA
gen_nb=200
sim_nb=10
N=100

def freq(N, p,gen_nb):
    
    pop = np.zeros((N,gen_nb), dtype=int) # Use pop to store the population structure of all each generation
   
    ind=rd.sample(range(N),int(p*N)) #follow up a 
    for i in ind:
        pop[i][0]=1
    
    #fréquence des allèles 1(A)
    freq=[p]
    proba_array = [1/N]*N
     
    for t in range(1,gen_nb):
        pop[:,t] = np.random.choice(pop[:,t-1], size=N, replace=True, p=proba_array)
    
        freq.append(sum(pop[:,t])/N)
        
    return pop,freq



for i in range(sim_nb):
    pop,f= freq(N, pA, gen_nb)
    title = 'allele A frequency over '+str(gen_nb)+' generations, N='+str(N)
    plt.title(title,fontsize = 15)
    plt.ylabel('allele A frequency',fontsize = 13)
    plt.xlabel('generations',fontsize = 13)
    plt.plot(np.arange(gen_nb), f)
    
plt.show()


# b/ Fixation probability


probaA_list = np.arange(0, 1.1, 0.1)
def fix_proba(N, p, gen_nb, sim_nb):    
    fix_num = 0
    
    for i in range(sim_nb):
        pop,freqs = freq(N, p, gen_nb)
    
        if 1 in freqs: # means the allele is fixed
            fix_num += 1
                
    return fix_num/sim_nb

fix=[]
#for the different proba values:
for p in range(len(probaA_list)):
    fix.append(fix_proba(N, probaA_list[p], 500, 50))

    
plt.plot(np.arange(0, 1.1, 0.1),fix)    
plt.title('Fixation probablility of allele A in dependence of intitial A frequency')
plt.ylabel('Probability of fixation of allele A')
plt.xlabel('p, Frequency of allele A ')

plt.show()




# 2/ COALESCENCE

p0=0.4
q0=1-p0
gen_nb=500
simulations_nb=10
N=100

def coalescence(N,n,p, simulations_nb):    
    #search the first generation at which the first coalsecent event appears
    #pick random individuals to be in the sample
    sample=rd.sample(range(N), n)
    
    #generation 
    somme_t=0
    for sim in range(simulations_nb):
        for t in range(gen_nb-1, -1, -1):
            ind_list=[]
            for indiv in sample:
                ind_list.append(np.random.randint(0, N))

            unique=list(set(ind_list))

            if len(unique)<len(ind_list):
                #print('y a eu coalescence, deux alleles ont un ancetre commun', t)
                somme_t+=(gen_nb-t)
                break
                #print('pas de coalescence')
        
    return(int(somme_t/simulations_nb))


n_sam=[2, 3,4,5]

x = n_sam
y1=[coalescence(N,n,p0, 30) for n in x]
y2 = [int((2*N)/(n*(n-1))) for n in x]

fig, ax = plt.subplots(1, figsize=(9, 7))
ax.plot(x,y1, label = "Simulated")
ax.plot(x,y2, label = "Expected")
plt.title("Estimated number of generations until the first coalescence", fontsize = 12)
plt.xlabel("Sample size", fontsize = 12)
plt.ylabel("Number of generations", fontsize = 12)
plt.legend(fontsize = 12)
plt.show()

for nb in n_sam:
    print('for', nb, ' individuals in sample, we have first coalescent event after',coalescence(N,nb,p0, 30) , 'generations')



#The model seen in class
n=[2, 3,4,5]
for nb in n:
    time=int((2*N)/(nb*(nb-1)))
    print(time)


