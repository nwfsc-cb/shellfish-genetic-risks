import argparse

import pandas as pd

import numpy as np

import itertools

import random

import time

import datetime

import os

import subprocess

import simuOpt

simuOpt.setOptions(alleleType='long', quiet=True)

import simuPOP as sim
    

def shellfishrisk(batch, reps, coreid, freq = None,pre_farm_years = 50, farm_years = 50, post_farm_years = 50,wild_N_init = 3000,rec_a = 750,sd_recruit = 0.01,
numWildOffspring_par = 10,wild_mig_rate_L = 0.035/12, wild_mig_rate_J = 0, wild_mig_rate_A = 0,seed_batch_size = 250,
numFarmOffspring_par = 20,
sd_seed = 0.01,
farm_nonB_recruit_factor = 0.25,
gamEsc_rec = 20,
source_idx = 0,
local_wild_idx = 0,
L_escape_rate = 0.07,
J_escape_rate = 0/12,
A_escape_rate = 0/12,
numGameteEscapeOffspring_par = 10
):
    '''An individual-based model of shellfish production, escape, and genetic impacts to wild populations'''
    
    batch_dir = os.path.join(os.getcwd(), 'results',batch)

    if not os.path.exists(batch_dir):
        os.makedirs(batch_dir, exist_ok=True)

    reps = int(reps)
    
    coreid = int(coreid)
    
    pre_farm_years = int(pre_farm_years)
    
    farm_years = int(farm_years)
    
    post_farm_years = int(post_farm_years)
    
     # wild popdy
    num_wild_pops = 3 # number wild pops (model structure currently only handles 3)
    wild_N_init = int(wild_N_init) # initial wild pop size
    rec_a = int(rec_a) # annual recruitment event size; density independent
    sd_recruit = 0.01 # standard deviation in wild recruitment
    numWildOffspring_par = int(numWildOffspring_par) # wild family size
    wild_mig_rate_L = 0.035/12 # larvae monthly mig rate, but currently only happens in one month
    wild_mig_rate_J = 0 # juvenile
    wild_mig_rate_A = 0 # adult
    
    # farm popdy
    seed_batch_size = int(seed_batch_size) # farm recruitment, i.e., number of seed planted on farm at once
    numFarmOffspring_par = int(numFarmOffspring_par) # farm family size
    sd_seed = 0.01 # standard deviation in farm recruitment / scale of seed production
    farm_nonB_recruit_factor = 0.25 # factor governing recruitment due to F2 adults on farm
    gamEsc_rec = int(gamEsc_rec) # monthly wild recruitment due to gamete escape
    source_idx = 0 # index of wild pop to collect broodstock from
    local_wild_idx = 0 # index of subpop that farm will escape to 
    # L_escape_rate = 0.07 # larvae monthly escape rate, but currently only happens in one month
    # J_escape_rate = 0/12 # juvenile
    # A_escape_rate = 0/12 # adult
    numGameteEscapeOffspring_par = int(numGameteEscapeOffspring_par) # hybrid / gamete escape family size
    
    
    progress_report = open(os.path.join(batch_dir,"progress.txt"), 'w')

    # print('Reading in parameter values')
    
    progress_report.write("Reading in parameter values\n")
    
    progress_report.flush()
    

    ################################################################################
#### ------------------------ Import modules ------------------------------ ####
################################################################################

 
    ################################################################################
    #### ----------------------- Define functions ----------------------------- ####
    ################################################################################

    
    ###########################################################################
    #### --- Genomic architecture & fitness fns

    def initEqAFs(loci_for_fst, eq_afs_path):
        '''Initialize FST loci allele frequencies  using equilibrium
        allele frequencies file'''
        eq_afs = pd.read_csv(eq_afs_path, delimiter = ",")
        random_year = np.random.choice(a=list(set(list(eq_afs['Year']))), size=1)[0]
        eq_afs_random_year = eq_afs[eq_afs['Year']==random_year]
        df = eq_afs_random_year[eq_afs_random_year['Adaptive']==False]
        for sp_idx, sp_name in enumerate(['wild1', 'wild2', 'wild3']):
            for l_idx in loci_for_fst:  
                this_sp_locus_afs = []
                for allele in [0,1]:
                    this_row = df.loc[(df['Subpop']==sp_name) & (df['Locus_index']==l_idx) & (df['Allele']==allele)]
                    this_sp_locus_afs.append(float(this_row['AlleleFrequency']))
                sim.initGenotype(pop, freq=this_sp_locus_afs, loci=l_idx, subPops=sp_idx)
    
    def make_fit_dict(stages, num_aloci_stage, sp_names, s, z, deg_add):
        '''Get dictionary of relative fitness of genotype, per locus, per life stage, per subpopulation.
        Parameter deg_add is a dictionary with each stage as keys and degree additive effects as values
        with 0 = complete dominance and 1 = complete additive effects. Simpler fitness landscape, with: 
        
        farm       wild1     wild2      wild3
        0     1    0   1     0   1      0   1
        1-z   1    1  1-s    1  1-2s    1  3-s
        '''
                 
        # create fitness dictionary, per locus, per life history stage, per subpopulation
        fit_dict = {}
        for sp_name in sp_names: # for each subpop
            fit_dict[sp_name] = {}
            for stage in stages: # for each life history stage: larva, juvenile, adult
                fit_dict[sp_name][stage] = {}
                deg_add_stage = deg_add[stage] # retrieve degree additive effects from dictionary
                for l_idx in range(num_aloci_stage): # for each locus
                    fit_dict[sp_name][stage][l_idx] = {}
                    allele_effects = {}                
                    if sp_name !='farm': # if advantage is first allele, 0
                        allele_effects[0] = 1/2
                        allele_effects[1] = (1-s)/2
                        fit_dict[sp_name][stage][l_idx]['0'] = allele_effects[0] + allele_effects[0]
                        fit_dict[sp_name][stage][l_idx]['2'] = allele_effects[0] + (1-s*deg_add_stage)/2
                        fit_dict[sp_name][stage][l_idx]['1'] = allele_effects[1] + allele_effects[1]
                    elif sp_name == 'farm': # if advantage is second allele, 1
                        allele_effects[0] = (1-z)/2
                        allele_effects[1] = 1/2
                        fit_dict[sp_name][stage][l_idx]['0'] = allele_effects[0] + allele_effects[0]
                        fit_dict[sp_name][stage][l_idx]['2'] = allele_effects[1] + (1-z*deg_add_stage)/2
                        fit_dict[sp_name][stage][l_idx]['1'] = allele_effects[1] + allele_effects[1]
        return(fit_dict)
        
    def init_geno_aloci(pop, sp_names, num_aloci, fit_dict, good_start_AF, bad_start_AF):
        '''Initialize starting allele frequencies depending on which 
        allele is beneficial per subpop and life stage.'''
        for sp_idx, sp_name in enumerate(sp_names): # for each subpop
            for l_idx in range(num_aloci): # for each adaptive locus
                stage = stages[l_idx//num_aloci_stage] # depending on locus, which stage affected
                within_stage_l_idx = l_idx % num_aloci_stage # get w/in stage locus index          
                if fit_dict[sp_name][stage][within_stage_l_idx]['0'] == 1: # order based on fit_dict
                    freq = [good_start_AF, bad_start_AF]
                else:
                    freq = [bad_start_AF, good_start_AF] 
                sim.initGenotype(pop, loci = l_idx, freq = freq, subPops = sp_idx) 
    
    def trans_geno(geno):
        '''Translate simupop genotye to 0, 1, or 2;
        0 = 00, 1 = 11, 2 = 01 or 10'''
        if 0 in geno and 1 in geno:
            return('2')
        elif 0 in geno and 1 not in geno:
            return('0')
        elif 1 in geno and 0 not in geno:
            return('1')
            
    def get_genos_aloci(ind):
        '''Return a list object of genotypes, in 0,1,2 format,
        in order of loci for adaptive loci.'''
        genos = []
        sim_genos = ind.genotype()
        num_sim_genos = int(len(sim_genos))
        first_alleles = sim_genos[0:num_aloci]
        second_alleles = sim_genos[int(num_sim_genos/2):int(num_sim_genos/2)+num_aloci]
        for zipped_geno in zip(first_alleles, second_alleles):
            genos.append(trans_geno(zipped_geno))     
        return(genos)
    
    def get_fitness(ind, pop):
        '''Get fitness for this individual, depending on genotype and environment'''
        stage = get_stage(ind) # get life history stage
        stage_index = stages.index(stage) # and index
        sp_names = pop.subPopNames()
        for sp_name in sp_names: # get subpop index of this individual
            if ind in pop.individuals(sp_name):
                this_sp_name = sp_name # get subpop name
        genos = get_genos_aloci(ind) 
        start_this_stage = stage_index*num_aloci_stage
        stop_this_stage = (stage_index+1)*num_aloci_stage
        genos_this_stage = genos[start_this_stage:stop_this_stage] # just genos affecting this stage
        fitness = 0
        for geno_count, geno_item in enumerate(genos_this_stage):
            locus_effect = fit_dict[this_sp_name][stage][geno_count][geno_item]
            fitness += locus_effect 
        fitness = fitness / len(genos_this_stage)    
        return(fitness)
    
    def assign_fitness_all(pop):
        '''Assign fitness for all individuals, depending on genotype and environment'''
        for ind in pop.individuals():
            stage = get_stage(ind) # get life history stage
            stage_index = stages.index(stage) # and index
            sp_names = pop.subPopNames()
            for sp_name in sp_names: # get subpop index of this individual
                if ind in pop.individuals(sp_name):
                    this_sp_name = sp_name # get subpop name
            genos = get_genos_aloci(ind) 
            start_this_stage = stage_index*num_aloci_stage
            stop_this_stage = (stage_index+1)*num_aloci_stage
            genos_this_stage = genos[start_this_stage:stop_this_stage] # get genos affecting this stage
            fitness = 0
            for geno_count, geno_item in enumerate(genos_this_stage): 
                locus_effect = fit_dict[this_sp_name][stage][geno_count][geno_item]
                fitness += locus_effect 
            fitness = fitness / len(genos_this_stage)
            ind.fitness = fitness
    
    ###########################################################################
    #### --- Demographic & life history fns
    
    def wo_farm_recruitment(pop):
        '''Provide next time step wild subpopulation sizes for before or after 
        farm phase; density-independent recruitment'''
        projected_sp_sizes = [int(Nt + rec_a*np.exp(sd_recruit*clim_dev[year_counter] - (sd_recruit*sd_recruit)/2)) \
             for idx, Nt in enumerate(list(pop.subPopSizes()))]
        return(projected_sp_sizes)
    
    def w_farm_recruitment(pop):
        '''Provide next time step wild subpopulation sizes for during farm phase; 
        density-independent recruitment'''
        current_sp_sizes = list(pop.subPopSizes())
        projected_sp_sizes = []
        for idx, Nt in enumerate(current_sp_sizes):      
            if idx !=farm_idx: # if not farm
                recruitment = rec_a*np.exp(sd_recruit*clim_dev[year_counter] - (sd_recruit*sd_recruit)/2)
                Nt1 = Nt + recruitment
                projected_sp_sizes.append(Nt1)    
            else: # if farm
                # see if at least 1 male and female nonBstock adults
                nonB_farm_adults = False
                farm_male, farm_female = False, False
                for ind in pop.individuals(farm_idx):
                    if ind.age >= repro_age and ind.broodstock == 0:
                        if ind.sex() == 1:
                            farm_male = True
                        else:
                            farm_female = True
                if farm_male and farm_female:
                    nonB_farm_adults = True  
                if nonB_farm_adults: # if are enough nonBstock adults
                    wild_recruitment = rec_a*np.exp(sd_recruit*clim_dev[year_counter] - (sd_recruit*sd_recruit)/2)
                    nonB_farm_recruitment = recruitment * farm_nonB_recruit_factor
                    Nt1 = Nt + nonB_farm_recruitment
                    projected_sp_sizes.append(Nt1)
                else: # if no, no recruitment on farm from nonB adults
                    projected_sp_sizes.append(Nt)
        projected_sp_sizes = [int(i) for i in projected_sp_sizes]    
        return(projected_sp_sizes)
            
    def get_mort_fitfactor(ind_fitness):
        '''Get fitness factor to multiply against mortality rate, making mortality selective'''
        output_start = 1 + var_in_mort
        output_end = 1 - var_in_mort
        input_start = 1 - z
        input_end = 1
        fitfactor = output_start + ((output_end - output_start) / (input_end - input_start)) * (ind_fitness - input_start)
        return(fitfactor)
        
    def mortality(pop, lar_annual_M, juv_annual_M, adult_annual_M, farm_reduced_mort):
        '''Kill individuals based on age based mortality rates & fitness,
        such that probability of death = life stage probability * fitness factor,
        where fitness factor =  fitness value constrained'''
        for sp_idx, sp_name in enumerate(pop.subPopNames()):
            if sp_name !='farm': # if not farm
                farm_wild_factor = 1
            else: # if farm
                farm_wild_factor = farm_reduced_mort   
            for ind in pop.individuals(sp_idx): 
                draw = random.random()
                fitfactor = get_mort_fitfactor(ind_fitness=ind.fitness)
                if ind.age < settle_age and draw < (lar_annual_M/12)*fitfactor*farm_wild_factor:
                    ind.kill = 1
                elif settle_age <= ind.age < repro_age and draw < (juv_annual_M/12)*fitfactor*farm_wild_factor: 
                    ind.kill = 1
                elif repro_age <= ind.age < max_age and draw < (adult_annual_M/12)*fitfactor*farm_wild_factor:
                    ind.kill = 1
                elif ind.age >= max_age:
                    ind.kill = 1   
            pop.removeSubPops([(sp_idx,vsp_n2i['to_kill'])]) # kill individuals
    
    def mark_broodstock(pop, source_idx, source_name, mat_vsp_name, num_broodstock, eq_sex):
        '''Identify and mark potential broodstock. Option of equal 
        numbers from each sex, or random and at least one from each sex.
        Marked as 1 if broodstock, otherwise info field remains 0.'''
            
        # store female and male idx that are in source pop & sexually mature
        mature_F_idx = []
        mature_M_idx = []
        for index in range(pop.subPopBegin(source_idx), pop.subPopEnd(source_idx)):
            ind = pop.individual(index)
            if ind.age >= repro_age:
                if ind.sex() == 2: 
                    mature_F_idx.append(index) # store index if mature and female
                else: 
                    mature_M_idx.append(index) # store index if mature and male
                    
        # option of equal sexes, otherwise, random and at least one male and one female
        idx_to_mark = []
        random.shuffle(mature_F_idx) # shuffle idx to get random sample below
        random.shuffle(mature_M_idx) 
        if eq_sex == True:
            sex_sample_sizes = [int(num_broodstock/2), num_broodstock-int(num_broodstock/2)]
            random.shuffle(sex_sample_sizes) # ensure not consistently more males or more females if odd num bstock  
            idx_to_mark +=  mature_F_idx[0:sex_sample_sizes[0]] # store female idx to mark
            idx_to_mark += mature_M_idx[0:sex_sample_sizes[1]]  # store male idx to mark
        else:
            idx_to_mark.append(mature_F_idx[0]) # get one female (make sure at least one female)
            idx_to_mark.append(mature_M_idx[0]) # get one male (make sure at least one male)
            combined = mature_F_idx[1:] + mature_M_idx[1:] # combine rest of indeces for random sample, regardless of sex
            combined_grab_idx = random.sample(combined,num_broodstock-2) # remainder of brstock, total-2
            idx_to_mark += combined_grab_idx
            
        # mark broodstock by switching infofield 'broodstock' to 1 for selected individuals
        for i in idx_to_mark:
            pop.individual(i).broodstock = 1
            
    def collect_bstock(pop, sp_idcs, farm_idx):
        '''Migrate marked broodstock from source subpopulation
        to farm subpopulation.'''
        for sp_idx in sp_idcs:
            for ind in pop.individuals([sp_idx]):
                if ind.broodstock == 1:
                    ind.migrate_to = farm_idx
                else:
                    ind.migrate_to = sp_idx
        sim.migrate(pop, mode=sim.BY_IND_INFO)
                
    def escape_J_A(pop, J_escape_rate, A_escape_rate, w_farm_wld_sp_idcs, local_wild_idx):
        '''Allow for escape from farm to local wild subpopulation for life stages expected to escape
        randomly across the year, J and A, such as due to storm damage'''
        for ind in pop.individuals('farm'):
            if ind.age < settle_age: # if L
                ind.migrate_to = farm_idx # mark individuals to remain on farm   
            elif settle_age <= ind.age < repro_age: # if J
                if np.random.choice(a=['escape','remain'], size=1, p=[J_escape_rate, 1-J_escape_rate]) == 'escape':
                    ind.migrate_to = local_wild_idx # mark individuals to migrate to local wild
                else:
                    ind.migrate_to = farm_idx # mark individuals to remain on farm   
            elif ind.age >= repro_age and ind.broodstock == 0: # if A and not broodstock
                if np.random.choice(a=['escape','remain'], size=1, p=[A_escape_rate, 1-A_escape_rate]) == 'escape':
                    ind.migrate_to = local_wild_idx # mark individuals to migrate to local wild
                else:
                    ind.migrate_to = farm_idx # mark individuals to remain on farm   
        for sp_idx in w_farm_wld_sp_idcs: 
            for ind in pop.individuals([sp_idx]):
                ind.migrate_to = sp_idx # mark wild individuals to remain
        sim.migrate(pop, mode=sim.BY_IND_INFO) # migrate farm animals to local wild
        
    def escape_L(pop, L_escape_rate, w_farm_wld_sp_idcs, local_wild_idx):
        '''Allow for escape from farm to local wild subpopulation for life stages expected to escape
        randomly across the year'''
        for ind in pop.individuals(farm_idx):
            if ind.age < settle_age and ind.F == 2: # if L and F2
                if np.random.choice(a=['escape','remain'], size=1, p=[L_escape_rate, 1-L_escape_rate]) == 'escape':
                    ind.migrate_to = local_wild_idx # mark individuals to migrate to local wild
                else:
                    ind.migrate_to = farm_idx # mark individuals to remain on farm   
            elif settle_age <= ind.age < repro_age: # if J
                ind.migrate_to = farm_idx # mark individuals to remain on farm   
            elif ind.age >= repro_age and ind.broodstock == 0: # if A and not broodstock
                ind.migrate_to = farm_idx # mark individuals to remain on farm   
        for sp_idx in w_farm_wld_sp_idcs: 
            for ind in pop.individuals([sp_idx]):
                ind.migrate_to = sp_idx # mark wild individuals to remain
        sim.migrate(pop, mode=sim.BY_IND_INFO) # migrate farm animals to local wild
        
    def gameteEscapeChooser(pop, subPop):
        '''Identify parent pairs that will hybridize across farm and wild1'''
        males = {}
        females = {}
        origin_indeces = [0, farm_idx]
        for origin_index in origin_indeces: # 0 = wild 1 always
            males[origin_index] = [x for x in pop.individuals(hybrid_index) \
                if x.sex() == 1 and x.return_to == origin_index]
            females[origin_index] = [x for x in pop.individuals(hybrid_index) \
                if x.sex() == 2 and x.return_to == origin_index]
        while True:
            random.shuffle(origin_indeces)                        
            yield males[origin_indeces[0]][random.randint(0, len(males[origin_indeces[0]]) - 1)], \
                females[origin_indeces[1]][random.randint(0, len(females[origin_indeces[1]]) - 1)]
                            
    def scaleGameteEscape(pop):
        '''Provide new subpop size for temp hybrid pop, to account for
        larve from gamete escape'''
        current_sp_sizes = list(pop.subPopSizes())
        projected_sp_sizes = []
        for idx, Nt in enumerate(current_sp_sizes):
            if idx == hybrid_index:
                Nt1 = Nt + gamEsc_rec
                projected_sp_sizes.append(int(Nt1))
            else:
                projected_sp_sizes.append(int(Nt))
        return(projected_sp_sizes)
                     
    def get_fitness_hybrid(ind, pop):
        '''Get fitness, but for temporary hybrid pop; treating hybrid pop same
        as wild1 environment.'''
        stage = get_stage(ind) # get life history stage
        stage_index = stages.index(stage) # and index
        this_sp_name = 'wild1' # get subpop name
        genos = get_genos_aloci(ind) 
        start_this_stage = stage_index*num_aloci_stage
        stop_this_stage = (stage_index+1)*num_aloci_stage
        genos_this_stage = genos[start_this_stage:stop_this_stage] # get genotypes that affect fitness in this stage
        fitness = 0
        for geno_count, geno_item in enumerate(genos_this_stage): # iterate through (character) 0, 1, 2s of string of genotypes
            mag = 1 / len(genos_this_stage)
            locus_effect = fit_dict[this_sp_name][stage][geno_count][geno_item]
            fitness += locus_effect # get magnitude of effect X locus fitness effect, add to running sum
        fitness = fitness / len(genos_this_stage)
        return(fitness)
        
    def wild_mig(pop, farm_boo, farm_idx, wild_names, wild_mig_rate_L, 
                 wild_mig_rate_J, wild_mig_rate_A, w_farm_wld_sp_idcs):
        '''Allow for migration among wild subpopulations. Wild migration rate is how many
        individuals leave a subpopulation, but some subpopulations can migrate to more than one
        subpopulation.'''  
        if farm_boo == False: # if not farm year
            sp_idcs = list(range(len(wild_names)))
        else: # if farm year
            sp_idcs = w_farm_wld_sp_idcs
            for ind in pop.individuals(farm_idx): # farm animals stay in place
                ind.migrate_to = farm_idx 
        for sp_idx in sp_idcs:
            for ind in pop.individuals(sp_idx):
                if ind.age < settle_age: # if L
                    if np.random.choice(a=['migrate','remain'], size=1, p=[wild_mig_rate_L, 1-wild_mig_rate_L]) == 'migrate': 
                        if sp_idx == sp_idcs[0]: # if first wild subpop, only migrates to next futher wild
                            ind.migrate_to = sp_idcs[1]
                        elif sp_idx == sp_idcs[-1]: # elif last wild subpop, only migrates to next nearest wild
                            ind.migrate_to = sp_idcs[-2]
                        else: # else, if a middle subpop, migrates to one of two neighboring subpops
                            if random.random() > 0.5:
                                ind.migrate_to = sp_idcs[sp_idcs.index(sp_idx)+1]
                            else:
                                ind.migrate_to = sp_idcs[sp_idcs.index(sp_idx)-1]
                    else: # doesn't migrate
                        ind.migrate_to = sp_idx 
                elif settle_age <= ind.age < repro_age: # if J
                    if np.random.choice(a=['migrate','remain'], size=1, p=[wild_mig_rate_J, 1-wild_mig_rate_J]) == 'migrate': 
                        if sp_idx == sp_idcs[0]: # if first wild subpop, only migrates to next futher wild
                            ind.migrate_to = sp_idcs[1]
                        elif sp_idx == sp_idcs[-1]: # elif last wild subpop, only migrates to next nearest wild
                            ind.migrate_to = sp_idcs[-2]
                        else: # else, if a middle subpop, migrates to one of two neighboring subpops
                            if random.random() > 0.5:
                                ind.migrate_to = sp_idcs[sp_idcs.index(sp_idx)+1]
                            else:
                                ind.migrate_to = sp_idcs[sp_idcs.index(sp_idx)-1]
                    else: # doesn't migrate
                        ind.migrate_to = sp_idx  
                elif ind.age >= repro_age: # if A
                    if np.random.choice(a=['migrate','remain'], size=1, p=[wild_mig_rate_A, 1-wild_mig_rate_A]) == 'migrate': 
                        if sp_idx == sp_idcs[0]: # if first wild subpop, only migrates to next futher wild
                            ind.migrate_to = sp_idcs[1]
                        elif sp_idx == sp_idcs[-1]: # elif last wild subpop, only migrates to next nearest wild
                            ind.migrate_to = sp_idcs[-2]
                        else: # else, if a middle subpop, migrates to one of two neighboring subpops
                            if random.random() > 0.5:
                                ind.migrate_to = sp_idcs[sp_idcs.index(sp_idx)+1]
                            else:
                                ind.migrate_to = sp_idcs[sp_idcs.index(sp_idx)-1]
                    else: # doesn't migrate
                        ind.migrate_to = sp_idx 
        sim.migrate(pop, mode=sim.BY_IND_INFO) # migrate farm animals to local wild
            
    def harvest(pop, farm_idx, harvest_rate, harvest_age, max_harvest_age, harvest_this_year):
        '''Identify a number of animals to harvest defined by harvest_rate * number_harvestable.
        Harvest all animals that are at least max_harvest_age. For the remaining animals to harvest,
        select randomly from animals above harvest_age. Remove from the population to simulate harvest.'''  
        num_harvestable = 0 # over harvest_age
        IDs_harvestable = []
        IDs_to_harvest = []
        for ind_idx in range(pop.subPopBegin(farm_idx), pop.subPopEnd(farm_idx)):
            ind = pop.individual(ind_idx)
            if ind.age >= harvest_age and ind.broodstock == 0:
                num_harvestable += 1   
                IDs_harvestable.append(ind.ind_id)
                if ind.age >= max_harvest_age:
                    IDs_to_harvest.append(ind.ind_id)   
        if len(IDs_to_harvest) < harvest_rate*num_harvestable: # harvest remaining over max age
            more_to_harvest = round(harvest_rate*num_harvestable - len(IDs_to_harvest))
            random.shuffle(IDs_harvestable)
            IDs_to_harvest += IDs_harvestable[0:more_to_harvest]    
        pop.removeIndividuals(IDs=IDs_to_harvest)
        harvest_this_year += len(IDs_to_harvest)
        
    def seed_production(pop, farm_idx, seed_batch_size, sd_seed):
        '''Provide next time step subpopulation size for farm,
        and give back current sp sizes for wild subpops'''
        current_sp_sizes = list(pop.subPopSizes())
        projected_sp_sizes = []
        for idx, Nt in enumerate(current_sp_sizes):
            if idx !=farm_idx: # if not farm
                projected_sp_sizes.append(Nt)
            else: # if farm
                Nt1 = Nt + int(seed_batch_size*np.exp(sd_seed-(sd_seed*sd_seed)/2))
                projected_sp_sizes.append(Nt1)
        return(projected_sp_sizes)
    
    def rem_bstock(pop):
        '''Remove broodstock after breeding within PyOperator'''
        pop.removeSubPops([(farm_idx,vsp_n2i['bstock'])]) # remove broodstock from farm population
        return(True)
        
    ###########################################################################
    #### --- Parameter deriving functions
            
    def get_w_farm_wld_sp_idcs(num_wild_pops, source_idx):
        '''Get list of wild pop subpop indeces, which
        depends on the broodstock source_idx. '''
        wld_sp_idcs = list(range(num_wild_pops+1))
        del wld_sp_idcs[source_idx+1]
        return(wld_sp_idcs)
    
    def get_sp_names(source_idx, farm_boo, wild_names):
        '''Get list of subpop names, which
        depends on whether the farm is included, 
        and the broodstock source_idx.'''
        if farm_boo == False: # of not farm year
            sp_names = wild_names
        elif farm_boo == True: # if farm year
            sp_names = []
            name_idx = 0
            for i in range(len(wild_names)):
                if i == source_idx:
                    sp_names.append(wild_names[name_idx])
                    sp_names.append('farm')
                    name_idx += 1
                else:
                    sp_names.append(wild_names[name_idx])
                    name_idx += 1
        return(sp_names)
    
    def make_vsp_name_to_index_dict(pop):
        '''Make a dictionary, with each key as 
        virtual subpopname, and each value as virtual
        subpop index. Assumes > 1 subPop'''
        vsp_n2i = {}
        for i in range(pop.numVirtualSubPop()):
            n = pop.subPopName([0,i]).split(' - ')[1].strip()
            vsp_n2i[n] = i
        return(vsp_n2i)
    
    def get_vsp_map(vsp_n2i):
        '''Get new vsp map and names; workaround because
        product splitter needed for nonB_adults but complicates
        other splitters'''
        
        raw_vsp_names = list(vsp_n2i.keys())
        desired_vsps = ['L', 'J', 'A',
         'leave', 'to_harvest',
         'not_bstock', 'bstock',
         'to_survive', 'to_kill',
         'nonB_adults']
        
        desired_vsp_dict = {}
        for desired_vsp in desired_vsps:
            this_vsp_idcs = []
            if desired_vsp == 'nonB_adults':
                for raw_vsp_idx, raw_vsp_name in enumerate(raw_vsp_names):
                    raw_as_list = [name.strip() for name in raw_vsp_name.split(',')]
                    if 'A' in raw_as_list and 'not_bstock' in raw_as_list:
                        this_vsp_idcs.append(raw_vsp_idx)
                desired_vsp_dict[desired_vsp] = this_vsp_idcs
            else:
                for raw_vsp_idx, raw_vsp_name in enumerate(raw_vsp_names):
                    raw_as_list = [name.strip() for name in raw_vsp_name.split(',')]
                    if desired_vsp in raw_as_list:
                        this_vsp_idcs.append(raw_vsp_idx)
                desired_vsp_dict[desired_vsp] = this_vsp_idcs
        return(desired_vsp_dict)
    
    def get_sps_to_mate_w_farm(wild_or_farm, farm_idx, wild_idcs):
        '''Get list of subpops indeces (as tuples) to mate.'''
        sps_to_mate = []    
        if wild_or_farm == 'wild': # just wild pops
            for wild_idx in wild_idcs:
                sps_to_mate.append((wild_idx, vsp_n2i['A']))
        elif wild_or_farm == 'farm': # just broodstock for seed production
            sps_to_mate.append((farm_idx, vsp_n2i['bstock']))
        elif wild_or_farm == 'wild_and_nonBfarm': # wild pops and non broodstock adults on farm
            for wild_idx in wild_idcs: 
                   sps_to_mate.append((wild_idx, vsp_n2i['A']))
            sps_to_mate.append((farm_idx, vsp_n2i['nonB_adults']))
        else:
            print('wild_or_farm parameter must be either "wild" or "farm" or "wild_and_nonBfarm"')
        return(sps_to_mate)
    
    def get_wo_farm_sps_to_mate(num_wild_pops):
        '''Get subpops to mate based on number of wild subpops and mature vsps'''
        return([(sp_idx, vsp_n2i['A']) for sp_idx in range(num_wild_pops)])
    
    def get_sps_to_clone(sp_idcs, age_vsp_idcs):
        '''Get list of subpops to clone'''
        return([x for x in itertools.product(sp_idcs, age_vsp_idcs)])
    
    def get_stage(ind):
        '''Get life history stage from ind age'''
        if ind.age < settle_age:
            stage = 'L'
        elif settle_age <= ind.age < repro_age:
            stage = 'J'     
        elif ind.age >= repro_age:
            stage = 'A'
        return(stage)
    
    ###########################################################################
    #### --- Response variable functions
            
    def store_afs(pop, rep, year_counter, farm_boo, farm_idx, afs_dict, loci_to_track):
        '''Stores allele frequences into dictionary'''
        if farm_boo == True:
            for herein_idx, sp_idx in enumerate(w_farm_wld_sp_idcs):  
                for locus in biallelic_loci:
                    for allele in [0,1]:
                        afs_dict[rep][year_counter][herein_idx][locus][allele] = \
                        pop.vars()['subPop'][sp_idx]['alleleFreq'][locus][allele]        
            for locus in biallelic_loci:            
                for allele in [0,1]:
                    afs_dict[rep][year_counter]['farm'][locus][allele] = \
                    pop.vars()['subPop'][farm_idx]['alleleFreq'][locus][allele] 
        else:
            for sp_idx in range(num_wild_pops):
                for locus in biallelic_loci:
                    for allele in [0,1]:
                        afs_dict[rep][year_counter][sp_idx][locus][allele] = \
                        pop.vars()['subPop'][sp_idx]['alleleFreq'][locus][allele]
                        
            
    def get_store_mfit(pop, sp_names, rvar_dict, rep, year_counter):
        '''Calculates the mean fitness of the subpopulation by using
        mean fitness = p^2(rel fit pp) + 2pq(rel fit pq) + q^2(rel fit qq)
        Stores the mean fitness by subpopulation in a response variable dictionary.''' 
        for sp_idx, sp_name in enumerate(sp_names):
            sum_mfit_this_subpop = 0
            for l_idx in range(num_aloci):
                stage = stages[l_idx//num_aloci_stage]
                within_stage_l_idx = l_idx % num_aloci_stage
                rfit_dict = {}
                rfit_dict['0'] = fit_dict[sp_name][stage][within_stage_l_idx]['0']
                rfit_dict['1'] = fit_dict[sp_name][stage][within_stage_l_idx]['1']
                rfit_dict['2'] = fit_dict[sp_name][stage][within_stage_l_idx]['2']            
                freq_pp = pop.dvars().subPop[sp_idx]['genoFreq'][l_idx][(0,0)]
                freq_pq = pop.dvars().subPop[sp_idx]['genoFreq'][l_idx][(1,0)] + \
                pop.dvars().subPop[sp_idx]['genoFreq'][l_idx][(0,1)] 
                freq_qq = pop.dvars().subPop[sp_idx]['genoFreq'][l_idx][(1,1)]
                mfit_this_locus = freq_pp*rfit_dict['0'] + freq_pq*rfit_dict['2'] \
                + freq_qq*rfit_dict['1'] # first sum to get mean
                sum_mfit_this_subpop += mfit_this_locus
            mfit_this_subpop = sum_mfit_this_subpop / num_aloci # then divide to get mean 
            if sp_name != 'farm':
                rvar_pop = wild_names.index(sp_name) 
            else:
                rvar_pop = sp_name
            rvar_dict[rep][year_counter][rvar_pop]['mfit'] =  mfit_this_subpop # store value
        
    def store_fst(pop, pair_rvar_dict, rep, year_counter, wild_pop_pairs, all_pop_pairs, farm_boo):
        '''Stores Fst from all loci for each population pair in a response variable dictionary.'''
        if farm_boo ==  False:
            pop_pairs = wild_pop_pairs
        else:
            pop_pairs = all_pop_pairs
        for pop_pair in pop_pairs:
            pair_rvar_dict[rep][year_counter][pop_pair]['Fst'] = pop.vars()[''.join(['F_st_', pop_pair])]  
            
    def store_fst_adaploci(pop, pair_rvar_dict, rep, year_counter, wild_pop_pairs, all_pop_pairs, farm_boo):
        '''Stores Fst at adaptive loci only for each population pair in a response variable dictionary.'''
        if farm_boo ==  False:
            pop_pairs = wild_pop_pairs
        else:
            pop_pairs = all_pop_pairs
        for pop_pair in pop_pairs:
            pair_rvar_dict[rep][year_counter][pop_pair]['Fst_aL'] = pop.vars()[''.join(['F_st_aL_', pop_pair])]   
            
    def store_fst_neuloci(pop, pair_rvar_dict, rep, year_counter, wild_pop_pairs, all_pop_pairs, farm_boo):
        '''Stores Fst at neutral loci only for each population pair in a response variable dictionary.'''
        if farm_boo ==  False:
            pop_pairs = wild_pop_pairs
        else:
            pop_pairs = all_pop_pairs
        for pop_pair in pop_pairs:
            pair_rvar_dict[rep][year_counter][pop_pair]['Fst_nL'] = pop.vars()[''.join(['F_st_nL_', pop_pair])]  
    
    def store_popsize(pop, rvar_dict, rep, year_counter, farm_boo, wild_subPopIndeces, farm_idx):
        '''Stores population size for each subPop in response variable dictionary.'''
        subPop_sizes = pop.vars()['subPopSize']
        if farm_boo == False: # if no farm, wild subpops in order
            for count, subPop_size in enumerate(subPop_sizes):
                rvar_dict[rep][year_counter][count]['popsize'] = subPop_size
        elif farm_boo == True:
            farm_subPopSize = subPop_sizes[farm_idx]
            wild_subPopSizes = subPop_sizes[0:farm_idx] + subPop_sizes[farm_idx+1:]
            for count, subPop_size in enumerate(wild_subPopSizes):
                rvar_dict[rep][year_counter][count]['popsize'] = subPop_size
            rvar_dict[rep][year_counter]['farm']['popsize'] = farm_subPopSize
            
    def get_and_store_het_nloci(pop, rvar_dict, rep, year_counter, w_farm_wld_sp_idcs, farm_boo, farm_idx):
        '''Get and store heterozygosity (proportion of heterozygous genotypes),
        at neutral loci.'''
        if farm_boo == False: # if not farm year
            for sp_idx in range(num_wild_pops):    
                rvar_dict[rep][year_counter][sp_idx]['het'] = \
                np.mean(list(pop.dvars().subPop[sp_idx]['heteroFreq'].values()))
        else: # if farm year
            for herein_idx, sp_idx in enumerate(w_farm_wld_sp_idcs):
                rvar_dict[rep][year_counter][herein_idx]['het'] = \
                np.mean(list(pop.dvars().subPop[sp_idx]['heteroFreq'].values()))
            rvar_dict[rep][year_counter]['farm']['het'] = \
                np.mean(list(pop.dvars().subPop[farm_idx]['heteroFreq'].values()))
                    
    def get_and_store_ar(pop, rvar_dict, rep, year_counter, farm_boo, w_farm_wld_sp_idcs, farm_idx): #, loci_for_ibd):
        '''Get and store allelic richness per subpop, number of unique alleles'''
        if farm_boo == False: # if not farm year
            for sp_idx in range(num_wild_pops):
                ar = 0
                for locus_idx in loci_for_ibd:
                    this_locus_allele_idcs = [i*num_total_loci + locus_idx for i in range(2*pop.subPopSize(sp_idx))]
                    alleles_this_locus = [pop.genotype(sp_idx)[i] for i in this_locus_allele_idcs]
                    ar += len(list(set(alleles_this_locus)))
                rvar_dict[rep][year_counter][sp_idx]['ar'] = ar    
        else: # if farm year
            for herein_idx, sp_idx in enumerate(w_farm_wld_sp_idcs):  
                ar = 0
                for locus_idx in loci_for_ibd:
                    this_locus_allele_idcs = [i*num_total_loci + locus_idx for i in range(2*pop.subPopSize(sp_idx))]
                    alleles_this_locus = [pop.genotype(sp_idx)[i] for i in this_locus_allele_idcs]
                    ar += len(list(set(alleles_this_locus)))
                rvar_dict[rep][year_counter][herein_idx]['ar'] = ar 
                
            ar = 0
            for locus_idx in loci_for_ibd:
                this_locus_allele_idcs = [i*num_total_loci + locus_idx for i in range(2*pop.subPopSize(farm_idx))]
                alleles_this_locus = [pop.genotype(farm_idx)[i] for i in this_locus_allele_idcs]
                ar += len(list(set(alleles_this_locus)))
            rvar_dict[rep][year_counter]['farm']['ar'] = ar
                       
    def store_rvars(pop, farm_boo, farm_idx, rep, year_counter, afs_dict, loci_to_track, 
                    rvar_dict, pair_rvar_dict, wild_subPopIndeces, num_aloci, sp_names):
        '''Store response variables'''
        get_store_mfit(pop=pop, sp_names=sp_names, rvar_dict=rvar_dict, rep=rep, year_counter=year_counter)
        store_fst(pop=pop, pair_rvar_dict=pair_rvar_dict, rep=rep, year_counter=year_counter, 
                  wild_pop_pairs=wild_pop_pairs, all_pop_pairs = all_pop_pairs, farm_boo=farm_boo)
        store_fst_adaploci(pop=pop, pair_rvar_dict=pair_rvar_dict, rep=rep, 
                           year_counter=year_counter, wild_pop_pairs=wild_pop_pairs, 
                           all_pop_pairs = all_pop_pairs, farm_boo=farm_boo)
        store_fst_neuloci(pop=pop, pair_rvar_dict=pair_rvar_dict, rep=rep, 
                          year_counter=year_counter, wild_pop_pairs=wild_pop_pairs, 
                          all_pop_pairs = all_pop_pairs, farm_boo=farm_boo)
        store_popsize(pop=pop, rvar_dict=rvar_dict, rep=rep, year_counter=year_counter, farm_boo=farm_boo,
                      wild_subPopIndeces=wild_subPopIndeces, farm_idx=farm_idx)
        store_afs(pop=pop, rep=rep, year_counter=year_counter, farm_boo=farm_boo, farm_idx=farm_idx,
                  afs_dict=afs_dict,loci_to_track=loci_to_track)
        get_and_store_het_nloci(pop=pop, rvar_dict=rvar_dict, rep=rep, year_counter=year_counter, 
                                w_farm_wld_sp_idcs=w_farm_wld_sp_idcs, farm_boo=farm_boo, farm_idx=farm_idx)
        get_and_store_ar(pop=pop, rvar_dict=rvar_dict, rep=rep, year_counter=year_counter,
                         farm_boo=farm_boo, w_farm_wld_sp_idcs=w_farm_wld_sp_idcs, farm_idx=farm_idx)
    
    ###########################################################################
    #### --- Reporting functions
    
    def write_to_log(logfilename, pop, farm_idx, farm_boo, farm_phase):
        '''Write line to log per year, with sp sizes and date, like
           wild1 wild2 wild3 farm (if applicable)'''
           
        sp_sizes = list(pop.subPopSizes())
        sp_sizes = [str(sp_size) for sp_size in sp_sizes]
        if farm_boo == False: # if not a farm year
            log_line =  datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
            log_line = '\t'.join([str(rep), farm_phase, str(year_counter), '\t'.join(sp_sizes), 'NA',
                               datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H:%M:%S')])
        else: # if farm year
            sp_sizes_ordered_w_farm = sp_sizes[0:farm_idx] + sp_sizes[farm_idx+1:] + [sp_sizes[farm_idx]]
            log_line = '\t'.join([str(rep), farm_phase, str(year_counter),'\t'.join(sp_sizes_ordered_w_farm),
                                datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H:%M:%S')]) 
        with open(logfilename, 'a') as logfile:
            logfile.write(log_line)
            logfile.write('\n')
                         
    def track_all_inds(sp_names, rep, year_counter, month):
        '''Report individual information, for producing life history table'''
        for sp_name in sp_names:
            for ind in pop.individuals(subPop=[sp_name]):            
                life_hist_report.write('\t'.join([str(rep),
                                                 str(year_counter),
                                                 str(month),
                                                 str(sp_name),
                                                 str(ind.ind_id),
                                                 str(ind.sex()),
                                                 str(ind.cohortYear),
                                                 str(round(ind.age, 3)),
                                                 str(ind.mother_id),
                                                 str(ind.father_id),
                                                 str(ind.meanParentAge)]))
                life_hist_report.write('\n')
                
    def track_adults(pop, sp_names, ageX, rep, year_counter):
        '''Report information for adults, inds having reached ageX'''
        for sp_name in sp_names:
            all_inds = {} # key = ind_id of all inds, value = [mom id, dad id]
            for ind in pop.individuals(subPop=[sp_name]):
                all_inds[ind.ind_id] = [ind.mother_id, ind.father_id]
                if ind.age >= ageX:
                    track_for_Ne.write('\t'.join([str(rep),
                                                  str(year_counter),
                                                  str(sp_name),
                                                  str(ind.ind_id),
                                                  str(ind.cohortYear),
                                                  str(round(ind.age, 3)),
                                                  str(all_inds[ind.ind_id][0]),
                                                  str(all_inds[ind.ind_id][1]),
                                                  str(ind.meanParentAge)]))
                    track_for_Ne.write('\n')
    
    
    # ################################################################################
    # #### ----------------------------- Parameters ----------------------------- ####
    # ################################################################################
    
    startTime = time.time()

    
    # manage argumetns from command line
    # parser = argparse.ArgumentParser(description="An individual-based model of shellfish production, escape, and genetic impacts to wild populations")
    # parser.add_argument("-b", "--batch", help="Name of batch to be applied to ouptut files, e.g., 'low_mig'",
    #                     type=str, required=True)
    # parser.add_argument("-r", "--reps", help="Number of replicates for simulation, e.g., 10", 
    #                     type=int, required=True)
    # parser.add_argument("-x", "--coreid", help="Identifier for combining results across cores, e.g., 1", 
    #                     type=int, required=True)
    # parser.add_argument("-f", "--freq", help="Path to equilibrium allele frequencies file; if not provided, will initialize biallelic neutral loci to 0.5 and 0.5.", 
    #                     type=str, required=False)
    # args = parser.parse_args()
    
    # years
    # pre_farm_years = 5 # num of years first (of three) section of the model
    # farm_years = 1 # num of years second section of the model
    # post_farm_years = 1 # num of years third section of the model
    
    # wild popdy
    # num_wild_pops = 3 # number wild pops (model structure currently only handles 3)
    # wild_N_init = 300 # initial wild pop size
    # rec_a = 750 # annual recruitment event size; density independent
    # sd_recruit = 0.01 # standard deviation in wild recruitment
    # numWildOffspring_par = 10 # wild family size
    # wild_mig_rate_L = 0.035/12 # larvae monthly mig rate, but currently only happens in one month
    # wild_mig_rate_J = 0 # juvenile
    # wild_mig_rate_A = 0 # adult
    # 
    # # farm popdy
    # seed_batch_size = 250 # farm recruitment, i.e., number of seed planted on farm at once
    # numFarmOffspring_par = 20 # farm family size
    # sd_seed = 0.01 # standard deviation in farm recruitment / scale of seed production
    # farm_nonB_recruit_factor = 0.25 # factor governing recruitment due to F2 adults on farm
    # gamEsc_rec = 20 # monthly wild recruitment due to gamete escape
    # source_idx = 0 # index of wild pop to collect broodstock from
    # local_wild_idx = 0 # index of subpop that farm will escape to 
    # L_escape_rate = 0.07 # larvae monthly escape rate, but currently only happens in one month
    # J_escape_rate = 0/12 # juvenile
    # A_escape_rate = 0/12 # adult
    # numGameteEscapeOffspring_par = 10 # hybrid / gamete escape family size
    # 
    # seasonal processes
    prob_repro_by_month = {'Jan':0, 'Feb':0, 'Mar':0, # for wild reproduction
                           'Apr':1, 'May': 0, 'Jun': 0, 
                           'Jul': 0, 'Aug':0, 'Sep':0, 
                           'Oct':0, 'Nov':0, 'Dec':0}
    prob_prodseed_by_month = {'Jan':1, 'Feb':0, 'Mar':0, # for production of seed in farm
                              'Apr':0, 'May': 0, 'Jun': 0, 
                              'Jul': 0, 'Aug':0, 'Sep':0, 
                              'Oct':0, 'Nov':0, 'Dec':0}
    prob_L_G_escape_by_month = {'Jan':0, 'Feb':0, 'Mar':0, # for larval and gamete
                            'Apr':1, 'May': 0, 'Jun': 0, 
                            'Jul': 0, 'Aug':0, 'Sep':0, 
                            'Oct':0, 'Nov':0, 'Dec':0}
    prob_J_A_escape_by_month = {'Jan':0, 'Feb':0, 'Mar':0, # for juvenile and adult escape
                            'Apr':0, 'May': 0, 'Jun': 0, 
                            'Jul': 0, 'Aug':0, 'Sep':0, 
                            'Oct':0, 'Nov':0, 'Dec':0}
    
    # loci & selection
    num_aloci_stage = 5 # number of adaptive loci per stage
    num_nloci_ibd = 10 # number of neutral loci for identity by descent
    num_nloci_fst = 10 # number of neutral loci for fst
    good_start_AF = 0.9 # starting allele freq for beneficial alleles
    s = 0.05 # selection coefficient, to differentiate among wild pops
    z = 0.5 # selection coefficient, to differentiate between farm and wild pops
    
    # model reports
    track_AFs = True # track biallelic allele freqs
    make_report_all_inds = True # report ind info per month, for all inds
    make_report_adults = False # report ind info per year, only for adults
    yr_interval_log = 1 # log to write each X years
    
    # constant
    wild_pop_pairs = ['wild1_wild2', 'wild2_wild3', 'wild1_wild3'] # name of pop pairs
    all_pop_pairs = ['wild1_wild2', 'wild2_wild3', 'wild1_wild3', 'farm_wild1', 'farm_wild2', 'farm_wild3'] # name of pop pairs
    
    # species life history and production: olympia oysters
    settle_age = 1/12 # age at settlement
    repro_age = 1 # age at maturity in years
    ageX = 1 # age to start enumerating adults for adult report
    max_age = 10 # max age in years
    harvest_age = 2 # minimum harvest age in years
    max_harvest_age = 3 # harvest all animals once they reach this age
    harvest_rate = 1/12 # proportion of animals to harvest on farm per month
    kill_used_bstock = False # True: kill when collecting new, False: return bstock to wild
    adult_annual_M = 0.3 # adult mortality rate per year
    lar_annual_M = 0.99 #adult_annual_M*10 # larval mortality rate per year
    juv_annual_M = 0.3 # juvenile mortality rate per year
    var_in_mort = 0.25 # proportion above and below, e.g., fit factor will range from 1-x to 1+x
    num_broodstock = 10 # number of broodstock to be collected at once
    eq_sex = True # boolean; T =  eq females & males in bstock; F = random sex ratio, at least one male & female
    stages = ['L', 'J', 'A'] # names for life history stages
    deg_add = {'L': 0, 'J': 1, 'A': 1} # degree additive effects for alleles per stage
    farm_reduced_mort = 0.1 # mortality on farm will be less than in wild by this factor
    xyrs_new_bstock = 1 # get fresh broodstock every x years
    
    # write parameter report (testing for now - idea is users can load and run from a parameter report instead of manually)
    
    # F = open(batch + '_shellfish_params.pkl', 'wb') 
    # 
    # params = {
    # 'prob_repro_by_month':prob_repro_by_month, 'prob_prodseed_by_month':prob_prodseed_by_month, 'prob_L_G_escape_by_month':prob_L_G_escape_by_month,
    # 'prob_L_G_escape_by_month':prob_L_G_escape_by_month, 'prob_J_A_escape_by_month':prob_J_A_escape_by_month, 'num_aloci_stage':num_aloci_stage,
    # 'num_nloci_ibd':num_nloci_ibd,
    # 'num_nloci_fst':num_nloci_fst,
    # 'good_start_AF':good_start_AF,
    # 's': s,
    # 'z': z
    # }
    # 
    # pickle.dump(params, F) 
    # 
    # F.close() 

    # inferred parameters (based on input parameters above)
    total_years = pre_farm_years + farm_years + post_farm_years
    total_years_farm_boo = [False]*pre_farm_years + [True]*farm_years + [False]*post_farm_years
    num_aloci = num_aloci_stage*len(stages) # number adaptive loci across stages
    num_total_loci = num_aloci + num_nloci_ibd + num_nloci_fst # total number of loci
    bad_start_AF = 1-good_start_AF
    wild_names = [''.join(['wild',str(i+1)]) for i in range(num_wild_pops)]
    age_vsp_idcs = list(range(len(stages)))
    source_name = wild_names[source_idx]
    farm_idx = source_idx+1
    wo_farm_sp_names = get_sp_names(source_idx=source_idx, farm_boo=False, wild_names=wild_names)
    wo_farm_sp_idcs = list(range(num_wild_pops))
    w_farm_sp_idcs = list(range(num_wild_pops+1))
    w_farm_sp_names = get_sp_names(source_idx=source_idx, farm_boo=True, wild_names=wild_names)
    w_farm_wld_sp_idcs = get_w_farm_wld_sp_idcs(source_idx=source_idx, num_wild_pops=3) # wild subpop indeces with farm
    w_farm_sps_to_clone = get_sps_to_clone(sp_idcs=w_farm_sp_idcs, age_vsp_idcs=age_vsp_idcs)
    loci_for_fst = range(num_aloci+num_nloci_ibd, num_total_loci)
    loci_for_ibd = range(num_aloci, num_aloci+num_nloci_ibd)
    biallelic_loci = [j for j in range(num_aloci)] + [k for k in range(num_aloci+num_nloci_ibd, num_total_loci)]
    w_farm_w1_w2_idcs = [w_farm_wld_sp_idcs[0], w_farm_wld_sp_idcs[1]]  
    w_farm_w2_w3_idcs = [w_farm_wld_sp_idcs[1], w_farm_wld_sp_idcs[2]]  
    w_farm_w1_w3_idcs = [w_farm_wld_sp_idcs[0], w_farm_wld_sp_idcs[2]]  
    batch = '_'.join([str(batch),str(coreid)])
    
    ################################################################################
    #### --- Initialize response var storage objects, per experimental test --- ####
    ################################################################################
    
    # print('batch', batch)
    # 
    progress_report.write("Batch " + str(batch) + "\n")
    
    progress_report.flush()
    
    logfilename = os.path.join(batch_dir, '_'.join(['batch',batch,'log.txt'])) # name log file with batch id
    if os.path.isfile(logfilename): # if exists, remove old to prevent appending
        os.remove(logfilename)
    
    # write header to log file
    with open(logfilename, 'a') as logfile:            
        logfile.write('rep\tfarm_phase\tyear\twild1_size\twild2_size\twild3_size\tfarm_size\ttime_stamp\n')
    
    # open report files and write header lines
    bstock_report = open(os.path.join(batch_dir,''.join(['batch_', batch, '_bstock_report.txt'])), 'w')
    bstock_report.write('Rep\tYear\tIndID_survivedAgeX\n')
    if make_report_adults:
        track_for_Ne = open(os.path.join(batch_dir,''.join(['batch_', batch, '_trackingForNE.txt'])), 'w')
        track_for_Ne.write('Rep\tYear\tSubpop\tIndID_survivedAgeX\tCohortYear\tAge\tMother_id\tFather_id\tmeanParentAge\n')
    if make_report_all_inds:
        life_hist_report = open(os.path.join(batch_dir,''.join(['batch_', batch, '_life_hist_report.txt'])), 'w')
        life_hist_report.write('Rep\tYear\tMonth\tSubpop\tIndID\tSex\tCohortYear\tAge\tMother_id\tFather_id\tmeanParentAge\n')
    
    # init dictionary to store rvars that occur per subpop, per year
    pop_rvar_dict = {}
    for rep in range(reps):
        pop_rvar_dict[rep] = {}
        for year in range(total_years):
            pop_rvar_dict[rep][year] = {}
            for wild_sp in range(num_wild_pops):
                pop_rvar_dict[rep][year][wild_sp] = {}
            if total_years_farm_boo[year] == True: # it's a farm year
                pop_rvar_dict[rep][year]['farm'] = {}
                
    # init dictionary to store allele frequencies
    afs_dict = {}
    for rep in range(reps):
        afs_dict[rep] = {}
        for year in range(total_years):
            afs_dict[rep][year] = {}
            for wild_sp in range(num_wild_pops):
                afs_dict[rep][year][wild_sp] = {}
                for l_idx in range(num_total_loci):
                    afs_dict[rep][year][wild_sp][l_idx] = {}
            if total_years_farm_boo[year] == True: # it's a farm year
                afs_dict[rep][year]['farm'] = {}
                for l_idx in range(num_total_loci):
                    afs_dict[rep][year]['farm'][l_idx] = {}   
                    
    # init dictionary to store rvars that occur per pair of subpops, per year
    pop_pair_rvar_dict = {}
    for rep in range(reps):
        pop_pair_rvar_dict[rep] = {}
        for year in range(total_years):
            pop_pair_rvar_dict[rep][year] = {}
            if total_years_farm_boo[year] == True: # it's a farm year
                for pop_pair in all_pop_pairs:
                    pop_pair_rvar_dict[rep][year][pop_pair] = {}
            else: # not a farm year
                for pop_pair in wild_pop_pairs:
                    pop_pair_rvar_dict[rep][year][pop_pair] = {}
    
    # init dictionary to store rvars that occur per subpop, across time
    temp_rvar_dict = {}
    for rep in range(reps):
        temp_rvar_dict[rep] = {}
        for sp in range(num_wild_pops):
            temp_rvar_dict[rep][sp] = {}
    
    # init dictionary to store harvest quantities
    harvest_dict = {}
    for rep in range(reps):
        harvest_dict[rep] = {}
        for year in range(total_years):
            harvest_dict[rep][year] = 0 # per rep per year, number inds harvested
            
    # repeat the model with same parameters to aggregate results
    for rep in range(reps):
        # print(' Rep', rep, datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))
        
        progress_report.write('Rep ' + str(rep) + " " + str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')))

        progress_report.write("\n")

        progress_report.flush()

        # assign climate deviation series per model rep, which will add error to recruitment
        clim_dev = np.zeros(total_years)
        for year in range(1, total_years):
            clim_dev[year] = 0.9*clim_dev[year-1] + np.random.normal(0,1,1)[0]
            
        # assign fitness dictionary to create fitness landscape per model rep
        fit_dict = make_fit_dict(stages=stages,
                                 num_aloci_stage=num_aloci_stage,
                                 sp_names=w_farm_sp_names, 
                                 s=s, z=z,
                                 deg_add=deg_add)
         
        year_counter = 0
        
        ################################################################################
        #### -------------------- Initialize wild subpopulations ------------------ ####
        ################################################################################
    
        # initialize a population with three wild subpopulations
        pop = sim.Population(size=[wild_N_init]*num_wild_pops, 
                            loci=num_total_loci,
                            subPopNames=wo_farm_sp_names,
                            infoFields=['fitness', 'migrate_to', 'ind_id','age',
                                        'meanParentAge', 'broodstock', 'kill',
                                        'nonB_adult', 'return_to', 'cohortYear', 
                                        'father_id', 'mother_id', 'F'])
    
        # initialize age in years, size in mm, sex, afs, individual IDs
        for ind in pop.individuals():
            ind.age = float(np.random.choice(a=list(np.arange(1,max_age,1/12)), size=1, replace=True, 
                             p=np.exp(-0.1*np.array(list(np.arange(1,max_age,1/12)))) \
                                             /sum(np.exp(-0.1*np.array(list(np.arange(1,max_age,1/12)))))))   
        sim.initSex(pop) 
        
        # loci sets in order: adaptive loci, neutral loci for IBD, then neutral loci for FST
        
        # initialize allele frequencies for adaptive loci differently than nuetral loci
        init_geno_aloci(pop=pop, sp_names=wo_farm_sp_names, num_aloci=num_aloci,
                        fit_dict=fit_dict, good_start_AF=good_start_AF, bad_start_AF=bad_start_AF)
        
        # assign only unique alleles so we can measure identity by descent for first set of neutral loci
        sim.initGenotype(pop, loci=loci_for_ibd, 
                         genotype=[i for i in range(2*pop.popSize()) for x in range(len(loci_for_ibd))])
        if freq is not None: # if provided input file for starting AFs at equilibrium
            initEqAFs(loci_for_fst= loci_for_fst, eq_afs_path=freq)
        else: # else, start at 50 / 50
            sim.initGenotype(pop, loci=loci_for_fst, 
                             freq=[0.5,0.5])
        sim.tagID(pop) # give each ind an indivdiual tag ID
    
        # set virtual splitters: for age stages, for whether broodstock, and for whether to kill in mortality or harvest
        pop.setVirtualSplitter(sim.CombinedSplitter([
            sim.ProductSplitter([sim.InfoSplitter(field='age', cutoff=[settle_age, repro_age], names=['L','J','A']), 
                                 sim.InfoSplitter(field='age', cutoff=[harvest_age], names=['leave','to_harvest']),
                                 sim.InfoSplitter(field='broodstock', values=[0,1], names=['not_bstock','bstock']),
                                 sim.InfoSplitter(field='kill', values=[0,1], names=['to_survive', 'to_kill'])])]))
        vsp_n2i = make_vsp_name_to_index_dict(pop) # get vsp name to index dict
        dict_for_vspMap = get_vsp_map(vsp_n2i)
        pop.setVirtualSplitter(sim.CombinedSplitter([
            sim.ProductSplitter([sim.InfoSplitter(field='age', cutoff=[settle_age, repro_age], names=['L','J','A']), 
                                 sim.InfoSplitter(field='age', cutoff=[harvest_age], names=['leave','to_harvest']),
                                 sim.InfoSplitter(field='broodstock', values=[0,1], names=['not_bstock','bstock']),
                                 sim.InfoSplitter(field='kill', values=[0,1], names=['to_survive', 'to_kill'])])],
                vspMap = list(dict_for_vspMap.values()),
                names = list(dict_for_vspMap.keys())))
        vsp_n2i = make_vsp_name_to_index_dict(pop) 
    
        # desginate subpops to clone and mate
        wo_farm_sps_to_clone = get_sps_to_clone(sp_idcs=w_farm_sp_idcs, age_vsp_idcs=age_vsp_idcs)
        wo_farm_sps_to_mate = get_wo_farm_sps_to_mate(num_wild_pops=num_wild_pops)
    
        ################################################################################
        #### --------------------------- Pre-farm period -------------------------- ####
        ################################################################################
    
        # print('  Pre-farm period')
        
        progress_report.write("Pre-farm period\n")

        progress_report.flush()
        
        for year in np.arange(pre_farm_years):
            
            for month in prob_repro_by_month:
                
                # track individual information for life history table
                if make_report_all_inds:
                    track_all_inds(sp_names=wo_farm_sp_names, rep=rep, 
                                   year_counter=year_counter, month=month)
                    
                # migrate, grow, age, and die
                wild_mig(pop=pop, farm_boo=False, farm_idx=farm_idx, wild_names=wild_names, wild_mig_rate_L=wild_mig_rate_L,
                         wild_mig_rate_J=wild_mig_rate_J, wild_mig_rate_A=wild_mig_rate_A, 
                         w_farm_wld_sp_idcs=w_farm_wld_sp_idcs) # migrate         
                pop.setIndInfo(values=list(np.array(pop.indInfo(field='age'))+1/12), field='age') # age
                assign_fitness_all(pop) # for selective mortality
                mortality(pop=pop, lar_annual_M=lar_annual_M, juv_annual_M=juv_annual_M, 
                          adult_annual_M=adult_annual_M, farm_reduced_mort=farm_reduced_mort) # die
    
                # reproduce based on probability to reproduce by month
                if random.random() < prob_repro_by_month[month]:
                    pop.evolve(preOps=[sim.PySelector(func=get_fitness)],
                            matingScheme=sim.HeteroMating([sim.CloneMating(subPops=wo_farm_sps_to_clone, weight=-1),
                                                           sim.RandomMating(subPops=wo_farm_sps_to_mate,
                                                                            numOffspring=numWildOffspring_par,
                                                                             ops=[sim.MendelianGenoTransmitter(),
                                                                                     sim.IdTagger(),
                                                                                  sim.PedigreeTagger(),
                                                                                 sim.SummaryTagger(mode=sim.MEAN, infoFields=['age', 'meanParentAge'])])],
                                                         subPopSize=wo_farm_recruitment(pop)),
                            postOps=[sim.Stat(alleleFreq=biallelic_loci, subPops=wo_farm_sp_idcs,
                                              vars=['alleleFreq','alleleFreq_sp']),
                                     sim.Stat(genoFreq=range(num_aloci), subPops=wo_farm_sp_idcs,
                                          vars=['genoFreq','genoFreq_sp']),
                                 sim.Stat(heteroFreq=loci_for_ibd, vars='heteroFreq_sp'),
                                 sim.Stat(structure=biallelic_loci, subPops=[0,1], 
                                          suffix='_wild1_wild2', vars=['F_st']),
                                 sim.Stat(structure=biallelic_loci, subPops=[1,2], 
                                          suffix='_wild2_wild3', vars=['F_st']),
                                 sim.Stat(structure=biallelic_loci, subPops=[0,2], 
                                          suffix='_wild1_wild3', vars=['F_st']),
                                 sim.Stat(structure=range(num_aloci), subPops=[0,1], 
                                          suffix='_aL_wild1_wild2', vars=['F_st']),
                                 sim.Stat(structure=range(num_aloci), subPops=[1,2], 
                                          suffix='_aL_wild2_wild3', vars=['F_st']),
                                 sim.Stat(structure=range(num_aloci), subPops=[0,2], 
                                          suffix='_aL_wild1_wild3', vars=['F_st']),
                                 sim.Stat(structure=loci_for_fst, subPops=[0,1], 
                                          suffix='_nL_wild1_wild2', vars=['F_st']),
                                 sim.Stat(structure=loci_for_fst, subPops=[1,2], 
                                          suffix='_nL_wild2_wild3', vars=['F_st']),
                                 sim.Stat(structure=loci_for_fst, subPops=[0,2], 
                                          suffix='_nL_wild1_wild3', vars=['F_st']),
                                 sim.Stat(popSize=True, subPops=wo_farm_sp_idcs)],
                            gen=1)
                    
                    for ind in pop.individuals():
                        if ind.age == 0:
                            ind.cohortYear = year_counter
                    
                if month == 'Dec': # pre farm
                    
                    # store response variables
                    store_rvars(pop=pop, farm_boo=False, farm_idx=farm_idx, rep=rep, year_counter=year_counter, 
                                afs_dict=afs_dict, loci_to_track=biallelic_loci,
                                rvar_dict=pop_rvar_dict, 
                                pair_rvar_dict=pop_pair_rvar_dict, wild_subPopIndeces=wo_farm_sp_idcs,
                                num_aloci=num_aloci, sp_names=wo_farm_sp_names)
                    
                    if make_report_adults:
                        track_adults(pop=pop, sp_names=wo_farm_sp_names, ageX=ageX, 
                                     rep=rep, year_counter=year_counter)
    
            if year_counter % yr_interval_log == 0:
                print('   Yr',year_counter, 
                      datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))
                write_to_log(logfilename=logfilename, pop=pop, farm_idx=farm_idx, 
                             farm_boo=False, farm_phase='pre-farm')
            year_counter += 1
            
            if year == pre_farm_years-1: # in last year, save population as object to get rvar temporal fst
                pop.save(os.path.join(batch_dir,''.join(['batch_', batch, '_t0.pop'])))
            
        ################################################################################
        #### ---------------------- During-farm period ---------------------------- ####
        ################################################################################
            
        # print('  During-farm period')
    
        progress_report.write("During-farm period\n")

        progress_report.flush()
        # define broodstock sp for breeding, define wild subpops for mating 
        created_farm = False
        broodstock_sp = get_sps_to_mate_w_farm(wild_or_farm='farm', farm_idx=farm_idx, wild_idcs= w_farm_wld_sp_idcs)
        w_farm_sps_to_mate_wild = get_sps_to_mate_w_farm(wild_or_farm='wild_and_nonBfarm', 
                                                         farm_idx=farm_idx, wild_idcs=w_farm_wld_sp_idcs)        
        
        for year in np.arange(farm_years):
            
            # track harvest per year
            harvest_this_year = 0
            
            for month in prob_repro_by_month:
                
                if created_farm == False: # arbitrarily create farm population at the very beginning of farm phase
                    mark_broodstock(pop, source_idx=source_idx, source_name=source_name, 
                                    mat_vsp_name='A', num_broodstock=num_broodstock, eq_sex=False)                   
                    sim.splitSubPops(pop, subPops=source_idx, infoFields='broodstock', names=[source_name, 'farm']) 
                    pop.evolve(preOps=[sim.PySelector(func=get_fitness)],
                            matingScheme=sim.HeteroMating([sim.CloneMating(subPops=w_farm_sps_to_clone, weight=-1),
                                                           sim.RandomMating(subPops=broodstock_sp,
                                                                            numOffspring=numFarmOffspring_par,
                                                                             ops=[sim.MendelianGenoTransmitter(),
                                                                                     sim.IdTagger(),
                                                                                    sim.PedigreeTagger(),
                                                                                 sim.SummaryTagger(mode=sim.MEAN, infoFields=['age', 'meanParentAge'])])],
                                                         subPopSize=seed_production(pop=pop, farm_idx=farm_idx, 
                                                                                    seed_batch_size=seed_batch_size, 
                                                                                    sd_seed=sd_seed)),
                            postOps=[sim.Stat(alleleFreq=biallelic_loci, subPops=w_farm_sp_idcs, 
                                          vars=['alleleFreq','alleleFreq_sp']),
                                     sim.Stat(genoFreq=range(num_aloci), subPops=w_farm_sp_idcs,
                                          vars=['genoFreq','genoFreq_sp']),
                                     sim.Stat(heteroFreq=loci_for_ibd, vars='heteroFreq_sp'),
                                 sim.Stat(structure=biallelic_loci, subPops=w_farm_w1_w2_idcs, 
                                          suffix='_wild1_wild2', vars=['F_st']),
                                 sim.Stat(structure=biallelic_loci, subPops=w_farm_w2_w3_idcs, 
                                          suffix='_wild2_wild3', vars=['F_st']),
                                 sim.Stat(structure=biallelic_loci, subPops=w_farm_w1_w3_idcs, 
                                          suffix='_wild1_wild3', vars=['F_st']),
                                 sim.Stat(structure=biallelic_loci, subPops=[farm_idx, w_farm_wld_sp_idcs[0]], 
                                          suffix='_farm_wild1', vars=['F_st']),
                                 sim.Stat(structure=biallelic_loci, subPops=[farm_idx, w_farm_wld_sp_idcs[1]], 
                                          suffix='_farm_wild2', vars=['F_st']),
                                 sim.Stat(structure=biallelic_loci, subPops=[farm_idx, w_farm_wld_sp_idcs[2]], 
                                          suffix='_farm_wild3', vars=['F_st']),   
                                 sim.Stat(structure=range(num_aloci), subPops=w_farm_w1_w2_idcs,
                                          suffix='_aL_wild1_wild2', vars=['F_st']),
                                 sim.Stat(structure=range(num_aloci), subPops=w_farm_w2_w3_idcs, 
                                          suffix='_aL_wild2_wild3', vars=['F_st']),
                                 sim.Stat(structure=range(num_aloci), subPops=w_farm_w1_w3_idcs, 
                                          suffix='_aL_wild1_wild3', vars=['F_st']),
                                 sim.Stat(structure=range(num_aloci), subPops=[farm_idx, w_farm_wld_sp_idcs[0]],
                                          suffix='_aL_farm_wild1', vars=['F_st']),
                                 sim.Stat(structure=range(num_aloci), subPops=[farm_idx, w_farm_wld_sp_idcs[1]], 
                                          suffix='_aL_farm_wild2', vars=['F_st']),
                                 sim.Stat(structure=range(num_aloci), subPops=[farm_idx, w_farm_wld_sp_idcs[2]], 
                                          suffix='_aL_farm_wild3', vars=['F_st']),
                                 sim.Stat(structure=loci_for_fst, subPops=w_farm_w1_w2_idcs, 
                                          suffix='_nL_wild1_wild2', vars=['F_st']),
                                 sim.Stat(structure=loci_for_fst, subPops=w_farm_w2_w3_idcs, 
                                          suffix='_nL_wild2_wild3', vars=['F_st']),
                                 sim.Stat(structure=loci_for_fst, subPops=w_farm_w1_w3_idcs, 
                                          suffix='_nL_wild1_wild3', vars=['F_st']), 
                                 sim.Stat(structure=loci_for_fst, subPops=[farm_idx, w_farm_wld_sp_idcs[0]], 
                                          suffix='_nL_farm_wild1', vars=['F_st']),
                                 sim.Stat(structure=loci_for_fst, subPops=[farm_idx, w_farm_wld_sp_idcs[1]], 
                                          suffix='_nL_farm_wild2', vars=['F_st']),
                                 sim.Stat(structure=loci_for_fst, subPops=[farm_idx, w_farm_wld_sp_idcs[2]], 
                                          suffix='_nL_farm_wild3', vars=['F_st']),
                                 sim.Stat(popSize=True, subPops=w_farm_sp_idcs)],
                            gen=1)  
                    created_farm = True # farm subpop made
                    
                for ind in pop.individuals():
                        if ind.age == 0:
                            ind.cohortYear = year_counter
                            
                for ind in pop.individuals(farm_idx):
                    if ind.age == 0:
                        ind.F = 1
                       
                # migrate, grow, age, and die
                wild_mig(pop=pop, farm_boo=True, farm_idx=farm_idx, wild_names=wild_names, 
                         wild_mig_rate_L=wild_mig_rate_L, wild_mig_rate_J=wild_mig_rate_J, 
                         wild_mig_rate_A=wild_mig_rate_A, w_farm_wld_sp_idcs=w_farm_wld_sp_idcs) # migrate   
                pop.setIndInfo(values=list(np.array(pop.indInfo(field='age'))+1/12), field='age') # age
                assign_fitness_all(pop) # for selective mortality
                mortality(pop=pop, lar_annual_M=lar_annual_M, juv_annual_M=juv_annual_M, 
                          adult_annual_M=adult_annual_M, farm_reduced_mort=farm_reduced_mort) # die
                
                # harvest (only if number of harvestable animals is above a threshold)
                if pop.subPopSize([farm_idx,vsp_n2i['to_harvest']]) > 20:
                    harvest(pop=pop, farm_idx=farm_idx, harvest_rate=harvest_rate, 
                            harvest_age=harvest_age, max_harvest_age=max_harvest_age,
                            harvest_this_year=harvest_this_year)  
    
                # broodstock make seed based on probability to produce seed by month
                if random.random() < prob_prodseed_by_month[month]:
                           
                    # if lacking 1 mature male or 1 mature female or this is a year for fresh bstock, get fresh bstock
                    min1male = False
                    min1female = False
                    for ind in pop.individuals(subPop=[farm_idx,vsp_n2i['bstock']]): # is it even iterating?
                        if ind.age >= repro_age: # if mature
                            if ind.sex() == 2: # if female
                                min1female = True
                            if ind.sex() == 1: # if male
                                min1male = True
    
                    if min1male == False or min1female == False or year_counter % xyrs_new_bstock == 0:
    
                        for ind in pop.individuals():
                            if ind.broodstock == 1:
                                for sp_name in pop.subPopNames(): # get subpop index of this individual
                                    if ind in pop.individuals(sp_name):
                                        this_sp_name = sp_name # get subpop name
                        
                        # either kill old bstock or return to wild
                        if kill_used_bstock == True:
                            pop.removeSubPops([(farm_idx,vsp_n2i['bstock'])])
                        else:
                            for sp_idx in range(len(pop.subPopNames())):
                                for ind in pop.individuals([sp_idx]):
                                    if ind.broodstock == 1:
                                        ind.migrate_to = source_idx
                                    else:
                                        ind.migrate_to = sp_idx
                            sim.migrate(pop, mode=sim.BY_IND_INFO)   
                        
                        # make sure no bstock, including recently migrated
                        pop.setIndInfo(values=0, field='broodstock') 
                        
                        # collect new bstock, by mark and collect
                        mark_broodstock(pop=pop, source_idx=source_idx, source_name=source_name, 
                                        mat_vsp_name='A', num_broodstock=num_broodstock, eq_sex=False)
                        collect_bstock(pop=pop, sp_idcs=w_farm_sp_idcs, farm_idx=farm_idx)
    
                    pop.evolve(preOps=[sim.PySelector(func=get_fitness)],
                            matingScheme=sim.HeteroMating([sim.CloneMating(subPops=w_farm_sps_to_clone, weight=-1),
                                                           sim.RandomMating(subPops=broodstock_sp,
                                                                            numOffspring=numFarmOffspring_par,
                                                                             ops=[sim.MendelianGenoTransmitter(),
                                                                                     sim.IdTagger(),
                                                                                  sim.PedigreeTagger(),
                                                                                 sim.SummaryTagger(mode=sim.MEAN, 
                                                                                                   infoFields=['age', 
                                                                                                              'meanParentAge'])])],
                                                         subPopSize=seed_production(pop=pop, farm_idx=farm_idx, 
                                                                                    seed_batch_size=seed_batch_size, 
                                                                                    sd_seed=sd_seed)),
                            postOps=[sim.Stat(alleleFreq=biallelic_loci, subPops=w_farm_sp_idcs, 
                                              vars=['alleleFreq','alleleFreq_sp']),
                                 sim.Stat(genoFreq=range(num_aloci), subPops=w_farm_sp_idcs,
                                          vars=['genoFreq','genoFreq_sp']),
                                 sim.Stat(heteroFreq=loci_for_ibd, vars='heteroFreq_sp'),                     
                                 sim.Stat(structure=biallelic_loci, subPops=w_farm_w1_w2_idcs, 
                                          suffix='_wild1_wild2', vars=['F_st']),
                                 sim.Stat(structure=biallelic_loci, subPops=w_farm_w2_w3_idcs, 
                                          suffix='_wild2_wild3', vars=['F_st']),
                                 sim.Stat(structure=biallelic_loci, subPops=w_farm_w1_w3_idcs, 
                                          suffix='_wild1_wild3', vars=['F_st']),
                                 sim.Stat(structure=biallelic_loci, subPops=[farm_idx, w_farm_wld_sp_idcs[0]], 
                                          suffix='_farm_wild1', vars=['F_st']),
                                 sim.Stat(structure=biallelic_loci, subPops=[farm_idx, w_farm_wld_sp_idcs[1]], 
                                          suffix='_farm_wild2', vars=['F_st']),
                                 sim.Stat(structure=biallelic_loci, subPops=[farm_idx, w_farm_wld_sp_idcs[2]], 
                                          suffix='_farm_wild3', vars=['F_st']),   
                                 sim.Stat(structure=range(num_aloci), subPops=w_farm_w1_w2_idcs,
                                          suffix='_aL_wild1_wild2', vars=['F_st']),
                                 sim.Stat(structure=range(num_aloci), subPops=w_farm_w2_w3_idcs, 
                                          suffix='_aL_wild2_wild3', vars=['F_st']),
                                 sim.Stat(structure=range(num_aloci), subPops=w_farm_w1_w3_idcs, 
                                          suffix='_aL_wild1_wild3', vars=['F_st']),
                                 sim.Stat(structure=range(num_aloci), subPops=[farm_idx, w_farm_wld_sp_idcs[0]],
                                          suffix='_aL_farm_wild1', vars=['F_st']),
                                 sim.Stat(structure=range(num_aloci), subPops=[farm_idx, w_farm_wld_sp_idcs[1]], 
                                          suffix='_aL_farm_wild2', vars=['F_st']),
                                 sim.Stat(structure=range(num_aloci), subPops=[farm_idx, w_farm_wld_sp_idcs[2]], 
                                          suffix='_aL_farm_wild3', vars=['F_st']),
                                 sim.Stat(structure=loci_for_fst, subPops=w_farm_w1_w2_idcs, 
                                          suffix='_nL_wild1_wild2', vars=['F_st']),
                                 sim.Stat(structure=loci_for_fst, subPops=w_farm_w2_w3_idcs, 
                                          suffix='_nL_wild2_wild3', vars=['F_st']),
                                 sim.Stat(structure=loci_for_fst, subPops=w_farm_w1_w3_idcs, 
                                          suffix='_nL_wild1_wild3', vars=['F_st']), 
                                 sim.Stat(structure=loci_for_fst, subPops=[farm_idx, w_farm_wld_sp_idcs[0]], 
                                          suffix='_nL_farm_wild1', vars=['F_st']),
                                 sim.Stat(structure=loci_for_fst, subPops=[farm_idx, w_farm_wld_sp_idcs[1]], 
                                          suffix='_nL_farm_wild2', vars=['F_st']),
                                 sim.Stat(structure=loci_for_fst, subPops=[farm_idx, w_farm_wld_sp_idcs[2]], 
                                          suffix='_nL_farm_wild3', vars=['F_st']),
                                 sim.Stat(popSize=True, subPops=w_farm_sp_idcs)],
                            gen=1)     
      
                    for ind in pop.individuals():
                        if ind.age == 0:
                            ind.cohortYear = year_counter
                        
                    for ind in pop.individuals(farm_idx):
                        if ind.age == 0:
                            ind.F = 1
                                
                # wild subpops reproduce based on probability to reproduce by month
                if random.random() < prob_repro_by_month[month]:                
                    pop.evolve(preOps=[sim.PySelector(func=get_fitness)],
                            matingScheme=sim.HeteroMating([sim.CloneMating(subPops=w_farm_sps_to_clone, weight=-1),
                                                           sim.RandomMating(subPops=w_farm_sps_to_mate_wild,
                                                                            numOffspring=numWildOffspring_par,
                                                                             ops=[sim.MendelianGenoTransmitter(),
                                                                                     sim.IdTagger(),
                                                                                  sim.PedigreeTagger(),
                                                                                 sim.SummaryTagger(mode=sim.MEAN, 
                                                                                                   infoFields=['age',
                                                                                                              'meanParentAge'])])],
                                                         subPopSize=w_farm_recruitment(pop)), 
                            postOps=[sim.Stat(alleleFreq=biallelic_loci, subPops=w_farm_sp_idcs, 
                                              vars=['alleleFreq','alleleFreq_sp']),
                                 sim.Stat(genoFreq=range(num_aloci), subPops=w_farm_sp_idcs,
                                          vars=['genoFreq','genoFreq_sp']),
                                 sim.Stat(heteroFreq=loci_for_ibd, vars='heteroFreq_sp'),   
                                 sim.Stat(structure=biallelic_loci, subPops=w_farm_w1_w2_idcs, 
                                          suffix='_wild1_wild2', vars=['F_st']),
                                 sim.Stat(structure=biallelic_loci, subPops=w_farm_w2_w3_idcs, 
                                          suffix='_wild2_wild3', vars=['F_st']),
                                 sim.Stat(structure=biallelic_loci, subPops=w_farm_w1_w3_idcs, 
                                          suffix='_wild1_wild3', vars=['F_st']),
                                 sim.Stat(structure=biallelic_loci, subPops=[farm_idx, w_farm_wld_sp_idcs[0]], 
                                          suffix='_farm_wild1', vars=['F_st']),
                                 sim.Stat(structure=biallelic_loci, subPops=[farm_idx, w_farm_wld_sp_idcs[1]], 
                                          suffix='_farm_wild2', vars=['F_st']),
                                 sim.Stat(structure=biallelic_loci, subPops=[farm_idx, w_farm_wld_sp_idcs[2]], 
                                          suffix='_farm_wild3', vars=['F_st']),   
                                 sim.Stat(structure=range(num_aloci), subPops=w_farm_w1_w2_idcs,
                                          suffix='_aL_wild1_wild2', vars=['F_st']),
                                 sim.Stat(structure=range(num_aloci), subPops=w_farm_w2_w3_idcs, 
                                          suffix='_aL_wild2_wild3', vars=['F_st']),
                                 sim.Stat(structure=range(num_aloci), subPops=w_farm_w1_w3_idcs, 
                                          suffix='_aL_wild1_wild3', vars=['F_st']),
                                 sim.Stat(structure=range(num_aloci), subPops=[farm_idx, w_farm_wld_sp_idcs[0]],
                                          suffix='_aL_farm_wild1', vars=['F_st']),
                                 sim.Stat(structure=range(num_aloci), subPops=[farm_idx, w_farm_wld_sp_idcs[1]], 
                                          suffix='_aL_farm_wild2', vars=['F_st']),
                                 sim.Stat(structure=range(num_aloci), subPops=[farm_idx, w_farm_wld_sp_idcs[2]], 
                                          suffix='_aL_farm_wild3', vars=['F_st']),
                                 sim.Stat(structure=loci_for_fst, subPops=w_farm_w1_w2_idcs, 
                                          suffix='_nL_wild1_wild2', vars=['F_st']),
                                 sim.Stat(structure=loci_for_fst, subPops=w_farm_w2_w3_idcs, 
                                          suffix='_nL_wild2_wild3', vars=['F_st']),
                                 sim.Stat(structure=loci_for_fst, subPops=w_farm_w1_w3_idcs, 
                                          suffix='_nL_wild1_wild3', vars=['F_st']), 
                                 sim.Stat(structure=loci_for_fst, subPops=[farm_idx, w_farm_wld_sp_idcs[0]], 
                                          suffix='_nL_farm_wild1', vars=['F_st']),
                                 sim.Stat(structure=loci_for_fst, subPops=[farm_idx, w_farm_wld_sp_idcs[1]], 
                                          suffix='_nL_farm_wild2', vars=['F_st']),
                                 sim.Stat(structure=loci_for_fst, subPops=[farm_idx, w_farm_wld_sp_idcs[2]], 
                                          suffix='_nL_farm_wild3', vars=['F_st']),
                                 sim.Stat(popSize=True, subPops=w_farm_sp_idcs)],
                            gen=1)
                    
                    for ind in pop.individuals():
                        if ind.age == 0:
                            ind.cohortYear = year_counter
                    
                    for ind in pop.individuals(farm_idx):
                        if ind.age == 0 and ind.F == 0:
                            ind.F = 2
                            
                ### escape
                
                # juvenile and adult escape 
                if random.random() < prob_J_A_escape_by_month[month]:
                    escape_J_A(pop=pop, 
                               J_escape_rate=J_escape_rate, 
                               A_escape_rate=A_escape_rate, 
                               w_farm_wld_sp_idcs=w_farm_wld_sp_idcs, 
                               local_wild_idx=local_wild_idx)
    
                 # larvae and gamete escape
                if random.random() < prob_L_G_escape_by_month[month]:
                    escape_L(pop=pop, 
                               L_escape_rate=L_escape_rate, 
                               w_farm_wld_sp_idcs=w_farm_wld_sp_idcs, 
                               local_wild_idx=local_wild_idx)
                    
                    ### gamete escape
                    
                    # ensure enough and track parents through hybrid sp
                    farm_male = False
                    farm_female = False
                    nonB_farm_adults = False
                    for ind in pop.individuals(farm_idx):
                        if ind.age >= repro_age and ind.broodstock == 0:
                            ind.nonB_adult = 1
                            ind.return_to = farm_idx
                            if ind.sex() == 1:
                                farm_male = True
                            else:
                                farm_female = True
                        elif ind.age >= repro_age and ind.broodstock == 1: 
                            ind.nonB_adult = 0
                    if farm_male and farm_female:
                        nonB_farm_adults = True
    
                    wild1_male = False
                    wild1_female = False
                    nonB_wild1_adults = False
                    for ind in pop.individuals(0): # local wild
                        if ind.age >= repro_age:
                            ind.nonB_adult = 1
                            ind.return_to = 0
                            if ind.sex() == 1:
                                wild1_male = True
                            else:
                                wild1_female = True
                               
                    if wild1_male and wild1_female:
                        nonB_wild1_adults = True
    
                    # if non broodstock adults present in farm (to provide gamete escape)
                    if nonB_farm_adults and nonB_wild1_adults:
    
                        # make temp hybrid pop 
                        sim.splitSubPops(pop, subPops=farm_idx, infoFields='nonB_adult', names=['farm', 'hybrid']) 
                        hybrid_index = farm_idx + 1
    
                        # migrate adults from wild1 to temp hybrid pop
                        for ind in pop.individuals(0):
                            if ind.nonB_adult:
                                ind.migrate_to = hybrid_index
                            else:
                                ind.migrate_to = 0       
                        all_but_wild1 = [x for x in range(len(wild_names)+2) if x != 0] # get indeces of subpops other than wild1
                        for subPop in all_but_wild1: 
                            for ind in pop.individuals([subPop]):
                                ind.migrate_to = subPop # mark wild individuals to remain    
                        sim.migrate(pop, mode=sim.BY_IND_INFO) # migrate local wild1 animals to hybrid
                        
                        pop.evolve(preOps=[sim.PySelector(func=get_fitness_hybrid)],
                                       matingScheme=sim.HeteroMating(
                                           [sim.CloneMating(weight=-1), 
                                            sim.HomoMating(sim.PyParentsChooser(gameteEscapeChooser),
                                            sim.OffspringGenerator(ops=[sim.MendelianGenoTransmitter(),
                                            sim.IdTagger(), sim.PedigreeTagger(), 
                                            sim.SummaryTagger(mode=sim.MEAN, 
                                                              infoFields=['age', 'meanParentAge'])],
                                                                   numOffspring=numGameteEscapeOffspring_par), 
                                                           subPops=hybrid_index)], subPopSize=scaleGameteEscape(pop)), 
                                   gen=1)
    
                        for ind in pop.individuals():
                            if ind.age == 0:
                                ind.cohortYear = year_counter
                                
                        # return adults back to respective subpops
                        for sp_idx in range(len(wild_names)+2):
                            if sp_idx == hybrid_index:
                                for ind in pop.individuals(sp_idx):
                                    ind.migrate_to = ind.return_to
                                    # new larvae should have this set to 0, which is wild1 where escape lands
                                    # adults should have a value either 0 or farm index, so they'll go back to theirs
                            else:
                                for ind in pop.individuals(sp_idx):
                                    ind.migrate_to = sp_idx
                        sim.migrate(pop, mode=sim.BY_IND_INFO) # migrate farm animals to local wild
    
                        # delete byrid farm-wild1 pop
                        pop.removeSubPops(hybrid_index)
    
                 # track individual information for life history table
                if make_report_all_inds:
                    track_all_inds(sp_names=w_farm_sp_names, rep=rep, 
                                   year_counter=year_counter, month=month)
    
                if month == 'Dec': # during farm
                                    
                    # store response variables
                    store_rvars(pop=pop, farm_boo=True, farm_idx=farm_idx, rep=rep, 
                                year_counter=year_counter, afs_dict=afs_dict,
                                loci_to_track=biallelic_loci, rvar_dict=pop_rvar_dict, 
                                pair_rvar_dict=pop_pair_rvar_dict, wild_subPopIndeces=wo_farm_sp_idcs,
                                num_aloci=num_aloci, sp_names=w_farm_sp_names)    
                    harvest_dict[rep][year_counter] = harvest_this_year
                    
                    if make_report_adults:
                        track_adults(pop=pop, sp_names=w_farm_sp_names, ageX=ageX, rep=rep, year_counter=year_counter)
              
                    # 20201222 - report which inds broodstock in each year
                    for ind in pop.individuals([farm_idx]):
                        if ind.broodstock == 1:
                            bstock_report.write('\t'.join([str(rep),str(year_counter),str(ind.ind_id)]))
                            bstock_report.write('\n')
                                    
            
            if year_counter % yr_interval_log == 0:
                print('   Yr', year_counter, datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))
                write_to_log(logfilename=logfilename, pop=pop, farm_idx=farm_idx, 
                             farm_boo=True, farm_phase='during-farm')
    
            year_counter += 1
        
        # remove farm subpop
        pop.removeSubPops(farm_idx)
        
        ################################################################################
        #### ------------------------ Post-farm period ---------------------------- ####
        ################################################################################
                
        # print('  Post-farm period')
        
        progress_report.write(" Post-farm period\n")

        progress_report.flush()
        for year in np.arange(post_farm_years):
        
            for month in prob_repro_by_month:
                
                # track individual information for life history table
                if make_report_all_inds:
                    track_all_inds(sp_names=wo_farm_sp_names, rep=rep, 
                                   year_counter=year_counter, month=month)
                
                # migrate, grow, age, and die
                wild_mig(pop=pop, farm_boo=False, farm_idx=farm_idx, wild_names=wild_names, 
                         wild_mig_rate_L=wild_mig_rate_L, wild_mig_rate_J=wild_mig_rate_J, 
                         wild_mig_rate_A=wild_mig_rate_A, w_farm_wld_sp_idcs=w_farm_wld_sp_idcs) # migrate               
                pop.setIndInfo(values=list(np.array(pop.indInfo(field='age'))+1/12), field='age') # age
                assign_fitness_all(pop) # assign fitness for selective mortality
                mortality(pop=pop, lar_annual_M=lar_annual_M, juv_annual_M=juv_annual_M, 
                          adult_annual_M=adult_annual_M, farm_reduced_mort=farm_reduced_mort) # die
                
                # reproduce based on probability to reproduce by month
                if random.random() < prob_repro_by_month[month]:
                    pop.evolve(preOps=[sim.PySelector(func=get_fitness)],
                            matingScheme=sim.HeteroMating([sim.CloneMating(subPops=wo_farm_sps_to_clone, weight=-1),
                                                           sim.RandomMating(subPops=wo_farm_sps_to_mate,
                                                                            numOffspring=numWildOffspring_par,
                                                                             ops=[sim.MendelianGenoTransmitter(),
                                                                                     sim.IdTagger(),
                                                                                  sim.PedigreeTagger(),
                                                                                 sim.SummaryTagger(mode=sim.MEAN, 
                                                                                                   infoFields=['age',
                                                                                                              'meanParentAge'])])],
                                                         subPopSize=wo_farm_recruitment(pop)),
                            postOps=[sim.Stat(alleleFreq=biallelic_loci, subPops=wo_farm_sp_idcs, 
                                              vars=['alleleFreq','alleleFreq_sp']),
                                 sim.Stat(genoFreq=range(num_aloci), subPops=wo_farm_sp_idcs,
                                          vars=['genoFreq','genoFreq_sp']),
                                 sim.Stat(heteroFreq=loci_for_ibd, vars='heteroFreq_sp'),
                                 sim.Stat(structure=biallelic_loci, subPops=[0,1], 
                                          suffix='_wild1_wild2', vars=['F_st']),
                                 sim.Stat(structure=biallelic_loci, subPops=[1,2], 
                                          suffix='_wild2_wild3', vars=['F_st']),
                                 sim.Stat(structure=biallelic_loci, subPops=[0,2], 
                                          suffix='_wild1_wild3', vars=['F_st']),
                                 sim.Stat(structure=range(num_aloci), subPops=[0,1], 
                                          suffix='_aL_wild1_wild2', vars=['F_st']),
                                 sim.Stat(structure=range(num_aloci), subPops=[1,2], 
                                          suffix='_aL_wild2_wild3', vars=['F_st']),
                                 sim.Stat(structure=range(num_aloci), subPops=[0,2], 
                                          suffix='_aL_wild1_wild3', vars=['F_st']),
                                 sim.Stat(structure=loci_for_fst, subPops=[0,1], 
                                          suffix='_nL_wild1_wild2', vars=['F_st']),
                                 sim.Stat(structure=loci_for_fst, subPops=[1,2], 
                                          suffix='_nL_wild2_wild3', vars=['F_st']),
                                 sim.Stat(structure=loci_for_fst, subPops=[0,2], 
                                          suffix='_nL_wild1_wild3', vars=['F_st']),
                                 sim.Stat(popSize=True, subPops=wo_farm_sp_idcs)],
                            gen=1)
                    
                                    
                    for ind in pop.individuals():
                        if ind.age == 0:
                            ind.cohortYear = year_counter
                    
                if month == 'Dec': # post farm
                    
                    # store response variables
                    store_rvars(pop=pop, farm_boo=False, farm_idx=farm_idx, rep=rep, year_counter=year_counter, 
                                afs_dict=afs_dict, loci_to_track=biallelic_loci, rvar_dict=pop_rvar_dict, 
                                pair_rvar_dict=pop_pair_rvar_dict, wild_subPopIndeces=wo_farm_sp_idcs,
                                num_aloci=num_aloci, sp_names=wo_farm_sp_names)
                    
                    # track adults for report
                    if make_report_adults:
                        track_adults(pop=pop, sp_names=wo_farm_sp_names, ageX=ageX, rep=rep, year_counter=year_counter)
                
            if year_counter % yr_interval_log == 0:
                print('   Yr', year_counter, datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))
                write_to_log(logfilename=logfilename, pop=pop, farm_idx=farm_idx, farm_boo=False, farm_phase='post-farm')
            year_counter += 1    
            
        ################################################################################
        #### ------------------- Temporal Fst, across years ----------------------- ####
        ################################################################################    
        
        # load time zero pop, rename the subpops, and combine with time one (end) pop
        pop_t0 = sim.loadPopulation(os.path.join(batch_dir,''.join(['batch_', batch, '_t0.pop'])))
        for name in pop_t0.subPopNames():
            # add _t0 to subpop names so distinct upon merge
            pop_t0.setSubPopName(''.join([''.join([name,'_t0'])]), pop.subPopNames().index(name)) 
        t1_names = pop.subPopNames()
        t0_names = [''.join([name, '_t0']) for name in t1_names]
    
        # merge population objects, now that time zero subpops renamed with t0 identifier
        pop.addIndFrom(pop_t0)
        
        # get pairwise Fst between t0 and now in each wild subpop
        for i in range(len(t1_names)):
            sim.stat(pop, structure=biallelic_loci, subPops=[i,i+len(t1_names)], 
                     suffix=''.join(['_temp_sp_', str(i)]), vars=['F_st'])
            temp_rvar_dict[rep][i]['Fst'] = pop.vars()[''.join(['F_st_temp_sp_', str(i)])]    
            
    ################################################################################
    #### --------------------- Write results to file -------------------------- ####
    ################################################################################
    
    # for pop_rvar_dict 
    pop_rvar_outfile = open(os.path.join(batch_dir,''.join(['batch_', batch, '_pop_rvars.txt'])), 'w')
    pop_rvar_outfile.write('Srep\tRep\tYear\tSubpop\tRvar\tValue\n')
    for rep in pop_rvar_dict:
        for year in pop_rvar_dict[rep]:
            for sp in pop_rvar_dict[rep][year]:
                for rvar in pop_rvar_dict[rep][year][sp]:
                    if sp !='farm':
                        sp_name = ''.join(['wild', str(sp+1)])
                    else:
                        sp_name = sp
                    pop_rvar_outfile.write('\t'.join([str(coreid), str(rep), str(year), sp_name,  
                                                      rvar, str(pop_rvar_dict[rep][year][sp][rvar])]))
                    pop_rvar_outfile.write('\n')
    pop_rvar_outfile.close()
    
    # for harvest
    harvest_outfile = open(os.path.join(batch_dir,''.join(['batch_', batch, '_harvest.txt'])), 'w')
    harvest_outfile.write('Srep\tRep\tYear\tIndsHarvested\n')
    for rep in harvest_dict:
        for year in harvest_dict[rep]:
            harvest_outfile.write('\t'.join([str(coreid), str(rep), str(year), str(harvest_dict[rep][year])]))
            harvest_outfile.write('\n')
    harvest_outfile.close()
    
    if track_AFs == True:
        afs_outfile = open(os.path.join(batch_dir,''.join(['batch_', batch, '_AFs.txt'])), 'w')
        afs_outfile.write('Srep\tRep\tYear\tSubpop\tLocus_index\tAllele\tAdaptive\tAdv\tAlleleFrequency\n')
        for rep in afs_dict:
            for year in afs_dict[rep]:
                for sp in afs_dict[rep][year]:
                    if sp !='farm':
                        sp_name = ''.join(['wild', str(sp+1)])
                    else:
                        sp_name = sp
                    for l_idx in biallelic_loci:
                        if l_idx in loci_for_fst:
                            adaptive_locus = False
                            for allele in [0,1]:
                                advantage = "NA"
                                afs_outfile.write('\t'.join([str(coreid), str(rep), str(year), sp_name, 
                                                             str(l_idx), str(allele), str(adaptive_locus),
                                                             str(advantage), 
                                                             str(afs_dict[rep][year][sp][l_idx][allele])]))
                                afs_outfile.write('\n')  
                                 
                        else:
                            adaptive_locus = True
                            for allele in [0,1]:
                                stage = stages[l_idx//num_aloci_stage]
                                within_stage_l_idx = l_idx % num_aloci_stage
                                if fit_dict[sp_name][stage][within_stage_l_idx][str(allele)] == 1:
                                    advantage = True
                                else:
                                    advantage = False
                                afs_outfile.write('\t'.join([str(coreid), str(rep), str(year), sp_name, 
                                                             str(l_idx), str(allele), str(adaptive_locus),
                                                             str(advantage), 
                                                             str(afs_dict[rep][year][sp][l_idx][allele])]))
                                afs_outfile.write('\n')    
        afs_outfile.close()
        
    # for pop_pair_rvar_dict
    pop_pair_rvar_outfile = open(os.path.join(batch_dir,''.join(['batch_', batch, '_pop_pair_rvars.txt'])), 'w')
    pop_pair_rvar_outfile.write('Srep\tRep\tYear\tPop_pair\tRvar\tValue\n')
    for rep in pop_pair_rvar_dict:
        for year in pop_pair_rvar_dict[rep]:
            for pop_pair in pop_pair_rvar_dict[rep][year]:
                for rvar in pop_pair_rvar_dict[rep][year][pop_pair]:
                    pop_pair_rvar_outfile.write('\t'.join([str(coreid), str(rep), str(year), pop_pair, rvar, 
                                                           str(pop_pair_rvar_dict[rep][year][pop_pair][rvar])]))
                    pop_pair_rvar_outfile.write('\n')
    pop_pair_rvar_outfile.close()
    
    # for temporal fst (and any other temporal rvars per subpop; none yet)
    temp_rvar_outfile = open(os.path.join(batch_dir,''.join(['batch_', batch, '_temp_rvars.txt'])), 'w')
    temp_rvar_outfile.write('Srep\tRep\tSubpop\tRvar\tValue\n')
    for rep in temp_rvar_dict:
        for sp in temp_rvar_dict[rep]:
            for rvar in temp_rvar_dict[rep][sp]:
                if sp !='farm':
                    sp_name = ''.join(['wild',str(sp+1)])
                else:
                    sp_name = sp
                temp_rvar_outfile.write('\t'.join([str(coreid), str(rep), sp_name, rvar, str(temp_rvar_dict[rep][sp][rvar])]))
                temp_rvar_outfile.write('\n')           
    temp_rvar_outfile.close()
    
    # close report files
    bstock_report.close() 
    if make_report_adults:
        track_for_Ne.close() 
    if make_report_all_inds:
        life_hist_report.close() 
    
    # print('Done! Run time:', round(time.time() - startTime, 2))
    
    # progress_report.write("Done! Run Time: ".join([str(round(time.time() - startTime, 2))]))

    progress_report.write("Done! Run Time: " + str(round(time.time() - startTime, 2)))

    progress_report.close()
