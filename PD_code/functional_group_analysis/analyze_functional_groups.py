#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 10:02:52 2019

@author: jglopez
"""

import pandas as pd
import subprocess

import statsmodels.api as sm1
import statsmodels as sm2
import os


exclude_steroids = False

#Load SMILES of drugs from screen
drug_file = 'drug_SMILES.tsv';
drug_SMILES = pd.read_csv(drug_file,delimiter = '\t')

if exclude_steroids:
    drug_SMILES = drug_SMILES.loc[~drug_SMILES['drug_class'].str.contains('steroid|Steroid'),:]
    drug_SMILES = drug_SMILES.reset_index(drop=True)


#Load groups for openbabel analysis
groups = pd.read_csv('SMARTS_groups.csv',sep = ',')

#Get rid of any duplicates in the groups
groups.drop_duplicates(subset ="name", keep = 'first', inplace = True) 
groups = groups.reset_index(drop=True)
drug_list = drug_SMILES['name'].tolist()
group_list = groups['name'].tolist()

subprocess.call('rm babel_results/*', shell=True)

valid_SMARTS = list()
for i in range(0,len(groups['name'])):
    
    #Search against drug list for functional group
    group_name = groups.loc[i,'name']
    group_SMARTS = groups.loc[i,'SMARTS']
    results_file = 'babel_results/' + group_name.replace(' ','_') + '_results.txt'
    babel_cmd = 'obgrep -n ' + '\'' + group_SMARTS+ '\' ' + drug_file + ' > ' + results_file
    subprocess.call(babel_cmd, shell=True)
    
    #Load the babel results and add to main df
    if os.path.isfile(results_file): 
        if os.stat(results_file).st_size != 0:
            temp_df = pd.read_csv(results_file,sep = '\t', header = None)    
            positive_drugs = list(set(temp_df[0].tolist()) & set(drug_list))
            positive_index = [drug_list.index(x) for x in positive_drugs]
            drug_SMILES.loc[:,group_name] = 0
            drug_SMILES.loc[positive_index,group_name] = 1

#Select out groups with more than one member, and more than one non-member
test_groups = drug_SMILES.iloc[:,4:].copy()
test_groups = test_groups.loc[:,test_groups.sum(axis=0) > 1] 
test_groups = test_groups.loc[:,test_groups.sum(axis=0) < (len(drug_list)-1)]  
test_groups = test_groups.columns

#Build final results dataframe
results_df = pd.DataFrame(columns = 
                             ['+group/+mod','+group','-group/+mod',
                              '-group','prop_group','prop_nongroup','direction',
                              'p_val','corrected_p_val','structural_class',
                              'SMARTS'],index = test_groups)

#Compute enrichment statistics and add in group metadata 
for i in range(0,len(test_groups)):
    group_name = test_groups[i]    
    group_df = drug_SMILES.loc[drug_SMILES[group_name] ==1 ,:].copy()
    nongroup_df = drug_SMILES.loc[drug_SMILES[group_name] ==0 ,:].copy()
    results_df.loc[group_name,'+group/+mod'] = group_df['modification'].sum()
    results_df.loc[group_name,'+group'] = len(group_df.index)
    results_df.loc[group_name,'-group/+mod'] = nongroup_df['modification'].sum()
    results_df.loc[group_name,'-group'] = len(nongroup_df.index)  
    results_df.loc[group_name,'prop_group'] = results_df.loc[group_name,'+group/+mod'] / results_df.loc[group_name,'+group']
    results_df.loc[group_name,'prop_nongroup'] = results_df.loc[group_name,'-group/+mod'] / results_df.loc[group_name,'-group']
    
    if results_df.loc[group_name,'prop_group'] > results_df.loc[group_name,'prop_nongroup']:
        results_df.loc[group_name,'direction'] = 'enriched'
    elif results_df.loc[group_name,'prop_group'] < results_df.loc[group_name,'prop_nongroup']: 
        results_df.loc[group_name,'direction'] = 'depleted'
        
    s_vector = [results_df.loc[group_name,'+group/+mod'],
                results_df.loc[group_name,'-group/+mod']]
    n_vector = [results_df.loc[group_name,'+group'],
                results_df.loc[group_name,'-group']]    
    z, results_df.loc[group_name,'p_val'] = sm1.stats.proportions_ztest(
            s_vector, n_vector, alternative = 'two-sided')                     
    
    matching_index = groups.name == group_name
    results_df.loc[group_name,'structural_class'] = groups['structural_class'][matching_index].values[0]
    results_df.loc[group_name,'SMARTS'] = groups['SMARTS'][matching_index].values[0]

#Perform multiple hypothesis testing
multitest_results = sm2.stats.multitest.multipletests(results_df.p_val.tolist(), alpha=0.01, method='fdr_bh', 
                                           is_sorted=False, returnsorted=False)      
results_df.corrected_p_val = multitest_results[1]        

#Sort by the corrected p value and save results
results_df.sort_values(by = 'corrected_p_val',inplace = True)
significant_results = results_df[results_df.corrected_p_val < 0.01]

if exclude_steroids:
    results_df.to_csv('functional_analysis_results_exclude_steroids.csv')

else:
    results_df.to_csv('functional_analysis_results_include_steroids.csv')
   
    
