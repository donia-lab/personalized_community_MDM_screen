#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 08:49:34 2019

@author: jglopez
"""

#This script tests whether different classes of drug origin (natural, derived,
# and synthetic) are enriched or depleted for MDM modifications

import pandas as pd
import statsmodels.api as sm
import statsmodels as sm2


#Define functions to be used in analysis
def get_mdm_frac(drug_df):
    frac = len(drug_df[drug_df['positive'] == 'Y'].index)/len(drug_df.index)
    return frac

def mdm_stats(df1,df2):
    n1 = len(df1.index)
    s1 = len(df1[df1['positive']=='Y'].index)
    n2 = len(df2.index)
    s2 = len(df2[df2['positive']=='Y'].index)
    z, p_value = sm.stats.proportions_ztest([s1, s2], [n1, n2],
                                            alternative = 'two-sided')
    return p_value

#Import data and get dataframes of different classes    
drug_classes = pd.read_csv('drug_classes.csv', names = 
                           ['code','drug','mw','positive','origin','class',
                            'growth'])

abx_drugs = drug_classes[drug_classes['class'] == 'antibiotic']

natural_drugs = drug_classes[drug_classes['origin'] == 'natural']
derivative_drugs = drug_classes[drug_classes['origin'] == 'derivative']
synthetic_drugs  = drug_classes[drug_classes['origin'] == 'synthetic']
nat_deriv_drugs = drug_classes[drug_classes['origin'] != 'synthetic']

steroid_drugs = drug_classes[drug_classes['class'].str.contains(pat = 'steroid')]
non_steroid_drugs = drug_classes[~drug_classes['class'].str.contains(pat = 'steroid')]
non_steroid_natderiv = non_steroid_drugs[non_steroid_drugs['origin'] != 'synthetic']
non_steroid_synth = synthetic_drugs[~synthetic_drugs['class'].str.contains(pat = 'steroid')]

nat_frac = get_mdm_frac(natural_drugs)
deriv_frac = get_mdm_frac(derivative_drugs)
synth_frac = get_mdm_frac(synthetic_drugs)
natderiv_frac = get_mdm_frac(nat_deriv_drugs)
steroid_frac = get_mdm_frac(steroid_drugs)
non_steroid_frac = get_mdm_frac(non_steroid_drugs)
non_s_natderiv_frac = get_mdm_frac(non_steroid_natderiv)
non_s_synth_frac = get_mdm_frac(non_steroid_synth)


nat_vs_deriv = mdm_stats(natural_drugs,derivative_drugs)
nat_vs_synth = mdm_stats(natural_drugs,synthetic_drugs)
deriv_vs_synth = mdm_stats(derivative_drugs,synthetic_drugs)
natderiv_vs_synth = mdm_stats(nat_deriv_drugs,synthetic_drugs)
natderiv_vs_synth_no_steroid = mdm_stats(non_steroid_natderiv,non_steroid_synth)
steroid_vs_non_steroid = mdm_stats(steroid_drugs,non_steroid_drugs)

#Test across all classes
class_comparison = pd.DataFrame(
        index = list(set(drug_classes['class'].tolist())),columns = ['n_total','frac','p'])
        
for drug_class in class_comparison.index:
    drug_df = drug_classes[drug_classes['class'] == drug_class]
    non_drug_df = drug_classes[drug_classes['class'] != drug_class]
    class_comparison.loc[drug_class,'n_total'] = len(drug_df.index)
    class_comparison.loc[drug_class,'frac'] = get_mdm_frac(drug_df)
    class_comparison.loc[drug_class,'p'] = mdm_stats(drug_df,non_drug_df)

multitest_results = sm2.stats.multitest.multipletests(class_comparison.p.tolist(), alpha=0.01, method='fdr_bh', 
                                           is_sorted=False, returnsorted=False)   
class_comparison['corrected_p'] = multitest_results[1]

    
#Make figure of nat/deriv vs synthetic MDM+ with all drugs
fontsize = 16
frac_df = pd.DataFrame({'frac':[natderiv_frac,synth_frac],
                       '': ['Natural/Derivative','Synthetic']})
ax = frac_df.plot.bar(x = '',y='frac',legend=None,width= 0.65,
                      color = '0.5',rot=0,fontsize=fontsize,
                      yticks = [0, 0.1, 0.2, 0.3])

ax.set_ylabel("Fraction MDM positive",rotation=90,fontsize=fontsize)
         
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

above1 = 0.04
above2 = 0.02
top_all = natderiv_frac + above1
top_nat = natderiv_frac + above2
top_synth = synth_frac + above2
ax.plot([0,1],[top_all, top_all],'k')
ax.plot([0,0],[top_nat,top_all],'k')
ax.plot([1,1],[top_synth,top_all],'k')
ax.set_ylim([0,0.3])
ax.text(0.5,top_all,'***',horizontalalignment =  'center',fontsize = 14)

axfig = ax.get_figure()
axfig.savefig('mdm_positive_fraction.eps',dpi=600,bbox_inches = 'tight',format = 'eps')


#Make figure of steroid vs non-steroid MDM+ fraction
frac_df = pd.DataFrame({'frac':[steroid_frac,non_steroid_frac],
                       '': ['Steroid','Non-steroid']})
ax1 = frac_df.plot.bar(x = '',y='frac',legend=None,width= 0.65,
                      color = '0.5',rot=0,fontsize=fontsize,
                      yticks = [0, 0.2, 0.4,0.6])

ax1.set_ylabel("Fraction MDM positive",rotation=90,fontsize=fontsize)
         
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

above1 = 0.04
above2 = 0.02
top_all = steroid_frac + above1
top_ster = steroid_frac + above2
top_non = non_steroid_frac + above2
ax1.plot([0,1],[top_all, top_all],'k')
ax1.plot([0,0],[top_ster,top_all],'k')
ax1.plot([1,1],[top_non,top_all],'k')
ax1.set_ylim([0,0.6])
ax1.text(0.5,top_all,'***',horizontalalignment =  'center',fontsize = 14)

axfig1 = ax1.get_figure()
axfig1.savefig('steroid_vs_non_fraction.eps',dpi=600,bbox_inches = 'tight',format = 'eps')


#Make figure of nat/deriv vs. synthetic with no steroids
frac_df = pd.DataFrame({'frac':[non_s_natderiv_frac,non_s_synth_frac],
                       '': ['Natural/Derivative','Synthetic']})
ax2 = frac_df.plot.bar(x = '',y='frac',legend=None,width= 0.65,
                      color = '0.5',rot=0,fontsize=fontsize,
                      yticks = [0, 0.1, 0.2,0.3])

ax2.set_ylabel("Fraction MDM positive",rotation=90,fontsize=fontsize)
         
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)

above1 = 0.04
above2 = 0.02
top_all = non_s_natderiv_frac + above1
top_ster = non_s_natderiv_frac + above2
top_non = non_s_synth_frac + above2
ax2.plot([0,1],[top_all, top_all],'k')
ax2.plot([0,0],[top_ster,top_all],'k')
ax2.plot([1,1],[top_non,top_all],'k')
ax2.set_ylim([0,0.3])
ax2.text(0.5,top_all+0.01,'n.s.',horizontalalignment =  'center',fontsize = 14)

axfig2 = ax2.get_figure()
axfig2.savefig('no_steroid_synthetic_vs_natderiv.eps',dpi=600,bbox_inches = 'tight',format = 'eps')

