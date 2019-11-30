#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 16:09:55 2018

@author: jglopez
"""

import pandas as pd
import numpy as np

n_point = 80
n_iter = 150

div_rare_unformatted = pd.read_table('shannon_80_40k.csv',sep = ',',
                        index_col=0,dtype=None)

new_index = [val.split('-')[1].split('_')[0] for val in div_rare_unformatted.columns]
div_rare = pd.DataFrame(columns = div_rare_unformatted.index, index = new_index, data = None)

for column in div_rare.columns:
    div_rare[column] = div_rare_unformatted.loc[column].values
    
div_rare = div_rare[div_rare.columns[div_rare.columns.str.contains(pat = 'HD1|seq')]]
div_rare = div_rare[div_rare.columns[div_rare.columns.str.contains(pat = '02|01|Feces|seq')]]
div_rare = div_rare[div_rare.columns[~div_rare.columns.str.contains(pat = 'PBSr|PSBr')]]    
#div_rare = div_rare[div_rare.columns[~div_rare.columns.str.contains(pat = 'BB|MRS|RCM')]]

for column in div_rare.columns:
        div_rare.rename(columns = {column: '.'.join(column.split('.')[1:3])},inplace=True)
div_rare.rename(columns = {'FecesextendedFrags': 'Feces'},inplace=True)
div_rare.rename(columns = {'BestMix.01': 'BMix.01'},inplace=True)
div_rare.rename(columns = {'BestMix.02': 'BMix.02'},inplace=True)
div_rare.rename(columns = {'GAM.01': 'mGAM.01'},inplace=True)
div_rare.rename(columns = {'GAM.02': 'mGAM.02'},inplace=True)

new_order = ['Feces', 'mGAM.02', 'mGAM.01', 'TYG.01', 'TYG.02', 'GMM.02', 'BB.01',
       'BHI.02', 'M17.02', 'Meat.02', 'BMix.02', 'GMM.01', 'BHI.01', 'LB.02',
       'RCM.01', 'BB.02', 'Meat.01', 'BMix.01', 'RCM.02', 'Liver.02',
       'Liver.01', 'BT.02', 'BT.01', 'LB.01', 'M17.01', 'MRS.01', 'TB.02',
       'TB.01', 'MRS.02']
div_rare = div_rare[new_order]

div_rare = div_rare.astype(float,errors='ignore')
div_rare = div_rare.replace('n/a',np.NaN)
av_rare = pd.DataFrame(columns=div_rare.columns,index =range(0,n_point))

#Get mean of rarefaction iterations
for column in div_rare.columns:
    for i in range(0,n_point):
        av_rare[column][i] = pd.to_numeric(div_rare[column][n_iter*i:(n_iter*i)+n_iter]).mean()
av_rare['reads_sampled'] = np.unique(div_rare.index.astype(float))

ax = av_rare.plot(x = 'reads_sampled', colormap='tab10',kind = 'line')        
ax.get_lines()[1].set_linestyle('--')
ax.get_lines()[1].set_color('black')
ax.get_lines()[2].set_color('black')

linestyles = ['-','--',':']

for i in range(0,29):
    ax.get_lines()[i].set_linestyle(linestyles[i % len(linestyles)])

ax.legend(bbox_to_anchor=(1.07, 1.025),loc = 'upper left',ncol=2)
ax.set_xlim([0,40000])
ax.set_ylabel('Shannon Diversity (ASV level)',rotation=90)
ax.set_xlabel('Sequences Sampled',rotation=0)

axfig = ax.get_figure()
axfig.savefig('PD_rarefaction_figure.png',dpi=600,bbox_inches = 'tight',format='png')