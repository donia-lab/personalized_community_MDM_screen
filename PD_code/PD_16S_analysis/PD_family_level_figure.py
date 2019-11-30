#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 10:30:00 2018

@author: jglopez
"""

import pandas as pd
import numpy as np
import matplotlib.patches as patches
import matplotlib.pyplot as plt

#Import table of family-level otus
otutable = pd.read_table('family_table_split_read.txt',sep = '\t',
                        index_col=0)

#Pull out desired samples from master table
times = ["01","02","Feces"]
desired_otus = otutable[otutable.columns[otutable.columns.str.contains(pat = 'HD1')]]
desired_otus = desired_otus[desired_otus.columns[desired_otus.columns.str.contains(pat = '|'.join(times))]]
desired_otus = desired_otus[desired_otus.columns[~desired_otus.columns.str.contains(pat = 'PBSr|PSBr|NYU|MD102|MD101')]]

#Filter to only have otus above 1%
abundant_otus = desired_otus.loc[(desired_otus > 0.01).any(axis=1)].copy()

#Rename samples to shorter names
for column in abundant_otus.columns:
        abundant_otus.rename(columns = {column: '.'.join(column.split('.')[1:3])},inplace=True)
abundant_otus.rename(columns = {'HD1.Feces.extendedFrags': 'HD1.Feces'},inplace=True)
abundant_otus.rename(columns = {'BestMix.01': 'BMix.01'},inplace=True)
abundant_otus.rename(columns = {'BestMix.02': 'BMix.02'},inplace=True)
abundant_otus.rename(columns = {'GAM.01': 'mGAM.01'},inplace=True)
abundant_otus.rename(columns = {'GAM.02': 'mGAM.02'},inplace=True)

#Rename families to only have family name
for index in abundant_otus.index:
    abundant_otus.rename(index = {index : index.split(';')[4][3:]}, inplace=True)

#Add "other" category and drop things not classified at the family level
abundant_otus.loc['Other'] = 1- abundant_otus.sum() + abundant_otus.loc['']
abundant_otus = abundant_otus.drop('')

#Prepare matrix form versions of the data
HD1 = desired_otus['HD1.Feces'].values
otu_array = desired_otus.values  

#Compute Jensen-Shannon divergence
JS = np.zeros(len(otu_array[1,:]))
for i in range(0,len(otu_array[1,:])):
    sum1 = 0
    sum2 = 0
    av_dist = (HD1 + otu_array[:,i])/2
    for j in range(0,len(HD1)):
        if HD1[j] == 0 or av_dist[j] == 0: 
            pass
        else:
            sum1 = sum1 + HD1[j]*np.log(av_dist[j]/HD1[j])

        if otu_array[j,i] == 0 or av_dist[j] == 0:
            pass
        else:
            sum2 = sum2 + otu_array[j,i]*np.log(av_dist[j]/otu_array[j,i])
    JS[i] = -0.5*sum1 - 0.5*sum2
    #JS[i] = np.linalg.norm(otu_array[:,i]-HD1,ord=1)


#Save unsorted version for reference
unsorted_otus = abundant_otus.copy()

#Sort df according to JS divergences
abundant_otus = abundant_otus.iloc[:,np.argsort(JS)]


#Plot stacked bar plot
ax = abundant_otus.T.plot.bar(stacked=True, 
                          colormap = "tab20",legend=False)
ax.legend(bbox_to_anchor=(1.05, 1.1))

axfig = ax.get_figure()

rect = patches.Rectangle((0.5,-0.02),2,1.04,linewidth=1.5,edgecolor='k',facecolor='none',clip_on=False,linestyle='--')
ax.add_patch(rect)

#Build JS inset axis
axins = axfig.add_axes([0.171, 0.92, 0.75, 0.1])
axins.spines['right'].set_visible(False)
axins.spines['bottom'].set_visible(False)
axins.spines['top'].set_visible(False)
axins.axes.get_xaxis().set_visible(False)
axins.set_ylabel('$D_{JS}$',rotation=0)
axins.yaxis.set_label_coords(-0.1,0.3,transform = None)
axins.set_ylim([0,1])
axins.set_yticks([0,1])
x = list(range(0,len(JS)-1))
wedgecolor = 'black'
lw = 3
plt.plot(x,JS[np.argsort(JS)][1:len(JS)],color = wedgecolor,linewidth = lw)
plt.plot([x[0],x[len(x)-1]],[0,0],color=wedgecolor,clip_on=False,linewidth = lw)
plt.plot([len(JS)-2, len(JS)-2],[0,JS[np.argsort(JS)][-1]],color=wedgecolor,linewidth = lw)
plt.plot([0, 0],[0,JS[np.argsort(JS)][1]],color=wedgecolor,linewidth = lw)

#Modified axis and label properties
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xticklabels(ax.get_xticklabels(),rotation=90)
ax.set_yticklabels(['0','1'])
ax.set_yticks([0,1])
ax.set_ylim([0,1])
foffset = -1.5
ax.set_xlim([foffset-0.6,ax.get_xlim()[1]])
ax.set_xticks([foffset+0.25] + list(range(1,len(JS))))

#Move reference composition away from main group
for container in ax.containers:
    plt.setp(container.get_children()[0],x=foffset)

axfig.savefig('PD_family_level_figure.eps',dpi=600,bbox_inches = 'tight',format = 'eps')
