#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 18:14:37 2018

@author: jglopez
"""

import pandas as pd
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt

#Load in asv table
otutable = pd.read_table('asv_table_split_read.txt',sep = '\t',
                        index_col=0)

#Sort out desired samples
times = ["02","Feces"]
desired_otus = otutable[otutable.columns[otutable.columns.str.contains(pat = '|'.join(times))]]
desired_otus = desired_otus[desired_otus.columns[~desired_otus.columns.str.contains(pat = 'PBSr|PSBr|NYU|MD102|MD101')]]

#Remove asvs not present in any sample
present_otus = desired_otus.loc[(desired_otus > 0).any(axis=1)]

#Create relative abundance table
rel_otus = present_otus.copy()
for column in rel_otus.columns:
    rel_otus[column] = present_otus[column]/present_otus[column].sum()

#Rename samples    
for column in rel_otus.columns:
        rel_otus.rename(columns = {column: '.'.join(column.split('.')[1:3])},inplace=True)
rel_otus.rename(columns = {'BestMix.02': 'BMix.02'},inplace=True)
rel_otus.rename(columns = {'GAM.02': 'mGAM.02'},inplace=True)


#Compute shared ASVs between GAM day 2 and HD-1
total_shared= 0;
no_above_1pc = 0;
no_above_1pc_shared = 0;
for index in rel_otus['Feces'].index:
    if rel_otus['mGAM.02'][index] > 0:
        total_shared = total_shared + rel_otus['Feces'][index]
    if rel_otus['Feces'][index] > 0.01:
        no_above_1pc = no_above_1pc + 1
        if rel_otus['mGAM.02'][index] > 0:
            no_above_1pc_shared = no_above_1pc_shared + 1

print('Of the ' + str(no_above_1pc) + ' ASVs present above 1% in the original fecal sample, '
      + str(no_above_1pc_shared) + ' of them are present in GAM day 2 condition.' 
      + ' The shared ASVs between the original fecal sample and GAM day 2 account for '
      + str(total_shared*100)[0:5] + '% of the original fecal sample.')



#Select only asvs aboe 1% in the feces
HD1_asv_selection = rel_otus.loc[(rel_otus> 0.01)['Feces']].copy()


#Sort all by feces abundance
HD1_asv_selection = HD1_asv_selection.sort_values(by = ['Feces'],ascending=False)

#Create "other" category
HD1_asv_selection.loc['Other'] = 1- HD1_asv_selection.sum()


#Compute shannon entropy with all asvs
asv_shannon = pd.DataFrame(index = rel_otus.columns,dtype=float,
                           data=None,columns=['H'])
for i in range(0,len(rel_otus.iloc[1,:])):
    sum = 0
    for j in range(0,len(rel_otus.iloc[:,1])):
        if rel_otus.iloc[j,i] > 0:
           sum = sum +  rel_otus.iloc[j,i]*np.log2(rel_otus.iloc[j,i])
        else:
            pass
    #asv_shannon.iloc[i] = math.pow(2,-sum)
    asv_shannon.iloc[i] = -sum

#Resort by shannon entropy
asv_shannon = asv_shannon.sort_values(by = ['H'],ascending=False)
HD1_asv_selection = HD1_asv_selection[asv_shannon.index]


fig, ax = plt.subplots()
cmap = matplotlib.cm.get_cmap('gist_rainbow')
pos_x = 100
pos_y = 100
otu_num = len(HD1_asv_selection.iloc[:,1])
sample_num = len(HD1_asv_selection.iloc[1,:])
len_x = 6.0
len_y = 6.0
boxes_y = 4
boxes_x = 4
box_width = 1
spacing_x = 1.6
spacing_y = 2.1
box_buffer = 1.5
def_size = 5
fig_pad_x = 2
fig_pad_y = 2
power = 1/3
size_list = []

#Make the collector boxes
for i in range(0,sample_num):
    box_pos_x = pos_x + (i % boxes_x)*spacing_x*len_x
    box_pos_y = pos_y - math.floor(i/boxes_x)*spacing_y*len_y
    t_l = [box_pos_x - box_buffer, box_pos_y + box_buffer]
    b_l =  [box_pos_x - box_buffer, box_pos_y - len_y + 1 - box_buffer]
    t_r = [box_pos_x + len_x + box_buffer - 1,box_pos_y + box_buffer]
    b_r = [box_pos_x + len_x + box_buffer -1 ,box_pos_y - len_y + 1 - box_buffer]
    
    ax.plot([t_l[0],t_r[0]],[t_l[1],t_r[1]],'k',linewidth = box_width)
    ax.plot([t_l[0],b_l[0]],[t_l[1],b_l[1]],'k',linewidth = box_width)
    ax.plot([b_l[0],b_r[0]],[b_l[1],b_r[1]],'k',linewidth = box_width)
    ax.plot([t_r[0],b_r[0]],[t_r[1],b_r[1]],'k',linewidth = box_width)
    
    plt.text(t_l[0]/2 + t_r[0]/2,t_l[1]+2,HD1_asv_selection.columns[i], 
             horizontalalignment = 'center',fontsize= 5.5)
    plt.text(t_l[0]/2 + t_r[0]/2,t_l[1]+0.5,'$H = ' 
             + str(round(asv_shannon.iloc[i].values[0],2)) + '$', 
             horizontalalignment = 'center',fontsize= 5.5)
    for j in range(0,otu_num):
        if HD1_asv_selection.iloc[j,i] > 0:
            if j == otu_num - 1:
                ASV_color = [0.7, 0.7, 0.7]
            else:
                ASV_color = cmap(j/otu_num)
            ASV_markersize = 1.05*def_size - math.pow(-np.log10(HD1_asv_selection.iloc[j,i]),power)*(def_size/math.pow(5,power)) 
            size_list.append(ASV_markersize)
            ax.plot(box_pos_x+(j % len_x),box_pos_y - math.floor(j/len_y),'o',
            color = ASV_color, markersize = ASV_markersize,
            fillstyle = 'full')
        else:
            pass

right_x = pos_x + (boxes_x-1)*spacing_x*len_x + len_x - 1
bottom_y = pos_y - (boxes_y-1)*spacing_y*len_y - len_y + 1

#Generate legend
leg_x = right_x - len_x + 0.5
leg_y = bottom_y + len_y + 1
leg_space = 2
text_offset = 1
ax.plot(leg_x,leg_y,'o',markersize = 5,color = cmap(0),fillstyle='left',markeredgewidth=0)
ax.plot(leg_x,leg_y,'o',markersize = 5,color='b',fillstyle='right',markeredgewidth=0)
ax.plot(leg_x, leg_y-leg_space, 'o', markersize = 5, color=[0.7, 0.7, 0.7], markeredgewidth=0)
ax.plot(leg_x, leg_y-2.5*leg_space, 'o', markersize = 0.05*def_size, color='k')
ax.plot(leg_x, leg_y-3.5*leg_space, 'o', markersize = def_size, color='k')

ax.text(leg_x + text_offset, leg_y, 'ASV 1-33',
        fontsize = 5.5,verticalalignment = 'center')
ax.text(leg_x + text_offset, leg_y-leg_space, 'Other',
        fontsize = 5.5,verticalalignment = 'center')
ax.text(leg_x + text_offset, leg_y-2.5*leg_space, '$= 10^{-5}$',
        fontsize = 5.5,verticalalignment = 'center')
ax.text(leg_x + text_offset, leg_y-3.5*leg_space, '$= 10^0 $',
        fontsize = 5.5,verticalalignment = 'center')
ax.text(leg_x-1, leg_y-4.5*leg_space, '(Rel. abundance)',
        fontsize = 4.5,verticalalignment = 'center')

ax.set_ylim([bottom_y-fig_pad_y, pos_y+3*fig_pad_y])
ax.set_xlim([pos_x -fig_pad_x, right_x + fig_pad_x ])
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')
ax.set_aspect('equal')
fig.savefig('PD_ASV_level_figure.eps',dpi=800,bbox_inches = 'tight',format='eps',pad_inches=0)

