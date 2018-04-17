## Employ pandas for visualizing changes in gene expression across timepoints between two groups.
## NOTE: No statistics is included in this analysis. This is purely for visualization purposes.

from matplotlib import pyplot as plt
import math
import csv
import numpy as np
import os
import scipy as sp
import pandas as pd

df = pd.read_csv('Normalized_counts.csv')

# This function will assign every gene a "trend category" number depending on whether it goes up or down over time
def quartiles(df, one, two, three, four):
    minimum = 300
    groups = []
    for index in range(len(df['sod_change'])):
        row = df.iloc[index]
        diff = abs(row[17] - row[14])
        if diff >= 0.6 and (row[12] >= minimum or row[13] >= minimum or row[15] >= minimum or row[16] >= minimum):
            item = row[-1]
            if item <= one:
                groups.append(1)
            elif item > one and item <= two:
                groups.append(2)
            elif item > two and item <= three:
                groups.append(3)
            elif item > three and item <= four:
                groups.append(4)
            elif item > four:
                groups.append(5)
            else:
                groups.append(0)
        else:
            groups.append(0)
    return(groups)

cell1 = df[['genes', 'cell1_wt_pre1', 'cell1_wt_pre2', 'cell1_wt_pre3', 'cell1_wt_post1', 'cell1_wt_post2', 'cell1_sod_pre1', 'cell1_sod_pre2', 'cell1_sod_pre3', 'cell1_sod_post1', 'cell1_sod_post2', 'cell1_sod_post3']]
cell2 = df[['genes', 'cell2_wt_pre2', 'cell2_wt_pre3', 'cell2_wt_post1', 'cell2_wt_post2', 'cell2_wt_post3',  'cell2_sod_pre1', 'cell2_sod_pre2', 'cell2_sod_pre3', 'cell2_sod_post1', 'cell2_sod_post2', 'cell2_sod_post3']]

cell1['wt_pre_avg'] = (cell1.cell1_wt_pre1 + cell1.cell1_wt_pre2 + cell1.cell1_wt_pre3)/3
cell1['wt_post_avg'] = (cell1.cell1_wt_post1 + cell1.cell1_wt_post2)/2
cell1['wt_change'] = cell1.wt_post_avg/cell1.wt_pre_avg

cell1['sod_pre_avg'] = (cell1.cell1_sod_pre1 + cell1.cell1_sod_pre2 + cell1.cell1_sod_pre3)/3
cell1['sod_post_avg'] = (cell1.cell1_sod_post1 + cell1.cell1_sod_post2 + cell1.cell1_sod_post3)/3
cell1['sod_change'] = cell1.sod_post_avg/cell1.sod_pre_avg
cell1['group'] = quartiles(cell1, 0.5, 0.9, 1.1, 1.4)

cell2['wt_pre_avg'] = (cell2.cell2_wt_pre2 + cell2.cell2_wt_pre3)/2
cell2['wt_post_avg'] = (cell2.cell2_wt_post1 + cell2.cell2_wt_post2 + cell2.cell2_wt_post3)/3
cell2['wt_change'] = cell2.wt_post_avg/cell2.wt_pre_avg

cell2['sod_pre_avg'] = (cell2.cell2_sod_pre1 + cell2.cell2_sod_pre2 + cell2.cell2_sod_pre3)/3
cell2['sod_post_avg'] = (cell2.cell2_sod_post1 + cell2.cell2_sod_post2 + cell2.cell2_sod_post3)/3
cell2['sod_change'] = cell2.sod_post_avg/cell2.sod_pre_avg
cell2['group'] = quartiles(cell2, 0.5, 0.9, 1.1, 1.4)

# Subdivide the genes into their respective "trend" categories
cell1_1 = cell1[cell1.group == 1]
cell1_1.reset_index(inplace=True, drop=True)
print('Cell_1 1: '+ str(len(cell1_1['genes'])))
cell1_2 = cell1[cell1.group == 2]
cell1_2.reset_index(inplace=True, drop=True)
print('Cell_1 2: '+ str(len(cell1_2['genes'])))
cell1_3 = cell1[cell1.group == 3]
cell1_3.reset_index(inplace=True, drop=True)
print('Cell_1 3: '+ str(len(cell1_3['genes'])))
cell1_4 = cell1[cell1.group == 4]
cell1_4.reset_index(inplace=True, drop=True)
print('Cell_1 4: '+ str(len(cell1_4['genes'])))
cell1_5 = cell1[cell1.group == 5]
cell1_5.reset_index(inplace=True, drop=True)
print('Cell_1 5: '+ str(len(cell1_5['genes'])))

cell2_1 = cell2[cell2.group == 1]
cell2_1.reset_index(inplace=True, drop=True)
print('\nCell_2 1: '+ str(len(cell2_1['genes'])))
cell2_2 = cell2[cell2.group == 2]
cell2_2.reset_index(inplace=True, drop=True)
print('Cell_2 2: '+ str(len(cell2_2['genes'])))
cell2_3 = cell2[cell2.group == 3]
cell2_3.reset_index(inplace=True, drop=True)
print('Cell_2 3: '+ str(len(cell2_3['genes'])))
cell2_4 = cell2[cell2.group == 4]
cell2_4.reset_index(inplace=True, drop=True)
print('Cell_2 4: '+ str(len(cell2_4['genes'])))
cell2_5 = cell2[cell2.group == 5]
cell2_5.reset_index(inplace=True, drop=True)
print('Cell_2 5: '+ str(len(cell2_5['genes'])))

# Plot the results:
#print(cell2.head())
plt.close('all')
group = cell2_2
#print(group['genes'])
index = 500
array = group.iloc[index]
timepoints = ['Pre', 'Post']

plt.plot(timepoints, array[12:14], color = 'black')
plt.plot(timepoints, array[15:17], color = 'red')

plt.plot(['Pre', 'Pre', 'Pre'], array[6:9], marker='o', color='red', linewidth=0)
plt.plot(['Post', 'Post', 'Post'], array[9:12], marker='o', color='red', linewidth=0)
plt.plot(['Pre', 'Pre', 'Pre'], array[1:4], marker='o', color='black', linewidth=0)
plt.plot(['Post', 'Post'], array[4:6], marker='o', color='black', linewidth=0)

plt.legend(['WT', 'SOD'])
plt.title(array[0] + ', Group ' + str(array[-1]))
