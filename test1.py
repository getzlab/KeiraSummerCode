import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import matplotlib.patches as mpatches
import matplotlib.colors as mpcolor
import matplotlib.cm as cm
from matplotlib.lines import Line2D

def treatment_process_centered(filename):

    df = pd.read_csv(filename, sep=',')
    print(df)
    patient_val = df['participant_id'].tolist()
    treatment_val = df['categories'].tolist()
    tx_start = df['start_date_dfd'].tolist()
    tx_end = df['stop_date_dfd'].tolist()
    stop_reason=df['stop_reason'].tolist()
    
    new_dict = {'patient': patient_val,
               'treatment': treatment_val,
               'tx_start': tx_start,
               'tx_end': tx_end,
               'stop':stop_reason}
    
    return new_dict

def biosp_process_centered(filename):
    df= pd.read_csv(filename, sep='\t')
    print(df)
    collection_date=df['collection_date_dfd'].tolist()
    patient_val=df['participant_id'].tolist()
    type=df['submitted_material_type'].tolist()

    new_dict={'patient':patient_val,
              'specimen':type,
              'coll_date':collection_date}
    
    return new_dict
    
def generate_patientYpos(patient_DF):
    """
    Return a dictionary with unique patients as keys, and 
    y positions as values, for plotting. Assumes the column with 
    the patient names is called "patient"
    """
    
    patient_ids=list(patient_DF['patient'])
    y_pos={}
    new_pos=0
    
    for patient_name in patient_ids:
        if patient_name not in y_pos:
            new_pos=new_pos+0.25
            y_pos[patient_name]=new_pos
    
    return y_pos

treats=pd.DataFrame(treatment_process_centered('Juric_Rapid_Autopsy_MASTER-treatments1.csv'))
treats=treats.head(400) #just for ease while coding

patient_pos=generate_patientYpos(treats)
treatment = treats['treatment'].unique()
patients = treats['patient'].unique()

cmap = plt.cm.Paired
colors = [cmap(i) for i in np.linspace(0, 1, len(treatment))]
colorpal=dict(zip(treatment,colors))
print(colorpal)


treat_len = len(treatment)

plt.rcParams.update({'font.size': 14})
plt.rcParams['svg.fonttype'] = 'none'
plt.figure(figsize=(15,15))

#f=plt.subplot(5,5,1)
ax=plt.gca()

bottom_length=0
top_length=15

plt.title('Drug Progression in Rapid Autopsy Data')
plt.xlabel('Time (Days)')
plt.ylabel('Patient')

for idx, row in treats.iterrows():
    y_pos=patient_pos[row["patient"]]
    #y_pos=1.0
    ax.hlines(y_pos, 
              float(row['tx_start']), 
              float(row['tx_end']), 
              color=colorpal[row['treatment']],
              linewidths=6, 
              zorder=-1, 
              label=str(row['treatment']))
    
for idx, row in treats.iterrows():
    y_pos=patient_pos[row["patient"]]
    if row['stop'] == 'Progression and relapse':
        ax.scatter(float(row['tx_end']), y_pos, marker='.', c='r', s=100, zorder=1, label=row['stop'])
    if row['stop'] == 'Completion of standard course':
        ax.scatter(float(row['tx_end']), y_pos, marker='v', c='k', s=50, zorder=1, label=row['stop'])
    if row['stop'] == 'Toxicity':
        ax.scatter(float(row['tx_end']), y_pos, marker='v', c='#008000', s=50, zorder=1, label=row['stop'])
    if row['stop'] == 'Other':
        ax.scatter(float(row['tx_end']), y_pos, marker='.', c='#808080', zorder=1, label=row['stop'])

#shading for each year
counter=-1
n=365
for i in range(0, 10000, n):
    counter=counter+1

    #Ignore every other chromosome for the patches
    if counter%2==0:continue

    newpatch=mpatches.Rectangle([i, bottom_length], n, top_length-bottom_length, ec="none",color='gray', alpha=0.1)
    ax.add_patch(newpatch)

ax.spines['top'].set_visible(True)
ax.spines['right'].set_visible(True)
ax.spines['left'].set_visible(True)    

ax.set_xlim(0, 10000)
ax.set_ylim(bottom_length, top_length)
ax.set_yticks(list(patient_pos.values()))
ax.set_yticklabels(list(patient_pos.keys()))

ax.vlines(0, bottom_length, top_length, color='black', linestyles='--')

legend_elements=[]

#Add the treatments
for farmac in colorpal:
    new_legend=Line2D([0], [0], color=colorpal[farmac], lw=5, label=farmac)
    legend_elements.append(new_legend)


#Add the stop dots
legend_elements.append(Line2D([0], [0], marker='.', color='w', label= 'Progression and Relapse', markerfacecolor='r', markersize=15))
legend_elements.append(Line2D([0], [0], marker='v', color='w', label='Completion of Course', markerfacecolor='k', markersize=10))
legend_elements.append(Line2D([0], [0], marker='v', color='w', label='Toxicity', markerfacecolor='#008000', markersize=10))
legend_elements.append(Line2D([0], [0], marker='.', color='w', label='Other', markerfacecolor='#808080', markersize=15))


ax.legend(handles=legend_elements, loc=1, borderaxespad=0.5,fontsize=8)
plt.tight_layout()
matplotlib.pyplot.show()