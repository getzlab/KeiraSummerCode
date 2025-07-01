import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.colors as mpcolor
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import math

def treatment_process_centered(filename, filename2):

    df = pd.read_csv(filename, sep=',')
    df2=pd.read_csv(filename2, sep='\t')
    df=pd.merge(df,df2, on='participant_id', how='inner')
    print(df)
    patient_val = df['participant_id'].tolist()
    treatment_val = df['categories'].tolist()
    tx_start = df['start_date_dfd'].tolist()
    tx_end = df['stop_date_dfd'].tolist()
    stop_reason=df['stop_reason'].tolist()
    tumor_morph=df['tumor_morphology'].tolist()
    drug_type=df['drugs'].tolist()
    
    new_dict = {'patient': patient_val,
               'treatment': treatment_val,
               'tx_start': tx_start,
               'tx_end': tx_end,
               'stop':stop_reason,
               'morph':tumor_morph,
               'drug':drug_type}
    
    return new_dict

    
def generate_patientYpos(patient_DF):
    """
    Return a dictionary with unique patients as keys, and 
    y positions as values, for plotting. Assumes the column with 
    the patient names is called "patient"
    """
    
    patient_ids=list(patient_DF['patient'])
    y_pos={}
    new_pos=0000
    
    for patient_name in patient_ids:
        if patient_name not in y_pos:
            new_pos=new_pos+1
            y_pos[patient_name]=new_pos
    
    return y_pos
treats=pd.DataFrame(treatment_process_centered('Juric_Rapid_Autopsy_MASTER-treatments1.csv', 'Juric_Rapid_Autopsy_MASTER-participants.txt'))
disease_name = input("Enter the tumor morphology to filter (e.g., CHOLANGIOCARCINOMA): ")
treats= treats[treats['morph'].str.contains(disease_name, case=False, na=False)]
treats['tx_end'] = pd.to_numeric(treats['tx_end'], errors='coerce')
xmax = math.ceil(treats['tx_end'].max())



patient_pos=generate_patientYpos(treats)
drugtype = treats['drug'].unique()
patients = treats['patient'].unique()

cmap = plt.cm.Paired
colors = [cmap(i) for i in np.linspace(0, 1, len(drugtype))]
colorpal=dict(zip(drugtype,colors))
print(colorpal)

plt.rcParams.update({'font.size': 14})
plt.rcParams['svg.fonttype'] = 'none'
plt.figure(figsize=(15,8))

#f=plt.subplot(5,5,1)
ax=plt.gca()

bottom_length=0
top_length = max(patient_pos.values()) + 1
plt.title('Drug Progression in Rapid Autopsy Data of '+ disease_name + ' Patients')
plt.xlabel('Time (Days)')
plt.ylabel('Patient')

for idx, row in treats.iterrows():
    y_pos=patient_pos[row["patient"]]
    #y_pos=1.0
    ax.hlines(y_pos, 
              float(row['tx_start']), 
              float(row['tx_end']), 
              color=colorpal[row['drug']],
              linewidths=6, 
              zorder=-1, 
              label=str(row['drug']))
    
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
for i in range(0, xmax, n):
    counter=counter+1

    #Ignore every other chromosome for the patches
    if counter%2==0:continue

    newpatch=mpatches.Rectangle([i, bottom_length], n, top_length-bottom_length, ec="none",color='gray', alpha=0.1)
    ax.add_patch(newpatch)

ax.spines['top'].set_visible(True)
ax.spines['right'].set_visible(True)
ax.spines['left'].set_visible(True)    

ax.set_xlim(0, xmax)
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

plt.savefig('carcinoma', dpi=300, bbox_inches='tight')

