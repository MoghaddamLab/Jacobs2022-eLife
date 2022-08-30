#! bin/usr/env python 
#
#
#Written by dave jacobs
#uses numpy to calculate AUC values based on z score averages per subject
#
#

import pandas as pd
import numpy as np
import csv
import os

#set directory where the zs-cored averages are
os.chdir('3hz_learn')
# pick region
region=input('Region?:')
# get the file. NOTE: change here for different file naming structers
file=region+'Action_allthedata.csv'


# read file
a=pd.read_csv(file,header=None)

#get all the headers (i.e. figure out how many columns need anlysis)
headers=a.columns[1:]
#a place to send data after analysis
masterlist=[]

# function to go through a column and pull the relevant data and get AUC. the first 4  rows are headers. the second row is if shock happened or not (used to adjust auc windows). Adjust these as needed
def getAUC(column,df=a):
	subject=df[column][3]#2
	session=df[column][2]#2
	block=df[column][0]#0
	shock=df[column][1]#1
#for a subject get the pre event, post event (action) AUC and their change score. aldo get the food change from after action execution
#these windows correapond to where the event took place on the zscore file (i.e. row 170) and assume 41 hz sampling. so 129:170 is 41 samples (1-sec) before action. these can be adjusted as needed
	if shock == '0' or shock =='0.0':
		preaction=np.trapz(df[column][129:170].to_numpy(dtype=float))
		postaction=np.trapz(df[column][170:211].to_numpy(dtype=float))
		changescore=postaction-preaction
		food=np.trapz(df[column][211:293].to_numpy(dtype=float))-postaction
		return (block,shock,session,subject,preaction,postaction,changescore,food)
	elif shock == '1' or shock =='1.0' or shock == 'shnorisk':
		print('sh')
		preaction=np.trapz(df[column][149:191].to_numpy(dtype=float))
		postaction=np.trapz(df[column][191:232].to_numpy(dtype=float))#shockoff + 1 sec
		changescore=postaction-preaction
		food=np.trapz(df[column][211:293].to_numpy(dtype=float))-postaction
		return (block,shock,session,subject,preaction,postaction,changescore,food)

# for each column get the auc dta
for name in headers:
	masterlist.append(getAUC(name))

#save AUC data as csv
AUCdf = pd.DataFrame.from_records(masterlist, columns =['Block','Shock', 'Session','Subject', 'preaction','postaction','actiondelta','food'])
os.chdir('../')
AUCdf.sort_values(by=['Shock','Session','Block','Subject']).transpose().to_csv(region+'AUClearnbyblock.csv',header=None)



