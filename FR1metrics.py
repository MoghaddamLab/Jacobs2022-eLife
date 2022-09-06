#! bin/usr/env python 
#
#
#Written by Dave Jacobs. Analysis code for basic fr1 data from raw coulbourn graphic state files. Cpould be repurposed for multiple programs by changing function inputs
#
#
#
#basic needed packages
import numpy as np
import csv
from glob import glob




#Function which uses the start and end states for trial begin and action to determine the the number of trials and RTs per trial. will return an average of all trials or a trial by trial list of times. 
#this function is run seperately for action and reward data because it is different start and end state combos. 
#prior is the state before epoch beginning
#shift is a way to adjusrt the time if there is some sort of delay imposed. (e.g. the 1 second delay for food drop means this is set to 1 for reward retrieval calculations)
def getlatencies(dataset,startrt,prior,endrt,shift):
    actionlist=[]
    pells=0
    mark=0
    IApokes=0
    actionmean='Nan'
    for sample in range(len(dataset)):
            if dataset[sample][9] == startrt and data2[sample][8] in prior:
                    mark=1
                    start=float(data2[sample][7])*20*.001
            elif dataset[sample][9] == endrt and mark == 1:
                    end=float(data2[sample][7])*20*.001
                    latency=end-start
                    actionlist.append(latency-shift)
                    actionmean=np.mean(actionlist)
                    pells=pells+1
                    mark=0
            # inactive nosepokes. Not relevant for elife research advance
            elif dataset[sample][8] == startrt and dataset[sample][12] == '1':
                IApokes=IApokes+1
    return (actionlist,actionmean,np.std(actionlist),pells,IApokes)




#main list to add data to
masterlist=[]
#get all raw files from folder. Since this is GS3 all exported files are txts
a=glob('TestData/*.txt')

#for a given file...
for sess in a:
    data=[]
    with open(sess,'r') as csvfile:
    	reader=csv.reader(csvfile,delimiter='\t')
    	for line in reader:
    		data.append(line)
    	csvfile.close()
    data2 =np.array(data)
    # get session and subject name from filename
    session=sess.split('\\')[1].split('.')[0].split('_')[1][1:]
    subject='S'+sess.split('\\')[1].split('.')[0].split('_')[0]
    #this one rat had a different transition numbers. in the end he wasn't included in the main dataset because he didnt get any PRT sessions
    if subject == 'S162':
        ActionRT=getlatencies(data2,'3',['6','3'],'4',0)
        FoodRT=getlatencies(data2,'4',['3'],'6',1)
    #get action and food rts and trial counts
    else:
        ActionRT=getlatencies(data2,'2','1','5',0)
        FoodRT=getlatencies(data2,'5',['2'],'1',1)

    #pulls relevant data based on function outputs 0- list of all trial times, 1 - average RT, 2 - std of RTs, 3 - number of trials (pellets), 4 - inactive pokes (irrelevant here)
    Actionmean=ActionRT[1]
    Actionstd=ActionRT[2]
    Foodmean=FoodRT[1]
    Foodstd=FoodRT[2]
    Pellets=ActionRT[3]
    InactivePokes=ActionRT[4]

    #append means for action and reward and trial counts for each subject and session to main dataframe
    masterlist.append([subject,session,Actionmean,Actionstd,Foodmean,Foodstd,Pellets,InactivePokes])

#write out the dataframe to csv after all  files are analyzed
with open('FR1Metrics.csv', 'w') as file:
    wr=csv.writer(file,lineterminator='\n')
    wr.writerow(['Subject','Session','Action_Mean','Action_std','Food_Mean','Food_Std','Pellets','InactivePokes'])
    for item in masterlist:
        wr.writerow(item)
