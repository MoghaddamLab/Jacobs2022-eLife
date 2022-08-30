#
#Writted by Dave Jacobs
#
#Code to get behavioral events from coulbourn files. Will get basic behavior data (trials, RTs) and get timestamps for each event. 
#Because the clock can get messed up after long periods in grahic state 3. these timestamps are interfaced the arduino timestamps for the feeder to adjust things when using for fiber photmetry peri-event analysis




import math
import numpy 
import csv
from glob import glob
from numpy import mean
from numpy import var
import re
import os as os

#THE GROUND TRUTH IS THE ONSET OF THE FIRST SEEK LIGHT

# gets all the behavior files which are tab delimeted txts in graphic state 3. run this file from the folder the raw data are stored in 
a=glob('*.txt')



t=0
n=0
csvseek=[]
csvpell=[]
csvfood=[]
csvseekvar=[]
csvfirstshocks=[]
csvfirsttrials=[]
csvIBIpokes=[]
Masteractionlat=[]
stds=[]
masterlist=[]

# for a given file.....
for row in a:
	#gets the data
	data=[]
	with open(row,'r') as csvfile:
		reader=csv.reader(csvfile,delimiter='\t')
		for line in reader:
			data.append(line)
		csvfile.close()
	data2=numpy.array(data)


	#makes some names based on the subject
	subj=data2[1][6].split('_')[0][1:]
	Subj=data2[1][6]+'_Stamps.csv'
	Subj2=data2[1][6]+'_Trials.csv'
	stop=len(data2)	

	# alot of open lists for storing data. Not all of these are used...
	IBIseek=0
	IBItake=0
	IBIfood=0
	DispBeh=0
	SeekActions=0
	FoodActions=0
	RRetrievals=0
	Cues=0
	pshock=0
	DispStamps=[]
	PMStamps=[]
	SeekStamps=[]
	FoodStamps=[]
	RetrievalStamps=[]
	CueStamps=[]
	BlockType=[]
	FBlockType=[]
	RBlockType=[]
	CueBlockType=[]
	OutcomeType=[]
	actionlatencies=[]
	foodlatencies=[]
	b1lat=[]
	fb1lat=[]
	postshocks=[]
	FoodDispBeh=0

	# these are the states based on graphic state that reflect given events and blocks
	IBIstates=['9','19','29']#['1','6','20'] BTW this is ITI states
	Seekstates=['3','13','23']
	seekactions=['4','20','30']
	pelletdrops=['5','15','25']
	ITIstarts=['4','12']
	CueOnsets=['2','12','22']
	shocks=['17','27']
	Outcomes=['5','14','17','24','27']
	ITIstates=[]

#uses a while loop to iterate through data. Was written when I didnt understand for looops correctly and didn't realize they are superioir but the code works
#Gets all the stamps and block type for each trial. normalizes times to the first cue onset
	while t+1<=stop:        
	    if data2[n][9] == '2' and data2[n][8] == '1':
	        GroundTruth=float(data2[n][7])*20*.001
	        Stamp=(float(data2[n][7])*20*.001)-GroundTruth
	        Start=(float(data2[n][7])*20*.001)-GroundTruth
	        Cues=Cues+1
	        CueStamps.append(Stamp)
	        CueBlockType.append(float(data2[n][9]))     
	        n=n+1
	        t=t+1
	    elif data2[n][9] in seekactions:
	    	SeekActions=SeekActions+1
	    	Stamp=(float(data2[n][7])*20*.001)-GroundTruth
	    	latency=Stamp-Start
	    	SeekStamps.append(Stamp)
	    	BlockType.append(float(data2[n][9]))
	    	postshocks.append(pshock)
	    	if float(data2[n][9]) == 4.0:
	    		actionlatencies.append(latency)
	    		b1lat.append(latency)
	    	else:
	    		actionlatencies.append(latency)
	    	n=n+1
	    	t=t+1
	    
	    elif data2[n][9] in Outcomes and data2[n][8] in seekactions:
	    	if data2[n][9] in shocks:
	    		OutcomeType.append(1)
	    		pshock=1
	    	elif data2[n][9] == '5':
	    		FoodActions=FoodActions+1
	    		Stamp=(float(data2[n][7])*20*.001)-GroundTruth
	    		fStart=(float(data2[n][7])*20*.001)-GroundTruth
	    		ftrig=1
	    		FoodStamps.append(Stamp)
	    		FBlockType.append(float(data2[n][9]))
	    		OutcomeType.append(0)
	    		pshock=0
	    	else:
	    		OutcomeType.append(0)
	    		pshock=0
	    	n=n+1
	    	t=t+1

	    elif data2[n][9] in pelletdrops:
	    	FoodActions=FoodActions+1
	    	Stamp=(float(data2[n][7])*20*.001)-GroundTruth
	    	fStart=(float(data2[n][7])*20*.001)-GroundTruth
	    	ftrig=0
	    	FoodStamps.append(Stamp)
	    	FBlockType.append(float(data2[n][9]))
	    	n=n+1
	    	t=t+1

	    elif data2[n][10] in ITIstarts:
	    	RRetrievals=RRetrievals+1
	    	Stamp=(float(data2[n][7])*20*.001)-GroundTruth
	    	flatency=Stamp-fStart
	    	if ftrig == 1:
	    		fb1lat.append(flatency)
	    		foodlatencies.append(flatency)
	    	else:
	    		foodlatencies.append(flatency)
	    	RetrievalStamps.append(Stamp)
	    	RBlockType.append(float(data2[n][8]))
	    	n=n+1
	    	t=t+1
	    elif data2[n][9] in CueOnsets:
	    	Cues=Cues+1
	    	Stamp=(float(data2[n][7])*20*.001)-GroundTruth
	    	Start=(float(data2[n][7])*20*.001)-GroundTruth
	    	CueStamps.append(Stamp)
	    	CueBlockType.append(float(data2[n][9]))
	    	n=n+1
	    	t=t+1
	    else:
	        t=t+1
	        n=n+1	        
	n=0
	t=0



	ForSeekcsv=list(zip(SeekStamps,BlockType,OutcomeType,actionlatencies2,actionlatencies,postshocks))
	ForFoodcsv=list(zip(FoodStamps,FBlockType))
	ForRetrievalcsv=list(zip(RetrievalStamps,RBlockType,foodlatencies2,foodlatencies))
	ForCuecsv=list(zip(CueStamps,CueBlockType))
	stds.append([subj,numpy.mean(b1lat),numpy.std(b1lat),numpy.mean(fb1lat),numpy.std(fb1lat)])

	SeekArray=numpy.array(ForSeekcsv)
	FoodArray=numpy.array(ForFoodcsv)
	RetrievalArray=numpy.array(ForFoodcsv)

	filestart=row.split('.')[0]+'_'

#note this code expects a folder name 'ProcessedBeh' to drop the analyzed data into
	os.chdir('../ProcessedBeh')

#save the data. These files need to be used with main photmetry script for it to run and parse everything properly
	with open(filestart+'Seeks.csv', 'w') as file:
		wr=csv.writer(file,lineterminator='\n')
		wr.writerow(['Times','BlockType','Shock?','Latencies','Znorm latency','PostShock?'])
		for item in ForSeekcsv:
			wr.writerow(item)

	with open(filestart+'Food.csv', 'w') as file:
		wr=csv.writer(file,lineterminator='\n')
		wr.writerow(['Times','BlockType'])
		for item in ForFoodcsv:
			wr.writerow(item)
	
	with open(filestart+'Retreivals.csv', 'w') as file:
		wr=csv.writer(file,lineterminator='\n')
		wr.writerow(['Times','BlockType','Latencies','Znorm latency'])
		for item in ForRetrievalcsv:
			wr.writerow(item)

	with open(filestart+'Cues.csv', 'w') as file:
		wr=csv.writer(file,lineterminator='\n')
		wr.writerow(['Times','BlockType'])
		for item in ForCuecsv:
			wr.writerow(item)
	os.chdir('../RawDAta')

######################################################

# code here will get the average per block per subject rather than trial by trial data as above
	mark=0
	block1counter=0
	block2counter=0
	block3counter=0
	fblock1counter=0
	fblock2counter=0
	fblock3counter=0
	rblock1counter=0
	rblock2counter=0
	rblock3counter=0
	seek1mean='Nan'
	seek2mean='Nan'
	seek3mean='Nan'
	food1mean='Nan'
	food2mean='Nan'
	food3mean='Nan'
	seek1first='Nan'
	seek2first='Nan'
	seek3first='Nan'
	seek1list=[]
	seek2list=[]
	seek3list=[]
	food1list=[]
	food2list=[]
	food3list=[]

# count number of trials, food drops, food entires
	for sample in range(len(ForSeekcsv)):
		if ForSeekcsv[sample][1] == 4.0:
			block1counter=block1counter+1
		elif ForSeekcsv[sample][1] == 20.0:
			block2counter=block2counter+1
		elif ForSeekcsv[sample][1] == 30.0:
			block3counter=block3counter+1


	for sample in range(len(ForFoodcsv)):
		if ForFoodcsv[sample][1] == 5.0:
			fblock1counter=fblock1counter+1
		elif ForFoodcsv[sample][1] == 15.0:
			fblock2counter=fblock2counter+1
		elif ForFoodcsv[sample][1] == 25.0:
			fblock3counter=fblock3counter+1

	for sample in range(len(ForRetrievalcsv)):
		if ForRetrievalcsv[sample][1] == 6.0:
			rblock1counter=rblock1counter+1
		elif ForRetrievalcsv[sample][1] == 16.0:
			rblock2counter=rblock2counter+1
		elif ForRetrievalcsv[sample][1] == 26.0:
			rblock3counter=rblock3counter+1


##Action latencies
	for sample in range(len(data2)):
		if data2[sample][9] == '3' and data2[sample][8] == '2':
			mark=1
			start=float(data2[sample][7])*20*.001
		
		elif data2[sample][9] == '4' and mark == 1:
			end=float(data2[sample][7])*20*.001
			latency=end-start
			seek1list.append(latency)
			seek1mean=mean(seek1list)
			seek1first=seek1list[0]
			mark=0
		elif data2[sample][9] == '13' and data2[sample][8] == '12':
			start=float(data2[sample][7])*20*.001
			mark=1
		elif data2[sample][9] == '20' and mark == 1:
			end=float(data2[sample][7])*20*.001
			latency=end-start
			seek2list.append(latency)
			seek2mean=mean(seek2list)
			seek2first=seek2list[0]
			mark=0
		elif data2[sample][9] == '23' and data2[sample][8] == '22':
			start=float(data2[sample][7])*20*.001
			mark=1
		elif data2[sample][9] == '30' and mark == 1:
			end=float(data2[sample][7])*20*.001
			latency=end-start
			seek3list.append(latency)
			seek3mean=mean(seek3list)
			seek3first=seek3list[0]
			mark=0
## Reward retrival latencies
	for sample in range(len(data2)):
		if data2[sample][9] == '5':
			mark=1
			start=float(data2[sample][7])*20*.001
		
		elif data2[sample][9] == '9' and mark == 1:
			end=float(data2[sample][7])*20*.001
			latency=end-start
			food1list.append(latency)
			food1mean=mean(food1list)
			mark=0
		elif data2[sample][9] == '15':
			start=float(data2[sample][7])*20*.001
			mark=1
		elif data2[sample][9] == '19' and mark == 1:
			end=float(data2[sample][7])*20*.001
			latency=end-start
			food2list.append(latency)
			food2mean=mean(food2list)
			mark=0
		elif data2[sample][9] == '25':
			start=float(data2[sample][7])*20*.001
			mark=1
		elif data2[sample][9] == '29' and mark == 1:
			end=float(data2[sample][7])*20*.001
			latency=end-start
			food3list.append(latency)
			food3mean=mean(food3list)
			mark=0

# save the subject averages
	subjs=[subj,subj,subj]
	Blocks=['0%','6%','10%']
	seektrials=[block1counter,block2counter,block3counter]
	foodtrials=[fblock1counter,fblock2counter,fblock3counter]
	Retrievals=[rblock1counter,rblock2counter,rblock3counter]

	seekmeans=[seek1mean,seek2mean,seek3mean]
	foodmeans=[food1mean,food2mean,food3mean]
#aggregate subject averages
	csvmetrics=list(zip(subjs,Blocks,seektrials,Retrievals,seekmeans,foodmeans))
	masterlist.append([subj,filestart.split('_')[1],block1counter,block2counter,block3counter,seek1mean,seek2mean,seek3mean,food1mean,food2mean,food3mean,seek1first,seek2first,seek3first])
#Write out subject averages
	with open(filestart+'Behavior.csv', 'w') as file:
		wr=csv.writer(file,lineterminator='\n')
		wr.writerow(['subject','Risk(%)','Seek Trials','Food Retrievals', 'Seek Latency','Reward Latency'])
		for item in csvmetrics:
			wr.writerow(item)

	with open('allBehavior2.csv', 'w') as file:
		wr=csv.writer(file,lineterminator='\n')
		wr.writerow(['subject','session','TrialsB1','TrialsB2','TrialsB3','ALatB1','ALatB2','ALatB3','f1','f2','f3','first1','first2','first3'])
		for item in masterlist:
			wr.writerow(item)


######################################################
