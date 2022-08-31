#! bin/usr/env python
#
#A plug in to do cross correlation analysis for each trial
# Written Dave Jacobs
#Does cross correlatoin for food and action epoch in junchols task. culd do it for the entire food and action epoch too but it is commented out currently
#
#	ASSUMES 41 HZ
#
import pandas as pd 
import numpy as np 
import os as os

#CROSS CORR FUNCTION ASSUMES 41 HTZ
#returns basic metrics [0] as well as a full correlation as a function of lag [1]
def crosscorr(st1,st2):
	cov=np.correlate(st1-st1.mean(),st2-st2.mean(),'full')   
	cor=cov/(len(st1)*st1.std()*st2.std())
	lag=np.arange(-len(st1)+1,len(st1))
	nolag=np.where(lag==0)[0][0]
	if abs(min(cor)) > abs(max(cor)):
		absmax=min(cor)
	else:
		absmax=max(cor)
	corrout=[lag[np.argmax(cor)],lag[np.argmax(cor)]/41, cor[np.argmax(cor)],lag[np.argmin(cor)],lag[np.argmin(cor)]/41, cor[np.argmin(cor)],cor[nolag],absmax]
	return corrout,cor.tolist()


#takes in the fial with z scored peri-event traces for each trial for each event. By default this raw file must alternate between regions 1 and region 1 for this particular script to run
Bothregions=pd.read_csv('AllActionTrials.csv',index_col=0,dtype=object).transpose()

#No filtering of excludable subjects is done here. That will be done later in the R markdown script that get the mean, sem, and n for the CCF by block (see CrossCorr_getavgs.Rmd)


trialstoget=np.arange(0,Bothregions.shape[0],2)

#open lists for data storing
Masterlist=[]
MasterlistFood=[]
MasterlistAlltrial=[]
Masterfx=[]
Masterfoodfx=[]

#does this for each column pair (mpfc-vta pair for a given trial)
for number in trialstoget:
	indexer=int(number)
	out=float(Bothregions.iloc[indexer][11])+float(Bothregions.iloc[indexer+1][11])
	Subject=Bothregions.iloc[indexer][0:11].tolist()+[out] #the out check for traces that had outliers (z score > 40)
	#will print an error if its not finding the right pair
	if Bothregions.iloc[indexer][7] == 'mPFC' and Bothregions.iloc[indexer+1][7] == 'VTA':
		pass
	else:
		print ('error: unmatched trial',Bothregions.iloc[indexer][7],Bothregions.iloc[indexer+1][7])
	# will print error if trials arent lining up
	if Bothregions.iloc[indexer][0] == Bothregions.iloc[indexer+1][0]:
		pass
	else:
		print ('error: skipped trial',Bothregions.iloc[indexer][0],Bothregions.iloc[indexer][0])


# set indexes based on where you want to look. 135:219 for example grabe +/- 1 sec around action 
	mPFCstream=Bothregions.iloc[indexer][135:219].astype(float).to_numpy()#135,219 for no shock but will depend on how much data you pulled around the event
	VTAstream=Bothregions.iloc[indexer+1][135:219].astype(float).to_numpy()#186:239 for shock
	FoodmPFCstream=Bothregions.iloc[indexer][220:].astype(float).to_numpy()
	FoodVTAstream=Bothregions.iloc[indexer+1][220:].astype(float).to_numpy()
	#AllmPFCstream=Bothregions.iloc[indexer][54:].astype(float).to_numpy()
	#AllVTAstream=Bothregions.iloc[indexer+1][54:].astype(float).to_numpy()
	Actioncor=crosscorr(mPFCstream,VTAstream)
	Actionfx=Actioncor[1]
	Actioncor=Actioncor[0]
	Foodcor=crosscorr(FoodmPFCstream,FoodVTAstream)[0]#1
	Foodfx=crosscorr(FoodmPFCstream,FoodVTAstream)[1]
	#Allcor=crosscorr(AllmPFCstream,AllVTAstream)[0]
	Finalout=Subject+Actioncor
	Fxout=Subject+Actionfx
	Finalfoodout=Subject+Foodcor
	FoodFxout=Subject+Foodfx
	#Finalallout=Subject+Allcor 
	Masterlist.append(Finalout)
	MasterlistFood.append(Finalfoodout)
	#MasterlistAlltrial.append(Finalallout)
	Masterfx.append(Fxout)
	Masterfoodfx.append(FoodFxout)

# write out summary numbers for each trial
Actdf=pd.DataFrame(Masterlist).sort_values(by=[5,10,8])
Actdf.columns=['Trial#','Food Norm Lat','Food Latency','Norm Lat','Latency','Session','Subject','R1','Block','Post-Shock','Shock','Outlier','Pos Lag','Pos Lag(sec)','r-value','Neg Lag','Neg Lag(sec)','Neg r-value','corr@nolag','Absolute Max Corr']
Actdf.transpose().drop(['R1']).transpose().to_csv('Action-CrossCorr.csv')#addheader=None
Fooddf=pd.DataFrame(MasterlistFood).sort_values(by=[5,10,8])
Fooddf.columns=['Trial#','Food Norm Lat','Food Latency','Norm Lat','Latency','Session','Subject','R1','Block','Post-Shock','Shock','Outlier','Pos Lag','Pos Lag(sec)','r-value','Neg Lag','Neg Lag(sec)','Neg r-value','corr@nolag','Absolute Max Corr']
Fooddf.transpose().drop(['R1']).transpose().to_csv('Food-CrossCorr.csv')#,header=None)
# AllTrialdf=pd.DataFrame(MasterlistAlltrial).sort_values(by=[5,9,8])
# AllTrialdf.columns=['Trial#','Food Norm Lat','Food Latency','Norm Lat','Latency','Session','Subject','R1','Block','Post-Shock','Shock','Outlier','Pos Lag','Pos Lag(sec)','r-value','Neg Lag','Neg Lag(sec)','Neg r-value','corr@nolag','Absolute Max Corr']
# AllTrialdf.transpose().drop(['R1']).to_csv('AllTrial-CrossCorr.csv',header=None)

# write out full cross correlation functions - read these into R code
Actiontrialsdf=pd.DataFrame(Masterfx).sort_values(by=[5,6,7])
Actiontrialsdf.rename(columns={0:'Trials',1:'FoodNormLatency',2:'FoodLatency',3:'NormLatency',4:'Latency',5:'Session',6:'Subject',7:'Region',8:'Block',9:'postshock',10:'shock',11:'outliertest'})
Actiontrialsdf.to_csv('Action-FullCrossCorr.csv',header=None)


Foodtrialsdf=pd.DataFrame(Masterfoodfx).sort_values(by=[5,6,7])
Foodtrialsdf.rename(columns={0:'Trials',1:'FoodNormLatency',2:'FoodLatency',3:'NormLatency',4:'Latency',5:'Session',6:'Subject',7:'Region',8:'Block',9:'postshock',10:'shock',11:'outliertest'})
Foodtrialsdf.to_csv('Reward-FullCrossCorr.csv',header=None)



