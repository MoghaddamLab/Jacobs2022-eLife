# code by Dave Jacobs
#
#Used to take individual subject avwerages and run permutation tests as a function of session or block
#
##take z scored photometry traces around events of interest
#returns p values for each timepoint and a p value list that will have '1' if the threshold is met


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


#this adjusts some things according to if the file was a shock average or not. That said its usually just easier to adjust them manually in lines 80-100 to suit the comparison one is trying to make 
shockorno=input('is shock? (y/n):  ')
if shockorno == 'y':
	shocktrig=1
	shocktrig2=3
	ttypeloc=0
	trialtype=['1','1.0']
else:
	shocktrig=2
	shocktrig2=0
	ttypeloc=1
	trialtype=['0','0.0'] # change to 0

data=pd.read_csv('3hz_dz/VTAAction_allthedata.csv',header=None).transpose()
data=data.drop([0])
#data2=data
data2=data[data[ttypeloc].isin(trialtype)]
#get shock control data (no shock trials) if needed
shctrdata=data[data[ttypeloc].isin(['0','0.0'])]


#exclude subjects here
#the fifth row can be used to exclude data by typing 'EX' however one could also just remove the irrelevent subjects from the dataframe
data2=data2[~data2[4].isin(['EX'])]


# function for permutation testing
def permutation(data,iterations,samps,m1,m2,plots=False,twotail=True):
	mean_diff=[]
	#gets a given timepoint and develops a null distribution of mean differences of shuffled data for n-iterations (I use 1000)
	for i in range(0,iterations):
		datas=np.array(data).astype(float)
		np.random.shuffle(datas)
		array_split= np.split(datas,[0,samps],axis=0)
		a1=array_split[1]
		a2=array_split[2]
		diff=a1.mean()-a2.mean()
		mean_diff.append(diff)
		#get mean
	#calculated the observed mean diff
	observed=m1-m2
	#counts for number of times observed is > the null counts
	counter=0
	if twotail== True: 
		for item in mean_diff:
			if abs(item) >= abs(observed): # make abs() and > for two tailed
				counter=counter+1
# depending on your directionality one may need to adjust the '<' here
	elif twotail== False: 
		for item in mean_diff:
			if item <= observed: # make abs() and > for two tailed
				counter=counter+1

#optional plotting of the null distribution and observed mean diff for sanity checks
	if plots == True:
		plt.clf()
		plt.hist(mean_diff,bins=20)
		plt.axvline(x=observed)
	else:
		pass
#returns p value (e.g. number of times observed > null/total number obs)
	return(counter/iterations)





#####



#list for aggreagating p values for each sample
pmaster=[]

#pull relevent sessions. adjust this so that the control sesssion is NOT in this list
sessions=data2[shocktrig].unique()[:1]#for no shock

#not required
ctrsess='1'

#go through each session and get comparisons 
for a in range(len(sessions)):
	data3=data2[data2[shocktrig]==sessions[a]]#2for no shock
	blocks=data3[shocktrig2].unique()#gets each block for iteration in for loop
	block=blocks#[blocks!='1']   # NOTE one may want to leave out block 1 and include each session of the goal is to compare across blocks per sesssion rather than across session per block. if so the 'shocktrig' value in the for loop on 107 should be adjusted
	for b in block:
		datafin=data3[data3[shocktrig2].isin([b])]#gesta  give session and block
		ctr=data2[data2[shocktrig2].isin([b])] # gets a block for all sessions and will filter out a control session on 104
		#ctr=data2[data2[shocktrig]==sessions[a]]#Additional options for choosing specific sessions or block as the control
		#ctr=ctr[ctr[shocktrig2]=='1']#
		ctr=ctr[ctr[shocktrig]=='Saline']#keep the saline session and block as a comparitor
		datafin=datafin.append(ctr,ignore_index = True)# compare saline to the other block and session
		plist=[]
		for i in range(5,373):#the first 6 rows shold be headers. the upper of the range will determine how far out you want to test. depends on epoch, sampling rate, etc
			forav=datafin.loc[0:,[shocktrig,i]]#split the data by what we are comparing
			forav[i]=forav[i].astype(float) 
			ms=forav.groupby([shocktrig]).mean() #get avgs
			firstav=ms.iloc[0].tolist()[0]
			secav=ms.iloc[1].tolist()[0]
			samples=datafin[shocktrig].value_counts()[0]#determine the sample size for one o fthe groups so we can shuffle the data properly for  generating null
			plist.append(permutation(datafin[i],1000,samples,firstav,secav,plots=False,twotail=True))# run permutationals NOTE two tail is set to true. append the p value to a list
		pmaster.append([sessions[a],b]+plist)# take all the p values and append to a df that has all p values for each timepoint per session or block

#writes out data
pd.DataFrame(pmaster).transpose().rename({0:'Session',1:'BlockCompared'}).to_csv('pvalues.csv',header=None)  


# A SEPERATE df TO DO THE CONSECUTIVE THRESHOLDING BELOW
pdf=pd.DataFrame(pmaster).transpose().rename({0:'Session',1:'BlockCompared'})



# create a run a function that clusters the p values (consecutive thresholding). I.e. this checks if a run of p values was less than .05 for n-consecutive samples. Default is 14 because this is the lowpass filter window for what we used. See Philip-BRessel(2017) for more
threshlist=[]
def clusterps(data, threshold=14,stepsize=1):
	for dset in data.columns:
		sig_spots=np.where(np.array(data[dset][2:].astype(float)<.05))
		p_locs=np.split(sig_spots[0], np.where(np.diff(sig_spots[0]) != stepsize)[0]+1)
		spots=[]
		for row in range(len(p_locs)):
			if len(p_locs[row]) >= threshold:
				spots.append(p_locs[row].tolist()[0:])
			else:
				pass
		base=['Nan']*(len(data)-2)
		for a in range(len(spots)):
			start=spots[a][0]
			end=spots[a][-1]
			base[start:end+1]=[1]*(end+1-start)
		threshlist.append(base)

clusterps(pdf)

# writeout consec threshold results, note no headers so but order will match the pdf dataframe
pd.DataFrame(threshlist).transpose().to_csv('consecutivethresholdlist.csv')

