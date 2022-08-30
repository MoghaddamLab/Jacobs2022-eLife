# written by Dave JAacobs, Moghaddam Lab

#Code to normalize ephys and photmetry data together and take cross correlation of the two
# 
#Take in traces of mean data in columns
#return correlation scores and cross correlation metircs as well as correlation scores for randmoly shuffled photometry data


import numpy as np
import scipy
from scipy import signal
import csv
import pandas as pd




# get covariance and nromalize to get r value for all time lags. Uses zero-padding
def crosscorr(st1,st2,samplerate=41):
	cov=np.correlate(st1-st1.mean(),st2-st2.mean(),'full')   
	cor=cov/(len(st1)*st1.std()*st2.std())
	lag=np.arange(-len(st1)+1,len(st1))
	nolag=np.where(lag==0)[0][0]
	if abs(min(cor)) > abs(max(cor)):
		absmax=min(cor)
	else:
		absmax=max(cor)
	corrout=[lag[np.argmax(cor)],lag[np.argmax(cor)]/samplerate, cor[np.argmax(cor)],lag[np.argmin(cor)],lag[np.argmin(cor)]/samplerate, cor[np.argmin(cor)],cor[nolag],absmax]
	return (corrout,cor) #,cor.tolist()

# lowpass filter function
def butter_lowpass_filter(data, cutoff=3, fs=63.5, order=2):
	nyq=0.5*fs
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = signal.butter(order, normal_cutoff, btype='low', analog=False)
    y = signal.filtfilt(b, a, data)
    return y

# min-max normalize function
def NormalizeData(data,fullset=None):
    return (data - np.min(data)) / (np.max(data) - np.min(data))


#read in data, add to numpy array
data=[]
with open('ephyscorrs.csv','r') as csvfile:
	reader=csv.reader(csvfile)
	for line in reader:
		data.append(line)
	csvfile.close()

data2=np.array(data)
stop=len(data2)	

#get length of data and set up dataframes fror lowpass filtered ephys data and correlation results
length=np.shape(data2)[1]
resamples=pd.DataFrame()
results=[['epoch','r val','p val','r val (shuff)','p val (shuff)']]
crosscorrdf=[]
shufflecrosscorrdf=[]
lowpassedephys=pd.DataFrame()

#get every other sample (i.e the ephys data)
for a in np.arange(0,24,2):

	if a < length-1: #and np.count_nonzero(data2[:,a][1:]) > 50:
		#ephys data
		s1=data2[:,a][1:]
		s1=s1[s1!=''].astype('float64')
		# lowpass ephys data
		s1=butter_lowpass_filter(s1, cutoff=3, fs=20, order=2)
		#min-max normalizeephys
		s1=NormalizeData(s1)
		
		# for the shorter periods of ephys data. adds in blank data so we can combine cue and action reward data into one dataframe in the end that isnt 'ragged'
		if len(s1) <50:
			resamples[data2[0][a]+'ephys']=s1.tolist()+['']*40
		else:
			resamples[data2[0][a]+'ephys']=s1


		resamp=len(s1)

		#Photometry data
		s2=data2[:,a+1][1:]
		s2=s2[s2!=''].astype('float64')
		#downsample to match ephys sample rate
		s2=scipy.signal.resample(s2,resamp)
		s2=NormalizeData(s2)
		s2_shuffle=s2.copy()
		np.random.shuffle(s2_shuffle)



		#r vals for zero-lag between phys and photmetry
		rval=scipy.stats.pearsonr(s1,s2)[0]
		pval=scipy.stats.pearsonr(s1,s2)[1]
		#r vals for zero-lag between phys and shuffled photmetry
		shuff_rval=scipy.stats.pearsonr(s1,s2_shuffle)[0]
		shuff_pval=scipy.stats.pearsonr(s1,s2_shuffle)[1]
		#save results
		results.append([data2[0][a],rval,pval,shuff_rval,shuff_pval])

		#see above
		if len(s2) <50:
			resamples[data2[0][a]]=s2.tolist()+['']*40
		else:
			resamples[data2[0][a]]=s2
		#runs cross correlation analysis
		crosscorrdf.append([data2[0][a]]+crosscorr(s1,s2,20)[0])
		shufflecrosscorrdf.append([data2[0][a]]+crosscorr(s1,s2_shuffle,20)[0])

	elif a >= length-1:
		pass 
	else:
		pass


#write out datsets
with open('results_corrs.csv','w') as file:
	wr=csv.writer(file,lineterminator='\n')
	for item in results:
		wr.writerow(item)

with open('results_crosscorrs.csv','w') as file:
	wr=csv.writer(file,lineterminator='\n')
	wr.writerow(['epoch','Pos Lag','Pos Lag(sec)','r-value','Neg Lag','Neg Lag(sec)','Neg r-value','corr@nolag','Absolute Max Corr'])
	for item in crosscorrdf:
		wr.writerow(item)

with open('shufffleresults_crosscorrs.csv','w') as file:
	wr=csv.writer(file,lineterminator='\n')
	wr.writerow(['epoch','Pos Lag','Pos Lag(sec)','r-value','Neg Lag','Neg Lag(sec)','Neg r-value','corr@nolag','Absolute Max Corr'])
	for item in shufflecrosscorrdf:
		wr.writerow(item)

#writes out the traces too for plotting
resamples.to_csv('downsampledfp.csv')
lowpassedephys.to_csv('lowpassephys.csv')