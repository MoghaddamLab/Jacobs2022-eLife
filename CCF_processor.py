#
#
#written by dave jacobs
#
#Quick code to take the cross correlation functions (averaged) and find the min max and location of max
#


import pandas as pd
import  numpy as np


# takes in an average cross correlation function file which has error (SD or SEM) and sample (n trials) in the next two columns respectivelt
DF=pd.read_csv('Food_crosscorr_forpy.csv',header=None)
data=[]
new=DF.transpose().to_numpy()
# find columns labeled 'mean' i.e. where the trace is
places=np.where(new == 'mean')[0]
#action or food data, since food data will have a different locatio nof the no-lag samples
ftype=input('action or food?:'  )

if ftype == 'action':
	zerospot=84
elif ftype== 'food':
	zerospot=163
else:
	raise ('value error')

#for each trace fine the minimum, maximum, the true peak (highest absolute value for max or min peak - used in paper), SD, n, and location of the true max. 
for a in places:
	maxima=max(new[a][3:].astype(np.float))
	minima=min(new[a][3:].astype(np.float))
	IDs=new[a][0:2].tolist()
	if abs(maxima) >= abs(minima):
		truemax=maxima
	else:
		truemax=minima
	location=np.where(new[a][3:].astype(np.float)==truemax)[0][0]
	Time=(location-zerospot)/41
	SD=float(new[a+1][location])
	N=float(new[a+2][location])
	data.append(IDs+[maxima,minima,truemax,SD,N,Time])

#return data as csv
out=pd.DataFrame(data)
out.columns=['session','risk','Max','min','True Max','SD','N','Time']
out.to_csv('CCFMetrics_'+ftype+'.csv')
   
