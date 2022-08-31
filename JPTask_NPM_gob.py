#! bin/usr/env python 
#written by dave Jacobs
#this works with Times_JPTask to get the trial data. Now writes out shock trials and blocks into the pandas dataframe for food iti and fr1 actions 3/17/20
#Main script for processing the Fiber Photometry data. Needs each subject to have a directory with sesssion as a subdirectory and within session will be recording data. behavior data etc
#outputs z scored perievent traces both for average over a subject and session as well as for each trial

# 4/14/21 will autodrop trials where a zscore > 40 happens for any time point. 

import scipy
import csv
import numpy as np
from numpy import mean
from scipy import signal
from statsmodels import robust
import matplotlib.pyplot as plt 
from pylab import *
from glob import glob
import os
import pandas
import shutil
from datetime import datetime
from scipy.stats import pearsonr
from scipy.signal import butter,filtfilt



# open dataframes for storing data
bigAction_mPFC=pandas.DataFrame()
bigAction_VTA=pandas.DataFrame()
bigCue_mPFC=pandas.DataFrame()
bigCue_VTA=pandas.DataFrame()
bigFood_mPFC=pandas.DataFrame()
bigFood_VTA=pandas.DataFrame()
AllActionTrials=pandas.DataFrame()
AllCueTrials=pandas.DataFrame()



#counters for big motion artifacts (zscore>40)
motionissues=0
motionissuesVTA=0
cuemotionissues=0
cuemotionissuesVTA=0
trialcounter=0

#some finctions for dropping a known bad trial andusing lowpass filter
def dropbadtrials (triallist,dframe,frametype):
    if frametype=='tpose':
        return(dframe.drop(columns=triallist))#for tpose
    elif frametype=='trialholder':
        return(dframe[~dframe.Trials.isin(triallist)])

#3 Hz cutoff
def butter_lowpass_filter(data, cutoff=3, fs=41, order=2):
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y
########################################################################################
#get all subjects
Subjects=glob('S*')#
Badtrial={}

for subject in Subjects:
    os.chdir(subject)
    print (subject)
    AVTA= 'Y'#if yes it will assume the VTA data are in the colums to the right of the mPFC columns

    sessions=glob('JPD*') #pull all sessions, note the JP prefix
    for sessarino in sessions:
    #get timestamp data and load files
        os.chdir(sessarino)
        a=glob('corrN*')[0]#nosepoke cue arduino stamp
        b=glob('feedercue*')[0]#feeder arduino stamps
        d=glob('*JPT*')[0]#times of things
        CueStamps=glob('*Cues*')[0]#couldbourne stamps for syncing
        SeekStamps=glob('*Seeks*')[0]
        FoodDrop=glob('*Food*')[0]
        FoodRet=glob('*Retreivals*')[0]

        # some metric things
        spd=sessarino.split('_')[0]
        

        ###########################################################################################3

        print('cuereads',datetime.now().time())
        data=[]


        with open(a,'r') as csvfile:
            reader=csv.reader(csvfile,delimiter=' ')
            for line in reader:
                data.append(line)
            csvfile.close()

        data2 =np.array(data)


        #get first cue as Ground Truth. offset by 1.4 sec due to arduino-bonsai delay
        StartCue=(float(data2[1][0])/1000)-(float(data2[0][0])/1000)+1.4
        offset=StartCue
        ###
        #find the trial specific ground truth 
        #coulbourn stamps
        data3=[]
        with open(FoodDrop,'r') as csvfile:
            reader=csv.reader(csvfile)
            for line in reader:
                data3.append(line)
            csvfile.close()
        data3=np.array(data3)
        Stamps=data3[1:][:,0]
        FoodDropList=[]
        for i in Stamps:
            FoodDropList.append(float(i)+offset)

#ground truth of pellet drop
        data3=[]
        with open(b,'r') as csvfile:
            reader=csv.reader(csvfile,delimiter=' ')
            for line in reader:
                data3.append(line)
            csvfile.close()

        data3=np.array(data3)
        Start=float(data3[0][0])/1000
        Foodtime=data3[1:][:,0]
        Foodtime2=[]
        for x in Foodtime:
            Foodtime2.append(float(x)/1000)
        Foodtime=Foodtime2[0::2]
        foodnorm=[]
        for x in range(len(Foodtime)):
            foodnorm.append((Foodtime[x]-Start)+1.4)
        both=np.array(list(zip(FoodDropList,foodnorm)))
        #find the trial by trial drift times to adjust other trials
        diffs=[]
        for x in range(len(both)):
            diffs.append(both[x][0]-both[x][1])

        #get coulbourn cue stamps
        data3=[]
        with open(CueStamps,'r') as csvfile:
            reader=csv.reader(csvfile)
            for line in reader:
                data3.append(line)
            csvfile.close()
        data3=np.array(data3)

        Stamps=data3[1:][:,0]
        CueList=[]
        for i in Stamps:
            CueList.append(float(i)+offset)

        # nomralize the cues
        allnums=list(zip(CueList,diffs))
        NormCueList=[]
        for num in range(len(allnums)):
            NormCueList.append(allnums[num][0]-allnums[num][1])


        # actions stamps
        datask=[]
        with open(SeekStamps,'r') as csvfile:
            reader=csv.reader(csvfile)
            for line in reader:
                datask.append(line)
            csvfile.close()
        datask=np.array(datask)
        Stamps=datask[1:][:,0]
        SeekList=[]
        for i in Stamps:
            SeekList.append(float(i)+offset)

        allnums=list(zip(SeekList,diffs))
        NormSeekList=[]
        for num in range(len(allnums)):
            NormSeekList.append(allnums[num][0]-allnums[num][1])

        shockorder=np.array(datask[1:][:,2]).tolist()
        aftershock=np.array(datask[1:][:,5]).tolist()
        actlatencies=np.array(datask[1:][:,3]).tolist()
        actlatenciesnorm=np.array(datask[1:][:,4]).tolist()
        blocksdecoded=[]
        blockd=np.array(datask[1:][:,1]).tolist()
        for a in range(len(blockd)):
            if blockd[a] == '4.0':
                blocksdecoded.append(1)
            elif blockd[a] == '20.0':
                blocksdecoded.append(2)
            elif blockd[a] == '30.0':
                blocksdecoded.append(3)
        #pellet stamps
        data3=[]
        with open(FoodDrop,'r') as csvfile:
            reader=csv.reader(csvfile)
            for line in reader:
                data3.append(line)
            csvfile.close()
        data3=np.array(data3)
        Stamps=data3[1:][:,0]
        FoodDropList=[]
        for i in Stamps:
            FoodDropList.append(float(i)+offset)

        allnums=list(zip(FoodDropList,diffs))
        NormFDList=[]
        for num in range(len(allnums)):
            NormFDList.append(allnums[num][0]-allnums[num][1])

        #Pellet Retrieval stamps
        data3=[]
        with open(FoodRet,'r') as csvfile:
            reader=csv.reader(csvfile)
            for line in reader:
                data3.append(line)
            csvfile.close()
        data3=np.array(data3)
        Stamps=data3[1:][:,0]
        RetrievalList=[]
        for i in Stamps:
            RetrievalList.append(float(i)+offset)

        allnums=list(zip(RetrievalList,diffs))
        NormRetList=[]
        for num in range(len(allnums)):
            NormRetList.append(allnums[num][0]-allnums[num][1])
        fdlatencies=np.array(data3[1:][:,2]).tolist()
        fdlatenciesnorm=np.array(data3[1:][:,3]).tolist()
        fdlatencies=fdlatencies+((len(actlatencies)-len(fdlatencies))*['NA']) 
        fdlatenciesnorm=fdlatenciesnorm+((len(actlatenciesnorm)-len(fdlatenciesnorm))*['NA']) 


        ####
        #this changes list names to work with rest of code

        cue=NormCueList
        seekresponse=NormSeekList
        foodcue=NormFDList
        foodresponse=NormRetList
        #################################################################################



        ##################################################
        # time decoder


        #Relevent Functions to convert weird flycapture time to seconds. also float some data
        def converttime(time):
            cycle1 = (int(time) >> 12) &  0x1FFF
            cycle2 = (int(time) >> 25) &  0x7f
            seconds = cycle2 + float(cycle1)/8000
            return seconds

        def rezero(x):
            x=x-cumetime[0]
            return x

        def cam1(x):
            camer1.append(float(x))

        def cam1red(x):
            camerred1.append(float(x))

        def cam2(x):
            camer2.append(float(x))

        def cam2red(x):
            camerred2.append(float(x))

        def floater(x):
            x=float(x)
            return x



        #########################################
        #Start reading in the data


        print('read fpdata',datetime.now().time())

        data=[]

        with open(d,'r') as csvfile:
            reader=csv.reader(csvfile,delimiter=' ')
            for line in reader:
                data.append(line)
            csvfile.close()
        data2=np.array(data)
        timestamps=data2[:,0]
        ints=timestamps.astype(np.int64)
#convert camera time to seconds
        timeconverted=list(map(converttime,ints))


        t=1
        n=0
        time=timeconverted[0]
        cumetime=[time]
#convert camera seconds to cumulative time
        stop=len(timeconverted)

        while t+1 <= stop:

            if (timeconverted[t]-timeconverted[t-1]) > 0:
                time=time + (timeconverted[t]-timeconverted[t-1])
                cumetime.append(time)
                t=t+1
            else:
                n=n+1
                time=timeconverted[t] +(n*128)
                cumetime.append(time)
                t=t+1

        cumetime2=list(map(rezero,cumetime))


        ############################################
        #get relevent camera data and rescale the control channel
        #assumes [1] and [2] are mpfc gcamp,iso  and [3][4] are vta gcamp,iso respectively


        camer1=[]
        camerred1=[]
        camer2=[]
        camerred2=[]


        list(map(cam1,data2[:,1]))
        list(map(cam1red,data2[:,2]))
        if data2[:,3][0] == '':
            pass
        else:    
            list(map(cam2,data2[:,3]))
            list(map(cam2red,data2[:,4]))
        ####rescale the red channel

        y=list(map(float,data2[:,1]))
        x=list(map(float,data2[:,2]))
        # 
        params=np.polyfit(x,y,1)
        # 
        slope=float(params[0])
        intercept=float(params[1])
        # 
        # 
        def rescale(x):
             x = (slope*x)+intercept
             return x
        # 
        scaledred=list(map(rescale,camerred1))


        #lowpass filter
        fs = 41.0       # sample rate, Hz      # desired cutoff frequency of the filter, Hz 
        nyq = 0.5 * fs  # Nyquist Frequency            
        camer1=butter_lowpass_filter(camer1).tolist()
        camerred1=butter_lowpass_filter(camerred1).tolist()
        if data2[:,3][0] == '':
            pass
        else:             
            camer2=butter_lowpass_filter(camer2).tolist()    
            camerred2=butter_lowpass_filter(camerred2).tolist()        




#saves the initial LP filter data
        newcsv=list(zip(cumetime2,camer1,camerred1,scaledred))
        newcsv2=list(zip(cumetime2,camer2,camerred2,scaledred))

        with open('newfile.csv', 'w') as file:
             wr=csv.writer(file,lineterminator='\n')
             wr.writerow(['time','cam1green','cam1red','scaled_red'])
             for item in newcsv:
                 wr.writerow(item)

        with open('newfile2.csv', 'w') as file:
             wr=csv.writer(file,lineterminator='\n')
             wr.writerow(['time','cam2green','cam2red','scaled_red(cam1)'])
             for item in newcsv2:
                 wr.writerow(item)
        ####################################################
        print('start finding data',datetime.now().time())
        #delta f/f


        deltadata=[]
        deltadata2=[]
        print ('Analyzing '+sessarino)

#get LP data in. now we will normalize the control to the gcamp channel for each trial
        with open('newfile.csv','r') as csvfile:
            reader=csv.reader(csvfile,delimiter=',')
            for line in reader:
                deltadata.append(line)
            csvfile.close()

        deltaarray=np.array(deltadata)
        stop=len(deltaarray)


        cam1deltaf=[]
        times2=[]
        green=[]
        red=[]

        #cam2
        cam2deltaf=[]
        green2=[]
        red2=[]

        ################

        ## timestamp seconds to sample #
        def tosample (x):
            x=round(x,2)
            return x
#get the samples needed to parse data i.e trial start and ends
        samplesneeded=list(map(tosample,cue))
        seekrespsamplesneeded=list(map(tosample,seekresponse))
        foodsamplesneeded=list(map(tosample,foodcue))
        foodrespsamplesneeded=list(map(tosample,foodresponse)) ###pellet drop not exit

        ##
#get all the trial samples for linear fit
        deltascorr=[]
        deltareds=[]
        deltasgreen=[]
        times=[]
        dg=[]
        dr=[]
        samples=[]
        t=1
        stop=len(samplesneeded)


        #####
        #get index
        t=0
        n=1
        stop=len(samplesneeded)
        indexed=[]

        while t+1<= stop:
            ti=float(samplesneeded[t])
             
            if ti <= (float(deltaarray[n][0])): #or ti > (float(deltaarray[n][0])):
                frame=n
                indexed.append(frame) 
                n=1
                t=t+1
            else:
                n=n+1

        print('start rescale')
                
                
        ##rescale by trial
        last=float(len(deltaarray))-2
        indexed.append(last)
        #indexed.append(195448) #cheeky add to get last trial if needed

#takes the prior 10 sec ITI up unbtil the start of next trial
        t=1
        stop=len(indexed)
        while t+1 <= stop:
            x=indexed[t-1]
            x=int(x)
            diff=indexed[t]

            t=t+1
            #41 htz
            x=x-410
# gets the gcamp and control data
            while x < diff:
                deltafgreen=float(deltaarray[x][1])#-gmean)/gmean
                dg.append(deltafgreen)
                samples.append(x)
                times.append(float(deltaarray[x][0]))

                deltafred=float(deltaarray[x][2])#-rmean)/rmean
                dr.append(deltafred)
                x=x+1

#get slope and rescale red to green channel;
            params=np.polyfit(dr,dg,1)
            #print (params)
            slope=float(params[0])
            intercept=float(params[1])

            def rescale2(x):
                tran = (slope*x)+intercept  
                return tran

            deltas2=list(map(rescale2,dr))

            for x in dg:
                deltasgreen.append(x)
            for x in deltas2:
                deltareds.append(x)

            
            dg=[]
            dr=[]


#calculate corrected deltaf

        deltascorr=[]
        stop=len(deltareds)

        t=0
        while t+1 <= stop:
            x=(deltasgreen[t]-deltareds[t])/deltareds[t] ##
           # x=(dg[t]-dr[t])/dr[t] ##    
            deltascorr.append(x)
            t=t+1


        print('start seek re',datetime.now().time())
        #########################
        # this holds the rescaled data df/f as well as the individual trace of gcamp and control
        newcsv=list(zip(times,samples,deltasgreen,deltareds,deltascorr))
        array=np.array(newcsv)
        ##################################
        #open lists for storing action data
        seeklist=[]
        seek2list=[]
        seekcuelist=[]
        peaklist=[]
        peakforcsv=[]

        #get index of cue onset and action completion
        t=0
        n=0
        stop=len(samplesneeded)#-1
        indexed=[]

        while t+1<= stop:
            ti=float(samplesneeded[t])
              
            
            if ti <= float(array[n][0]): #or ti > float(array[n][0]):
                frame=n
                indexed.append(frame) 
                n=0
                t=t+1
            else:
                n=n+1
        t=0
        n=0
        stop=len(seekrespsamplesneeded)#-1
        srindexed=[]

        while t+1<= stop:
            ti=seekrespsamplesneeded[t]
             
            if ti <= float(array[n][0]):# and ti <= float(array[n][0]+.02):
                frame=n
                srindexed.append(frame) 
                n=0
                t=t+1
            else:
                n=n+1
                

        print('start action data')
        ######
        #get action datafor each trial
        t=0
        stop=len(srindexed) #srindexed
        while t+1<=stop:
            x=srindexed[t] #srindexed
           # print (x)
            y=indexed[t]
            x=int(x)
            y=int(y)
#NOTE z score baseline is set here. default is -4 to -2 before action
            BL=mean(list((map(float,array[x-164:x-82][:,4]))))
            BLs=np.std(list((map(float,array[x-164:x-82][:,4])))) #[y-4100:y-2][:,3]
            maxchange=float(max(array[x:x+60][:,4]))
            

            #window is 2-sec prior until 3 sec after response
            a=list(map(float,array[x-164:x+205][:,4])) #194
            seeklist.append(list(map(float,a)))

            t=t+1
            

            hold=[]
#get the z score            
            for x in a:
                b=(x-BL)/BLs
                hold.append(b)
            hold=np.array(hold)

            hold=hold.tolist()
# check for outliers
            outliercheck=len(np.where(abs(np.array(hold))>40)[0])  
            if outliercheck >=1:
                motionissues=motionissues+1
                trialcounter=trialcounter+1
            else:
                trialcounter=trialcounter+1 
            seek2list.append([outliercheck,(np.trapz(hold[164:205])-np.trapz(hold[122:164])),np.trapz(hold[205:])]+hold)

        ##
# do the same for the cue onset       
        t=0
        stop=len(indexed) #srindexed
        while t+1<=stop:
            x=indexed[t] #srindexed
           # print (x)
            y=indexed[t]
            x=int(x)
            y=int(y)
            BL=mean(list((map(float,array[x-164:x-82][:,4]))))
            BLs=np.std(list((map(float,array[x-164:x-82][:,4])))) #[y-4100:y-2][:,3]
            maxchange=float(max(array[x:x+41][:,4]))
            

            ##
            
            peakchange=(maxchange)
            
            #window is 2-sec prior until 3 sec after response
            a=list(map(float,array[x-164:x+205][:,4])) #194
            t=t+1
            

            hold=[]
            for x in a:
                b=(x-BL)/BLs
                hold.append(b)
            hold=np.array(hold)

            hold=hold.tolist()
            outliercheck=len(np.where(abs(np.array(hold))>40)[0])   
            if outliercheck >=1:
                cuemotionissues=cuemotionissues+1
            else:
                pass                     
            seekcuelist.append([outliercheck,(np.trapz(hold[164:205])-np.trapz(hold[122:164])),np.trapz(hold[205:])]+hold) 

        print('action over',datetime.now().time())

        ###############################################################################################################3

#makes a plot to look at the overal;l signal over the session
        print('plot',datetime.now().time())
        for a in seekresponse:
            plt.axvline(x=a,color='b')
        plt.plot(times,deltasgreen,color='g')
        plt.plot(times,deltareds,color='r')
        plt.title('mPFC Streams '+spd)
        plt.xlabel('Seconds')
        plt.ylabel('Intensity')
        plt.legend("GT",loc='best')
        plt.savefig('mPFCStreams_'+spd+'.png',dpi=400)
        plt.clf()
        plt.close()



        ####
#makes directory to save the data into for mPFC

        print('startcsv',datetime.now().time())
        name=CueStamps.split('_')
        rnum=str(round(rand(),4))
        name=name[0]+name[1]+'mPFC_Analysis_'

        if os.path.exists(name):
            shutil.rmtree(name)
        os.mkdir(name)
        os.chdir(name)


#save files
        with open('seekres.csv', 'w') as file:
            wr=csv.writer(file,lineterminator='\n')
            for item in seeklist:
                wr.writerow(item)
                
        with open('seekreszscore.csv', 'w') as file:
            wr=csv.writer(file,lineterminator='\n')
            for item in seek2list:
                wr.writerow(item)


        with open('seekcuezscore.csv', 'w') as file:
            wr=csv.writer(file,lineterminator='\n')
            for item in seekcuelist:
                wr.writerow(item)



        ####
#this takes the save files and uses pandas to integrate behavior data (i.e. trial type, risk block) into the dataframes and get the averages for a subject
# its ugly and should have been a function to save space but it works fine

        df=pandas.read_csv('seekres.csv')
        tpose=df.transpose()
        tpose.to_csv('seekres.csv')

        ######### for the action
        df=pandas.read_csv('seekreszscore.csv')
        tpose=df.transpose().to_csv('seekreszscore.csv')
        tpose=pandas.read_csv('seekreszscore.csv')
        tpose=tpose.transpose()
        trialsdf=arange(1,(len(seekresponse)+1)).tolist()
        tpose.insert(0,'Block',blocksdecoded)
        tpose.insert(1,'shock',shockorder)
        

        ##
        trialholder=tpose.copy()
        trialholder.insert(1,'postshock',aftershock)
        trialholder.insert(0,'Subject',subject)
        trialholder.insert(0,'Session',spd)
        trialholder.insert(2,'Region','mPFC')
        trialholder.insert(0,'Trials',trialsdf)
        if subject+sessarino in Badtrial:
            trialholder=dropbadtrials(Badtrial[subject+sessarino],trialholder,'trialholder')
        else:
            pass   
        trialholder.insert(1,'Latency',actlatencies)
        trialholder.insert(1,'NormLatency',actlatenciesnorm)
        trialholder.insert(1,'FoodLatency',fdlatencies)
        trialholder.insert(1,'Food NormLatency',fdlatenciesnorm)
        trialholder=trialholder.rename(columns={0:'outliertest'})
        #trialholder=trialholder[trialholder['outliertest']==0] 
        AllActionTrials=AllActionTrials.append(trialholder)
        ##
        tpose=tpose.transpose()
        tpose.columns=trialsdf
        tpose=tpose.transpose().rename(columns={0:'outliertest'})
        tpose=tpose[tpose['outliertest']==0]
        tpose=tpose.transpose()
        tpose=tpose.transpose().astype(float).groupby(['Block','shock']).mean()
        tpose.insert(0,'Subject',subject)
        tpose.insert(0,'Session',spd)
        bigAction_mPFC=bigAction_mPFC.append(tpose)
        tpose.transpose().to_csv('ActionMeans.csv')

# for the cue
        df=pandas.read_csv('seekcuezscore.csv')
        tpose=df.transpose().to_csv('seekcuezscore.csv')
        tpose=pandas.read_csv('seekcuezscore.csv')
        tpose=tpose.transpose()
        tpose.insert(0,'Block',blocksdecoded)
        tpose.insert(1,'shock',shockorder)
###
        trialholder=tpose.copy()
        trialholder.insert(1,'postshock',aftershock)
        trialholder.insert(0,'Subject',subject)
        trialholder.insert(0,'Session',spd)
        trialholder.insert(2,'Region','mPFC')
        trialholder.insert(0,'Trials',trialsdf)
        trialholder.insert(1,'Latency',actlatencies)
        trialholder.insert(1,'NormLatency',actlatenciesnorm)
        trialholder.insert(1,'FoodLatency',fdlatencies)
        trialholder.insert(1,'Food NormLatency',fdlatenciesnorm)    
        trialholder=trialholder.rename(columns={0:'outliertest'})
        #trialholder=trialholder[trialholder['outliertest']==0]            
        AllCueTrials=AllCueTrials.append(trialholder)

####
        tpose=tpose.transpose()
        trialsdf=arange(1,(len(seekresponse)+1)).tolist()
        tpose.columns=trialsdf
        tpose=tpose.transpose().rename(columns={0:'outliertest'})
        tpose=tpose[tpose['outliertest']==0]
        tpose=tpose.transpose()
        tpose=tpose.transpose().astype(float).groupby(['Block','shock']).mean()
        tpose.insert(0,'Subject',subject)
        tpose.insert(0,'Session',spd)
        bigCue_mPFC=bigCue_mPFC.append(tpose)        
        tpose.transpose().to_csv('SeekCueMeans.csv')


        os.chdir('../')


        ##################33
# this will do the same as above but for the VTA data. We dont need to reload all the arduino coulbourn data however since its the same session

        if AVTA == 'Y':

            #cam2
            with open('newfile2.csv','r') as csvfile:
                reader=csv.reader(csvfile,delimiter=',')
                for line in reader:
                    deltadata2.append(line)
                csvfile.close()
            deltaarray=np.array(deltadata2)
            stop=len(deltaarray)

            cam1deltaf=[]
            times2=[]
            green=[]
            red=[]

            #cam2
            cam2deltaf=[]
            green2=[]
            red2=[]

            ################

            ## timestamp seconds to sample #


            def tosample (x):
                x=round(x,2)
                return x

            samplesneeded=list(map(tosample,cue))
            seekrespsamplesneeded=list(map(tosample,seekresponse))
            foodsamplesneeded=list(map(tosample,foodcue))
            foodrespsamplesneeded=list(map(tosample,foodresponse)) ###pellet drop not exit

            ##
            #get all the trial samples for regression
            deltascorr=[]
            deltareds=[]
            deltasgreen=[]
            times=[]
            dg=[]
            dr=[]
            samples=[]
            t=1
            stop=len(samplesneeded)


            #####
            #get index
            t=0
            n=1
            stop=len(samplesneeded)
            indexed=[]

            while t+1<= stop:
                ti=float(samplesneeded[t])
                 
                if ti <= (float(deltaarray[n][0])): 
                    frame=n
                    indexed.append(frame) 
                    n=1
                    t=t+1
                else:
                    n=n+1


                    
                    
            ##rescale by trial
            last=float(len(deltaarray))-2
            indexed.append(last)
            #indexed.append(195448) #cheeky add to get last trial if needed
                
            t=1
            stop=len(indexed)

# normalize and get df/f
            while t+1 <= stop:
                x=indexed[t-1]
                x=int(x)
                diff=indexed[t]
                t=t+1
                #41 htz
                x=x-410
                while x < diff:
                    deltafgreen=float(deltaarray[x][1])#-gmean)/gmean
                    dg.append(deltafgreen)
                    samples.append(x)
                    times.append(float(deltaarray[x][0]))

                    deltafred=float(deltaarray[x][2])#-rmean)/rmean
                    dr.append(deltafred)
                    x=x+1

                # #get slope and rescale red to green channel;
                params=np.polyfit(dr,dg,1)
                #print (params)
                slope=float(params[0])
                intercept=float(params[1])

                def rescale2(x):
                    tran = (slope*x)+intercept  # change to "- intercept"
                    return tran

                deltas2=list(map(rescale2,dr))

                for x in dg:
                    deltasgreen.append(x)
                for x in deltas2:
                    deltareds.append(x)

                
                dg=[]
                dr=[]


            #corrected deltaf

            deltascorr=[]
            stop=len(deltareds)

            t=0
            while t+1 <= stop:
                x=(deltasgreen[t]-deltareds[t])/deltareds[t] ##
               # x=(dg[t]-dr[t])/dr[t] ##    
                deltascorr.append(x)
                t=t+1



            #########################
#VTA gcamp, control, and df/f data
            newcsv=list(zip(times,samples,deltasgreen,deltareds,deltascorr))
            array=np.array(newcsv)
            ##################################
            #get all the trial samples for regression
            seeklist=[]
            seek2list=[]
            seekcuelist=[]


            #####
            #get index
            t=0
            n=0
            stop=len(samplesneeded)#-1
            indexed=[]

            while t+1<= stop:
                ti=float(samplesneeded[t])
                  
                
                if ti <= float(array[n][0]): #or ti > float(array[n][0]):
                    frame=n
                    indexed.append(frame) 
                    n=0
                    t=t+1
                else:
                    n=n+1
            t=0
            n=0
            stop=len(seekrespsamplesneeded)#-1
            srindexed=[]

            while t+1<= stop:
                ti=seekrespsamplesneeded[t]
                 
                if ti <= float(array[n][0]):# and ti <= float(array[n][0]+.02):
                    frame=n
                    srindexed.append(frame) 
                    n=0
                    t=t+1
                else:
                    n=n+1
                    


            ######
            #VTA action data
            t=0
            stop=len(srindexed) #srindexed
            while t+1<=stop:
                x=srindexed[t] #srindexed
                y=indexed[t]
                x=int(x)
                y=int(y)
                #Again the z score baseline can be set here
                BL=mean(list((map(float,array[x-164:x-82][:,4]))))
                BLs=np.std(list((map(float,array[x-164:x-82][:,4])))) 
                maxchange=float(max(array[x:x+60][:,4]))
                #window is 2-sec prior until 3 sec after response
                a=list(map(float,array[x-164:x+205][:,4])) #194
                seeklist.append(list(map(float,a)))
                t=t+1
                
                hold=[]
                for x in a:
                    b=(x-BL)/BLs
                    hold.append(b)
                hold=np.array(hold)
                hold=hold.tolist()
# check for outliers
                outliercheck=len(np.where(abs(np.array(hold))>40)[0])
                if outliercheck >=1:
                    motionissuesVTA=motionissuesVTA+1
                else:
                    pass                          
                seek2list.append([outliercheck,(np.trapz(hold[164:205])-np.trapz(hold[122:164])),np.trapz(hold[205:])]+hold)    
            ##
            t=0
            stop=len(indexed) #srindexed
            while t+1<=stop:
                x=indexed[t] #srindexed
                y=indexed[t]
                x=int(x)
                y=int(y)
                BL=mean(list((map(float,array[x-164:x-82][:,4]))))
                BLs=np.std(list((map(float,array[x-164:x-82][:,4])))) 
                maxchange=float(max(array[x:x+41][:,4]))
                

                
                peakchange=(maxchange)
                
                #window is 2-sec prior until 3 sec after response
                a=list(map(float,array[x-164:x+205][:,4])) #194
                t=t+1
                
                #its at 63 htz and there is a .089 sec time lag. So the seek response comes in at sample 132
                hold=[]
                for x in a:
                    b=(x-BL)/BLs
                    hold.append(b)
                hold=np.array(hold)
                hold=hold.tolist()
                outliercheck=len(np.where(abs(np.array(hold))>40)[0])
                if outliercheck >=1:
                    cuemotionissuesVTA=cuemotionissuesVTA+1
                else:
                    pass                                              
                seekcuelist.append([outliercheck,(np.trapz(hold[164:205])-np.trapz(hold[122:164])),np.trapz(hold[205:])]+hold)    
###############################################################################################################33

# VTA plot of whole session trace
            for a in seekresponse:
                plt.axvline(x=a,color='b')
            plt.plot(times,deltasgreen,color='g')
            plt.plot(times,deltareds,color='r')
            plt.title('VTA Streams '+ spd)
            plt.xlabel('Seconds')
            plt.ylabel('Intensity')
            plt.legend("GT",loc='best')
            plt.savefig('VTAStreams'+spd+'.png',dpi=400)
            plt.clf()
            plt.close()



#make a directory for the VTA
            name=CueStamps.split('_')
            rnum=str(round(rand(),4))
            name=name[0]+name[1]+'VTA_Analysis_'
            if os.path.exists(name):
                shutil.rmtree(name)

            os.mkdir(name)
            os.chdir(name)


#save the data
            with open('seekres.csv', 'w') as file:
                wr=csv.writer(file,lineterminator='\n')
                for item in seeklist:
                    wr.writerow(item)
                    
            with open('seekreszscore.csv', 'w') as file:
                wr=csv.writer(file,lineterminator='\n')
                for item in seek2list:
                    wr.writerow(item)


            with open('seekcuezscore.csv', 'w') as file:
                wr=csv.writer(file,lineterminator='\n')
                for item in seekcuelist:
                    wr.writerow(item)


            ####
#like the mPFC goes through the same process of integrating behavior infor in photomrety info to get averages
            import pandas

            df=pandas.read_csv('seekres.csv')

            tpose=df.transpose()
            tpose.to_csv('seekres.csv')



            #########
            df=pandas.read_csv('seekreszscore.csv')
            tpose=df.transpose().to_csv('seekreszscore.csv')
            tpose=pandas.read_csv('seekreszscore.csv')
            tpose=tpose.transpose()
            trialsdf=arange(1,(len(seekresponse)+1)).tolist()            
            tpose.insert(0,'Block',blocksdecoded)
            tpose.insert(1,'shock',shockorder)
            ##
            trialholder=tpose.copy()
            trialholder.insert(1,'postshock',aftershock)
            trialholder.insert(0,'Subject',subject)
            trialholder.insert(0,'Session',spd)
            trialholder.insert(2,'Region','VTA')
            trialholder.insert(0,'Trials',trialsdf)     
            trialholder.insert(1,'Latency',actlatencies)
            trialholder.insert(1,'NormLatency',actlatenciesnorm)     
            trialholder.insert(1,'FoodLatency',fdlatencies)
            trialholder.insert(1,'Food NormLatency',fdlatenciesnorm)
            trialholder=trialholder.rename(columns={0:'outliertest'})
            #trialholder=trialholder[trialholder['outliertest']==0] #we will exclude these outliers later in the R code
            AllActionTrials=AllActionTrials.append(trialholder)
            ##            
            tpose=tpose.transpose()
            tpose.columns=trialsdf
            tpose=tpose.transpose().rename(columns={0:'outliertest'})
            tpose=tpose[tpose['outliertest']==0]
            tpose=tpose.transpose()
            tpose=tpose.transpose().astype(float).groupby(['Block','shock']).mean()
            tpose.insert(0,'Subject',subject)
            tpose.insert(0,'Session',spd)
            bigAction_VTA=bigAction_VTA.append(tpose)        
            tpose.transpose().to_csv('ActionMeans.csv')
        ###

            df=pandas.read_csv('seekcuezscore.csv')
            tpose=df.transpose().to_csv('seekcuezscore.csv')
            tpose=pandas.read_csv('seekcuezscore.csv')
            tpose=tpose.transpose()
            tpose.insert(0,'Block',blocksdecoded)
            tpose.insert(1,'shock',shockorder)
        
            trialholder=tpose.copy()
            trialholder.insert(1,'postshock',aftershock)
            trialholder.insert(0,'Subject',subject)
            trialholder.insert(0,'Session',spd)
            trialholder.insert(2,'Region','VTA')
            trialholder.insert(0,'Trials',trialsdf)
            trialholder.insert(1,'Latency',actlatencies)
            trialholder.insert(1,'NormLatency',actlatenciesnorm)
            trialholder.insert(1,'FoodLatency',fdlatencies)
            trialholder.insert(1,'Food NormLatency',fdlatenciesnorm)
            trialholder=trialholder.rename(columns={0:'outliertest'})
            #trialholder=trialholder[trialholder['outliertest']==0]          
            AllCueTrials=AllCueTrials.append(trialholder)


            tpose=tpose.transpose()
            trialsdf=arange(1,(len(seekresponse)+1)).tolist()
            tpose.columns=trialsdf
            tpose=tpose.transpose().rename(columns={0:'outliertest'})
            tpose=tpose[tpose['outliertest']==0]
            tpose=tpose.transpose()
            tpose=tpose.transpose().astype(float).groupby(['Block','shock']).mean()
            tpose.insert(0,'Subject',subject)
            tpose.insert(0,'Session',spd)
            bigCue_VTA=bigCue_VTA.append(tpose)
            tpose.transpose().to_csv('SeekCueMeans.csv')

            ############



            os.chdir('../../')


        else:
            os.chdir('../')
    os.chdir('../')


#now we can save out all the data to big csvs once we have run all subjects and sesssions
#NOTE thhat rows 14 and 15 are AUCs and not traces themselves. It will be very obvious based on the difference in magnitude if you forget
bigAction_mPFC.sort_values(by=['shock','Session','Block','Subject']).transpose().to_csv('mPFCAction_allthedata.csv')
bigCue_mPFC.sort_values(by=['shock','Session','Block','Subject']).transpose().to_csv('mPFCCue_allthedata.csv')

if AVTA == 'Y':
    bigAction_VTA.sort_values(by=['shock','Session','Block','Subject']).transpose().to_csv('VTAAction_allthedata.csv')
    bigCue_VTA.sort_values(by=['shock','Session','Block','Subject']).transpose().to_csv('VTACue_allthedata.csv')

#all action trials for cross correlation this also holds the food response data
AllActionTrials=AllActionTrials.replace(to_replace=[r'JPDD.*',r'JPDS.*'], value=['Diazepam','Saline'], regex=True)
AllActionTrials.sort_values(by=['shock','Session','Block','Subject','Trials']).transpose().to_csv('AllActionTrials.csv')

AllCueTrials=AllCueTrials.replace(to_replace=[r'JPDD.*',r'JPDS.*'], value=['Diazepam','Saline'], regex=True)
AllCueTrials.sort_values(by=['shock','Session','Block','Subject','Trials']).transpose().to_csv('AllCueTrials.csv')



bigAction_mPFC.replace(to_replace=[r'JPDD.*',r'JPDS.*'], value=['Diazepam','Saline'], regex=True).sort_values(by=['shock','Session','Block','Subject']).transpose().to_csv('mPFCAction_allthedata.csv')
bigCue_mPFC.replace(to_replace=[r'JPDD.*',r'JPDS.*'], value=['Diazepam','Saline'], regex=True).sort_values(by=['shock','Session','Block','Subject']).transpose().to_csv('mPFCCue_allthedata.csv')

if AVTA == 'Y':
    bigAction_VTA.replace(to_replace=[r'JPDD.*',r'JPDS.*'], value=['Diazepam','Saline'], regex=True).sort_values(by=['shock','Session','Block','Subject']).transpose().to_csv('VTAAction_allthedata.csv')
    bigCue_VTA.replace(to_replace=[r'JPDD.*',r'JPDS.*'], value=['Diazepam','Saline'], regex=True).sort_values(by=['shock','Session','Block','Subject']).transpose().to_csv('VTACue_allthedata.csv')


#let us know how many outliers we found...
print ('cue extreme artifacts:'+str(cuemotionissues)+' and '+str(cuemotionissuesVTA))
print ('action+ food extreme artifacts:'+str(motionissues)+' and '+str(motionissuesVTA)+' of '+ str(trialcounter) + ' trials')
