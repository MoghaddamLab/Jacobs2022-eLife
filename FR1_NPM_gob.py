#! bin/usr/env python 
#
#Written by dave jacob. a repourposed photometry script to process basic FR1 data. added to functionalisty to split fr1 trials by early 1-30, middle 31-60, and late 61-90
#
#this works with bonsai stanmps to get the trial data. writes out food iti and fr1 actions 3/17/20
# double checked action data 9/26/2020
#
#
#package import
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
from datetime import datetime
from scipy.signal import butter,filtfilt



# function to cluster tone timestamps by the first tone
def clustertones(stamps,thresh=6):
    tones=[]
    run=1
    start=0
    while run == 1:
        tones.append(stamps[start])
        nexttrial=np.where(stamps>stamps[start]+thresh)
        if np.shape(nexttrial)[1] > 0:
            start=nexttrial[0][0]
        else:
            run = 0
    return tones
#lowpass filter function. assumes 41 Hz and default is 3 Hz cutoff
def butter_lowpass_filter(data, cutoff=3, fs=41, order=2):
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y


#data frames for mass aggreagation
bigAction_mPFC=pandas.DataFrame()
mPFC_ActionStrat=pandas.DataFrame()
bigAction_VTA=pandas.DataFrame()
bigCue_mPFC=pandas.DataFrame()
bigCue_VTA=pandas.DataFrame()
bigFood_mPFC=pandas.DataFrame()
bigFood_VTA=pandas.DataFrame()
VTA_ActionStrat=pandas.DataFrame()
AllActionTrials=pandas.DataFrame()

#keep track of outliers (zscore>40)
motionissues=0
motionissuesVTA=0
cuemotionissues=0
cuemotionissuesVTA=0
trialcounter=0
######################################################################################33
Subjects=glob('S*')#get all subjects

for subject in Subjects:
    os.chdir(subject)
    print (subject)
    #get all subjects. we only use 251,252, 256, 283, 318,320,321 because the other ones didnt do any PRT sessions due to lost fibers in the FR1 training
    vtasubjs=['S251','S252','S256','S282','S283','S284','S318','S319','S320','S321']
#get days 1-5
    sessions=glob('D[1,2,3,4,5]') #?*

    for sessarino in sessions:
    #get timestamp data and load files
        os.chdir(sessarino)
        a=glob('corrN*')[0]
        b=glob('feedercue*')[0]
        d=glob('*FR1*')[0]#photometry raw signal file
        c=glob('Tone*')[0]
        #CueStamps=glob('*Cues*')[0] # dont need these because we can use the NP light for action stamps based on way the behavior program was writtem
        #SeekStamps=glob('*Seeks*')[0]
        #FoodDrop=glob('*Food*')[0]
        #FoodRet=glob('*Retreivals*')[0]

        # get FR1 session
        spd=sessarino.split('_')[0]


        ###########################################################################################3

        print('cuereads',datetime.now().time())
        data=[]

# reads in action timees
        with open(a,'r') as csvfile:
            reader=csv.reader(csvfile,delimiter=' ')
            for line in reader:
                data.append(line)
            csvfile.close()

        data2 =np.array(data)[:,0].astype(float)
        Actions=(data2[2::2]-data2[0])/1000
        

        data=[]

# reads in feeder onset times
        with open(b,'r') as csvfile:
            reader=csv.reader(csvfile,delimiter=' ')
            for line in reader:
                data.append(line)
            csvfile.close()

        data2 =np.array(data)[:,0].astype(float)
        FoodRet=(data2[2::2]-data2[0])/1000

        data=[]

#reads in cue onset times and clusters them (since each cue is a pulse it is stamped several times, this just takes the first onset)
        with open(c,'r') as csvfile:
            reader=csv.reader(csvfile,delimiter=' ')
            for line in reader:
                data.append(line)
            csvfile.close()

        data2 =np.array(data)[:,0].astype(float)
        TrialCue=(data2[1::2]-data2[0])/1000

        cue=clustertones(TrialCue)



        ####
        #Change list names to work with rest of code. makes stuff into lists

        #cue=Actions[:89].tolist() 
        seekresponse=Actions[:89].tolist()
        foodcue=FoodRet.tolist() 
        foodresponse=FoodRet.tolist()


        # cause of the weird stamping an extra stamp is occasionally added when the file closes
        if len(seekresponse) > len(cue):
            del seekresponse [-1]
        else:
            pass 




        ##################################################
        # time decoder for flycapture cyclic timekeeping


        #Relevent Functions
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
        #Start reading in FP data


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

        timeconverted=list(map(converttime,ints))

#get converted linear time and normalize to first timestamp to make cumulative time in session
        t=1
        n=0
        time=timeconverted[0]
        cumetime=[time]

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
        #get relevent camera data and rescale the control channel NOTE: this rescale uses the whole session and is not used in the 2022 elife paper. We use the trial by trial data to rescale  in the paper (see later in code)


        camer1=[]
        camerred1=[]
        camer2=[]
        camerred2=[]


        list(map(cam1,data2[:,1]))
        list(map(cam1red,data2[:,2]))
        #if there is no VTA data, i.e. column 4 is empty
        if data2[:,3][0] == '':
            pass
        else:    
            list(map(cam2,data2[:,3]))
            list(map(cam2red,data2[:,4]))
        ####rescale the red channel (not used in paper)
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

        #apply owpass filter
        fs = 41.0       # sample rate, Hz     
        nyq = 0.5 * fs  # Nyquist Frequency            
        camer1=butter_lowpass_filter(camer1).tolist()
        camerred1=butter_lowpass_filter(camerred1).tolist()
        if data2[:,3][0] == '':
            pass
        else:             
            camer2=butter_lowpass_filter(camer2).tolist()    
            camerred2=butter_lowpass_filter(camerred2).tolist()        


#save out low passed data
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

#start parsing by tril        
        print('start finding data',datetime.now().time())
        #delta f/f


        deltadata=[]
        deltadata2=[]
        print ('Analyzing '+sessarino)


        with open('newfile.csv','r') as csvfile:
            reader=csv.reader(csvfile,delimiter=',')
            for line in reader:
                deltadata.append(line)
            csvfile.close()

        deltaarray=np.array(deltadata)
        stop=len(deltaarray)

#dmpfc data
        cam1deltaf=[]
        times2=[]
        green=[]
        red=[]
# vta data
        #cam2
        cam2deltaf=[]
        green2=[]
        red2=[]

        ################

        ## timestamp seconds to sample #


        def tosample (x):
            x=round(x,2)
            return x
#get the time points needed andn index locations
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
             
            if ti <= (float(deltaarray[n][0])): #or ti > (float(deltaarray[n][0])):
                frame=n
                indexed.append(frame) 
                n=1
                t=t+1
            else:
                n=n+1

        print('start rescale',datetime.now().time())
                
                
        ##rescale by trial
        last=float(len(deltaarray))-2
        indexed.append(last)
        #indexed.append(195448) #cheeky add to get last trial if needed
            
        t=1
        stop=len(indexed)

#get 10 sec before and after a trial and rescaleall the data
        while t+1 <= stop:
            x=indexed[t-1]
            x=int(x)
            diff=indexed[t]
            t=t+1
            #10 secs @ 41 htz 
            x=x-410

            if x < 1:
                x=1
                print ('early start')
            else:
                pass

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


        #get corrected corrected deltaf after rescaling the red channel

        deltascorr=[]
        stop=len(deltareds)

        t=0
        while t+1 <= stop:
            x=(deltasgreen[t]-deltareds[t])/deltareds[t] ##
           # x=(dg[t]-dr[t])/dr[t] ##    
            deltascorr.append(x)
            t=t+1


        print('start action re',datetime.now().time())
        #########################
        #get all the delta f and times data


        newcsv=list(zip(times,samples,deltasgreen,deltareds,deltascorr))
        array=np.array(newcsv)



        ##################################
        #get all the trial samples for each action
        seeklist=[]
        seek2list=[]
        seekcuelist=[]
        peaklist=[]
        peakforcsv=[]

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
                

        print('start zseek',datetime.now().time())
        ######
        #action data
        t=0
        stop=len(srindexed) #srindexed
        while t+1<=stop:
            x=srindexed[t] #srindexed
           # print (x)
            y=indexed[t]
            x=int(x)
            y=int(y)
            BL=mean(list((map(float,array[x-164:x-82][:,4]))))
            BLs=np.std(list((map(float,array[x-164:x-82][:,4])))) #BL z score metrics period (default is -4 to -2 sec)
            maxchange=float(max(array[x:x+60][:,4]))
            
            peakchange=(maxchange)
            
            #window is 4-sec prior until 5 sec after response
            a=list(map(float,array[x-164:x+205][:,4])) #194
            seeklist.append(list(map(float,a)))

            t=t+1

            hold=[]
            #normalize to make z score
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
  #save data for later
            seek2list.append([outliercheck,(np.trapz(hold[164:205])-np.trapz(hold[122:164])),np.trapz(hold[205:])]+hold)            

        ##
        # do the same for the cue
        t=0
        stop=len(indexed) #srindexed
        while t+1<=stop:
            x=indexed[t] #srindexed
           # print (x)
            y=indexed[t]
            x=int(x)
            y=int(y)
            BL=mean(list((map(float,array[x-164:x-82][:,4]))))
            BLs=np.std(list((map(float,array[x-164:x-82][:,4])))) 
            maxchange=float(max(array[x:x+41][:,4]))
                        
            peakchange=(maxchange)
            
            a=list(map(float,array[x-164:x+205][:,4])) 

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
            

        print('actionover',datetime.now().time())

        ###############################################################################################################33
# this is an old way to do this. We didnt use this analysis because we go off of the food drop time which is always one second after action for unpunished trials or one sec after shock for punished
        foodlist=[]
        food2list=[]
        foodcuelist=[]
        foodcuelistnonz=[]
        Foodpeakcsv=[]

        #####
        #get index
        t=0
        n=0
        stop=len(foodsamplesneeded)#-1
        indexed=[]

        while t+1<= stop:
            ti=float(foodsamplesneeded[t])
              
            
            if ti <= float(array[n][0]): #or ti > float(array[n][0]):
                frame=n
                indexed.append(frame) 
                n=0
                t=t+1
            else:
                n=n+1
        t=0
        n=0
        stop=len(foodrespsamplesneeded)#-1
        frindexed=[]

        while t+1<= stop:
            ti=foodrespsamplesneeded[t]
             
            if ti <= float(array[n][0]):# and ti <= float(array[n][0]+.02):
                frame=n
                frindexed.append(frame) 
                n=0
                t=t+1
            else:
                n=n+1
                


        ######
        #action data
        t=0
        stop=len(frindexed) #srindexed
        while t+1<=stop:
            x=frindexed[t] #srindexed
           # print (x)
            y=indexed[t]
            x=int(x)
            y=int(y)
            BL=mean(list((map(float,array[x-164:x-82][:,4]))))
            BLs=np.std(list((map(float,array[x-164:x-82][:,4])))) #[y-4100:y-2][:,3]
            maxchange=float(max(array[x:x+60][:,4]))
            
            peakchange=(maxchange)

            a=list(map(float,array[x-82:x+328][:,4])) #194
            foodlist.append(list(map(float,a)))

            t=t+1
            
            hold=[]
            for x in a:
                b=(x-BL)/BLs
                hold.append(b)
            hold=np.array(hold)

            hold=hold.tolist()
            food2list.append(hold)
        ##
        t=0
        stop=len(indexed) 
        while t+1<=stop:
            x=indexed[t] 
           # print (x)
            y=indexed[t]
            x=int(x)
            y=int(y)
            BL=mean(list((map(float,array[x-164:x-82][:,4]))))
            BLs=np.std(list((map(float,array[x-164:x-82][:,4])))) #[y-4100:y-2][:,3]
            maxchange=float(max(array[x:x+60][:,4]))
            
            ###latency
            #latmaxchange=float(max(array[y:x][:,4]))    
            ##
            
            peakchange=(maxchange)
            
            #window is 2-sec prior until 3 sec after response
            a=list(map(float,array[x-82:x+205][:,4])) #194
            foodcuelistnonz.append(list(map(float,a)))
            #peaklist.append(peakchange)
            #latpeaklist.append(latmaxchange)
            t=t+1
            
            #its at 63 htz and there is a .089 sec time lag. So the seek response comes in at sample 132
            hold=[]
            for x in a:
                b=(x-BL)/BLs
                hold.append(b)
            hold=np.array(hold)
            # hold=hold.reshape(-1,2).mean(axis=1)
            # hold=hold.reshape(-1,2).mean(axis=1)
            hold=hold.tolist()
            foodcuelist.append(hold) 




#make a plot to see the overall trace to check for any crazy changes in global signal
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
        #save out the data

        print('startcsv',datetime.now().time())
        #name=CueStamps.split('_')
        rnum=str(round(rand(),4))
        #rnum allows for making a uniquely named directory for storage. This was updated with shutil later on to make it cleaner
        name=subject+sessarino+'mPFC_Analysis_'+rnum


        os.mkdir(name)
        os.chdir(name)


#write out the trial processed data for the subject
        with open('Actions.csv', 'w') as file:
            wr=csv.writer(file,lineterminator='\n')
            for item in seeklist:
                wr.writerow(item)
                
        with open('Actionsszscore.csv', 'w') as file:
            wr=csv.writer(file,lineterminator='\n')
            for item in seek2list:
                wr.writerow(item)


        with open('cuezscore.csv', 'w') as file:
            wr=csv.writer(file,lineterminator='\n')
            for item in seekcuelist:
                wr.writerow(item)

        with open('foodcuezscore.csv', 'w') as file:
            wr=csv.writer(file,lineterminator='\n')
            for item in foodcuelist:
                wr.writerow(item)
        with open('foodITIz.csv', 'w') as file:
            wr=csv.writer(file,lineterminator='\n')
            for item in food2list:
                wr.writerow(item)



        ####
        #IDK maybe use PANDAS

        df=pandas.read_csv('Actions.csv')

        tpose=df.transpose()
        tpose.to_csv('Actions.csv')



        #########
        df=pandas.read_csv('Actionsszscore.csv')
        tpose=df.transpose().to_csv('Actionsszscore.csv')
        tpose=pandas.read_csv('Actionsszscore.csv')
        tpose=tpose.transpose()
        trialsdf=arange(1,(len(seekresponse)+1)).tolist()
        #tpose.insert(0,'Block',blocksdecoded)
        #patch added to catogorize data based on trial location
        trialt=[]
        for a in trialsdf:
            if a < 31:
                cat='early'
                trialt.append(cat)
            elif a > 30 and a <61:
                cat='mid'
                trialt.append(cat)
            elif a > 60:
                cat= 'mlate'
                trialt.append(cat)
        

        ##
        trialholder=tpose.copy()
        stratifiedtrials=trialholder.copy()
        stratifiedtrials.insert(0,'Type',trialt)
        stratifiedtrials=stratifiedtrials.groupby(['Type']).mean()
        stratifiedtrials.insert(0,'Subject',subject)
        stratifiedtrials.insert(0,'Session',spd)
        #stratified is the same data but split by trial grouping set in line 759-767
        mPFC_ActionStrat=mPFC_ActionStrat.append(stratifiedtrials)
        #stratifiedtrials.transpose().to_csv('meansbytype.csv')

        #trialholder.insert(1,'postshock',aftershock)
        trialholder.insert(0,'Subject',subject)
        trialholder.insert(0,'Session',spd)
        trialholder.insert(2,'Region','mPFC')
        trialholder.insert(0,'Trials',trialsdf)
        #trialholder.insert(1,'Latency',actlatencies)
        #trialholder.insert(1,'NormLatency',actlatenciesnorm)
        #trialholder.insert(1,'FoodLatency',fdlatencies)
        #trialholder.insert(1,'Food NormLatency',fdlatenciesnorm)        
        AllActionTrials=AllActionTrials.append(trialholder)
        ##
        tpose=tpose.transpose()
        tpose.columns=trialsdf
        tpose=tpose.transpose().astype(float).mean()
        ids=pandas.Series([subject,spd])
        tpose=ids.append(tpose)
        #this store the average for each subject. used for the plotting
        bigAction_mPFC[subject,spd]=tpose
        tpose.transpose().to_csv('ActionMeans.csv')

        ############


        df=pandas.read_csv('cuezscore.csv')
        tpose=df.transpose().to_csv('cuezscore.csv')
        tpose=pandas.read_csv('cuezscore.csv')
        tpose=tpose.transpose()
        # tpose.insert(0,'Block',blocksdecoded)
        # tpose.insert(1,'shock',shockorder)
        tpose=tpose.transpose()
        trialsdf=arange(1,(len(cue)+1)).tolist()
        tpose.columns=trialsdf
        tpose=tpose.transpose().astype(float).mean()
        ids=pandas.Series([subject,spd])
        tpose=ids.append(tpose)
        bigCue_mPFC[subject,spd]=tpose       
        tpose.transpose().to_csv('CueMeans.csv')





        ############

        print('endcsv',datetime.now().time())


        os.chdir('../')


        ##################33
        #OK run the VTA Stuff too using the same process

        if subject in vtasubjs:

            #new file 2 holds data for cam2 which is the VTA streams
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
                 
                if ti <= (float(deltaarray[n][0])): #or ti > (float(deltaarray[n][0])):
                    frame=n
                    indexed.append(frame) 
                    n=1
                    t=t+1
                else:
                    n=n+1
            ##rescale by trial again
            last=float(len(deltaarray))-2
            indexed.append(last)
            #indexed.append(195448) #cheeky add to get last trial if needed
                
            t=1
            stop=len(indexed)


            while t+1 <= stop:
                x=indexed[t-1]
                x=int(x)
                diff=indexed[t]

                t=t+1
                #41 htz
                x=x-410

                if x < 1:
                    x=1
                    print ('early start')
                else:
                    pass


                while x < diff:
                    deltafgreen=float(deltaarray[x][1])
                    dg.append(deltafgreen)
                    samples.append(x)
                    times.append(float(deltaarray[x][0]))

                    deltafred=float(deltaarray[x][2])
                    dr.append(deltafred)
                    x=x+1

                # #get slope and rescale red to green channel;
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
            #array with corrected DF/f


            newcsv=list(zip(times,samples,deltasgreen,deltareds,deltascorr))
            array=np.array(newcsv)



            ##################################
#open lists for storing action data
            seeklist=[]
            seek2list=[]
            seekcuelist=[]
            peaklist=[]
            peakforcsv=[]


            #####
            #get index of actions
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
            #action data
            t=0
            stop=len(srindexed) #srindexed
            while t+1<=stop:
                x=srindexed[t] #srindexed
               # print (x)
                y=indexed[t]
                x=int(x)
                y=int(y)
                BL=mean(list((map(float,array[x-164:x-82][:,4]))))
                BLs=np.std(list((map(float,array[x-164:x-82][:,4])))) #z window is -4 to -2
                maxchange=float(max(array[x:x+60][:,4]))
                
                peakchange=(maxchange)

                a=list(map(float,array[x-164:x+205][:,4])) #194
                seeklist.append(list(map(float,a)))

                t=t+1

                hold=[]
                for x in a:
                    b=(x-BL)/BLs
                    hold.append(b)
                hold=np.array(hold)

                hold=hold.tolist()
                outliercheck=len(np.where(abs(np.array(hold))>40)[0])
                if outliercheck >=1:
                    motionissuesVTA=motionissuesVTA+1
                else:
                    pass                          
                seek2list.append([outliercheck,(np.trapz(hold[164:205])-np.trapz(hold[122:164])),np.trapz(hold[205:])]+hold)    
            ##For the cue period
            t=0
            stop=len(indexed) #srindexed
            while t+1<=stop:
                x=indexed[t] #srindexed
               # print (x)
                y=indexed[t]
                x=int(x)
                y=int(y)
                BL=mean(list((map(float,array[x-164:x-82][:,4]))))
                BLs=np.std(list((map(float,array[x-164:x-82][:,4])))) 
                maxchange=float(max(array[x:x+41][:,4]))
                
                
                peakchange=(maxchange)

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
                    cuemotionissuesVTA=cuemotionissuesVTA+1
                else:
                    pass                                              
                seekcuelist.append([outliercheck,(np.trapz(hold[164:205])-np.trapz(hold[122:164])),np.trapz(hold[205:])]+hold)    
                



            ###############################################################################################################33
# again this is not used in the paper since it is always 1 sec after action or shock. 
            foodlist=[]
            food2list=[]
            foodcuelist=[]
            foodcuelistnonz=[]
            Foodpeakcsv=[]

            #####
            #get index
            t=0
            n=0
            stop=len(foodsamplesneeded)#-1
            indexed=[]

            while t+1<= stop:
                ti=float(foodsamplesneeded[t])
                  
                
                if ti <= float(array[n][0]): #or ti > float(array[n][0]):
                    frame=n
                    indexed.append(frame) 
                    n=0
                    t=t+1
                else:
                    n=n+1
            t=0
            n=0
            stop=len(foodrespsamplesneeded)#-1
            frindexed=[]

            while t+1<= stop:
                ti=foodrespsamplesneeded[t]
                 
                if ti <= float(array[n][0]):# and ti <= float(array[n][0]+.02):
                    frame=n
                    frindexed.append(frame) 
                    n=0
                    t=t+1
                else:
                    n=n+1
                    


            ######
            #action data
            t=0
            stop=len(frindexed) #srindexed
            while t+1<=stop:
                x=frindexed[t] #srindexed
               # print (x)
                y=indexed[t]
                x=int(x)
                y=int(y)
                BL=mean(list((map(float,array[x-164:x-82][:,4]))))
                BLs=np.std(list((map(float,array[x-164:x-82][:,4]))))
                maxchange=float(max(array[x:x+60][:,4]))
                
                peakchange=(maxchange)
                

                a=list(map(float,array[x-82:x+328][:,4])) #194
                foodlist.append(list(map(float,a)))
                #peaklist.append(peakchange)
                #latpeaklist.append(latmaxchange)
                t=t+1
                
                #its at 63 htz and there is a .089 sec time lag. So the seek response comes in at sample 132
                hold=[]
                for x in a:
                    b=(x-BL)/BLs
                    hold.append(b)
                hold=np.array(hold)
                hold=hold.tolist()
                food2list.append(hold)
            ##
            t=0
            stop=len(indexed) 
            while t+1<=stop:
                x=indexed[t] 
               # print (x)
                y=indexed[t]
                x=int(x)
                y=int(y)
                BL=mean(list((map(float,array[x-164:x-82][:,4]))))
                BLs=np.std(list((map(float,array[x-164:x-82][:,4])))) #z window
                maxchange=float(max(array[x:x+60][:,4]))

                
                peakchange=(maxchange)
                

                a=list(map(float,array[x-82:x+205][:,4])) 
                foodcuelistnonz.append(list(map(float,a)))
                t=t+1
                
                hold=[]
                for x in a:
                    b=(x-BL)/BLs
                    hold.append(b)
                hold=np.array(hold)
                hold=hold.tolist()
                foodcuelist.append(hold) 


#VTA plot of whole session trace for QC
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



            #cam2
            
            rnum=str(round(rand(),4))
            name=subject+sessarino+'VTA_Analysis_'+rnum


            os.mkdir(name)
            os.chdir(name)



            with open('Actions.csv', 'w') as file:
                wr=csv.writer(file,lineterminator='\n')
                for item in seeklist:
                    wr.writerow(item)
                    
            with open('Actionsszscore.csv', 'w') as file:
                wr=csv.writer(file,lineterminator='\n')
                for item in seek2list:
                    wr.writerow(item)

            with open('cuezscore.csv', 'w') as file:
                wr=csv.writer(file,lineterminator='\n')
                for item in seekcuelist:
                    wr.writerow(item)

            with open('foodcuezscore.csv', 'w') as file:
                wr=csv.writer(file,lineterminator='\n')
                for item in foodcuelist:
                    wr.writerow(item)
            with open('foodITIz.csv', 'w') as file:
                wr=csv.writer(file,lineterminator='\n')
                for item in food2list:
                    wr.writerow(item)




            ####
            #IDK maybe use PANDAS
            import pandas

            df=pandas.read_csv('Actions.csv')

            tpose=df.transpose()
            tpose.to_csv('Actions.csv')



            #########
            df=pandas.read_csv('Actionsszscore.csv')
            tpose=df.transpose().to_csv('Actionsszscore.csv')
            tpose=pandas.read_csv('Actionsszscore.csv')
            tpose=tpose.transpose()
            trialsdf=arange(1,(len(seekresponse)+1)).tolist()            
            # add row to code trials by #
            trialt=[]
            for a in trialsdf:
                if a < 31:
                    cat='early'
                    trialt.append(cat)
                elif a > 30 and a <61:
                    cat='mid'
                    trialt.append(cat)
                elif a > 60:
                    cat= 'mlate'
                    trialt.append(cat)
            

            trialholder=tpose.copy()
            stratifiedtrials=trialholder.copy()
            stratifiedtrials.insert(0,'Type',trialt)
            stratifiedtrials=stratifiedtrials.groupby(['Type']).mean()
            stratifiedtrials.insert(0,'Subject',subject)
            stratifiedtrials.insert(0,'Session',spd)
            #stratified trial average goes here
            VTA_ActionStrat=VTA_ActionStrat.append(stratifiedtrials)
            # trialholder.insert(1,'postshock',aftershock)
            trialholder.insert(0,'Subject',subject)
            trialholder.insert(0,'Session',spd)
            trialholder.insert(2,'Region','VTA')
            trialholder.insert(0,'Trials',trialsdf)     
#all action trials not the averages               
            AllActionTrials=AllActionTrials.append(trialholder)
            ##            
            tpose=tpose.transpose()
            tpose.columns=trialsdf
            tpose=tpose.transpose().astype(float).mean()
            ids=pandas.Series([subject,spd])
            tpose=ids.append(tpose)
            #average for all trials a subject goes here
            bigAction_VTA[subject,spd]=tpose       
            tpose.transpose().to_csv('ActionMeans.csv')

            ############


            df=pandas.read_csv('cuezscore.csv')
            tpose=df.transpose().to_csv('cuezscore.csv')
            tpose=pandas.read_csv('cuezscore.csv')
            tpose=tpose.transpose()
            tpose=tpose.transpose()
            trialsdf=arange(1,(len(cue)+1)).tolist()
            tpose.columns=trialsdf
            tpose=tpose.transpose().astype(float).mean()
            ids=pandas.Series([subject,spd])
            tpose=ids.append(tpose)
            bigCue_VTA[subject,spd]=tpose  
            tpose.transpose().to_csv('CueMeans.csv')

            os.chdir('../../')


        else:
            os.chdir('../')
    os.chdir('../')


#write out the data for dmPFC, resort option provided after the #

bigAction_mPFC.transpose().to_csv('mPFCAction_allthedata.csv')#sort_values(by=['Session','Subject']).transpose().to_csv('mPFCAction_allthedata.csv')
bigFood_mPFC.transpose().to_csv('mPFCFood_allthedata.csv')#sort_values(by=['shock','Session','Block','Subject']).transpose().to_csv('mPFCFood_allthedata.csv')
bigCue_mPFC.transpose().to_csv('mPFCCue_allthedata.csv')#sort_values(by=['shock','Session','Block','Subject']).transpose().to_csv('mPFCCue_allthedata.csv')

#write out VTA data, resort option provided after the #

if AVTA == 'Y':
    bigAction_VTA.to_csv('VTAAction_allthedata.csv')#sort_values(by=['Session','Subject']).transpose().to_csv('VTAAction_allthedata.csv')
    bigCue_VTA.to_csv('VTACue_allthedata.csv')#sort_values(by=['shock','Session','Block','Subject']).transpose().to_csv('VTACue_allthedata.csv')
    bigFood_VTA.to_csv('VTAFood_allthedata.csv')#sort_values(by=['shock','Session','Block','Subject']).transpose().to_csv('VTAFood_allthedata.csv')

#save out all trials if wanted. I commented it out since its a large file and not done in the paper
#AllActionTrials.to_csv('AllActionTrials.csv')#sort_values(by=['shock','Session','Block','Subject','Trials']).transpose().to_csv('AllActionTrials.csv')


