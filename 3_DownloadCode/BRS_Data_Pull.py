# Installing Packages
from gwpy.timeseries import TimeSeries , TimeSeriesDict
import numpy as np
from gwpy.plot import Plot
import pandas as pd
from gwpy.time import tconvert,to_gps
import matplotlib.pyplot as plt
import calendar
import sys
from inputimeout import inputimeout , TimeoutOccurred
import os
from pathlib import Path  

#Setting up THIS IS WHERE YOU CHANGE THINGS, this is the manual input, will change this so you can upload a excel sheet?

month=int(input("Enter Month Number: " ))
day=input("Enter Day (DD): ") #ENTER DAY NUMBER
year=input("Enter Year (YYYY): ") #ENTER YEAR
time=input("Enter time (HH:MM:SS): ")
dt = int(input('Enter duration in minutes: '))*60 #ENTER TIME IN MINUTES, CONVERTS IT TO SECONDS
site=input("Which sites do you want data from? (H/L/B):  ")

print(f'Is this correct? \n Month={month} Day={day} Year={year} Time={time} Duration={dt/60} Site(s)={site}') #CHECKING IS TIME ENTERED IS CORRECT

if input('[y/n]?  ')=='n': #SETTING UP SO USER CONFIRMS THAT THIS IS CORRECT
    print("Program will now end (Shoshana Fix this)")
    sys.exit()

StartingTimeString=f"{calendar.month_abbr[month]} {day} {year} {time}.000" #getting starting time to form that gwpy can read
starttime=to_gps(StartingTimeString) #GETTING TIME TO GPS

#CONVERTING TIME TO FILE FORMAT
StartingTimeString=StartingTimeString.replace(" ", "_")
StartingTimeString=StartingTimeString.replace(":", "-")
StartingTimeString=StartingTimeString.replace(".", "-")


# Loading in all channel names, seperated out by seismometer, BRS's and Vault channels to prevent issues with differing frequencies
L1_SEI=['L1:ISI-GND_STS_ITMX_X_DQ','L1:ISI-GND_STS_ITMX_Y_DQ','L1:ISI-GND_STS_ITMX_Z_DQ','L1:ISI-GND_STS_ITMY_X_DQ','L1:ISI-GND_STS_ITMY_Y_DQ','L1:ISI-GND_STS_ITMY_Z_DQ','L1:ISI-GND_STS_ETMX_X_DQ','L1:ISI-GND_STS_ETMX_Y_DQ','L1:ISI-GND_STS_ETMX_Z_DQ','L1:ISI-GND_STS_ETMY_X_DQ','L1:ISI-GND_STS_ETMY_Y_DQ','L1:ISI-GND_STS_ETMY_Z_DQ'] #all L1  seismometers
L1_BRS=['L1:ISI-GND_BRS_ETMX_RY_OUT_DQ','L1:ISI-GND_BRS_ETMY_RX_OUT_DQ','L1:ISI-GND_BRS_ITMX_RY_OUT_DQ','L1:ISI-GND_BRS_ITMY_RX_OUT_DQ']

H1_SEI=['H1:ISI-GND_STS_ETMX_X_DQ','H1:ISI-GND_STS_ETMX_Y_DQ','H1:ISI-GND_STS_ETMX_Z_DQ','H1:ISI-GND_STS_ETMY_X_DQ','H1:ISI-GND_STS_ETMY_Y_DQ','H1:ISI-GND_STS_ETMY_Z_DQ','H1:ISI-GND_STS_ITMY_X_DQ','H1:ISI-GND_STS_ITMY_Y_DQ','H1:ISI-GND_STS_ITMY_Z_DQ']
H1_BRS=['H1:ISI-GND_BRS_ETMX_RY_OUT_DQ','H1:ISI-GND_BRS_ETMY_RX_OUT_DQ']
VAULT=['H1:PEM-VAULT_SEIS_1030X195Y_STS2_X_DQ','H1:PEM-VAULT_SEIS_1030X195Y_STS2_Y_DQ']



if site=="H" or site=="B": #Getting data for Hanford site
    FolderName=f'H1_{StartingTimeString}' #MAKING NAME FOR NEW FOLDER TO KEEP MY COMPUTER SLIGHTLY ORGANIZED
    if not os.path.exists(FolderName): #CHECKING TO SEE IF A FOLDER WITH THIS NAME ALREADY EXISTIS 
        os.makedirs(FolderName) #MAKES NEW FOLDER
    TypeNames=[VAULT,H1_BRS,H1_SEI]
    for j in range(3): #THIS IS ITERATING OVER ALL OF THE TYPES OF SENSOR HOPEFULLY
        H1_Data=TimeSeriesDict.fetch(TypeNames[j],starttime, starttime+dt) #FETCHING ALL OF THE CHANNELS WILL BE BIG OBJECT
        df_H1 = (pd.DataFrame(data=H1_Data)) #CONVERTING ALL OF THE CHANNELS TO DF
        H1_Timing = (pd.DataFrame(data=((list(H1_Data.values())[0]).times))) #THE NIGHTMARE THAT IS GETTING THE TIMING DATA. THIS IS THE SAME FOR ALL OF THE CHANNELS IN EACH GROUP
    
        for i in range(len(TypeNames[j])): #data for each channe
            df_Channel=df_H1.iloc[:,i].copy() #SEPERATING OUT EACH CHANNEL FOR EACH GROUP, THIS SHOULD BE A TEMPERARY VARIABLE
            df=H1_Timing.merge(df_Channel,left_index=True,right_index=True) #COMBINING THE TIMING DATA AND THE CHANNEL DATA
            FileName=f"{TypeNames[j][i].replace(':', '_')}_{StartingTimeString}_{str((list(H1_Data.values())[0]).sample_rate).replace('.0 ','')}.csv" #INSANE SET UP TO SAVE FILE NAMES
            path=Path(FolderName+'/'+FileName) #MAKE PATH FOR NEW FILE
            df.to_csv(path,index=False,header=False) #SAVING DATA TO CSV FILE
    i=0 #I HATE NOT KNOWING EXACTLY WHAT VARIABLE IS WHAT 
    j=0
    try:
        plotask = inputimeout(prompt='Do you want to see a plot? [y/n]:  ', timeout=10)
    except TimeoutOccurred:
        plotask ='n'
    if plotask=='y':
        plot=Plot(H1_Data)
        plt.show()
    print('All done with Hanford data :)')


#NOW LIVINGSTON DATA
if site=="L" or site=="B": #Getting data for Livingston site
    FolderName=f'L1_{StartingTimeString}' #MAKING NAME FOR NEW FOLDER TO KEEP MY COMPUTER SLIGHTLY ORGANIZED
    if not os.path.exists(FolderName): #CHECKING TO SEE IF A FOLDER WITH THIS NAME ALREADY EXISTIS 
        os.makedirs(FolderName) #MAKES NEW FOLDER
    TypeNames=[L1_BRS,L1_SEI]
    for j in range(2): #THIS IS ITERATING OVER ALL OF THE TYPES OF SENSOR
        L1_Data=TimeSeriesDict.fetch(TypeNames[j],starttime, starttime+dt) #FETCHING ALL OF THE CHANNELS WILL BE BIG OBJECT
        df_L1 = (pd.DataFrame(data=L1_Data)) #CONVERTING ALL OF THE CHANNELS TO DF
        L1_Timing = (pd.DataFrame(data=((list(L1_Data.values())[0]).times))) #THE NIGHTMARE THAT IS GETTING THE TIMING DATA. THIS IS THE SAME FOR ALL OF THE CHANNELS IN EACH GROUP
    
        for i in range(len(TypeNames[j])): #data for each channe
            df_Channel=df_L1.iloc[:,i].copy() #SEPERATING OUT EACH CHANNEL FOR EACH GROUP, THIS SHOULD BE A TEMPERARY VARIABLE
            df=L1_Timing.merge(df_Channel,left_index=True,right_index=True) #COMBINING THE TIMING DATA AND THE CHANNEL DATA
            FileName=f"{TypeNames[j][i].replace(':', '_')}_{StartingTimeString}_{str((list(L1_Data.values())[0]).sample_rate).replace('.0 ','')}.csv" #INSANE SET UP TO SAVE FILE NAMES
            path=Path(FolderName+'/'+FileName) #MAKE PATH FOR NEW FILE
            df.to_csv(path,index=False,header=False) #SAVING DATA TO CSV FILE

    try: #ASKING IF YOU WANT TO LOOK AT A PLOT
        plotask = inputimeout(prompt='Do you want to see a plot? [y/n]:  ', timeout=10) #IF NO RESPONCE IN 10 SECONDS ASSUME NO
    except TimeoutOccurred:
        plotask ='n'
    if plotask=='y':
        plot=Plot(L1_Data) #THIS SHOULD BE THE SEISMOMETER DATA
        plt.show()
    print('All done with Livingston data :)')


