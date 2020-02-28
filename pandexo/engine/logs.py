import os 
import pandas as pd
from sqlalchemy import *
import datetime

def jwst_log(final):
    '''Logging for Website;Used for tracking usage;Only storing instrument data'''
    try:
        logs = os.environ.get('pandexo_logs')
        engine = create_engine('sqlite:///' + os.path.join(logs))
    except: 
        pass 
    instrument = final['pandeia_input']['configuration']['instrument']['instrument']
    mode = final['pandeia_input']['configuration']['instrument']['mode']
    ff = final['pandeia_input']['configuration']['instrument']['filter']
    if ff == None: ff = 'None'
    aperture = final['pandeia_input']['configuration']['instrument']['aperture']
    disperser = final['pandeia_input']['configuration']['instrument']['disperser']
    subarray = final['pandeia_input']['configuration']['detector']['subarray']
    if subarray == None: subarray = "None"
    calc_type= final['pandexo_input']['planet']['type']
    time =  datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    pddb = pd.DataFrame({'instrument':[instrument], 'calc_type':[calc_type],
    	'mode':[mode],'ff':[ff], 'aperture':[aperture],'disperser':[disperser],'subarray':[subarray]},index=[time])
    #print(pddb)
    pddb.to_sql('jwst_cycle1',engine,if_exists='append')

    return


def hst_log(final):
    '''Logging for Website;Used for tracking usage;Used for tracking usage;Only storing instrument data'''
    try:
        logs = os.environ.get('pandexo_logs')
        engine = create_engine('sqlite:///' + os.path.join(logs))
    except: 
        pass 
    instrument = final['pandeia_input']['configuration']['instrument']['instrument']
    subarray = final['pandeia_input']['configuration']['detector']['subarray']
    nsamp = final['pandeia_input']['configuration']['detector']['nsamp']
    samp_seq = final['pandeia_input']['configuration']['detector']['samp_seq']
    disperser = final['pandeia_input']["configuration"]['instrument']['disperser']
    scanDirection = final['pandeia_input']["strategy"]["scanDirection"]
    schedulability = final['pandeia_input']["strategy"]["schedulability"]

    calc_type= final['pandexo_input']['planet']['type']


    time =  datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    pddb = pd.DataFrame({'instrument':[instrument], 'calc_type':[calc_type],
    	'nsamp':[nsamp], 'samp_seq':[samp_seq],'disperser':[disperser],'subarray':[subarray],'scanDirection':[scanDirection],
    	'schedulability':[schedulability]},index=[time])
    #print(pddb)
    pddb.to_sql('hst_cycle26',engine,if_exists='append')
    return
