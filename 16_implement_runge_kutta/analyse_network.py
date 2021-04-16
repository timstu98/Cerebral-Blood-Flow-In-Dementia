import pandas as pd
import numpy as np
import pickle
from pathlib import Path
from tqdm import tqdm
print('start)')
folder = Path.cwd().parent
from pydub import AudioSegment
from pydub.playback import play

number = '002'

print('1')

main_path = folder / '_storage' / 'main' / number

print('2')

network = pd.read_csv(main_path / 'network.csv').drop('Unnamed: 0',axis=1)

# pickle_in = open(main_path / 'network.pickle'  , 'rb' )
# network = pickle.load(pickle_in)
# pickle_in.close()

print('3')

pickle_in2 = open(main_path / 'changes.pickle' , 'rb' )
RK4 = pickle.load(pickle_in2)
pickle_in2.close()

print('4')

pickle_in3 = open(main_path / 'vessels.pickle' , 'rb' )
vessels_alltime = pickle.load(pickle_in3)
pickle_in3.close()

print('5')

def combine_network_changes(network,RK4):
    network['dphidt'] = np.nan
    network['dApdt'] = np.nan
    network['dDpdt'] = np.nan
    network['dAndt'] = np.nan
    network['dDndt'] = np.nan
    network['dcdt'] = np.nan
    for i in tqdm(range(len(network)-5)):
        network.at[i,'dphidt'] = RK4[i][4]
        network.at[i,'dApdt'] = RK4[i][0]
        network.at[i,'dDpdt'] = RK4[i][1]
        network.at[i,'dAndt'] = RK4[i][2]
        network.at[i,'dDndt'] = RK4[i][3]
        network.at[i,'dcdt'] = RK4[i][5]


    network.to_csv(main_path / 'combined_network.csv')

def combine_vessels_changes(vessels_alltime,RK4):
    for i in tqdm(range(len(vessels_alltime)-5)):
        vessels_alltime[i]['dptdt'] = RK4[i][6]['dptdt']
        vessels_alltime[i]['dSoutdt'] = RK4[i][7]['dSoutdt']
        vessels_alltime[i]['dSindt'] = RK4[i][8]['dSindt']

    save = pd.concat(vessels_alltime)

    save.to_csv(main_path / 'combined_vessels.csv')

combine_network_changes(network,RK4)
combine_vessels_changes(vessels_alltime,RK4)



sound = AudioSegment.from_mp3('minecraft.mp3')
play(sound)