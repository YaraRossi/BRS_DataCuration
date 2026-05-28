from obspy import read_inventory
from functions import rot_vec
import numpy
from numpy import pi
'''root_xml = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/metadata'
inv = read_inventory(f"{root_xml}/LIGO_HAN_LIV.xml")

for net in inv:
    for sta in net:
        print(net.code, sta.code, sta.latitude, sta.longitude)
        for cha in sta:
            print(f"  {cha.location_code}.{cha.code}  az={cha.azimuth:.3f}  dip={cha.dip:.1f}  sr={cha.sample_rate}")

'''

xyz_vector = numpy.asarray([1,0,4])
'''az_clockqisepos=pi/4 # azimuth of y
angle_anticlockwisepos = -az_clockqisepos
psi= -angle_anticlockwisepos'''
psi= (225)/180*pi
enz_vector = rot_vec(0,0,psi,xyz_vector)

print(enz_vector)

root_xml = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/metadata'
inv = read_inventory(f"{root_xml}/LIGO_HAN_LIV.xml")
#inv = response  # whatever you pass as metadata
print(inv[0][0][0].response.instrument_sensitivity)
for st in inv[0][0][0].response.response_stages:
    print(st.stage_sequence_number, getattr(st, "input_units", None), getattr(st, "output_units", None))
