import numpy
from obspy import read, UTCDateTime, read_inventory
from functions import rot_vec

root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'
root_xml = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/metadata'
folder = '2024-04-02'

## 1. READ DATA
# BRS
LIV_IY_DJX = read('%s/%s/LI_LIV_IY_DJX.ms' % (root_data, folder))
LIV_IX_DJY = read('%s/%s/LI_LIV_IX_DJY.ms' % (root_data, folder))
# SEIS
LIV_IX_DHX = read('%s/%s/LI_LIV_IX_DHX.ms' % (root_data, folder))
LIV_IX_DHY = read('%s/%s/LI_LIV_IX_DHY.ms' % (root_data, folder))
LIV_IX_DHZ = read('%s/%s/LI_LIV_IX_DHZ.ms' % (root_data, folder))

## 2. FILTER DATA
# BRS
LIV_IY_DJX = LIV_IY_DJX.filter('bandpass', freqmin=0.008, freqmax = 100, zerophase=True).taper(0.02)
LIV_IX_DJY = LIV_IX_DJY.filter('bandpass', freqmin=0.008, freqmax = 100, zerophase=True).taper(0.02)
# SEIS
LIV_IX_DHX = LIV_IX_DHX.filter('bandpass', freqmin=0.08, freqmax = 100, zerophase=True).taper(0.02)
LIV_IX_DHY = LIV_IX_DHY.filter('bandpass', freqmin=0.08, freqmax = 100, zerophase=True).taper(0.02)
LIV_IX_DHZ = LIV_IX_DHZ.filter('bandpass', freqmin=0.08, freqmax = 100, zerophase=True).taper(0.02)

response = read_inventory('%s/LIGO_HAN_LIV.xml' %root_xml)

## 3. REMOVE RESPONSE
# BRS
LIV_IY_DJX_rem = LIV_IY_DJX.copy()
LIV_IX_DJY_rem = LIV_IX_DJY.copy()
LIV_IY_DJX_rem[0].remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
LIV_IX_DJY_rem[0].remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)

# SEIS
LIV_IX_DHX_rem = LIV_IX_DHX.copy()
LIV_IX_DHY_rem = LIV_IX_DHY.copy()
LIV_IX_DHZ_rem = LIV_IX_DHZ.copy()
LIV_IX_DHX_rem[0].remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
LIV_IX_DHY_rem[0].remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
LIV_IX_DHZ_rem[0].remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)

## 4. TURN FROM XYZ to ENU
# BRS
channel = response.select(network='LI', station='LIV', location='IX', channel='DJY')[0][0][0]
az = channel.azimuth/180*numpy.pi
LIV_I_DJ = numpy.array([LIV_IY_DJX[0].data, LIV_IX_DJY[0].data, LIV_IX_DJY[0].data * 0])
LIV_I_DJ_baz = rot_vec(0, 0, az, LIV_I_DJ)
LIV_I_E_orig = LIV_IY_DJX.copy()
LIV_I_N_orig = LIV_IX_DJY.copy()
LIV_I_E_orig[0].data = LIV_I_DJ_baz[0, :]
LIV_I_N_orig[0].data = LIV_I_DJ_baz[1, :]
LIV_I_DJ  = LIV_I_E_orig
LIV_I_DJ += LIV_I_N_orig
print(LIV_I_DJ)

# SEIS
LIV_I_DH  = LIV_IX_DHX_rem
LIV_I_DH += LIV_IX_DHY_rem
LIV_I_DH += LIV_IX_DHZ_rem
LIV_I_DH.rotate(method="->ZNE", inventory=response)

print(LIV_I_DH)