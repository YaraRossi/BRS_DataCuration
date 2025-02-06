from Yara_trialcodes.source_conversion_functions import param2M0
from pyrocko import moment_tensor as pmt

# Usage example
param_file = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/3_SW4/Mendocino_2024/basic_inversion.param'
rupt_file = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/3_SW4/Mendocino_2024/basic_inversion.rupt'
log_file = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/3_SW4/Mendocino_2024/log_file.txt'
srf_file = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/3_SW4/Mendocino_2024/basic_inversion.srf'
'''param2rupt(param_file, rupt_file)
mudpy2srf(rupt_file, log_file, stf_dt=0.1, stf_type='triangle', inv_or_rupt='rupt', Ndip=None,hypocenter=None,
          time_pad=5.0, minSTFpoints=16, integrate=False, Mw=None)
'''

M0 = param2M0(param_file)
print(M0)

mo = 4.028e19
m = pmt.MomentTensor(strike= 280, dip=84, rake=-179, scalar_moment=mo)

print(m)


'''
log_file = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/3_SW4/test_david/ibaraki_3D.000002.log'
rupt_file = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/3_SW4/test_david/ibaraki_3D.000002.rupt'
srf_file = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/3_SW4/test_david/ibaraki_3D.000002.srf'
mudpy2srf(rupt_file, log_file, stf_dt=0.1, stf_type='triangle', inv_or_rupt='rupt', Ndip=None,hypocenter=None,
          time_pad=5.0, minSTFpoints=16, integrate=False, Mw=None)'''
