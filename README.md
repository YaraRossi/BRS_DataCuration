# BRS_DataCuration


## 1_conversion_codes:
This folder includes codes that convert data from LIGO to miniseed

The file "make_xml_BRS_STS2_HAN_LIV.py" creates an xml file "LIGO_HAN_LIV.xml" using two generic instrument responses 
for the BRS and STS2 saved in the files "LIGO_BRS_generic_digitizer.xml" and "LIGO_STS2_generic_digitizer.xml". This xml 
file only corrects for the phase error in the STS2 and converts the STS2 from nm/s to m/s units. 

### TODO's:
1. The BRS data we get  from LIGO is already corrected for the instrument response. In the future this would need to be 
changed, to also include the instrument response from the BRS.
2. the miniseed is saved as floats and therefore the filesize is large. We need to change that to only having integers 
which will then include an additional instrument response correction.
3. Include correct coordinates and elevation.



## 2_analysis_codes:
This folder contains different codes to perform seismic analysis of the LIGO data.
