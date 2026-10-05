# BRS_DataCuration


## 0_Conversions:
Users of our data do not need this folder.
This folder contains code that convert data from LIGO to miniseed format. The file "make_xml_BRS_STS2_HAN_LIV.py" creates an xml file "LIGO_HAN_LIV.xml" using two generic instrument responses for the BRS and STS2 saved in the files "LIGO_BRS_generic_digitizer.xml" and "LIGO_STS2_generic_digitizer.xml". This xml file only corrects for the phase error in the STS2 and converts the STS2 from nm/s to m/s units. 

## 0_DownloadCode:
Users of our data do not need this folder.
This folder contains the code to download the data from the LIGO facilities.

## 1_Basic_Analysis:
This folder provides basic analysis codes and plots to visualize data. It also includes two files showing how to remove the response of the seismometers at Handford and Livingston.

## 2_Paper_Table:
This folder provides the code for the event table shown in the manuscript.

## *_Paper_Figure_*
These three folder provide the code for the BRS noise sources Section, the rotatin error analysis, Maps and any spectra plot in the manuscript.

