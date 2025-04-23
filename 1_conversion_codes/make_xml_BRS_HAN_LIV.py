import obspy
from obspy import read_inventory
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
from obspy.clients.nrl import NRL

root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data'
root_xml = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/metadata'

# link to station xml documentation
# http://docs.fdsn.org/projects/stationxml/en/latest/reference.html#

# We'll first create all the various objects. These strongly follow the
# hierarchy of StationXML files.
inv = Inventory(
    # We'll add networks later.
    networks=[],
    # The source should be the id whoever create the file.
    source="ObsPy-Tutorial")

net = Network(
    # This is the network code according to the SEED standard.
    code="LI",
    # A list of stations. We'll add one later.
    stations=[],
    description="All stations from the LIGO sites in the USA, i.e. Hanford and Livingston. Hanford houses two BRS"+
                "sensors collocated with STS2 at the external ends of the lIGO arms, plus a single STS2 at the "+
                "internal corner. Livingston houses 4 BRS sensors with collocated STS-2, at both ends (external &" +
                "internal) of each LIGO arm.",
    # Start-and end dates are optional.
    start_date=obspy.UTCDateTime(2016, 1, 2))

########################################################################################################################
################################################ LIGO Hanford Site #####################################################
########################################################################################################################
sta_han = Station(
    # This is the station code according to the SEED standard.
    code="HAN",
    latitude=46.45510,
    longitude=-119.40758,
    elevation=162.0,
    creation_date=obspy.UTCDateTime(2016, 1, 2),
    site=Site(name="Hanford LIGO Site"))

##############################################################
#### BRS data channels DJX and DJY

cha_han_DJX_EY = Channel(
    # This is the channel code according to the SEED standard.
    code="DJX",
    # This is the location code according to the SEED standard.
    location_code="EY",
    # Note that these coordinates can differ from the station coordinates.
    latitude=46.433908,
    longitude=-119.44963,
    elevation=162,
    depth=0.0, # done
    azimuth=324.0007, # done
    dip=0.0, # done
    sample_rate=256) # done
cha_han_DJY_EX = Channel(
    code="DJY",
    location_code="EX",
    latitude=46.484156,
    longitude=-119.438288,
    elevation=162,
    depth=0.0, # done
    azimuth=234.0007, # done
    dip=0.0,  # done
    sample_rate=256) # done

##############################################################
#### Seis LIGO data channels DH1, DHY and DHZ for locations ES and EX

# Location EX
cha_han_DHX_EX = Channel(
    code="DHX",
    location_code="EX",
    latitude=46.484156,
    longitude=-119.438288,
    elevation=162,
    depth=0.0, # done
    azimuth=324.0007, # done
    dip=0, # done
    sample_rate=512) # done
cha_han_DHY_EX = Channel(
    code="DHY",
    location_code="EX",
    latitude=46.484156,
    longitude=-119.438288,
    elevation=162,
    depth=0.0, # done
    azimuth=234.0007, # done
    dip=0, # done
    sample_rate=512) # done
cha_han_DHZ_EX = Channel(
    code="DHZ",
    location_code="EX",
    latitude=46.484156,
    longitude=-119.438288,
    elevation=162,
    depth=0.0, # done
    azimuth=0, # done
    dip=-90, # done
    sample_rate=512) # done

# Location EY
cha_han_DHX_EY = Channel(
    code="DHX",
    location_code="EY",
    latitude=46.433908,
    longitude=-119.44963,
    elevation=162,
    depth=0.0, # done
    azimuth=324.0007, # done
    dip=0, # done
    sample_rate=512) # done
cha_han_DHY_EY = Channel(
    code="DHY",
    location_code="EY",
    latitude=46.433908,
    longitude=-119.44963,
    elevation=162,
    depth=0.0, # done
    azimuth=234.0007, # done
    dip=0, # done
    sample_rate=512) # done
cha_han_DHZ_EY = Channel(
    code="DHZ",
    location_code="EY",
    latitude=46.433908,
    longitude=-119.44963,
    elevation=162,
    depth=0.0, # done
    azimuth=0, # done
    dip=-90, # done
    sample_rate=512) # done

##############################################################
#### Seis HAN Vault data channels DHX and DHY
# Location VA
cha_han_DHX_VA = Channel(
    code="DHX",
    location_code="VA",
    latitude=46.461704,
    longitude=-119.417683,
    elevation=162,
    depth=0,
    azimuth=323.7681,
    dip=0,
    sample_rate=256)
cha_han_DHY_VA = Channel(
    code="DHY",
    location_code="VA",
    latitude=46.461704,
    longitude=-119.417683,
    elevation=162,
    depth=0,
    azimuth=233.6271,
    dip=0,
    sample_rate=256)

########################################################################################################################
################################################ LIGO Livingston Site #####################################################
########################################################################################################################
sta_liv = Station(
    code="LIV",
    latitude=30.56287,
    longitude=-90.77422,
    elevation=19.0,
    creation_date=obspy.UTCDateTime(2018, 1, 2),
    site=Site(name="Livingston LIGO Site, "))

##############################################################
#### BRS data channels DJX and DJY

cha_liv_DJX_IY = Channel(
    code="DJX",
    location_code="IY",
    latitude=0,
    longitude=0,
    elevation=19,
    depth=0.0, # done
    azimuth=252.2835, # done
    dip=0.0, # done
    sample_rate=256) # done
cha_liv_DJY_IX = Channel(
    code="DJY",
    location_code="IX",
    latitude=0,
    longitude=0,
    elevation=19,
    depth=0.0, # done
    azimuth=162.2835, # done
    dip=0.0, # done
    sample_rate=256) # done
cha_liv_DJX_EY = Channel(
    code="DJX",
    location_code="EY",
    latitude=0,
    longitude=0,
    elevation=19,
    depth=0.0, # done
    azimuth=252.2835, # done
    dip=0.0, # done
    sample_rate=256) # done
cha_liv_DJY_EX = Channel(
    code="DJY",
    location_code="EX",
    latitude=0,
    longitude=0,
    elevation=19,
    depth=0.0, # done # done
    azimuth=162.2835, # done
    dip=0.0, # done
    sample_rate=256) # done

##############################################################
#### Seis LIGO data channels DHX, DHY and DHZ for locations IW, IS, EW and ES

# Location IW
cha_liv_DHX_IX = Channel(
    code="DHX",
    location_code="IX",
    latitude=0,
    longitude=0,
    elevation=19,
    depth=0.0, # done
    azimuth=252.2835, # done
    dip=0, # done
    sample_rate=0) # done
cha_liv_DHY_IX = Channel(
    code="DHY",
    location_code="IX",
    latitude=0,
    longitude=0,
    elevation=19,
    depth=0.0, # done
    azimuth=162.2835, # done
    dip=0, # done
    sample_rate=0) # done
cha_liv_DHZ_IX = Channel(
    code="DHZ",
    location_code="IX",
    latitude=0,
    longitude=0,
    elevation=19,
    depth=0.0, # done
    azimuth=0, # done
    dip=-90, # done
    sample_rate=0) # done

# Location IY
cha_liv_DHX_IY = Channel(
    code="DHX",
    location_code="IY",
    latitude=0,
    longitude=0,
    elevation=19,
    depth=0.0, # done
    azimuth=252.2835, # done
    dip=0, # done
    sample_rate=0) # done
cha_liv_DHY_IY = Channel(
    code="DHY",
    location_code="IY",
    latitude=0,
    longitude=0,
    elevation=19,
    depth=0.0, # done
    azimuth=162.2835, # done
    dip=0, # done
    sample_rate=0) # done
cha_liv_DHZ_IY = Channel(
    code="DHZ",
    location_code="IY",
    latitude=0,
    longitude=0,
    elevation=19,
    depth=0.0, # done
    azimuth=0, # done
    dip=-90, # done
    sample_rate=0) # done

# Location EX
cha_liv_DHX_EX = Channel(
    code="DHX",
    location_code="EX",
    latitude=0,
    longitude=0,
    elevation=19,
    depth=0.0, # done
    azimuth=252.2835, # done
    dip=0, # done
    sample_rate=0) # done
cha_liv_DHY_EX = Channel(
    code="DHY",
    location_code="EX",
    latitude=0,
    longitude=0,
    elevation=19,
    depth=0.0, # done
    azimuth=162.2835, # done
    dip=0, # done
    sample_rate=0) # done
cha_liv_DHZ_EX = Channel(
    code="DHZ",
    location_code="EX",
    latitude=0,
    longitude=0,
    elevation=19,
    depth=0.0, # done
    azimuth=0, # done
    dip=-90, # done
    sample_rate=0) # done

# Location EY
cha_liv_DHX_EY = Channel(
    code="DHX",
    location_code="EY",
    latitude=0,
    longitude=0,
    elevation=19,
    depth=0.0, # done
    azimuth=252.2835, # done
    dip=0, # done
    sample_rate=0) # done
cha_liv_DHY_EY = Channel(
    code="DHY",
    location_code="EY",
    latitude=0,
    longitude=0,
    elevation=19,
    depth=0.0, # done
    azimuth=162.2835, # done
    dip=0, # done
    sample_rate=0) # done
cha_liv_DHZ_EY = Channel(
    code="DHZ",
    location_code="EY",
    latitude=0,
    longitude=0,
    elevation=19,
    depth=0.0, # done
    azimuth=0, # done
    dip=-90, # done
    sample_rate=0) # done

'''
Does this station even exist?
##############################################################
#### Seis LIV Vault data channels DHX and DHY
# Location VA
cha_liv_DHX_VA = Channel(
    code="DHX",
    location_code="VA",
    latitude=???,
    longitude=???,
    elevation=???,
    depth=???,
    azimuth=???,
    dip=???,
    sample_rate=???)
cha_liv_DHY_VA = Channel(
    code="DHY",
    location_code="VA",
    latitude=???,
    longitude=???,
    elevation=???,
    depth=???,
    azimuth=???,
    dip=???,
    sample_rate=???)
cha_liv_DHY_VA = Channel(
    code="DHY",
    location_code="VA",
    latitude=???,
    longitude=???,
    elevation=???,
    depth=???,
    azimuth=???,
    dip=???,
    sample_rate=???)'''

# By default this accesses the NRL online. Offline copies of the NRL can
# also be used instead
'''nrl = NRL()
nrl.sensors()'''
# The contents of the NRL can be explored interactively in a Python prompt,
# see API documentation of NRL submodule:
# http://docs.obspy.org/packages/obspy.clients.nrl.html
# Here we assume that the end point of data logger and sensor are already
# known:
# make my own response files. using the one from blueSeis-3A.

brs_response = read_inventory('%s/LIGO_BRS_generic_digitizer.xml' %root_xml)
network = brs_response[0]  # First network
station = network.stations[0]  # First station
channel = station.channels[0]  # First channel
response_BRS = channel.response

sts_response = read_inventory('%s/LIGO_STS2_generic_digitizer.xml' %root_xml)
# this response file will turn NM/S into M/S
network = sts_response[0]  # First network
station = network.stations[0]  # First station
channel = station.channels[0]  # First channel
response_seisligo = channel.response


'''response_BRS = nrl.get_response( # doctest: +SKIP
    sensor_keys=['UW', 'BRS', '360 seconds'],
    datalogger_keys=['REF TEK', 'RT 130 & 130-SMA', '1', '200'])

response_seisvault = nrl.get_response( # doctest: +SKIP
    sensor_keys=['Streckeisen', 'STS-2', '360 seconds'],
    datalogger_keys=['REF TEK', 'RT 130 & 130-SMA', '1', '200'])

response_seisligo = nrl.get_response( # doctest: +SKIP
    sensor_keys=['Streckeisen', 'STS-2', '1500','3 - installed 04/97 to present'],
    datalogger_keys=['REF TEK', 'RT 130 & 130-SMA', '1', '200'])'''


# Now tie it all together.
# add response to channels HAN
cha_han_DJX_EY.response = response_BRS
cha_han_DJY_EX.response = response_BRS

cha_han_DHX_EX.response = response_seisligo
cha_han_DHY_EX.response = response_seisligo
cha_han_DHZ_EX.response = response_seisligo
cha_han_DHX_EY.response = response_seisligo
cha_han_DHY_EY.response = response_seisligo
cha_han_DHZ_EY.response = response_seisligo
cha_han_DHX_VA.response = response_seisligo
cha_han_DHY_VA.response = response_seisligo

# add channels to station Hanford
sta_han.channels.append(cha_han_DJX_EY)
sta_han.channels.append(cha_han_DJY_EX)
sta_han.channels.append(cha_han_DHX_EX)
sta_han.channels.append(cha_han_DHY_EX)
sta_han.channels.append(cha_han_DHZ_EX)
sta_han.channels.append(cha_han_DHX_EY)
sta_han.channels.append(cha_han_DHY_EY)
sta_han.channels.append(cha_han_DHZ_EY)
sta_han.channels.append(cha_han_DHX_VA)
sta_han.channels.append(cha_han_DHY_VA)

# add response to channels LIV
cha_liv_DJX_IY.response = response_BRS
cha_liv_DJY_IX.response = response_BRS
cha_liv_DJX_EY.response = response_BRS
cha_liv_DJY_EX.response = response_BRS
cha_liv_DHX_IY.response = response_seisligo
cha_liv_DHY_IY.response = response_seisligo
cha_liv_DHZ_IY.response = response_seisligo
cha_liv_DHX_IX.response = response_seisligo
cha_liv_DHY_IX.response = response_seisligo
cha_liv_DHZ_IX.response = response_seisligo
cha_liv_DHX_EY.response = response_seisligo
cha_liv_DHY_EY.response = response_seisligo
cha_liv_DHZ_EY.response = response_seisligo
cha_liv_DHX_EX.response = response_seisligo
cha_liv_DHY_EX.response = response_seisligo
cha_liv_DHZ_EX.response = response_seisligo

# add channels to station Hanford
sta_liv.channels.append(cha_liv_DJX_IY)
sta_liv.channels.append(cha_liv_DJY_IX)
sta_liv.channels.append(cha_liv_DJX_EY)
sta_liv.channels.append(cha_liv_DJY_EX)
sta_liv.channels.append(cha_liv_DHX_IY)
sta_liv.channels.append(cha_liv_DHY_IY)
sta_liv.channels.append(cha_liv_DHZ_IY)
sta_liv.channels.append(cha_liv_DHX_IX)
sta_liv.channels.append(cha_liv_DHY_IX)
sta_liv.channels.append(cha_liv_DHZ_IX)
sta_liv.channels.append(cha_liv_DHX_EY)
sta_liv.channels.append(cha_liv_DHY_EY)
sta_liv.channels.append(cha_liv_DHZ_EY)
sta_liv.channels.append(cha_liv_DHX_EX)
sta_liv.channels.append(cha_liv_DHY_EX)
sta_liv.channels.append(cha_liv_DHZ_EX)

# add station Hanford to Network LI
net.stations.append(sta_han)
# add station Hanford to Network LI
net.stations.append(sta_liv)






inv.networks.append(net)

# And finally write it to a StationXML file. We also force a validation against
# the StationXML schema to ensure it produces a valid StationXML file.
#
# Note that it is also possible to serialize to any of the other inventory
# output formats ObsPy supports.
inv.write("%s/LIGO_HAN_LIV.xml" %root_xml, format="stationxml", validate=True)