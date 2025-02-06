import obspy
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
from obspy.clients.nrl import NRL

root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data'


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
    description="All stations from the LIGO sites in the USA, i.e. Hanford and Livingston.",
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
    depth=0.0,
    azimuth=323.7681,
    dip=0.0,
    sample_rate=256)

cha_han_DJY_EX = Channel(
    # This is the channel code according to the SEED standard.
    code="DJY",
    # This is the location code according to the SEED standard.
    location_code="EX",
    # Note that these coordinates can differ from the station coordinates.
    latitude=46.484156,
    longitude=-119.438288,
    elevation=162,
    depth=0.0,
    azimuth=233.6271,
    dip=0.0,
    sample_rate=256)

##############################################################
#### Seis LIGO data channels DH1, DHY and DHZ for locations ES and EX

# Location EX
cha_han_DHX_EX = Channel(
    # This is the channel code according to the SEED standard.
    code="DHX",
    # This is the location code according to the SEED standard.
    location_code="EX",
    # Note that these coordinates can differ from the station coordinates.
    latitude=46.484156,
    longitude=-119.438288,
    elevation=162,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=512)
cha_han_DHY_EX = Channel(
    # This is the channel code according to the SEED standard.
    code="DHY",
    # This is the location code according to the SEED standard.
    location_code="EX",
    # Note that these coordinates can differ from the station coordinates.
    latitude=46.484156,
    longitude=-119.438288,
    elevation=162,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=512)
cha_han_DHZ_EX = Channel(
    # This is the channel code according to the SEED standard.
    code="DHZ",
    # This is the location code according to the SEED standard.
    location_code="EX",
    # Note that these coordinates can differ from the station coordinates.
    latitude=46.484156,
    longitude=-119.438288,
    elevation=162,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=512)

# Location EY
cha_han_DHX_EY = Channel(
    # This is the channel code according to the SEED standard.
    code="DHX",
    # This is the location code according to the SEED standard.
    location_code="EY",
    # Note that these coordinates can differ from the station coordinates.
    latitude=46.433908,
    longitude=-119.44963,
    elevation=162,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=512)
cha_han_DHY_EY = Channel(
    # This is the channel code according to the SEED standard.
    code="DHY",
    # This is the location code according to the SEED standard.
    location_code="EY",
    # Note that these coordinates can differ from the station coordinates.
    latitude=46.433908,
    longitude=-119.44963,
    elevation=162,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=512)
cha_han_DHZ_EY = Channel(
    # This is the channel code according to the SEED standard.
    code="DHZ",
    # This is the location code according to the SEED standard.
    location_code="EY",
    # Note that these coordinates can differ from the station coordinates.
    latitude=46.433908,
    longitude=-119.44963,
    elevation=162,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=512)

##############################################################
#### Seis HAN Vault data channels DHX and DHY
# Location VA
cha_han_DHX_VA = Channel(
    # This is the channel code according to the SEED standard.
    code="DHX",
    # This is the location code according to the SEED standard.
    location_code="VA",
    # Note that these coordinates can differ from the station coordinates.
    latitude=46.461704,
    longitude=-119.417683,
    elevation=162,
    depth=???,
    azimuth=???,
    dip=0,
    sample_rate=256)
cha_han_DHY_VA = Channel(
    # This is the channel code according to the SEED standard.
    code="DHY",
    # This is the location code according to the SEED standard.
    location_code="VA",
    # Note that these coordinates can differ from the station coordinates.
    latitude=46.461704,
    longitude=-119.417683,
    elevation=162,
    depth=???,
    azimuth=???,
    dip=0,
    sample_rate=256)

########################################################################################################################
################################################ LIGO Livingston Site #####################################################
########################################################################################################################
sta_liv = Station(
    # This is the station code according to the SEED standard.
    code="LIV",
    latitude=30.56287,
    longitude=-90.77422,
    elevation=19.0,
    creation_date=obspy.UTCDateTime(2016, 1, 2),
    site=Site(name="Livingston LIGO Site"))

##############################################################
#### BRS data channels DJX and DJY

cha_liv_DJX_IY = Channel(
    # This is the channel code according to the SEED standard.
    code="DJX",
    # This is the location code according to the SEED standard.
    location_code="IY",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=0.0,
    sample_rate=256)

cha_liv_DJY_IX = Channel(
    # This is the channel code according to the SEED standard.
    code="DJY",
    # This is the location code according to the SEED standard.
    location_code="IX",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=0.0,
    sample_rate=256)

cha_liv_DJX_EY = Channel(
    # This is the channel code according to the SEED standard.
    code="DJX",
    # This is the location code according to the SEED standard.
    location_code="EY",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=0.0,
    sample_rate=256)

cha_liv_DJY_EX = Channel(
    # This is the channel code according to the SEED standard.
    code="DJY",
    # This is the location code according to the SEED standard.
    location_code="EX",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=0.0,
    sample_rate=256)

##############################################################
#### Seis LIGO data channels DHX, DHY and DHZ for locations IW, IS, EW and ES

# Location IW
cha_liv_DHX_IX = Channel(
    # This is the channel code according to the SEED standard.
    code="DHX",
    # This is the location code according to the SEED standard.
    location_code="IX",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=???)
cha_liv_DHY_IX = Channel(
    # This is the channel code according to the SEED standard.
    code="DHY",
    # This is the location code according to the SEED standard.
    location_code="IX",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=???)
cha_liv_DHZ_IX = Channel(
    # This is the channel code according to the SEED standard.
    code="DHZ",
    # This is the location code according to the SEED standard.
    location_code="IX",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=???)

# Location IY
cha_liv_DHX_IY = Channel(
    # This is the channel code according to the SEED standard.
    code="DHX",
    # This is the location code according to the SEED standard.
    location_code="IY",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=???)
cha_liv_DHY_IY = Channel(
    # This is the channel code according to the SEED standard.
    code="DHY",
    # This is the location code according to the SEED standard.
    location_code="IY",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=???)
cha_liv_DHZ_IY = Channel(
    # This is the channel code according to the SEED standard.
    code="DHZ",
    # This is the location code according to the SEED standard.
    location_code="IY",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=???)

# Location EX
cha_liv_DHX_EX = Channel(
    # This is the channel code according to the SEED standard.
    code="DHX",
    # This is the location code according to the SEED standard.
    location_code="EX",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=???)
cha_liv_DHY_EX = Channel(
    # This is the channel code according to the SEED standard.
    code="DHY",
    # This is the location code according to the SEED standard.
    location_code="EX",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=???)
cha_liv_DHZ_EX = Channel(
    # This is the channel code according to the SEED standard.
    code="DHZ",
    # This is the location code according to the SEED standard.
    location_code="EX",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=???)

# Location EY
cha_liv_DHX_EY = Channel(
    # This is the channel code according to the SEED standard.
    code="DHX",
    # This is the location code according to the SEED standard.
    location_code="EY",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=???)
cha_liv_DHY_EY = Channel(
    # This is the channel code according to the SEED standard.
    code="DHY",
    # This is the location code according to the SEED standard.
    location_code="EY",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=???)
cha_liv_DHZ_EY = Channel(
    # This is the channel code according to the SEED standard.
    code="DHZ",
    # This is the location code according to the SEED standard.
    location_code="EY",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=???)

'''
Does this station even exist?
##############################################################
#### Seis LIV Vault data channels DHX and DHY
# Location VA
cha_liv_DHX_VA = Channel(
    # This is the channel code according to the SEED standard.
    code="DHX",
    # This is the location code according to the SEED standard.
    location_code="VA",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=???,
    depth=???,
    azimuth=???,
    dip=???,
    sample_rate=???)
cha_liv_DHY_VA = Channel(
    # This is the channel code according to the SEED standard.
    code="DHY",
    # This is the location code according to the SEED standard.
    location_code="VA",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=???,
    depth=???,
    azimuth=???,
    dip=???,
    sample_rate=???)
cha_liv_DHY_VA = Channel(
    # This is the channel code according to the SEED standard.
    code="DHY",
    # This is the location code according to the SEED standard.
    location_code="VA",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=???,
    depth=???,
    azimuth=???,
    dip=???,
    sample_rate=???)'''

# By default this accesses the NRL online. Offline copies of the NRL can
# also be used instead
nrl = NRL()
# The contents of the NRL can be explored interactively in a Python prompt,
# see API documentation of NRL submodule:
# http://docs.obspy.org/packages/obspy.clients.nrl.html
# Here we assume that the end point of data logger and sensor are already
# known:
response_BRS = nrl.get_response( # doctest: +SKIP
    sensor_keys=['UW', 'BRS', '360 seconds'],
    datalogger_keys=['REF TEK', 'RT 130 & 130-SMA', '1', '200'])


response_seisvault = nrl.get_response( # doctest: +SKIP
    sensor_keys=['Streckeisen', 'STS-1', '360 seconds'],
    datalogger_keys=['REF TEK', 'RT 130 & 130-SMA', '1', '200'])

response_seisligo = nrl.get_response( # doctest: +SKIP
    sensor_keys=['Streckeisen', 'STS-1', '360 seconds'],
    datalogger_keys=['REF TEK', 'RT 130 & 130-SMA', '1', '200'])


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
cha_han_DHX_VA.response = response_seisvault
cha_han_DHY_VA.response = response_seisvault

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
'''does this instrument exist?
cha_liv_DHX_VA.response = response_seisvault_liv
cha_liv_DHY_VA.response = response_seisvault_liv
cha_liv_DHZ_VA.response = response_seisvault_liv'''

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
'''does this instrument exist?
sta_liv.channels.append(cha_liv_DHX_VA)
sta_liv.channels.append(cha_liv_DHY_VA)
sta_liv.channels.append(cha_liv_DHZ_VA)'''

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
inv.write("%s/LIGO_HAN_LIV.xml" %root_data, format="stationxml", validate=True)