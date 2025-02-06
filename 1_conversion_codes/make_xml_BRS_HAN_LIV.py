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
#### BRS data channels DJ1 and DJ2

cha_han_DJ1_ES = Channel(
    # This is the channel code according to the SEED standard.
    code="DJ1",
    # This is the location code according to the SEED standard.
    location_code="ES",
    # Note that these coordinates can differ from the station coordinates.
    latitude=46.433908,
    longitude=-119.44963,
    elevation=162,
    depth=0.0,
    azimuth=323.7681,
    dip=0.0,
    sample_rate=256)

cha_han_DJ2_EN = Channel(
    # This is the channel code according to the SEED standard.
    code="DJ2",
    # This is the location code according to the SEED standard.
    location_code="EN",
    # Note that these coordinates can differ from the station coordinates.
    latitude=46.484156,
    longitude=-119.438288,
    elevation=162,
    depth=0.0,
    azimuth=233.6271,
    dip=0.0,
    sample_rate=256)

##############################################################
#### Seis LIGO data channels DH1, DH2 and DH3 for locations ES and EN

# Location EN
cha_han_DH1_EN = Channel(
    # This is the channel code according to the SEED standard.
    code="DH1",
    # This is the location code according to the SEED standard.
    location_code="EN",
    # Note that these coordinates can differ from the station coordinates.
    latitude=46.484156,
    longitude=-119.438288,
    elevation=162,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=512)
cha_han_DH2_EN = Channel(
    # This is the channel code according to the SEED standard.
    code="DH2",
    # This is the location code according to the SEED standard.
    location_code="EN",
    # Note that these coordinates can differ from the station coordinates.
    latitude=46.484156,
    longitude=-119.438288,
    elevation=162,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=512)
cha_han_DH3_EN = Channel(
    # This is the channel code according to the SEED standard.
    code="DH3",
    # This is the location code according to the SEED standard.
    location_code="EN",
    # Note that these coordinates can differ from the station coordinates.
    latitude=46.484156,
    longitude=-119.438288,
    elevation=162,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=512)

# Location ES
cha_han_DH1_ES = Channel(
    # This is the channel code according to the SEED standard.
    code="DH1",
    # This is the location code according to the SEED standard.
    location_code="ES",
    # Note that these coordinates can differ from the station coordinates.
    latitude=46.433908,
    longitude=-119.44963,
    elevation=162,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=512)
cha_han_DH2_ES = Channel(
    # This is the channel code according to the SEED standard.
    code="DH2",
    # This is the location code according to the SEED standard.
    location_code="ES",
    # Note that these coordinates can differ from the station coordinates.
    latitude=46.433908,
    longitude=-119.44963,
    elevation=162,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=512)
cha_han_DH3_ES = Channel(
    # This is the channel code according to the SEED standard.
    code="DH3",
    # This is the location code according to the SEED standard.
    location_code="ES",
    # Note that these coordinates can differ from the station coordinates.
    latitude=46.433908,
    longitude=-119.44963,
    elevation=162,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=512)

##############################################################
#### Seis HAN Vault data channels DH1 and DH2
# Location VA
cha_han_DH1_VA = Channel(
    # This is the channel code according to the SEED standard.
    code="DH1",
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
cha_han_DH2_VA = Channel(
    # This is the channel code according to the SEED standard.
    code="DH2",
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
#### BRS data channels DJ1 and DJ2

cha_liv_DJ1_IS = Channel(
    # This is the channel code according to the SEED standard.
    code="DJ1",
    # This is the location code according to the SEED standard.
    location_code="IS",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=0.0,
    sample_rate=256)

cha_liv_DJ2_IW = Channel(
    # This is the channel code according to the SEED standard.
    code="DJ2",
    # This is the location code according to the SEED standard.
    location_code="IW",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=0.0,
    sample_rate=256)

cha_liv_DJ1_ES = Channel(
    # This is the channel code according to the SEED standard.
    code="DJ1",
    # This is the location code according to the SEED standard.
    location_code="ES",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=0.0,
    sample_rate=256)

cha_liv_DJ2_EW = Channel(
    # This is the channel code according to the SEED standard.
    code="DJ2",
    # This is the location code according to the SEED standard.
    location_code="EW",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=0.0,
    sample_rate=256)

##############################################################
#### Seis LIGO data channels DH1, DH2 and DH3 for locations IW, IS, EW and ES

# Location IW
cha_liv_DH1_IW = Channel(
    # This is the channel code according to the SEED standard.
    code="DH1",
    # This is the location code according to the SEED standard.
    location_code="IW",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=???)
cha_liv_DH2_IW = Channel(
    # This is the channel code according to the SEED standard.
    code="DH2",
    # This is the location code according to the SEED standard.
    location_code="IW",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=???)
cha_liv_DH3_IW = Channel(
    # This is the channel code according to the SEED standard.
    code="DH3",
    # This is the location code according to the SEED standard.
    location_code="IW",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=???)

# Location IS
cha_liv_DH1_IS = Channel(
    # This is the channel code according to the SEED standard.
    code="DH1",
    # This is the location code according to the SEED standard.
    location_code="IS",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=???)
cha_liv_DH2_IS = Channel(
    # This is the channel code according to the SEED standard.
    code="DH2",
    # This is the location code according to the SEED standard.
    location_code="IS",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=???)
cha_liv_DH3_IS = Channel(
    # This is the channel code according to the SEED standard.
    code="DH3",
    # This is the location code according to the SEED standard.
    location_code="IS",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=???)

# Location EW
cha_liv_DH1_EW = Channel(
    # This is the channel code according to the SEED standard.
    code="DH1",
    # This is the location code according to the SEED standard.
    location_code="EW",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=???)
cha_liv_DH2_EW = Channel(
    # This is the channel code according to the SEED standard.
    code="DH2",
    # This is the location code according to the SEED standard.
    location_code="EW",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=???)
cha_liv_DH3_EW = Channel(
    # This is the channel code according to the SEED standard.
    code="DH3",
    # This is the location code according to the SEED standard.
    location_code="EW",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=???)

# Location ES
cha_liv_DH1_ES = Channel(
    # This is the channel code according to the SEED standard.
    code="DH1",
    # This is the location code according to the SEED standard.
    location_code="ES",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=???)
cha_liv_DH2_ES = Channel(
    # This is the channel code according to the SEED standard.
    code="DH2",
    # This is the location code according to the SEED standard.
    location_code="ES",
    # Note that these coordinates can differ from the station coordinates.
    latitude=???,
    longitude=???,
    elevation=19,
    depth=0.0,
    azimuth=???,
    dip=???,
    sample_rate=???)
cha_liv_DH3_ES = Channel(
    # This is the channel code according to the SEED standard.
    code="DH3",
    # This is the location code according to the SEED standard.
    location_code="ES",
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
#### Seis LIV Vault data channels DH1 and DH2
# Location VA
cha_liv_DH1_VA = Channel(
    # This is the channel code according to the SEED standard.
    code="DH1",
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
cha_liv_DH2_VA = Channel(
    # This is the channel code according to the SEED standard.
    code="DH2",
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
cha_liv_DH2_VA = Channel(
    # This is the channel code according to the SEED standard.
    code="DH2",
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
    sensor_keys=['Streckeisen', 'STS-1', '360 seconds'],
    datalogger_keys=['REF TEK', 'RT 130 & 130-SMA', '1', '200'])


response_seisvault = nrl.get_response( # doctest: +SKIP
    sensor_keys=['Streckeisen', 'STS-1', '360 seconds'],
    datalogger_keys=['REF TEK', 'RT 130 & 130-SMA', '1', '200'])

response_seisligo = nrl.get_response( # doctest: +SKIP
    sensor_keys=['Streckeisen', 'STS-1', '360 seconds'],
    datalogger_keys=['REF TEK', 'RT 130 & 130-SMA', '1', '200'])


# Now tie it all together.
# add response to channels HAN
cha_han_DJ1_ES.response = response_BRS
cha_han_DJ2_EN.response = response_BRS
cha_han_DH1_EN.response = response_seisligo
cha_han_DH2_EN.response = response_seisligo
cha_han_DH3_EN.response = response_seisligo
cha_han_DH1_ES.response = response_seisligo
cha_han_DH2_ES.response = response_seisligo
cha_han_DH3_ES.response = response_seisligo
cha_han_DH1_VA.response = response_seisvault
cha_han_DH2_VA.response = response_seisvault

# add channels to station Hanford
sta_han.channels.append(cha_han_DJ1_ES)
sta_han.channels.append(cha_han_DJ2_EN)
sta_han.channels.append(cha_han_DH1_EN)
sta_han.channels.append(cha_han_DH2_EN)
sta_han.channels.append(cha_han_DH3_EN)
sta_han.channels.append(cha_han_DH1_ES)
sta_han.channels.append(cha_han_DH2_ES)
sta_han.channels.append(cha_han_DH3_ES)
sta_han.channels.append(cha_han_DH1_VA)
sta_han.channels.append(cha_han_DH2_VA)

# add response to channels LIV
cha_liv_DJ1_IS.response = response_BRS
cha_liv_DJ2_IW.response = response_BRS
cha_liv_DJ1_ES.response = response_BRS
cha_liv_DJ2_EW.response = response_BRS
cha_liv_DH1_IS.response = response_seisligo
cha_liv_DH2_IS.response = response_seisligo
cha_liv_DH3_IS.response = response_seisligo
cha_liv_DH1_IW.response = response_seisligo
cha_liv_DH2_IW.response = response_seisligo
cha_liv_DH3_IW.response = response_seisligo
cha_liv_DH1_ES.response = response_seisligo
cha_liv_DH2_ES.response = response_seisligo
cha_liv_DH3_ES.response = response_seisligo
cha_liv_DH1_EW.response = response_seisligo
cha_liv_DH2_EW.response = response_seisligo
cha_liv_DH3_EW.response = response_seisligo
'''does this instrument exist?
cha_liv_DH1_VA.response = response_seisvault_liv
cha_liv_DH2_VA.response = response_seisvault_liv
cha_liv_DH3_VA.response = response_seisvault_liv'''

# add channels to station Hanford
sta_liv.channels.append(cha_liv_DJ1_IS)
sta_liv.channels.append(cha_liv_DJ2_IW)
sta_liv.channels.append(cha_liv_DJ1_ES)
sta_liv.channels.append(cha_liv_DJ2_EW)
sta_liv.channels.append(cha_liv_DH1_IS)
sta_liv.channels.append(cha_liv_DH2_IS)
sta_liv.channels.append(cha_liv_DH3_IS)
sta_liv.channels.append(cha_liv_DH1_IW)
sta_liv.channels.append(cha_liv_DH2_IW)
sta_liv.channels.append(cha_liv_DH3_IW)
sta_liv.channels.append(cha_liv_DH1_ES)
sta_liv.channels.append(cha_liv_DH2_ES)
sta_liv.channels.append(cha_liv_DH3_ES)
sta_liv.channels.append(cha_liv_DH1_EW)
sta_liv.channels.append(cha_liv_DH2_EW)
sta_liv.channels.append(cha_liv_DH3_EW)
'''does this instrument exist?
sta_liv.channels.append(cha_liv_DH1_VA)
sta_liv.channels.append(cha_liv_DH2_VA)
sta_liv.channels.append(cha_liv_DH3_VA)'''

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