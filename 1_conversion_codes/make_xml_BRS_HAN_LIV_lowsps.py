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
                "internal) of each LIGO arm."+
                "[elevation, HAN, LIV] Althouse, W.E., Hand, S.D., Jones, L.K., Lazzarini, A., Weiss, R.. 2001, Precision " +
                "alignment of the LIGO 4km arms using the dual-frequency differential global positioning. Review" +
                "of scientific instruments." +
                "[lat, lon of HAN] estimated using Shoshana's Phone in 2025" +
                "[azimuth, HAN, LIV] from Shoshana, where exactly???",
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
    elevation=142.555, # done
    creation_date=obspy.UTCDateTime(2016, 1, 2),
    site=Site(name="Hanford LIGO Site"))

##############################################################
#### BRS data channels DJX and DJY

cha_han_BJX_EY = Channel(
    # This is the channel code according to the SEED standard.
    code="BJX",
    # This is the location code according to the SEED standard.
    location_code="EY",
    # Note that these coordinates can differ from the station coordinates.
    latitude=46.43388683, # Shoshana Phone, error 6.1m
    longitude=-119.449623, # Shoshana Phone, error= 6.1m
    elevation=142.555, # done
    depth=0.0, # done
    azimuth=324.0007, # done
    dip=0.0, # done
    sample_rate=10) # done
cha_han_BJY_EX = Channel(
    code="BJY",
    location_code="EX",
    latitude=46.48433894, # Shoshana Phone, error= 2.97m
    longitude=-119.438079, # Shoshana Phone, error= 2.97m
    elevation=142.555, # done
    depth=0.0, # done
    azimuth=234.0007, # done
    dip=0.0,  # done
    sample_rate=10) # done

##############################################################
#### Seis LIGO data channels DHX, DHY and DHZ for locations EX and EY, IY

# Location EX
cha_han_BHX_EX = Channel(
    code="BHX",
    location_code="EX",
    latitude=46.48433894, # Shoshana Phone, error= 2.97m
    longitude=-119.438079, # Shoshana Phone, error= 2.97m
    elevation=142.555, # done
    depth=0.0, # done
    azimuth=324.0007, # done
    dip=0, # done
    sample_rate=10) # done
cha_han_BHY_EX = Channel(
    code="BHY",
    location_code="EX",
    latitude=46.48433894, # Shoshana Phone, error= 2.97m
    longitude=-119.438079, # Shoshana Phone, error= 2.97m
    elevation=142.555, # done
    depth=0.0, # done
    azimuth=234.0007, # done
    dip=0, # done
    sample_rate=10) # done
cha_han_BHZ_EX = Channel(
    code="BHZ",
    location_code="EX",
    latitude=46.48433894, # Shoshana Phone, error= 2.97m
    longitude=-119.438079, # Shoshana Phone, error= 2.97m
    elevation=142.555, # done
    depth=0.0, # done
    azimuth=0, # done
    dip=-90, # done
    sample_rate=10) # done

# Location EY
cha_han_BHX_EY = Channel(
    code="BHX",
    location_code="EY",
    latitude=46.43388683, # Shoshana Phone, error 6.1m
    longitude=-119.449623, # Shoshana Phone, error= 6.1m
    elevation=142.555, # done
    depth=0.0, # done
    azimuth=324.0007, # done
    dip=0, # done
    sample_rate=10) # done
cha_han_BHY_EY = Channel(
    code="BHY",
    location_code="EY",
    latitude=46.43388683, # Shoshana Phone, error 6.1m
    longitude=-119.449623, # Shoshana Phone, error= 6.1m
    elevation=142.555, # done
    depth=0.0, # done
    azimuth=234.0007, # done
    dip=0, # done
    sample_rate=10) # done
cha_han_BHZ_EY = Channel(
    code="BHZ",
    location_code="EY",
    latitude=46.43388683, # Shoshana Phone, error 6.1m
    longitude=-119.449623, # Shoshana Phone, error= 6.1m
    elevation=142.555, # done
    depth=0.0, # done
    azimuth=0, # done
    dip=-90, # done
    sample_rate=10) # done

# Paul is estimating these.
'''# Location IY
cha_han_BHX_IY = Channel(
    code="BHX",
    location_code="IY",
    latitude=???, # Shoshana Phone, error 6.1m
    longitude=???, # Shoshana Phone, error= 6.1m
    elevation=142.555, # done
    depth=0.0, # done
    azimuth=324.0007, # done
    dip=0, # done
    sample_rate=10) # done
cha_han_BHY_IY = Channel(
    code="BHY",
    location_code="IY",
    latitude=???, # Shoshana Phone, error 6.1m
    longitude=???, # Shoshana Phone, error= 6.1m
    elevation=142.555, # done
    depth=0.0, # done
    azimuth=234.0007, # done
    dip=0, # done
    sample_rate=10) # done
cha_han_BHZ_IY = Channel(
    code="BHZ",
    location_code="EY",
    latitude=???, #
    longitude=???, # 
    elevation=142.555, # done
    depth=0.0, # done
    azimuth=0, # done
    dip=-90, # done
    sample_rate=10) # done'''

##############################################################
#### Seis HAN Vault data channels DHX and DHY
# Location VA
cha_han_BHX_VA = Channel(
    code="BHX",
    location_code="VA",
    latitude=46.461704,
    longitude=-119.417683,
    elevation=157, # done
    depth=0,
    azimuth=324.0007,  # done
    dip=0,  # done
    sample_rate=10) # done
cha_han_BHY_VA = Channel(
    code="BHY",
    location_code="VA",
    latitude=46.461704,
    longitude=-119.417683,
    elevation=157, # done
    depth=0,
    azimuth=234.0007,  # done
    dip=0,  # done
    sample_rate=10) # done

########################################################################################################################
################################################ LIGO Livingston Site #####################################################
########################################################################################################################
sta_liv = Station(
    code="LIV",
    latitude=30.56287,
    longitude=-90.77422,
    elevation=-6.574,
    creation_date=obspy.UTCDateTime(2018, 1, 2),
    site=Site(name="Livingston LIGO Site, "))

##############################################################
#### BRS data channels DJX and DJY

cha_liv_BJX_IY = Channel(
    code="BJX",
    location_code="IY",
    latitude=30.5630716, # done, Michael
    longitude=-90.7743621, # done, Michael
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=252.2835, # done
    dip=0.0, # done
    sample_rate=10) # done
cha_liv_BJY_IX = Channel(
    code="BJY",
    location_code="IX",
    latitude=30.5630722, # done, Michael
    longitude=-90.7743672, # done, Michael
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=162.2835, # done
    dip=0.0, # done
    sample_rate=10) # done
cha_liv_BJX_EY = Channel(
    code="BJX",
    location_code="EY",
    latitude=30.5286229, # done, Michael
    longitude=-90.7613708, # done, Michael
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=252.2835, # done
    dip=0.0, # done
    sample_rate=10) # done
cha_liv_BJY_EX = Channel(
    code="BJY",
    location_code="EX",
    latitude=30.5521257, # done, Michael
    longitude=-90.8140179, # done, Michael
    elevation=-6.574, # done
    depth=0.0, # done # done
    azimuth=162.2835, # done
    dip=0.0, # done
    sample_rate=10) # done

##############################################################
#### Seis LIGO data channels DHX, DHY and DHZ for locations IW, IS, EW and ES

# Location IW
cha_liv_BHX_IX = Channel(
    code="BHX",
    location_code="IX",
    latitude=30.5630722, # done, Michael
    longitude=-90.7743672, # done, Michael
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=252.2835, # done
    dip=0, # done
    sample_rate=10) # done
cha_liv_BHY_IX = Channel(
    code="BHY",
    location_code="IX",
    latitude=30.5630722, # done, Michael
    longitude=-90.7743672, # done, Michael
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=162.2835, # done
    dip=0, # done
    sample_rate=10) # done
cha_liv_BHZ_IX = Channel(
    code="BHZ",
    location_code="IX",
    latitude=30.5630722, # done, Michael
    longitude=-90.7743672, # done, Michael
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=0, # done
    dip=-90, # done
    sample_rate=10) # done

# Location IY
cha_liv_BHX_IY = Channel(
    code="BHX",
    location_code="IY",
    latitude=30.5630716, # done, Michael
    longitude=-90.7743621, # done, Michael
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=252.2835, # done
    dip=0, # done
    sample_rate=10) # done
cha_liv_BHY_IY = Channel(
    code="BHY",
    location_code="IY",
    latitude=30.5630716, # done, Michael
    longitude=-90.7743621, # done, Michael
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=162.2835, # done
    dip=0, # done
    sample_rate=10) # done
cha_liv_BHZ_IY = Channel(
    code="BHZ",
    location_code="IY",
    latitude=30.5630716, # done, Michael
    longitude=-90.7743621, # done, Michael
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=0, # done
    dip=-90, # done
    sample_rate=10) # done

# Location EX
cha_liv_BHX_EX = Channel(
    code="BHX",
    location_code="EX",
    latitude=30.5521257, # done, Michael
    longitude=-90.8140179, # done, Michael
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=252.2835, # done
    dip=0, # done
    sample_rate=10) # done
cha_liv_BHY_EX = Channel(
    code="BHY",
    location_code="EX",
    latitude=30.5521257, # done, Michael
    longitude=-90.8140179, # done, Michael
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=162.2835, # done
    dip=0, # done
    sample_rate=10) # done
cha_liv_BHZ_EX = Channel(
    code="BHZ",
    location_code="EX",
    latitude=30.5521257, # done, Michael
    longitude=-90.8140179, # done, Michael
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=0, # done
    dip=-90, # done
    sample_rate=10) # done

# Location EY
cha_liv_BHX_EY = Channel(
    code="BHX",
    location_code="EY",
    latitude=30.5286229, # done, Michael
    longitude=-90.7613708, # done, Michael
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=252.2835, # done
    dip=0, # done
    sample_rate=10) # done
cha_liv_BHY_EY = Channel(
    code="BHY",
    location_code="EY",
    latitude=30.5286229, # done, Michael
    longitude=-90.7613708, # done, Michael
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=162.2835, # done
    dip=0, # done
    sample_rate=10) # done
cha_liv_BHZ_EY = Channel(
    code="BHZ",
    location_code="EY",
    latitude=30.5286229, # done, Michael
    longitude=-90.7613708, # done, Michael
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=0, # done
    dip=-90, # done
    sample_rate=10) # done

'''
Does this station even exist?
##############################################################
#### Seis LIV Vault data channels DHX and DHY
# Location VA
cha_liv_BHX_VA = Channel(
    code="BHX",
    location_code="VA",
    latitude=???,
    longitude=???,
    elevation=???,
    depth=???,
    azimuth=???,
    dip=???,
    sample_rate=???)
cha_liv_BHY_VA = Channel(
    code="BHY",
    location_code="VA",
    latitude=???,
    longitude=???,
    elevation=???,
    depth=???,
    azimuth=???,
    dip=???,
    sample_rate=???)
cha_liv_BHY_VA = Channel(
    code="BHY",
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
cha_han_BJX_EY.response = response_BRS
cha_han_BJY_EX.response = response_BRS

cha_han_BHX_EX.response = response_seisligo
cha_han_BHY_EX.response = response_seisligo
cha_han_BHZ_EX.response = response_seisligo
cha_han_BHX_EY.response = response_seisligo
cha_han_BHY_EY.response = response_seisligo
cha_han_BHZ_EY.response = response_seisligo
cha_han_BHX_VA.response = response_seisligo
cha_han_BHY_VA.response = response_seisligo

# add channels to station Hanford
sta_han.channels.append(cha_han_BJX_EY)
sta_han.channels.append(cha_han_BJY_EX)
sta_han.channels.append(cha_han_BHX_EX)
sta_han.channels.append(cha_han_BHY_EX)
sta_han.channels.append(cha_han_BHZ_EX)
sta_han.channels.append(cha_han_BHX_EY)
sta_han.channels.append(cha_han_BHY_EY)
sta_han.channels.append(cha_han_BHZ_EY)
sta_han.channels.append(cha_han_BHX_VA)
sta_han.channels.append(cha_han_BHY_VA)

# add response to channels LIV
cha_liv_BJX_IY.response = response_BRS
cha_liv_BJY_IX.response = response_BRS
cha_liv_BJX_EY.response = response_BRS
cha_liv_BJY_EX.response = response_BRS
cha_liv_BHX_IY.response = response_seisligo
cha_liv_BHY_IY.response = response_seisligo
cha_liv_BHZ_IY.response = response_seisligo
cha_liv_BHX_IX.response = response_seisligo
cha_liv_BHY_IX.response = response_seisligo
cha_liv_BHZ_IX.response = response_seisligo
cha_liv_BHX_EY.response = response_seisligo
cha_liv_BHY_EY.response = response_seisligo
cha_liv_BHZ_EY.response = response_seisligo
cha_liv_BHX_EX.response = response_seisligo
cha_liv_BHY_EX.response = response_seisligo
cha_liv_BHZ_EX.response = response_seisligo

# add channels to station Hanford
sta_liv.channels.append(cha_liv_BJX_IY)
sta_liv.channels.append(cha_liv_BJY_IX)
sta_liv.channels.append(cha_liv_BJX_EY)
sta_liv.channels.append(cha_liv_BJY_EX)
sta_liv.channels.append(cha_liv_BHX_IY)
sta_liv.channels.append(cha_liv_BHY_IY)
sta_liv.channels.append(cha_liv_BHZ_IY)
sta_liv.channels.append(cha_liv_BHX_IX)
sta_liv.channels.append(cha_liv_BHY_IX)
sta_liv.channels.append(cha_liv_BHZ_IX)
sta_liv.channels.append(cha_liv_BHX_EY)
sta_liv.channels.append(cha_liv_BHY_EY)
sta_liv.channels.append(cha_liv_BHZ_EY)
sta_liv.channels.append(cha_liv_BHX_EX)
sta_liv.channels.append(cha_liv_BHY_EX)
sta_liv.channels.append(cha_liv_BHZ_EX)

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
inv.write("%s/LIGO_HAN_LIV_lowsps.xml" %root_xml, format="stationxml", validate=True)
inv.write("LIGO_HAN_LIV_lowsps.xml", format="stationxml", validate=True)