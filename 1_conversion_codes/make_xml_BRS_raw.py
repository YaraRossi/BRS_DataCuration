import obspy
from obspy import read_inventory
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
from obspy.clients.nrl import NRL
from obspy.core.inventory.response import PolesZerosResponseStage, Response

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

cha_han_DJX_EY = Channel(
    # This is the channel code according to the SEED standard.
    code="DJX",
    # This is the location code according to the SEED standard.
    location_code="EY",
    # Note that these coordinates can differ from the station coordinates.
    latitude=46.43388683, # Shoshana Phone, error 6.1m
    longitude=-119.449623, # Shoshana Phone, error= 6.1m
    elevation=142.555, # done
    depth=0.0, # done
    azimuth=324.0007, # done
    dip=0.0, # done
    sample_rate=256) # done
cha_han_DJY_EX = Channel(
    code="DJY",
    location_code="EX",
    latitude=46.48433894, # Shoshana Phone, error= 2.97m
    longitude=-119.438079, # Shoshana Phone, error= 2.97m
    elevation=142.555, # done
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
    latitude=46.48433894, # Shoshana Phone, error= 2.97m
    longitude=-119.438079, # Shoshana Phone, error= 2.97m
    elevation=142.555, # done
    depth=0.0, # done
    azimuth=324.0007, # done
    dip=0, # done
    sample_rate=512) # done
cha_han_DHY_EX = Channel(
    code="DHY",
    location_code="EX",
    latitude=46.48433894, # Shoshana Phone, error= 2.97m
    longitude=-119.438079, # Shoshana Phone, error= 2.97m
    elevation=142.555, # done
    depth=0.0, # done
    azimuth=234.0007, # done
    dip=0, # done
    sample_rate=512) # done
cha_han_DHZ_EX = Channel(
    code="DHZ",
    location_code="EX",
    latitude=46.48433894, # Shoshana Phone, error= 2.97m
    longitude=-119.438079, # Shoshana Phone, error= 2.97m
    elevation=142.555, # done
    depth=0.0, # done
    azimuth=0, # done
    dip=-90, # done
    sample_rate=512) # done

# Location EY
cha_han_DHX_EY = Channel(
    code="DHX",
    location_code="EY",
    latitude=46.43388683, # Shoshana Phone, error 6.1m
    longitude=-119.449623, # Shoshana Phone, error= 6.1m
    elevation=142.555, # done
    depth=0.0, # done
    azimuth=324.0007, # done
    dip=0, # done
    sample_rate=512) # done
cha_han_DHY_EY = Channel(
    code="DHY",
    location_code="EY",
    latitude=46.43388683, # Shoshana Phone, error 6.1m
    longitude=-119.449623, # Shoshana Phone, error= 6.1m
    elevation=142.555, # done
    depth=0.0, # done
    azimuth=234.0007, # done
    dip=0, # done
    sample_rate=512) # done
cha_han_DHZ_EY = Channel(
    code="DHZ",
    location_code="EY",
    latitude=46.43388683, # Shoshana Phone, error 6.1m
    longitude=-119.449623, # Shoshana Phone, error= 6.1m
    elevation=142.555, # done
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
    elevation=142.555, # done
    depth=0,
    azimuth=324.0007,  # done
    dip=0,  # done
    sample_rate=256)
cha_han_DHY_VA = Channel(
    code="DHY",
    location_code="VA",
    latitude=46.461704,
    longitude=-119.417683,
    elevation=142.555, # done
    depth=0,
    azimuth=234.0007,  # done
    dip=0,  # done
    sample_rate=256)

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

cha_liv_DJX_IY = Channel(
    code="DJX",
    location_code="IY",
    latitude=30.5630716, # done, Michael
    longitude=-90.7743621, # done, Michael
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=252.2835, # done
    dip=0.0, # done
    sample_rate=256) # done
cha_liv_DJY_IX = Channel(
    code="DJY",
    location_code="IX",
    latitude=30.5630722, # done, Michael
    longitude=-90.7743672, # done, Michael
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=162.2835, # done
    dip=0.0, # done
    sample_rate=256) # done
cha_liv_DJX_EY = Channel(
    code="DJX",
    location_code="EY",
    latitude=0,
    longitude=0,
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=252.2835, # done
    dip=0.0, # done
    sample_rate=256) # done
cha_liv_DJY_EX = Channel(
    code="DJY",
    location_code="EX",
    latitude=0,
    longitude=0,
    elevation=-6.574, # done
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
    latitude=30.5630722, # done, Michael
    longitude=-90.7743672, # done, Michael
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=252.2835, # done
    dip=0, # done
    sample_rate=0) # done
cha_liv_DHY_IX = Channel(
    code="DHY",
    location_code="IX",
    latitude=30.5630722, # done, Michael
    longitude=-90.7743672, # done, Michael
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=162.2835, # done
    dip=0, # done
    sample_rate=0) # done
cha_liv_DHZ_IX = Channel(
    code="DHZ",
    location_code="IX",
    latitude=30.5630722, # done, Michael
    longitude=-90.7743672, # done, Michael
    elevation=-6.574, # done
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
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=252.2835, # done
    dip=0, # done
    sample_rate=0) # done
cha_liv_DHY_IY = Channel(
    code="DHY",
    location_code="IY",
    latitude=0,
    longitude=0,
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=162.2835, # done
    dip=0, # done
    sample_rate=0) # done
cha_liv_DHZ_IY = Channel(
    code="DHZ",
    location_code="IY",
    latitude=0,
    longitude=0,
    elevation=-6.574, # done
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
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=252.2835, # done
    dip=0, # done
    sample_rate=0) # done
cha_liv_DHY_EX = Channel(
    code="DHY",
    location_code="EX",
    latitude=0,
    longitude=0,
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=162.2835, # done
    dip=0, # done
    sample_rate=0) # done
cha_liv_DHZ_EX = Channel(
    code="DHZ",
    location_code="EX",
    latitude=0,
    longitude=0,
    elevation=-6.574, # done
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
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=252.2835, # done
    dip=0, # done
    sample_rate=0) # done
cha_liv_DHY_EY = Channel(
    code="DHY",
    location_code="EY",
    latitude=0,
    longitude=0,
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=162.2835, # done
    dip=0, # done
    sample_rate=0) # done
cha_liv_DHZ_EY = Channel(
    code="DHZ",
    location_code="EY",
    latitude=0,
    longitude=0,
    elevation=-6.574, # done
    depth=0.0, # done
    azimuth=0, # done
    dip=-90, # done
    sample_rate=0) # done

# By default this accesses the NRL online. Offline copies of the NRL can
# also be used instead

# The contents of the NRL can be explored interactively in a Python prompt,
# see API documentation of NRL submodule:
# http://docs.obspy.org/packages/obspy.clients.nrl.html
# Here we assume that the end point of data logger and sensor are already
# known:

brs_response = read_inventory('%s/LIGO_BRS_generic_digitizer.xml' %root_xml)
network = brs_response[0]  # First network
station = network.stations[0]  # First station
channel = station.channels[0]  # First channel
response_BRS = channel.response

# add poles and zeros to response
zeros = [complex(0.000121284, 0.00358795), complex(0.000121284, -0.00358795)]
poles = [complex(0.0005, 0.000866025), complex(0.0005, -0.000866025)]
gain = 1.54362

pz_stage = PolesZerosResponseStage(
    stage_sequence_number=1,
    stage_gain=gain,
    stage_gain_frequency=0.02,
    pz_transfer_function_type='LAPLACE (RADIANS/SECOND)',
    normalization_factor=gain,
    normalization_frequency=1.0,
    zeros=zeros,
    poles=poles,
    input_units="m/s",       # or your correct input units
    output_units="V"         # or your correct output units
)
digitizer_stages = [stage for stage in channel.response.response_stages if not isinstance(stage, PolesZerosResponseStage)]
channel.response = Response(response_stages=[pz_stage] + digitizer_stages)
response_BRS = channel.response


################## Seismometer response #########################
sts_response = read_inventory('%s/LIGO_STS2_generic_digitizer.xml' %root_xml)
# this response file will turn NM/S into M/S
network = sts_response[0]  # First network
station = network.stations[0]  # First station
channel = station.channels[0]  # First channel
response_seisligo = channel.response


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
inv.write("LIGO_HAN_LIV.xml", format="stationxml", validate=True)