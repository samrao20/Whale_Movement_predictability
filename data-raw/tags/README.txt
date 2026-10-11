README

This data package is published by the Movebank Data Repository (www.datarepository.movebank.org). As of the time of publication, a version of this published animal tracking dataset can be viewed on Movebank (www.movebank.org) in the study "Movements of Australia's east coast humpback whales" (Movebank Study ID 3030068329). Individual attributes in the data files are defined below, in the NERC Vocabulary Server at http://vocab.nerc.ac.uk/collection/MVB and in the Movebank Attribute Dictionary at www.movebank.org/cms/movebank-content/movebank-attribute-dictionary. Metadata describing this data package are maintained at https://datacite.org.

This data package includes the following data files:
Movements of Australia's east coast humpback whales.csv
Movements of Australia's east coast humpback whales-reference-data.csv

Data package citation:
Andrews-Goff V, Gales N, Childerhouse S, Laverick S, Polanowski A, Double M. 2023. Data from: Australia’s east coast humpback whales: satellite tag derived movements on breeding grounds, feeding grounds and along the northern and southern migration. Movebank Data Repository. https://doi.org/10.5441/001/1.294

These data are described in the following written publication:
Andrews-Goff V, Gales N, Childerhouse SJ, Laverick SM, Polanowski AM, Double MC. 2023. Australia’s east coast humpback whales: satellite tag derived movements on breeding grounds, feeding grounds and along the northern and southern migration. Biodivers Data J. 11:e114729. https://doi.org/10.3897/BDJ.11.e114729

-----------

Terms of Use
This data file is licensed by the Creative Commons Zero (CC0 1.0) license. The intent of this license is to facilitate the re-use of works. The Creative Commons Zero license is a "no rights reserved" license that allows copyright holders to opt out of copyright protections automatically extended by copyright and other laws, thus placing works in the public domain with as little legal restriction as possible. However, works published with this license must still be appropriately cited following professional and ethical standards for academic citation.

We highly recommend that you contact the data creator if possible if you will be re-using or re-analyzing data in this file. Researchers will likely be interested in learning about new uses of their data, might also have important insights about how to properly analyze and interpret their data, and/or might have additional data they would be willing to contribute to your project. Feel free to contact us at support@movebank.org if you need assistance contacting data owners.

For additional information, see

License description: http://creativecommons.org/publicdomain/zero/1.0

General Movebank Terms of Use: www.movebank.org/cms/movebank-content/general-movebank-terms-of-use

Movebank user agreement: https://www.movebank.org/cms/movebank-content/data-policy#user_agreement

-----------

Data Attributes
These definitions come from the Movebank Attribute Dictionary, available at http://vocab.nerc.ac.uk/collection/MVB and www.movebank.org/node/2381.

algorithm marked outlier: Identifies events marked as outliers using a user-selected filter algorithm in Movebank. Outliers have the value TRUE. Information about how outliers were defined can be provided in 'outlier comments'. Units: none
Entity described: event
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000011

animal group ID: The name or identifier for an associated group, for example, the herd, pack, nest or family group. If the group changes over time, this can be recorded using the event term 'event group ID'. Units: none
Entity described: individual
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000340

animal ID: An individual identifier for the animal, provided by the data owner. Values are unique within the study. If the data owner does not provide an Animal ID, an internal Movebank animal identifier is sometimes shown.
Example: TUSC_CV5
Units: none
Entity described: individual
Same as: individual local identifier
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000016

animal life stage: The age class or life stage of the animal at the beginning of the deployment. Can be years or months of age or terms such as 'adult', 'subadult' and 'juvenile'. Best practice is to define units in the values if needed (e.g. '2 years').
Example: juvenile, adult
Units: none
Entity described: deployment
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000018

animal sex: The sex of the animal. Allowed values are
m = male;
f = female;
u = unknown.
Format: controlled list
Entity described: individual
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000023

animal taxon: The scientific name of the taxon on which the tag was deployed, as defined by the Integrated Taxonomic Information System (ITIS). If the species name can not be provided, this should be the lowest level taxonomic rank that can be determined and that is used in the ITIS taxonomy. Additional information can be provided using the term 'taxon detail'. The values 'test' and 'calibration' identify events relevant to animal tracking studies that should not be associated with a taxon.
Format: controlled list
Entity described: individual
Same as: individual taxon canonical name
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000024

Argos lat1: Solution 1. platform latitude in degrees and thousandths of degrees, Argos diagnostic data (definition from Argos User's Manual 2011).
Example: 19.493
Units: decimal degrees, WGS84 reference system
Entity described: event
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000032

Argos lat2: Solution 2. platform latitude in degrees and thousandths of degrees, Argos diagnostic data (definition from Argos User's Manual 2011).
Example: 14.773
Units: decimal degrees, WGS84 reference system
Entity described: event
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000033

Argos LC: The location class retrieved from Argos, Argos diagnostic data. Classes are based on the type of location (Argos Doppler Shift or GPS) and the number of messages received during the satellite pass. Location classes in order of decreasing accuracy are G (GPS), 3, 2, 1, 0, A, B and Z (definition from Argos User's Manual V1.6.6, 2016).
Format: controlled list
Entity described: event
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000034

Argos location algorithm: The processing algorithm used by Argos to estimate locations using Doppler shift. If the location data represent model output rather than the original estimates from Argos, also use 'modelled'. Values are chosen from a controlled list:
least squares = locations were calculated by Argos using a least-squares analysis;
Kalman = locations were calculated by Argos using Kalman filtering.
Format: controlled list
Entity described: event
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000278

Argos lon1: Solution 1. platform longitude in degrees and thousandths of degrees, Argos diagnostic data (definition from Argos User's Manual 2011).
Example: 99.712
Units: decimal degrees, WGS84 reference system
Entity described: event
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000035

Argos lon2: Solution 2. platform longitude in degrees and thousandths of degrees, Argos diagnostic data (definition from Argos User's Manual 2011).
Example: 120.286
Units: decimal degrees, WGS84 reference system
Entity described: event
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000036

attachment type: The way a tag is attached to an animal. Details can be provided in 'attachment comments'. Values are chosen from a controlled list:
backpack-harness = The tag is attached to the animal using a backpack-style harness;
collar = The tag is attached by a collar around the animal's neck;
ear-tag = The tag is attached to the animal's ear;
fin mount = The tag is attached to the animal's fin;
glue = The tag is attached to the animal using glue;
harness = The tag is attached to the animal using a harness;
implant = The tag is placed under the skin of the animal;
leg-band = The tag is attached as a leg band or ring;
leg-loop-harness = The tag is attached to the animal using a leg-loop-style harness;
none = No tag was attached, e.g., for observations using natural markings;
other = The tag is attached using another method;
subcutaneous-anchor = The tag is attached using one or more anchors attached underneath the animal's skin;
suction-cup = The tag is attached using one or more suction cups;
sutures = The tag is attached by one or more sutures;
tape = The tag is attached to the animal using tape.
Format: controlled list
Entity described: deployment
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000052

comments: Additional information about events that is not described by other event data terms.
Example: we observed the animal foraging (see photo BT25)
Units: none
Entity described: event
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000069

deploy off timestamp: The timestamp when the tag deployment ended. Data records recorded after this day and time are not associated with the animal related to the deployment. Values are typically defined by the data owner, and in some cases are created automatically during data import. Further information can be provided in 'deployment end type' and 'deployment end comments'.
Example: 2009-10-01 12:00:00.000
Format: yyyy-MM-dd HH:mm:ss.SSS
Units: UTC or GPS time
Entity described: deployment
Same as: deploy off date
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000077

deploy on latitude: The geographic latitude of the location where the animal was released. Intended primarily for cases in which the animal release location has higher accuracy than that derived from sensor data.
Example: 27.3516
Units: decimal degrees, WGS84 reference system
Entity described: deployment
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000078

deploy on longitude: The geographic longitude of the location where the animal was released. Intended primarily for cases in which the animal release location has higher accuracy than that derived from sensor data.
Example: -97.3321
Units: decimal degrees, WGS84 reference system
Entity described: deployment
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000079

deploy on timestamp: The timestamp when the tag deployment started. Data records recorded before this day and time are not associated with the animal related to the deployment. Values are typically defined by the data owner, and in some cases are created automatically during data import.
Example: 2008-08-30 18:00:00.000
Format: yyyy-MM-dd HH:mm:ss.SSS
Units: UTC or GPS time
Entity described: deployment
Same as: deploy on date
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000081

deployment ID: A unique identifier for the deployment of a tag on animal, provided by the data owner. If the data owner does not provide a Deployment ID, an internal Movebank deployment identifier may sometimes be shown.
Example: Jane_42818
Units: none
Entity described: deployment
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000085

duty cycle: Remarks associated with the duty cycle of a tag during the deployment, describing the times it is on/off and the frequency at which it transmits or records data. Units and time zones should be defined in the remarks.
Example: 15-min fixes from 8:00-18:00 local time (0:00-10:00 UTC)
Units: none
Entity described: deployment
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000086

event ID: An identifier for the set of values associated with each event, i.e. sensor measurement. A unique event ID is assigned to every time-location or other time-measurement record in Movebank. If multiple measurements are included within a single row of a data file, they will share an event ID. If users import the same sensor measurement to Movebank multiple times, a separate event ID will be assigned to each.
Example: 14328243575
Units: none
Entity described: event
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000103

location lat: The geographic latitude of the location as estimated by the sensor.
Example: -41.0982423
Units: decimal degrees, WGS84 reference system
Entity described: event
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000145

location long: The geographic longitude of the location as estimated by the sensor. Positive values are east of the Greenwich Meridian, negative values are west of it.
Example: -121.1761111
Units: decimal degrees, WGS84 reference system
Entity described: event
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000146

manipulation type: The way in which the animal was manipulated during the deployment. Additional information can be provided using 'manipulation comments'. Changes in manipulation status during deployment can be identified using 'manipulation status'. Values are chosen from a controlled list:
confined = The animal's movement was restricted to within a defined area;
domesticated = The animal is domesticated, for example, is a house pet or part of a managed herd;
manipulated-other = The animal was manipulated in some other way, such as a physiological manipulation;
none = The animal received no treatment other than tag attachment and related measurements and sampling (if applicable);
reintroduction = The animal has been reintroduced as part of wildlife conservation or management efforts;
relocated = The animal was released from a site other than the one at which it was captured.
Format: controlled list
Entity described: deployment
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000155

modelled: Can be used to identify locations or light level values that are modelled, interpolated, or otherwise added or changed based on original sensor measurements. These may include movement model outputs, locations added to create a dataset with equal time intervals, etc. Details about the procedures used can be provided in 'data processing software'. Allowed values are TRUE or FALSE.
Units: none
Entity described: event
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000162

outlier comments: A description or reference for methods used to define outliers in 'manually marked outlier', 'algorithm marked outlier', 'import marked outlier' and/or 'manually marked valid'. If outliers were removed prior to import, this field can be used to describe which data were excluded.
Example: Movebank value range filter set to flag records with "GPS HDOP" => 20, retaining records with null values.
Units: none
Entity described: deployment
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000374

sensor type: The type of sensor with which data were collected. All sensors are associated with a tag id, and tags can contain multiple sensor types. Each event record in Movebank is assigned one sensor type. If values from multiple sensors are reported in a single event, the primary sensor is used. Values are chosen from a controlled list:
acceleration = The sensor collects acceleration data;
accessory-measurements = The sensor collects accessory measurements, such as battery voltage;
acoustic-telemetry = The sensor transmits an acoustic signal that is detected by receivers to determine location;
argos-doppler-shift = The sensor location is estimated by Argos using Doppler shift;
barometer = The sensor records air or water pressure;
bird-ring = The animal is identified by a band or ring that has a unique identifier;
gps = The sensor uses GPS to determine location;
gyroscope = The sensor records angular velocity;
heart-rate = The sensor records or is used to calculate heart rate;
magnetometer = The sensor records the magnetic field;
natural-mark = The animal is identified by a unique natural marking;
orientation = Quaternion components describing the orientation of the tag are derived from accelerometer and gyroscope measurements;
proximity = The sensor identifies proximity to other tags;
radio-transmitter = The sensor transmits a radio signal that is detected by receivers to determine location;
sigfox-geolocation = The sensor location is determined by Sigfox using the received signal strength indicator;
solar-geolocator = The sensor collects light levels, which are used to determine position (for processed locations);
solar-geolocator-raw = The sensor collects light levels, which are used to determine position (for raw light-level measurements);
solar-geolocator-twilight = The sensor collects light levels, which are used to determine position (for twilights calculated from light-level measurements).
Format: controlled list
Entity described: event
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000170

study name: The name of the study in Movebank.
Example: Coyotes, Kays and Bogan, Albany NY
Units: none
Entity described: study
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000173

tag ID: An identifier for the tag, provided by the data owner. Values are unique within the study. If the data owner does not provide a tag ID, an internal Movebank tag identifier may sometimes be shown.
Example: 2342
Units: none
Entity described: tag
Same as: tag local identifier
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000181

tag manufacturer name: The company or person that produced the tag.
Example: Holohil
Units: none
Entity described: tag
Same as: manufacturer name
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000183

tag model: The model of the tag.
Example: T61
Units: none
Entity described: tag
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000185

tag readout method: The way the data are received from the tag. Values are chosen from a controlled list:
ISS = Data are transferred via the International Space Station;
LPWAN = Data are transferred through a low-power wide-area network, such as LoRa or Sigfox;
multiple = Data are acquired using multiple methods;
none = Data are obtained without use of an animal-borne tag, such as by observing a unique marking;
other-wireless = Data are transferred via another form of wireless data transfer, such as a VHF transmitter/receiver;
phone-network = Data are transferred via a phone network, such as GSM or AMPS;
satellite = Data are transferred via satellite;
tag-retrieval = The tag must be physically retrieved in order to obtain the data;
telemetry-network = Data are obtained through a radio or acoustic telemetry network;
Wi-Fi/Bluetooth = Data are transferred via a local Wi-Fi or Bluetooth system.
Format: controlled list
Entity described: deployment
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000188

timestamp: The date and time corresponding to a sensor measurement or an estimate derived from sensor measurements.
Example: 2008-08-14 18:31:00.000
Format: yyyy-MM-dd HH:mm:ss.SSS
Units: UTC or GPS time
Entity described: event
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000200

visible: Determines whether an event is visible on the Movebank map. Allowed values are TRUE or FALSE. Values are calculated automatically, with TRUE indicating the event has not been flagged as an outlier by 'algorithm marked outlier', 'import marked outlier' or 'manually marked outlier', or that the user has overridden the results of these outlier attributes using 'manually marked valid' = TRUE.
Units: none
Entity described: event
URI: http://vocab.nerc.ac.uk/collection/MVB/current/MVB000209

-----------

More Information
For more information about this repository, see www.movebank.org/cms/movebank-content/data-repository, the FAQ at https://www.datarepository.movebank.org/faq, or contact us at support@movebank.org.