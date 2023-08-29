
# Introduction

A customizable script to acquire data [waveforms, earthquake catalogs, phase picks] using FDSN Web Services and the Client provided by [obspy](https://docs.obspy.org/).

# Requirements

* Python 3.X
* Pandas 
* Obspy

# How to

## Introduction

The [script](https://github.com/mmesim/fdsnws-download/blob/main/fdsnws-download.py) is structured in a way that it is highly customizable in order to provide maximum flexibility to the user. 

**FDSNWS URL**: the FDSNWS address from which the data will be requested can be defined with credentials( http://user:password@myfdsnws.somewhere:8080) or without it (http://myfdsnws.somewhere:8080).

**Other parameters**: In the different functions the user can define criteria to filter the earthquake catalog, and to define the window length of the requested waveforms (*extratime*). The format type of the data is also flexible and can be defined in the script (e.g **SAC** instead of **MSEED**). 

## Earthquake catalog

To download the events occurred during the defined time period and store them in **csv format** to catalog.csv, run the following:

<pre>
python fdsnws-download.py 'http://myfdsnws:8080' 2023-04-19T12:00:00 2023-04-19T12:03:00 \
    > catalog.csv
</pre>

## Phase picks

While a csv file is easy to handle, some details cannot be stored in that format (e.g. phase picks). For this reason it is possible to specify a **destination folder**, where each **event** is stored **in [QUAKEML](https://quake.ethz.ch/quakeml/)** format along with the station inventory in the same format:

<pre>
python fdsnws-download.py 'http://myfdsnws:8080' 2023-04-19T12:00:00 2023-04-19T12:03:00 \
    output-catalog-dir > catalog.csv
</pre>

*Note*: Replace **output-catalog-dir**  with a folder name that makes sense to you. 

## Waveform data

 It is possible to **download the waveforms** too. The script loads the previously downloaded catalog data (csv and QUAKEML files) and then it downloads the waveforms for each event. By default it fetches the waveforms of the stations associated to the event picks and the latest pick time is used to determine the waveform length (oring time ~ latest pick time):

<pre>
 python fdsnws-download.py 'http://user:password@myfdsnws:8080' --waveforms catalog-dir catalog.csv
</pre>

*Note*: Replace **catalog-dir** and **catalog.csv** with the folder name and the csv file downloaded previously.

Additionally it is possible to manually specify the length of the waveforms to download and the list of stations to use:

<pre>
 python fdsnws-download.py --waveforms catalog-dir catalog.csv [length] [stations_list]
</pre>

*Note*:
- **length** is in seconds **e.g.** 3.5, 0.030
- **stations_list** is a comma separated list of station filters where each filter can be: "net", "net.sta", "net.sta.loc" or "net.sta.loc.cha" and it supports wildcards such as *,?,(,),| **e.g.** "CH,GR.(STA01|STA02).\*.HH?" - download all stations of the network CH (identical to CH.\*.\*.\*) plus the stations GT.STA01.\*.HH? and GT.STA02.\*.HH?

# Post-processing

Here we demonstrate how to use the csv catalog and the events stored in QUAKEML format. In this example each event is loaded and the corresponding arrivals are printed. The waveforms are loaded too and plotted.

```python
import pandas as pd
import obspy as ob
from pathlib import Path

csv_catalog = "catalog.csv"
xml_folder  = "mydirectory"  # folder where the event XML files were downloaded

### Load the csv catalog
cat = pd.read_csv(csv_catalog, na_filter=False)

#
# One interesting thing of using pandas to read the csv file is that you can
# now filter the events as you please, e.g.
#
#   cat = cat[ cat["evaluation_mode"] == "manual" ] # select only manual events
#
#   cat = cat[ cat["num_phase"] >= 8 ] # at least 8 picks
#
#   cat = cat[ cat["depth"] > -1400 & cat["depth"] < -1200 ] # events between 1200~1400 meters
#

#
# Loop through the catalog by event
#
for row in cat.itertuples():
  #
  # you can fetch any csv column with 'row.column'
  #
  ev_id = row.id
  print(f"Processing event {ev_id}")

  #
  # We want to access all event information, so we load its XML file
  # with obspy
  #
  cat = ob.core.event.read_events( Path(xml_folder, f"ev{ev_id}.xml"))
  ev= cat[0] # we know we stored only one event in this catalog

  #
  # We can now have access to all events information e.g. access picks
  #

  # get preferred origin
  o = ev.preferred_origin()
  print(f"Origin {o.longitude} {o.latitude} {o.depth} {o.time}")

  # loop trough origin arrivals
  for a in o.arrivals:
    #
    # find the pick associated with the current arrival
    #
    for p in ev.picks:
      if p.resource_id == a.pick_id:
        wfid = p.waveform_id
        print(f"Found pick {p.time} @ {wfid.network_code}.{wfid.station_code}.{wfid.location_code}.{wfid.channel_code}")
        break

  #
  # We can also load the waveforms for this event
  #
  trace = ob.read( Path(xml_folder, f"ev{ev_id}.mseed"))
  trace.plot()

```

Similarly it is possible to load the inventory XML with obspy and use the [Inventory class](https://docs.obspy.org/packages/autogen/obspy.core.inventory.inventory.Inventory.html) to fetch the required information.

```python
import obspy as ob
from pathlib import Path

xml_folder  = "mydirectory"  # folder where the event XML files were downloaded

### Load the previously downloaded inventory
inv = ob.core.inventory.inventory.read_inventory( Path(xml_folder, "inventory.xml") )

### Print the stations active at a specific point in time
ref_time = ob.UTCDateTime("2022-10-01 12:00:00")

stations = {}
for net in inv.networks:
    if ref_time < net.start_date or ref_time > net.end_date:
        continue
    for sta in net.stations:
        if ref_time < sta.start_date or ref_time > sta.end_date:
            continue
        for cha in sta.channels:
            if ref_time < cha.start_date or ref_time > cha.end_date:
                continue
            wfid = f"{net.code}.{sta.code}.{cha.location_code}.{cha.code}"
            stations[wfid] = (cha.latitude, cha.longitude, cha.elevation, cha.azimuth, cha.dip, cha.sample_rate, cha.description, cha.comments)

for (sta, info) in stations.items():
    print(f"{sta}: {info}")

#
# sometimes it might be easier to fetch the coordinates like this
#
inv.get_coordinates("NET.STA.LOC.CHA", ob.UTCDateTime("2022-10-01"))
```







