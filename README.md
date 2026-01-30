
# Introduction

This script utilizes FDSN Web Services and the ObsPy library to acquire seismological data, including waveforms, earthquake catalogs, and phase picks.

# Requirements

* Python 3.X
* Pandas 
* Obspy

# How to

## How to Use

### Getting Started

To use this script, simply download the `fdsnws-download.py` file. You do not need to install it as a Python package. The script is designed for direct execution and can be customized by editing its source code to extend its functionalities.

**FDSNWS URL**: This is the address of the FDSN Web Service. It can include credentials (e.g., `http://user:password@myfdsnws.somewhere:8080`) or be without them (e.g., `http://myfdsnws.somewhere:8080`). Replace `[FDSNWS_URL]` in the examples below with your specific service URL.

## Earthquake catalog

To download events within a defined time period and store them in **CSV format** to `catalog.csv`, run the following command:

```bash
python fdsnws-download.py '[FDSNWS_URL]' '2023-04-19T12:00:00' '2023-04-19T12:03:00' \
    > catalog.csv
```

The `> catalog.csv` part of the command redirects the script's standard output, which contains the catalog data, into a file named `catalog.csv`.

## Phase picks

While a csv file is easy to handle, some details cannot be stored in that format (e.g. phase picks). For this reason it is possible to specify a **destination folder**, where each **event** is stored **in [QUAKEML](https://quake.ethz.ch/quakeml/)** format along with the station inventory in the same format:

```bash
python fdsnws-download.py '[FDSNWS_URL]' '2023-04-19T12:00:00' '2023-04-19T12:03:00' \
    output-catalog-dir
```

*Note*: Replace `output-catalog-dir` with a descriptive folder name. This directory will store each event as a QUAKEML file and its corresponding station inventory. 

## Event Waveform data

 It is possible to **download event waveforms**. The script loads the previously downloaded catalog data (CSV and QUAKEML files) and then downloads waveforms for each event. By default, it fetches waveforms from stations associated with event picks, and the latest pick time determines the waveform length (origin time ~ latest pick time):

```bash
python fdsnws-download.py '[FDSNWS_URL]' --event-waveforms catalog-dir catalog.csv
```

*Note*: Replace `catalog-dir` with the name of the directory containing the QUAKEML files (e.g., `output-catalog-dir` from the previous step) and `catalog.csv` with the path to the CSV catalog file downloaded earlier.

Additionally it is possible to manually specify the length of the waveforms to download and the list of stations to use:

```bash
python fdsnws-download.py '[FDSNWS_URL]' --waveforms catalog-dir catalog.csv \
    [length-before-event:length-after-event] [station-filter]
```

*Note*:
- **length-before-event** waveform length in seconds (e.g. 3.5, 0.030) to download before the event time
- **length-after-event** waveform length in seconds (e.g. 3.5, 0.030) to download after the event time
- **stations-filter** is an optional comma separated list of station filters where each filter can be: "net", "net.sta", "net.sta.loc" or "net.sta.loc.cha" and it supports wildcards such as *,?,(,),| 

e.g. Download 3 seconds of waveforms before the event and 10 after and download all stations of the network CH (identical to CH.\*.\*.\*) plus the stations GT.STA01.\*.HH? and GT.STA02.\*.HH?

```bash
python fdsnws-download.py '[FDSNWS_URL]' --waveforms catalog-dir catalog.csv \
    "3:10" "CH,GR.(STA01|STA02).*.HH?"
```

## Waveform data at any point in time

To download waveforms at specific points in time, you can use the following syntax

```bash
python fdsnws-download.py '[FDSNWS_URL]' --waveforms inventory.xml [start-time] [length] [station-filter]
```

*Note*:
- **inventory.xml**: This refers to an FDSN StationXML file, which can be generated when downloading phase picks (e.g., in `output-catalog-dir`).
- **start-time**: waveform start time to download
- **length** waveform length in seconds
- **stations-filter** is an optional comma separated list of station filters where each filter can be: "net", "net.sta", "net.sta.loc" or "net.sta.loc.cha" and it supports wildcards such as *,?,(,),|

E.g. Download 20 seconds of data of any station at 2024-04-30T00:00:00:

```bash
python fdsnws-download.py '[FDSNWS_URL]' --waveforms catalog-dir/inventory.xml "2024-04-30T00:00:00" 20
```

Only station NN.STA1 and NN.STA2

```bash
python fdsnws-download.py '[FDSNWS_URL]' --waveforms catalog-dir/inventory.xml "2024-04-30T00:00:00" 20 "NN.STA1,NN.STA2"
```

# Post-processing

Here we demonstrate how to use the csv catalog and the events stored in QUAKEML format. In this example each event is loaded and the corresponding arrivals are printed. The waveforms are loaded too and plotted.

```python
import pandas as pd
import obspy as ob
from pathlib import Path
from obspy.core import UTCDateTime

csv_catalog = "catalog.csv"
xml_folder  = "output-catalog-dir"  # Folder where event QUAKEML files and inventory were downloaded

#
# Load the csv catalog
#
cat = pd.read_csv(csv_catalog, na_filter=False, parse_dates=False, converters={"isotime": lambda d: UTCDateTime(d)})

# Check the column types
print(cat.dtypes)

#
# One interesting thing of using pandas to read the csv file is that you can
# now filter the events as you please, e.g.
#
#   cat = cat[ cat["evaluation_mode"] == "manual" ] # select only manual events
#
#   cat = cat[ cat["num_phase"] >= 8 ] # at least 8 picks
#
#   cat = cat[ (cat["depth"] > -1.500) & (cat["depth"] < -1.200) ] # events between 1200~1400 meters
#
#   cat[ (cat.isotime > UTCDateTime("2024-04-30T05:00:00.369")) & (cat.isotime < UTCDateTime("2024-04-30T05:22:00.2369")) ]
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
  print(f"Origin {o.time} {o.longitude} {o.latitude} {o.depth} {o.evaluation_mode}")

  if o.quality:
      print(f" associated_phase_count {o.quality.associated_phase_count}")
      print(f" used_phase_count {o.quality.used_phase_count}")
      print(f" associated_station_count {o.quality.associated_station_count}")
      print(f" used_station_count {o.quality.used_station_count}")
      print(f" standard_error {o.quality.standard_error}")
      print(f" azimuthal_gap {o.quality.azimuthal_gap}")
      print(f" minimum_distance {o.quality.minimum_distance}")
      print(f" median_distance {o.quality.median_distance}")
      print(f" maximum_distance {o.quality.maximum_distance}")


  # Loop through origin arrivals and find associated picks
  for a in o.arrivals:
    #
    # find the pick associated with the current arrival
    #
    for p in ev.picks:
      if p.resource_id == a.pick_id:
        wfid = p.waveform_id
        print(f" Pick {p.evaluation_mode} {a.phase} @ {wfid.network_code}.{wfid.station_code}.{wfid.location_code}.{wfid.channel_code} residual {a.time_residual} distance {a.distance} deg {p.time}")
        break

  # Load and plot waveforms for this event
  trace = ob.read( Path(xml_folder, f"ev{ev_id}.mseed"))
  trace.plot()

```

Similarly it is possible to load the inventory XML with obspy and use the [Inventory class](https://docs.obspy.org/packages/autogen/obspy.core.inventory.inventory.Inventory.html) to fetch the required information.

```python
import obspy as ob
from pathlib import Path

xml_folder  = "output-catalog-dir"  # Folder where event QUAKEML files and inventory were downloaded

# Load the previously downloaded inventory
inv = ob.core.inventory.inventory.read_inventory( Path(xml_folder, "inventory.xml") )

# Print stations active at a specific point in time
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

# Alternatively, fetch coordinates like this
inv.get_coordinates("NET.STA.LOC.CHA", ob.UTCDateTime("2022-10-01"))
```







