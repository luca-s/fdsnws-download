
# Introduction

A customizable script to acquire data [waveforms, earthquake catalogs, phase picks] using FDSN Web Services and the Client provided by [obspy](https://docs.obspy.org/).

# Requirements

* Python 3.X
* Pandas 
* Obspy

# How to

## Introduction

The [script](https://github.com/mmesim/fdsnws-download/blob/main/fdsnws-download.py) is structured in a way that it is highly customizable in order to provide maximum flexibility to the user. 

**Define Client**: The first important step is to define the [Client](https://docs.obspy.org/packages/autogen/obspy.clients.fdsn.client.Client.html) from which the data will be requested. In lines 33-38 the user can define either a Client from the Obspy list [e.g. ETH], or alternatively provide a password protected FDSN web service or a custom FDSN web service that doesn't require credentials. 

**Other parameters**: In the different functions the user can define criteria to filter the earthquake catalog, and to define the window length of the requested waveforms (*extratime*). The format type of the data is also flexible and can be defined in the script (e.g **SAC** instead of **MSEED**). 

## Earthquake catalog

To download the events occurred during the defined time period and store them in **csv format** to catalog.csv, run the following:

<pre>
python fdsnws-download.py "2023-04-19T12:00:00" "2023-04-19T12:03:00" \
    > catalog.csv
</pre>

## Phase picks

While a csv file is easy to handle, some details cannot be stored in that format (e.g. phase picks). For this reason it is possible to specify a **destination folder**, where each **event** is stored **in [QUAKEML](https://quake.ethz.ch/quakeml/)** format along with the station inventory in the same format:

<pre>
python fdsnws-download.py "2023-04-19T12:00:00" "2023-04-19T12:03:00" \
    output-catalog-dir > catalog.csv
</pre>

*Note*: Replace **output-catalog-dir**  with a folder name that makes sense to you. 

## Waveform data

 It is possible to **download the waveforms** too. The script uses the previously downloaded catalog data (csv and QUAKEML files) and will download the waveforms for each event:

<pre>
 python fdsnws-download.py --waveforms catalog-dir catalog.csv
</pre>

Or, if you want to specify the lenght [sec] of the waveforms to download:

<pre>
 python fdsnws-download.py --waveforms catalog-dir catalog.csv 1.5
</pre>

*Note*: Replace **catalog-dir** and **catalog.csv** with the folder name and the csv file downloaded previously.

By default only the 

# Post-processing

Here we demonstrate how to use the csv catalog and the events stored in QUAKEML format. In this example each event is loaded and the corresponding arrivals are printed. The waveforms are loaded too and plotted.

```python
import pandas as pd
import obspy as ob
from pathlib import Path

csv_catalog = "catalog.csv"
xml_folder  = "mydirectory"  # folder where the event XML files were downloaded

### Load the csv catalog
cat = pd.read_csv(csv_catalog, dtype=str, na_filter=False)

#Loop through the catalog by event
for row in cat.itertuples():
  #
  # you can fetch any csv column with 'row.column'
  #
  print(f"Processing event {row.id}")

  #
  # We want to access all event information, so we load its XML file
  # with obspy
  #
  cat = ob.core.event.read_events( Path(xml_folder, f"ev{row.id}.xml"))
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
  trace = ob.read( Path(xml_folder, f"ev{row.id}.mseed"))
  trace.plot()

```

Similarly it is possible to load the inventory XML with obspy and use the [Inventory class](https://docs.obspy.org/packages/autogen/obspy.core.inventory.inventory.Inventory.html) to fetch the required information.

```python
import obspy as ob

inv = ob.core.inventory.inventory.read_inventory("inventory.xml")
print(inv.get_contents())
inv.get_coordinates("NET.STA.LOC.CHA", ob.UTCDateTime("2022-10-01"))
```







