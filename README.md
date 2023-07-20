# fdsnws-download

A customizable example script on how to use a FDSN Web Service programmatically with obspy to download an earthquake catalog.

## Running the script

**First of all edit the script and set the FDSNWS URL you want**.

The script can be run like this (download events between the passed dates and store it to catalog.csv):

<pre>
$ python fdsnws-download.py "2023-04-19T12:00:00" "2023-04-19T12:03:00" \
  > catalog.csv
</pre>

While a csv file is easy to handle, some details cannot be stored in that format (e.g. the picks). For this reason it is possible to specify a **destination folder** where each **event** is stored **in QUAKEML** (that file can later be loaded with obspy to gain full event information):

<pre>
$ python fdsnws-download.py "2023-04-19T12:00:00" "2023-04-19T12:03:00" \
  output-catalog-folder > catalog.csv
</pre>

Finally It is possible to **download the waveforms** too. That will be slow though, so fetching the catalog without waveforms first is recommended:

<pre>
$ python fdsnws-download.py "2023-04-19T12:00:00" "2023-04-19T12:03:00" \
  output-catalog-folder --waveforms > catalog.csv
</pre>

## Using the downloaded data

This is an example on how to use the csv catalog and the xml events. In this example each event is loaded and its arrivals are printed.

<pre>
import pandas as pd
import obspy as ob
from pathlib import Path

csv_catalog = "catalog.csv"
xml_folder  = "somewhere"  # folder where the event XML files were downloaded

#
# Load the csv catalog
#
cat = pd.read_csv(csv_catalog, dtype=str, na_filter=False)

#
# Loop through the catalog by event
#
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
    # We can now have access to all events information e.g.
    #

    # get preferred origin
    o = ev.preferred_origin()
    print(f"Origin {o.longitude} {o.latitude} {o.depth} {o.time}")

    # loop trough origin arrivals
    for a in o.arrivals:
      # find the pick associated with the current arrival
      for p in ev:
        if p.resource_id == a.pick_id:
          wfid = p.waveform_id
          print(f"Found pick {p.time} @ {wfid.network_code}.{wfid.station_code}.{wfid.location_code}.{wfid.channel_code}")
          break
</pre>

Similarly it is possible to load the inventory XML with obspy and use the [Inventory class](https://docs.obspy.org/packages/autogen/obspy.core.inventory.inventory.Inventory.html) to fetch the required information.

<pre>
import obspy as ob

inv = ob.core.inventory.inventory.read_inventory("inventory.xml")

print(inv.get_contents())

inv.get_coordinates("8R.V0102..JDD", ob.UTCDateTime("2022-10-01"))

</pre>

