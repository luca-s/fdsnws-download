#!/usr/bin/env python3

import sys
from pathlib import Path
import pandas as pd
import obspy as ob
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.core.event import Catalog
from obspy.core.stream import Stream


def magToSize(scaler, mag):
  sigma = 1e6 # [Pa] stress drop assumption, could also be higher than 1MPa
  M0 = 10**((mag + 6.03)*3/2) # [Nm] Scalar seismic moment from magnitude (Kanamori)
  SR = (M0 * 7/16 * 1/sigma)**(1/3) # [m]  Source radii from seismic moment and stress drop (Keilis-Borok, 1959)
  return SR * scaler  # scale the physical unit by something sensible for the plot


def main():

  if len(sys.argv) < 3:
    print("Usage:")
    print(f" {sys.argv[0]} start-date end-date")
    print(f" {sys.argv[0]} start-date end-date output-catalog-dir")
    print(f" {sys.argv[0]} start-date end-date --waveforms catalog-dir catalog.csv")
    sys.exit(0)

  starttime = UTCDateTime(sys.argv[1])
  endtime = UTCDateTime(sys.argv[2])

  #
  # Define a client that connects to a FDSN Web Service
  # https://docs.obspy.org/packages/autogen/obspy.clients.fdsn.client.Client.html
  #
  client = Client('ClientName') # When using a client  without credentials 
  #client = Client(base_url="http://myfdsnws.somewhere:8080", user="user", password="pass") # When using password protected clients
  #client = Client(base_url="http://myfdsnws.somewhere:8080" ) # When using custom client without credentials


  if len(sys.argv) == 3:
    download_catalog(client, None, starttime, endtime)
  elif len(sys.argv) == 4 and sys.argv[3] != "--waveforms":
    catdir = sys.argv[3]
    download_catalog(client, catdir, starttime, endtime)
  elif len(sys.argv) == 6 and sys.argv[3] == "--waveforms":
    catdir = sys.argv[4]
    catfile = sys.argv[5]
    download_waveform(client, catdir, catfile)
  else:
    print("Wrong syntax", file=sys.stderr)


def download_catalog(client, catdir, starttime, endtime):

  if catdir:
    # create folder if it doesn't exist
    Path(catdir).mkdir(parents=True, exist_ok=True)

  #
  # Download the inventory
  # https://docs.obspy.org/packages/autogen/obspy.core.inventory.inventory.Inventory.html
  #
  if catdir:
    inventory = client.get_stations(network="*", station="*", location="*", channel="*",
                                    starttime=starttime,  endtime=endtime,
                                    level="channel")
    inventory.write(Path(catdir, "inventory.xml"), format="STATIONXML")

  #
  # Load catalog:
  #  includeallorigins=False: we only want preferred origins, too heavy to load all orgins
  #  includearrivals=False: too slow to load, we can load the picks later in the loop
  #  includeallmagnitudes=False: optional, it depends on what magnitude we want to work with
  cat = client.get_events(starttime=starttime, endtime=endtime, 
                          includeallorigins=False,
                          includearrivals=False,
                          includeallmagnitudes=False)

  # csv file header
  print("id,time,latitude,longitude,depth,mag_type,mag,mag_plot_size,method_id,evaluation_mode,author,rms,az_gap,num_phase,num_station")

  #
  # Loop through the catalog and extract the information we need
  #
  id = 0
  for ev in cat.events:
    id += 1

    # ev = https://docs.obspy.org/packages/autogen/obspy.core.event.event.Event.html

    #
    # Since `cat` doesn't contain arrivals (for performance reason), we need to load them now
    #
    ev_id = ev.resource_id.id.removeprefix(
        'smi:org.gfz-potsdam.de/geofon/')  # this prefix is added by obspy
    origin_cat = client.get_events(
        eventid=ev_id, includeallorigins=False, includearrivals=True)
    if len(origin_cat.events) != 1:
      raise Exception("Something went wrong")
    ev_with_picks = origin_cat.events[0]
    if ev.resource_id != ev_with_picks.resource_id:
      raise Exception("Something went wrong")

    #
    # get preferred origin from event (an events might contain multiple origins: e.g. different locators or
    # multiple updates for the same origin, we only care about the preferred one)
    # https://docs.obspy.org/packages/autogen/obspy.core.event.origin.Origin.html
    #
    o = ev_with_picks.preferred_origin()
    if o is None:
        print(f"No preferred origin for event {ev_id}: skip", file=sys.stderr)
        continue

    #
    # filter out non manual events ?
    #
    #if o.evaluation_mode != "manual":
    #    continue

    #
    # get preferred magnitude from event (multiple magnitude might be present, we only care about the preferred one)
    # https://docs.obspy.org/packages/autogen/obspy.core.event.magnitude.Magnitude.html
    #
    m = ev_with_picks.preferred_magnitude()
    mag = -99  # default value in case magnitude is not computed for this event
    mag_type = "?"
    mag_size = 0  # convert magnitude to a size that can be used for plotting
    if m is not None:
      mag = m.mag
      mag_type = m.magnitude_type
      mag_size = magToSize(15, mag)

    #
    # loop trough origin arrivals
    #
    used_stations = set()
    for a in o.arrivals:
      #
      # find the pick associated with the current arrival
      # https://docs.obspy.org/packages/autogen/obspy.core.event.origin.Arrival.html
      # https://docs.obspy.org/packages/autogen/obspy.core.event.origin.Pick.html
      #
      for p in ev_with_picks.picks:
        if p.resource_id == a.pick_id:
          wfid = p.waveform_id
          used_stations.add(wfid.network_code + "." + wfid.station_code +
              "." + (wfid.location_code if wfid.location_code else ""))
          break

    #
    # Here it is possible to filter the catalog by certain criteria:
    # See also:
    #    https://docs.obspy.org/packages/autogen/obspy.core.event.origin.OriginQuality.html
    #
    if True:  # example of filtering: len(o.arrivals) > 8 and len(used_stations) > 4 and o.quality.azimuthal_gap < 180 and mag > 2.0 

      # Write csv entry for this event
      print(f"{id},{o.time},{o.latitude},{o.longitude},{o.depth},{mag_type},{mag},{mag_size},{o.method_id},{o.evaluation_mode},{o.creation_info.author},{o.quality.standard_error},{o.quality.azimuthal_gap},{len(o.arrivals)},{len(used_stations)}")

      #
      # Write extended event information as xml and waveform data
      #
      if catdir:
        cat_out = Catalog()
        cat_out.append(ev_with_picks)
        cat_out.write(Path(catdir, f"ev{id}.xml"), format="QUAKEML")


def download_waveform(client, catdir, catfile):
  #
  # Load the csv catalog
  #
  cat = pd.read_csv(catfile, dtype=str, na_filter=False)
  print(f"Loaded catalog {catfile}: found {len(cat)} events", file=sys.stderr)

  #
  # Loop through the catalog by event
  #
  for row in cat.itertuples():

    ev_id = row.id

    #
    # you can fetch any csv column with 'row.column'
    #
    print(f"Processing event {ev_id}", file=sys.stderr)

    dest = Path(catdir, f"ev{ev_id}.mseed")
    if dest.is_file():
      print(f"File {dest} exists: skip event", file=sys.stderr)
      continue

    #
    # We want to access all event information, so we load its XML file
    # with obspy
    #
    evfile = Path(catdir, f"ev{ev_id}.xml")
    print(f"Reading {evfile}", file=sys.stderr)
    tmpcat = ob.core.event.read_events(evfile)
    ev = tmpcat[0] # we know we stored only one event in this catalog

    #
    # get preferred origin from event
    #
    o = ev.preferred_origin()
    if o is None:
        print(f"No preferred origin for event {ev_id}: skip", file=sys.stderr)
        continue

    #
    # loop trough origin arrivals
    #
    bulk = []
    for a in o.arrivals:
      #
      # find the pick associated with the current arrival
      #
      for p in ev.picks:
        if p.resource_id == a.pick_id:
          wfid = p.waveform_id
          #
          # Keep track of the pick waveform to download
          #
          extratime = 1.0  # [sec] how much waveform to download after the pick
          starttime = o.time
          endtime = p.time+extratime
          bulk.append( (wfid.network_code, wfid.station_code, wfid.location_code, wfid.channel_code, starttime, endtime) )
          break

    #
    # Download the waveforms
    #
    waveforms = None
    if bulk:
      print(f"Downloading {len(bulk)} waveforms", file=sys.stderr)
      try:
        # https://docs.obspy.org/packages/autogen/obspy.clients.fdsn.client.Client.get_waveforms_bulk.html
        # https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.html
        waveforms = client.get_waveforms_bulk(bulk)
      except Exception as e:
        print(f"Cannot fetch waveforms for event {id}: {e}", file=sys.stderr)

    if waveforms:
      waveforms.trim(starttime=o.time)
      print(f"Saving {dest}", file=sys.stderr)
      waveforms.write(dest, format="MSEED")

if __name__ == '__main__':
  main()
