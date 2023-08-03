#!/usr/bin/env python3

import sys
import re
import pandas as pd
import obspy as ob
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.core.event import Catalog
from obspy.core.stream import Stream
from pathlib import Path


#
# Class to parse and verify a station filter. A station filter is
# something in the form:
#   net, net.sta, net.sta.loc, net.sta.loc.cha
# Each of net,sta,loc,cha supports wildcards such as *,?,(,),|
#
class StationNameFilter:

    def __init__(self):
        self.rules = {}

    def set_rules(self, rules):
        self.rules = {}
        return self.add_rules(rules)

    def add_rules(self, rules):
        #
        # split multiple comma separated rules
        #
        for singleRule in rules.split(","):
            # split the rule in NET STA LOC CHA
            tokens = singleRule.split(".")
            if len(tokens) == 1: # only network
                tokens.append("*")
                tokens.append("*")
                tokens.append("*")
            elif len(tokens) == 2: # only network.station
                tokens.append("*")
                tokens.append("*")
            elif len(tokens) == 3: # only network.station.location
                tokens.append("*")
            elif len(tokens) == 4: # all network.station.location.channel
                pass
            else:
                print(f"Error: check station filter syntax ({rules})",
                        file=sys.stderr)
                return False
            #
            # Check for valid characters
            #
            valid_str = re.compile("[A-Z|a-z|0-9|\?|\*|\||\(|\)]*")
            for tok in tokens:
                if valid_str.fullmatch(tok) is None:
                  print(f"Error: check station filter syntax ({rules})",
                          file=sys.stderr)
                  return False

            fullRule = f"{tokens[0]}.{tokens[1]}.{tokens[2]}.{tokens[3]}"
            #
            # convert user special characters (* ? .) to regex equivalent
            #
            reRule = re.sub(r"\.", r"\.", fullRule)  # . becomes \.
            reRule = re.sub(r"\?", ".",   reRule)    # ? becomes .
            reRule = re.sub(r"\*", ".*",  reRule)   # * becomes .*
            reRule = re.compile(reRule)

            self.rules[fullRule] = reRule
        return True

    def match(self, waveform_id):
        for (fullRule, reRule) in self.rules.items():
            if reRule.fullmatch(waveform_id) is not None:
                return True
        return False

    def print(self):
        for (fullRule, reRule) in self.rules.items():
            print(f"{fullRule} (regular expression: {reRule.pattern})")


#
# Utility function that returns the stations active at a specific point in time
#
def get_stations(inventory, ref_time, sta_filters):
    stations = []
    for net in inventory.networks:
        if ref_time < net.start_date or ref_time > net.end_date:
            continue
        for sta in net.stations:
            if ref_time < sta.start_date or ref_time > sta.end_date:
                continue
            for cha in sta.channels:
                if ref_time < cha.start_date or ref_time > cha.end_date:
                    continue
                wfid = f"{net.code}.{sta.code}.{cha.location_code}.{cha.code}"
                if sta_filters.match(wfid):
                    stations.append( (net.code, sta.code, cha.location_code, cha.code) )
    return stations


#
# Utility function that transforms magnitude to a phisical size, so that it can
# be used as the size of an event when plotting it
#
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
    print(f" {sys.argv[0]} --waveforms catalog-dir catalog.csv [wf-length] [wf-filter]")
    sys.exit(0)

  #
  # Define a client that connects to a FDSN Web Service
  # https://docs.obspy.org/packages/autogen/obspy.clients.fdsn.client.Client.html
  #
  client = Client('ClientName') # When using a client  without credentials 
  #client = Client(base_url="http://myfdsnws.somewhere:8080", user="user", password="pass") # When using password protected clients
  #client = Client(base_url="http://myfdsnws.somewhere:8080" ) # When using custom client without credentials


  if sys.argv[1] == "--waveforms" and len(sys.argv) >= 4:
    catdir = sys.argv[2]
    catfile = sys.argv[3]
    wflength = None
    if len(sys.argv) >= 5:
      wflength = float(sys.argv[4])  # [sec] how much waveform to download
    sta_filters = None
    if len(sys.argv) >= 6:
      sta_filters = StationNameFilter()
      if not sta_filters.set_rules(sys.argv[5]):
          return
    download_waveform(client, catdir, catfile, wflength, sta_filters)
  elif len(sys.argv) == 3 or len(sys.argv) == 4:
    starttime = UTCDateTime(sys.argv[1])
    endtime = UTCDateTime(sys.argv[2])
    catdir = None
    if len(sys.argv) == 4:
      catdir = sys.argv[3]
    download_catalog(client, catdir, starttime, endtime)
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
                                    level="response")
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


def download_waveform(client, catdir, catfile, wflength=None, sta_filters=None):
  #
  # Load the csv catalog
  #
  print(f"Loading catalog {catfile}...", file=sys.stderr)
  cat = pd.read_csv(catfile, dtype=str, na_filter=False)
  print(f"Found {len(cat)} events", file=sys.stderr)

  #
  # Load the inventory
  #
  print(f"Loading inventory {Path(catdir, 'inventory.xml')}...", file=sys.stderr)
  inv = ob.core.inventory.inventory.read_inventory( Path(catdir, "inventory.xml") )

  if sta_filters:
    print(f"Using the following station filters...", file=sys.stderr)
    sta_filters.print()

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
    # loop trough origin arrivals to find the used stations and the
    # latest pick time
    #
    last_pick_time = None
    auto_sta_filters = StationNameFilter()
    for a in o.arrivals:
      #
      # find the pick associated with the current arrival
      #
      for p in ev.picks:
        if p.resource_id == a.pick_id:
          #
          # Keep track of the last pick time
          #
          if last_pick_time is None or p.time > last_pick_time:
              last_pick_time = p.time
          #
          # Keep track of the stations used by the picks
          #
          wfid = p.waveform_id
          auto_sta_filters.add_rules(
            f"{wfid.network_code}.{wfid.station_code}.{wfid.location_code if wfid.location_code else ''}.{wfid.channel_code}"
          )
          break
    #
    # Find the stations list that were active at this event time
    #
    stations = get_stations(
      inventory=inv, ref_time=o.time,
      sta_filters=(auto_sta_filters if sta_filters is None else sta_filters)
    )

    #
    # Now that we know the last pick time we can fix the time window for all waveforms
    #
    starttime = o.time
    if wflength is None:
      endtime = last_pick_time
      endtime += (endtime - starttime) * 0.3 # add extra 30% 
    else:
      endtime = starttime + wflength # set user defined length

    bulk = []
    for (network_code, station_code, location_code, channel_code) in stations:
      bulk.append( (network_code, station_code, location_code, channel_code, starttime, endtime) )

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
      waveforms.trim(starttime=starttime, endtime=endtime)
      print(f"Saving {dest}", file=sys.stderr)
      waveforms.write(dest, format="MSEED")

if __name__ == '__main__':
  main()
