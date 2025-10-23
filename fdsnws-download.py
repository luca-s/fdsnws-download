#!/usr/bin/env python3

import sys
import re
from pathlib import Path
from urllib.parse import urlparse
import pandas as pd
import obspy as ob
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.core.event import Catalog
from obspy.core.stream import Stream
from obspy.clients.fdsn.header import FDSNTimeoutException
from obspy.clients.fdsn.header import FDSNRequestTooLargeException
from obspy.clients.fdsn.header import FDSNNoDataException

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
            valid_str = re.compile(r"[A-Z|a-z|0-9|\?|\*|\||\(|\)]*")
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
# Define a client that connects to a FDSN Web Service
# https://docs.obspy.org/packages/autogen/obspy.clients.fdsn.client.Client.html
#
def create_client(url, timeout):
  parsed = urlparse(url)
  base_url = parsed._replace(netloc=f"{parsed.hostname}:{parsed.port}").geturl()
  if parsed.username and parsed.password:
      print(f"Connecting to {base_url} with credentials", file=sys.stderr)
      return Client(base_url=base_url, user=parsed.username, password=parsed.password, timeout=timeout)
  else:
      print(f"Connecting to {base_url}", file=sys.stderr)
      return Client(base_url=base_url, timeout=timeout)

def usage():
    print("Usage:")
    print(f" {sys.argv[0]} fdsnws start-date end-date [--timeout sec]")
    print(f" {sys.argv[0]} fdsnws start-date end-date output-catalog-dir [--timeout sec]")
    print(f" {sys.argv[0]} fdsnws --waveforms catalog-dir catalog.csv [length-before-event:length-after-event] [station-filter] [--timeout sec]")
    print("\nNotes:")
    print(" Dates should be given in any format supported by obspy.UTCDateTime e.g. 2023-04-19T12:00:00 (UTC)")
    print(" --timeout is optional and specify the fdsnws connection timeout in seconds. Useful to deal with large requests that require a long time to be fulfilled")

def main():

  if len(sys.argv) < 4:
    usage()
    sys.exit(0)

  timeout = 300
  if sys.argv[-2] == "--timeout":
      timeout = float(sys.argv[-1])
      sys.argv = sys.argv[:-2]

  if sys.argv[2] == "--waveforms" and len(sys.argv) >= 5:
    client = create_client(sys.argv[1], timeout)
    catdir = sys.argv[3]
    catfile = sys.argv[4]
    length_before = None
    length_after = None
    if len(sys.argv) >= 6:
      tokens = sys.argv[5].split(":")
      if len(tokens) != 2:
        print(f"Wrong format for waveform length option: {sys.argv[5]}", file=sys.stderr)
        return
      length_before = float(tokens[0]) # [sec]
      length_after  = float(tokens[1]) # [sec]
    sta_filters = None
    if len(sys.argv) >= 7:
      sta_filters = StationNameFilter()
      if not sta_filters.set_rules(sys.argv[6]):
          return

    download_waveform(client, catdir, catfile, length_before, length_after, sta_filters)

  elif len(sys.argv) == 4 or len(sys.argv) == 5:
    client = create_client(sys.argv[1], timeout)
    starttime = UTCDateTime(sys.argv[2])
    endtime = UTCDateTime(sys.argv[3])
    catdir = None
    if len(sys.argv) == 5:
      catdir = sys.argv[4]

    download_catalog(client, catdir, starttime, endtime)

  else:
    print("Wrong syntax", file=sys.stderr)
    usage()


def download_catalog(client, catdir, starttime, endtime):

  if catdir:
    # create folder if it doesn't exist
    Path(catdir).mkdir(parents=True, exist_ok=True)

  #
  # Download the inventory
  # https://docs.obspy.org/packages/autogen/obspy.core.inventory.inventory.Inventory.html
  #
  if catdir:
    inv_file = Path(catdir, "inventory.xml")
    if not inv_file.is_file():
      print(f"Downloading inventory {inv_file}...", file=sys.stderr)
      try:
        inventory = client.get_stations(network="*", station="*", location="*", channel="*",
                                        starttime=starttime,  endtime=endtime,
                                        level="response")
        inventory.write(inv_file, format="STATIONXML")
      except Exception as e:
        print(f"Cannot download inventory: skip it ({e})", file=sys.stderr)
    else:
      print(f"Inventory file {inv_file} exists: do not download it again", file=sys.stderr)

  # csv file header
  print("id,event,origin,isotime,latitude,longitude,depth,event_type,mag_type,magnitude,method_id,evaluation_mode,author,rms,az_gap,num_phase")

  chunkstart = starttime
  chunkend   = endtime
  id = 0
  while chunkstart < endtime:

    if chunkend > endtime:
        chunkend = endtime

    print(f"Downloading {chunkstart} ~ {chunkend}...", file=sys.stderr)

    #
    # Load catalog:
    #  includeallorigins=False: we only want preferred origins, too heavy to load all orgins
    #  includearrivals=False: too slow to load, we can load the picks later in the loop
    #  includeallmagnitudes=False: optional, it depends on what magnitude we want to work with
    cat = None
    try:
        cat = client.get_events(starttime=chunkstart, endtime=chunkend, orderby="time-asc",
                                includeallorigins=False,
                                includearrivals=False,
                                includeallmagnitudes=False)
    except (FDSNRequestTooLargeException, FDSNTimeoutException) as e:
        if isinstance(e, FDSNTimeoutException):
            print("FDSNTimeoutException: splitting request in smaller ones (consider --timeout option)...", file=sys.stderr)
        else:
            print("FDSNRequestTooLargeException: splitting request in smaller ones...", file=sys.stderr)
        chunklen   = (chunkend - chunkstart) / 2.
        chunkend   = chunkstart + chunklen
        continue
    except FDSNNoDataException as e:
        print(f"FDSNNoDataException. No data between {chunkstart} ~ {chunkend}...", file=sys.stderr)

    #
    # Prepare next chunk to load
    #
    chunklen   = (chunkend - chunkstart)
    chunkstart = chunkend
    chunkend   += chunklen
    if (endtime - chunkend) < chunklen: # join leftover
        chunkend = endtime

    if cat is None:
        continue

    print(f"{len(cat.events)} events found", file=sys.stderr)

    #
    # Loop through the catalog and extract the information we need
    #
    for ev in cat.events:
      # ev: https://docs.obspy.org/packages/autogen/obspy.core.event.event.Event.html

      id += 1

      #
      # get preferred origin from event (an events might contain multiple origins: e.g. different locators or
      # multiple updates for the same origin, we only care about the preferred one)
      # https://docs.obspy.org/packages/autogen/obspy.core.event.origin.Origin.html
      #
      o = ev.preferred_origin()
      if o is None:
          print(f"No preferred origin for event {ev.resource_id}: skip it", file=sys.stderr)
          continue

      #
      # get preferred magnitude from event (multiple magnitude might be present, we only care about the preferred one)
      # https://docs.obspy.org/packages/autogen/obspy.core.event.magnitude.Magnitude.html
      #
      m = ev.preferred_magnitude()
      mag = None
      mag_type = "?"
      if m is not None:
        mag = m.mag
        mag_type = m.magnitude_type

      #
      # Write csv entry for this event
      #
      print(f"{id},{ev.resource_id},{o.resource_id},{o.time},{o.latitude},{o.longitude},{o.depth},"
            f"{ev.event_type if ev.event_type else ''},"
            f"{mag_type},{mag},"
            f"{o.method_id},{o.evaluation_mode},"
            f"{o.creation_info.author if o.creation_info else ''},"
            f"{o.quality.standard_error if o.quality else ''},"
            f"{o.quality.azimuthal_gap if o.quality else ''},"
            f"{o.quality.used_phase_count if o.quality else ''}")

      if catdir is None:
          continue

      ev_file = Path(catdir, f"ev{id}.xml")
      if ev_file.is_file():
        print(f"File {ev_file} exists: skip event", file=sys.stderr)
        continue

      #
      # Since `cat` doesn't contain arrivals and multiple origins/magnitudes (for performance reason),
      # we need to load them now
      #
      ev_id = str(ev.resource_id)
      try:
        client.get_events(eventid=ev_id,
                          includeallorigins=True,
                          includearrivals=True,
                          includeallmagnitudes=True,
                          filename=ev_file)
      except FDSNTimeoutException as e:
        print(f"FDSNTimeoutException while retrieving event {ev_id} (consider --timeout option)", file=sys.stderr)
      except Exception as e:
        print(f"While retrieving event {ev_id} an exception occurred: {e}", file=sys.stderr)


def download_waveform(client, catdir, catfile, length_before=None, length_after=None, sta_filters=None):

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
    evtime = UTCDateTime(row.isotime)

    #
    # you can fetch any csv column with 'row.column'
    #
    print(f"Processing event {ev_id} {evtime}", file=sys.stderr)

    wf_file = Path(catdir, f"ev{ev_id}.mseed")
    if wf_file.is_file():
      print(f"File {wf_file} exists: skip event", file=sys.stderr)
      continue

    #
    # loop trough origin arrivals to find the used stations and the
    # latest pick time
    #
    last_pick_time = None
    auto_sta_filters = StationNameFilter()

    if sta_filters is None or length_before is None or length_after is None:

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

      if o.time != evtime: # this should never happen
        print(f"Event {ev_id} prferred origin time is not the same as event time: skip event", file=sys.stderr)
        continue

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
            auto_sta_filters.add_rules(p.waveform_id.id)
            break

    #
    # Find the stations list that were active at this event time
    #
    stations = get_stations(
      inventory=inv, ref_time=evtime,
      sta_filters=(auto_sta_filters if sta_filters is None else sta_filters)
    )

    #
    # Fix the time window for all waveforms
    #
    if length_before is None:
        starttime = evtime
    else:
        starttime = evtime - length_before

    if length_after is None: 
      endtime = last_pick_time
      endtime += (endtime - starttime) * 0.3 # add extra 30% 
    else:
      endtime = evtime + length_after

    bulk = []
    for (network_code, station_code, location_code, channel_code) in stations:
      bulk.append( (network_code, station_code, location_code, channel_code, starttime, endtime) )

    #
    # Download the waveforms
    #
    waveforms = None
    if bulk:
      print(f"Downloading {len(bulk)} waveforms between {starttime} ~ {endtime} ", file=sys.stderr)
      try:
        # https://docs.obspy.org/packages/autogen/obspy.clients.fdsn.client.Client.get_waveforms_bulk.html
        # https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.html
        waveforms = client.get_waveforms_bulk(bulk)
      except FDSNTimeoutException as e:
        print(f"FDSNTimeoutException: Cannot fetch waveforms for event {id} (consider --timeout option)", file=sys.stderr)
      except Exception as e:
        print(f"Cannot fetch waveforms for event {id}: {e}", file=sys.stderr)

    if waveforms:
      waveforms.trim(starttime=starttime, endtime=endtime)
      print(f"Saving {wf_file}", file=sys.stderr)
      waveforms.write(wf_file, format="MSEED")

if __name__ == '__main__':
  main()
