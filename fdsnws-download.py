#!/usr/bin/env python3

import sys
from pathlib import Path
import obspy as ob
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.core.event import Catalog
from obspy.core.stream import Stream


def magToSize(scaler, mag):
    # [Pa] stress drop assumption, could also be higher than 1MPa, but go for that for nowâ€¦
    sigma = 1e6
    # [Nm] Scalar seismic moment from magnitude (Kanamori)
    M0 = 10**((mag + 6.03)*3/2)
    # [m]  Source radii from seismic moment and stress drop (Keilis-Borok, 1959)
    SR = (M0 * 7/16 * 1/sigma)**(1/3)
    return SR * scaler  # scale the physical unit by something sensible for the plot


def main():

    if len(sys.argv) < 3:
        print("Usage:")
        print(f" {sys.argv[0]} start-date end-date [output-catalog-folder] [--waveforms]")
        sys.exit(0)

    starttime = UTCDateTime(sys.argv[1])
    endtime = UTCDateTime(sys.argv[2])

    catdir = None
    if len(sys.argv) >= 4:
        catdir = sys.argv[3]
        # make out folder if it doesn't exitst
        Path(catdir).mkdir(parents=True, exist_ok=True)

    include_wf = False
    if len(sys.argv) >= 5 and sys.argv[4] == "--waveforms":
        include_wf = True

    #
    # Create a client that connects to a FDSN Web Service
    # https://docs.obspy.org/packages/autogen/obspy.clients.fdsn.client.Client.html
    #
    client = Client(base_url="http://myfdsnws.somewhere:8080") #, user="user", password="pass") optional

    #
    # Download the inventory
    # https://docs.obspy.org/packages/autogen/obspy.core.inventory.inventory.Inventory.html
    #
    if catdir:
        inventory = client.get_stations(network="*", station="*",  starttime=starttime,  endtime=endtime)
        inventory.write(Path(catdir, "inventory.xml"), format="STATIONXML")

    #
    # Load catalog:
    #  includeallorigins=False: we only want preferred origins, too heavy to load all orgins
    #  includearrivals=False: too slow to load, we can load the picks later in the loop
    #  includeallmagnitudes=False: optional, it depends on what magnitude we want to work with
    cat = client.get_events(starttime=starttime, endtime=endtime,
                            includeallorigins=False, includearrivals=False, includeallmagnitudes=False)

    # csv file header
    print("id,time,latitude,longitude,depth,mag,mag_plot_size,rms,az_gap,num_phase,num_station")

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
        # https://docs.obspy.org/packages/autogen/obspy.core.event.origin.Origin.html#obspy.core.event.origin.Origin
        #
        o = ev_with_picks.preferred_origin()

        #
        # get preferred magnitude from event (multiple magnitude might be present, we only care about the preferred one)
        # https://docs.obspy.org/packages/autogen/obspy.core.event.magnitude.Magnitude.html#obspy.core.event.magnitude.Magnitude
        #
        m = ev_with_picks.preferred_magnitude()
        mag = -99  # default value in case magnitude is not computed for this event
        mag_size = 0  # convert magnitude to a size that can be used for plotting
        if m is not None:
            mag = m.mag
            mag_size = magToSize(15, mag)

        #
        # loop trough origin arrivals
        #
        used_stations = set()
        waveforms = Stream()
        for a in o.arrivals:
            #
            # find the pick associated with the current arrival
            # https://docs.obspy.org/packages/autogen/obspy.core.event.origin.Arrival.html#obspy.core.event.origin.Arrival
            # https://docs.obspy.org/packages/autogen/obspy.core.event.origin.Pick.html#obspy.core.event.origin.Pick
            #
            for p in ev_with_picks.picks:
                if p.resource_id == a.pick_id:
                    wfid = p.waveform_id
                    used_stations.add(wfid.network_code + "." + wfid.station_code +
                                      "." + (wfid.location_code if wfid.location_code else ""))
                    #
                    # Download the waveform relative to the pick time
                    # https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.html
                    #
                    if include_wf:
                        try:
                            timewin = 10.000  # [sec] how much waveform to download
                            starttime = p.time-timewin/2
                            endtime = p.time+timewin/2
                            st = client.get_waveforms(
                                wfid.network_code, wfid.station_code, wfid.location_code, wfid.channel_code, starttime=starttime, endtime=endtime)
                            waveforms.extend(st)
                        except Exception as e:
                            print(
                                f"Cannot load trace {wfid.network_code}.{wfid.station_code}.{wfid.location_code}.{wfid.channel_code} {starttime}-{endtime}: {e}", file=sys.stderr)
                    break

        #
        # Here it is possible to filter the catalog by certain criteria:
        #
        # example of filtering: len(o.arrivals) > 8 and len(used_stations) > 4 and o.quality.azimuthal_gap < 180 and mag > 2.0
        if True:

            #
            # Write csv entry for this event
            #
            print(f"{id},{o.time},{o.latitude},{o.longitude},{o.depth},{mag},{mag_size},{o.quality.standard_error},{o.quality.azimuthal_gap},{len(o.arrivals)},{len(used_stations)}")

            #
            # Write extended event information as xml and waveform data
            #
            if catdir:
                cat_out = Catalog()
                cat_out.append(ev)
                cat_out.write(Path(catdir, f"ev{id}.xml"), format="QUAKEML")
                if waveforms:
                    waveforms.write(
                        Path(catdir, f"ev{id}.mseed"), format="MSEED")


if __name__ == '__main__':
    main()
