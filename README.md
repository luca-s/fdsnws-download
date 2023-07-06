# fdsnws-download

A customizable example script on how to use a FDSN Web Service programmatically with obspy to download an earthquake catalog.

The script can be run like this (download events between the passed dates and store it to catalog.csv):

<pre>
$ python fdsnws-download.py "2023-04-19T12:00:00" "2023-04-19T12:03:00" > catalog.csv
</pre>

While a csv file is easy to handle, some details cannot be stored in that format (e.g. the picks). For this reason it is possible to specify a destination folder where each event is stored in QUAKEML (that file can later be loaded with obspy to gain full event information):

<pre>
$ python fdsnws-download.py "2023-04-19T12:00:00" "2023-04-19T12:03:00" output-catalog-folder > catalog.csv
</pre>

Finally It is possible to download the waveforms too. That will be slow though, so fetching the catalog without waveforms first is recommended:

<pre>
$ python fdsnws-download.py "2023-04-19T12:00:00" "2023-04-19T12:03:00" output-catalog-folder --waveforms > catalog.csv
</pre>

