# Nuclease Assay % Cut and graphing traces

analyse-nuclease.R

Get percent cut from nuclease assay sequencer traces

The script will filter the peaks of a size you choose for the cut fragment and uncut fragment.
Then sum areas at each size for each trace.
Then calculate a percent cut.

Exported as csv to graph - I prefer graphing in GraphPad Prism - easier to manipulate and alter graph features.


traces-script.R

This script will align ladder to FAM trace, convert scan size to bp, and export sizing data for graphing.
This script exports more data points than you can get in Genemapper so your graphs are more precise.

