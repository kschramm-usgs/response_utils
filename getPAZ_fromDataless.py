#!/bin/env python
''' method to get poles and zeros from dataless'''

# start with the method from the obspy documentation.
import obspy 
from obspy.io.xseed import Parser
from obspy import UTCDateTime

st = obspy.read('/msd/IU_ANMO/2017/123/00_BH1.512.seed')
parser = Parser('/dcc/metadata/dataless/DATALESS.IU.ANMO.seed')

time=UTCDateTime(2017,123)

paz = parser.getPAZ('IU.ANMO.00.BHZ',time)
print(paz)

