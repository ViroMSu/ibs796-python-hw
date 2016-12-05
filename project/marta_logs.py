#!/usr/bin/python

import requests
import json
import time
import datetime
import argparse
from collections import Counter
from pandas import Series
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser("commandline arguments")
parser.add_argument("--getlogs", help="option to begin pulling api data, state interval in seconds", type=int)
parser.add_argument("--stats", help="given api data logfile name, get stats on adherence")
parser.add_argument("--bus", help="bus to get stats for, default is bus 16")
args = parser.parse_args()

if args.getlogs:
	logs = True
else:
	logs = False

def getBuses(route=''):
	#Base URL for MARTA API
	base = 'http://developer.itsmarta.com/BRDRestService/BRDRestService.svc/'
	# If user does not input a value for route number, use 'GetAllBus' API call
	if route == '':
		query = 'GetAllBus'
	# Else, use 'GetBusByRoute' API call with user-defined route number
	else:
		query = 'GetBusByRoute/' + route
	# Formulate URL request and format response as json object
	response = requests.get(base + query, timeout=30)
	#automatically decoded
	buses = response.json()
	output = ''
	
	# For each bus in response, print a few pieces of data.
	for bus in buses:
		output = bus['ROUTE'] + '  LAT:' + bus['LATITUDE'] + '  LON:' + bus['LONGITUDE'] + '  ADHER:' + bus['ADHERENCE'] + '  DIRECTION:' + bus['DIRECTION']
		print(output)

while logs:
	now = datetime.datetime.now()
	hour = now.hour
	if hour >= 5 & hour <= 23:
		print 'Time:', str(now.time())
		getBuses()
		time.sleep(args.getlogs)
	else:
		time.sleep(180000)

if args.stats:
	lfile = open(args.stats, 'r')

if args.bus:
	bus = args.bus + ' '
else:
	bus = '16 '

time = []
adher = []
direc = []
not_eof = True
while not_eof:
	line = lfile.readline().strip()
	if line.startswith('Time'):
		timept = line.split()[1]
	if line.startswith(bus):
		adherpt = line.split()[3][6:]
		direcpt = line.split()[4][10:]
		time.append(timept)
		adher.append(int(adherpt))
		direc.append(direcpt)
	if line == "":
		not_eof = False

direc2 = Series(direc)
uniq_direc = direc2.unique()
direc1 = uniq_direc[0]
direc2 = uniq_direc[1]

time_freq =  Counter(time)
time2 = Series(time)
uniq_time = time2.unique()

avg_adher1 = []
avg_adher2 = []

time_index = 0
for x in uniq_time:
	num_buses = time_freq[x]
	bus_adher = adher[time_index:time_index+num_buses]
	bus_direc = direc[time_index:time_index+num_buses]
	direc1_adher = []
	direc2_adher = [] 
	for i,x2 in enumerate(bus_direc):
		if x2 == direc1:
			direc1_adher.append(bus_adher[i])
		else:
			direc2_adher.append(bus_adher[i])
	adher1 = np.mean(direc1_adher)
	avg_adher1.append(adher1)
	adher2 = np.mean(direc2_adher)
	avg_adher2.append(adher2)
	time_index += num_buses

tick_freq = range(0, len(avg_adher1), 100)
all_labels = list(uniq_time)
labels = []
for x in all_labels:
	labels.append(x[0:5])

f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
ax1.plot(range(len(avg_adher1)), avg_adher1, 'b-')
ax1.set_title('Average adherence of ' + direc1 + ' bus ' + bus)
ax2.plot(range(len(avg_adher2)), avg_adher2, 'r-')
ax2.set_title('Average adherence of ' + direc2 + ' bus ' + bus)
plt.ylabel('Adherence')
plt.xlabel('Time')
plt.xticks(tick_freq, labels[0::100], rotation='vertical')
plt.show()

lfile.close()
