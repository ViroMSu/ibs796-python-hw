#!/usr/bin/python

import requests
import json
import datetime
import sys

# By default, function searches for all MARTA routes.
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
	return(buses)

def getTrains(station=''):
	#URL for MARTA API for all trains
	base = 'http://developer.itsmarta.com/RealtimeTrain/RestServiceNextTrain/GetRealtimeArrivals?apikey=af04335a-0366-40c1-9e26-061af015c462'
	prediction = False
	# Formulate URL request and format response as json object
	response = requests.get(base, timeout=30)
	trains = response.json()
	output = ''
	
	# For each train in response, print a few pieces of data.
	for train in trains:
		if prediction:
			output = 'dir: ' + train['Direction'] + '  Waiting_Seconds:' + str(train['Waiting_Seconds']) + '  LOCATION:' + train['Location'] + '  TRACK:' + train['Track'] + '\n' + '  ID:' + train['Train_ID'] + '\n'
		else:
			output = train['Direction'] + '  LAT:' + train['Latitude'] + '  LON:' + train['Longitude'] + '  LOCATION:' + train['Location'] + '  TRACK:' + train['Track'] + '\n' + '  ID:' + train['Train_ID'] + '\n'
		print(output)

user_location = raw_input("Please enter location: ")
user_route = raw_input("Please enter bus route number: ")
user_direction = raw_input("N, S, E, W? ") 
user_destination = raw_input("Please enter destination: ")
other_route = raw_input('Is your destination a stop for another bus? If so, enter route number. Otherwise, press enter ')
if other_route != '':
	other_direc = raw_input('Direction of other route? N, S, E, W? ')

import googlemaps
#google api
gmaps = googlemaps.Client(key='AIzaSyD7me1LkBkQly7kd8hR1s4CDHbNwhr9HAE')

#Find long and lat of user location
from googlemaps import convert
location = gmaps.places(user_location, location = 'Atlanta')
#Check for more than one result
loc_results = location['results']
addresses = []
for result in loc_results:
	addresses.append(str(result['formatted_address']))

if len(addresses) > 1:
	for i in addresses:
		print i
	location_index = raw_input("Starting at an index of 0, which address did you mean? ")
	address = addresses[int(location_index)]
elif len(addresses) == 0:
	sys.exit("Initial location yielded no results. Please retry")
else:
	address = addresses[0]	

#Redo for destination
destination = gmaps.places(user_destination, location = 'Atlanta')
dest_results = destination['results']
addresses2 = []
for result in dest_results:
	addresses2.append(str(result['formatted_address']))
if len(addresses2) > 1:
	for i in addresses2:
		print i
	location_index2 = raw_input("Starting at an index of 0, which address did you mean? ")
	address2 = addresses2[int(location_index2)]
elif len(addresses2) == 0:
	sys.exit("Destination location yielded no results. Please retry")
else:
	address2 = addresses2[0]	


#Can just use the intersection
gmaps_distance = gmaps.distance_matrix(address, address2, mode='transit', transit_mode='bus', units='imperial', traffic_model='pessimistic')
time_est1 = str(gmaps_distance['rows'][0]['elements'][0]['duration']['text'])
time_est = int(time_est1.split()[0])

#find closest bus
bus_data = getBuses(user_route)
busdist = []
bustime = []
for bus in bus_data:
	if bus['DIRECTION'].startswith(user_direction.upper()):
		loc = bus['LATITUDE'] + ', ' + bus['LONGITUDE']
		bustime.append(int(bus['ADHERENCE']))
		gquery = gmaps.distance_matrix(address, loc, units='imperial')
		distance = str(gquery['rows'][0]['elements'][0]['distance']['text'])
		busdist.append(float(distance.split()[0]))
			
busindex = busdist.index(min(busdist))
timeliness = bustime[busindex]

if other_route != '':
	bus_data2 = getBuses(other_route)
	busdist2 = []
	bustime2 = []
	for bus in bus_data2:
		if bus['DIRECTION'].startswith(other_direc.upper()):
			loc2 = bus['LATITUDE'] + ', ' + bus['LONGITUDE']
			bustime2.append(int(bus['ADHERENCE']))
			gquery2 = gmaps.distance_matrix(address2, loc2, units='imperial')
			distance2 = str(gquery2['rows'][0]['elements'][0]['distance']['text'])
			busdist2.append(float(distance2.split()[0]))
	busindex2 = busdist2.index(min(busdist2))
	timeliness2 = bustime2[busindex2]

if timeliness == 0:
	print 'Bus', user_route, 'is on time, estimated arrival time to destination is', time_est, 'min'
elif timeliness < 0:
	print 'Bus', user_route, 'is %d minutes late, estimated arrival time to destination is' % (abs(timeliness)), time_est-timeliness, 'min'
elif timeliness > 0:
	print 'Bus', user_route, 'is %d minutes early, estimated arrival time to destination is' % (abs(timeliness)), time_est-timeliness, 'min'

if other_route != '':
	if timeliness2 == 0:
		print 'Bus', other_route, 'is on time'
	elif timeliness2 < 0:
		print 'Bus', other_route, 'is', timeliness2, 'min late'
	elif timeliness2 > 0:
		print 'Bus', other_route, 'is', timeliness2, 'min early'
