#!/usr/bin/env python

#import urllib2
import urllib.request
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
	response = urllib.request.urlopen(base + query, timeout=30)
	str_response = response.readall().decode('utf-8')
	buses = json.loads(str_response)
	
	
	
	# Prints entirety of json response
	#print(buses)
	
	output = ''
	
	# For each bus in response, print a few pieces of data.
	for bus in buses:
		output = bus['ROUTE'] + '  LAT:' + bus['LATITUDE'] + '  LON:' + bus['LONGITUDE'] + '  ADHER:' + bus['ADHERENCE'] + '  VEHICLE:' + bus['VEHICLE'] + '\n' + '  TRIP:' + bus['TRIPID']
		
		print(output)
		
	checkTripIds(buses)
	
def getTrains(station=''):
	
	#Base URL for MARTA API
	base = 'http://gispd/RTTService/RTTService.svc/'
	prediction = False
	
	# If user does not input a value for station name, use 'GetAllTrains' API call
	if station == '':
		query = 'GetAllTrains'
	
	# Else, use 'GetTrainPrediction' API call with user-defined station name
	else:
		query = 'GetTrainPrediction/' + station
		prediction = True
	
	# Formulate URL request and format response as json object
	response = urllib.request.urlopen(base + query, timeout=30)
	str_response = response.readall().decode('utf-8')
	trains = json.loads(str_response)
	
	# Prints entirety of json response
	#print(trains)
	
	output = ''
	
	# For each bus in response, print a few pieces of data.
	for train in trains:
		
		if prediction:
			output = 'dir: ' + train['Direction'] + '  Waiting_Seconds:' + str(train['Waiting_Seconds']) + '  LOCATION:' + train['Location'] + '  TRACK:' + train['Track'] + '\n' + '  ID:' + train['Train_ID'] + '\n'
		else:
			output = train['Direction'] + '  LAT:' + train['Latitude'] + '  LON:' + train['Longitude'] + '  LOCATION:' + train['Location'] + '  TRACK:' + train['Track'] + '\n' + '  ID:' + train['Train_ID'] + '\n'
		print(output)

# Determines how many of the buses have TripIDs
def checkTripIds(buses):
	count = 0
	total = len(buses)
	for bus in buses:
		tripId = bus['TRIPID']
		# print("Current Trip ID: " + tripId)
		
		if (tripId == ''):
			count += 1
	print(str(count) + " of " + str(total) + " do not have Trip IDs\n\n")
