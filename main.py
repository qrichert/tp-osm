#!/usr/bin/env python
# coding: utf-8

# # TP Open Street Maps
# 
# Some useful OSMnx resources:
# 
# * [OSMnx GitHub repo](https://github.com/gboeing/osmnx)
# * [OSMnx doc](https://osmnx.readthedocs.io/)
# * [Camp sites](https://wiki.openstreetmap.org/wiki/Tag:tourism%3Dcamp_site)
# * [Tut](https://automating-gis-processes.github.io/2018/notebooks/L6/network-analysis.html#Saving-shortest-paths-to-disk)
# * [Other tut](https://automating-gis-processes.github.io/2018/notebooks/L6/retrieve_osm_data.html)
# 

# First, let's import the shit we need and set things up

# In[69]:


import osmnx as ox
import osmnxalt as oxalt
import matplotlib.pyplot as plt
import networkx as nx
from descartes import PolygonPatch
from shapely.geometry import Polygon, MultiPolygon
from matplotlib.collections import LineCollection
from pprint import pprint
import os
import math
import time

get_ipython().run_line_magic('matplotlib', 'inline')
ox.config(log_console=True, use_cache=True)

if not os.path.isdir('render'):
	os.mkdir('render')


# How far around the places do we look for camp sites?

# In[70]:


# Max dist in meters
# We don't want to drive more than 17 km for camping
# Also with more than 17 km it makes map generation too slow
K = 17 * 1000

print('Search for camp sites within a {} km bounding box'.format(K/1000))


# Here's a list of places

# In[71]:


placesToVisit = [
	'Brumath, France',
	# 'Krautwiller, France',
	'Strasbourg, France',
	# 'Weitbruch, France',
	# 'Gries, France',
	# 'Niederschaeffolsheim, France',
	# 'Haguenau, France',
	# 'Colmar, France',
	# 'Mulhouse, France',
	'Valras-Plage, France',
	'Paris, France',
	'New York, USA',
	# 'Los Angeles, USA',
]

print(placesToVisit)


# To use them more easily, in a standardized way, we convert them to a dictionary.
# 
# In the end, the dictionary will look like this:
# 
# ```python
# {
# 	'City 1, Country': {
# 		'coordinates': (48.7309406, 7.708107), # (lat, lon)
# 		'camp_sites': {
# 			'osmid': {
# 				'name': 'Camp site name',
# 				'phone': '03 88 00 00 00'
# 			},
# 			'osmid': {
# 				'name': 'Camp site name',
# 				'phone': '03 88 00 00 00'
# 			}
# 		},
# 		'selected_camp_site': 'Selected camp site ID'
# 	}
# }
# ```

# In[72]:


places = {}

for place in placesToVisit:
	# We use them as keys with an empty dictionary as value
	places[place] = {}

pprint(places)


# Extracting the geo coordinates from the places, because you can't give a distance
# with ox.pois_from_place(), only with ox.pois_from_point()

# In[73]:


for place in places.keys():
	# ox.geocode() returns a tuple of (lat, lon)
	places[place]['coordinates'] = ox.geocode(place)

pprint(places)


# The API sets its tags to 'nan' when the tag is there but not set.
# This is just a quick function to test whether a tag is NaN or not

# In[74]:


def isNaN(value):
	# math.isnan() fails on strings
	try:
		# math.isnan() works on numbers
		return math.isnan(value)
	except TypeError:
		# If it fails 'value' is not a number, so it can't be of NaN type
		return False


# A function to get the distance between two points using haversine.
# Returns the distance in meters
# Points are (lat, lon) tuples

# In[75]:


def distanceBetween(point1, point2):
	# Tuples to Lists
	point1 = list(point1)
	point2 = list(point2)

	return ox.great_circle_vec(lat1=point1[0],
	                           lng1=point1[1],
	                           lat2=point2[0],
	                           lng2=point2[1])


# Getting those camp sites

# In[76]:


for place, coordinates in places.items():
	print('Searching camp sites for {}...'.format(place))

	places[place]['camp_sites'] = {}

	nearbyCampSites = oxalt.pois_from_point(places[place]['coordinates'],
	                                       distance=K,
	                                       tags={'tourism':'camp_site'})

	# Converting GeoDataFrame to usable list of camp sites
	# See DataFrame.to_dict() (parent class of GeoDataFrame)
	nearbyCampSites = nearbyCampSites.to_dict(orient='records')

	# Here we go
	for camp in nearbyCampSites:
		# pprint(camp)

		# Now we can extract the informations we care about
		campInfos = {}

		# They always have coordinates
		# Almost every time there is a name
		# Sometimes there is a website and/or a phone number
		# Rarely there are start
		# Once there was camp_site
		# Never seen the others
		infosWeWant = [
			'name',  # string Name of this campsite
			'phone',  # string
			'website',  # string
			'addr',  # string
			'swimming_pool',  # (yes|no) Whether camp site has a swimming pool
			'bbq',  # (yes|no) Whether a bbq (grill) is available
			'camp_site',  # string (basic|standard|serviced|deluxe) A classification scheme.
			'stars',  # string (1-5) A classification scheme additional or alternatively to the camp_site-classification. The 5-star-system is widely used in France
			'booking',  # (yes|no|recommended|required) first come, first serve only?
			'power_supply'  # (yes|no) Whether the site has electricity (e.g. for caravans)
		]

		for info in infosWeWant:
			if info in camp:
				# I'd rather have None than those weird 'nan's. Nones are easier to test
				if not isNaN(camp[info]):
					campInfos[info] = camp[info]
				else:
					campInfos[info] = None

		# Coordinates tuple (lat, lon)
		campInfos['coordinates'] = (camp['geometry'].centroid.y,
		                            camp['geometry'].centroid.x)

		# We save the camp under its original OSM id
		places[place]['camp_sites'][camp['osmid']] = campInfos


# Now we have a list of places with coordinates and each place has a list
# of nearby camp sites, each having coordinates too.
# 
# What we want to do is select the best camp site for each place.
# 
# To do this, we start with the "extras". Because if certain camp sites
# have a swimming pool, we don't want to look at those which don't ;)
# 
# But if no camp site has a swimming pool, we'll look at something else.
# 
# If no camp site has a particular advantage, we go with the nearest one.

# In[77]:


# An ordered list of our selection criteria
selectionCriteria = [
	'swimming_pool',  # (yes|no) Whether camp site has a swimming pool
	'bbq',  # (yes|no) Whether a bbq (grill) is available
	'stars',  # string (1-5) A classification scheme additional or alternatively to the camp_site-classification. The 5-star-system is widely used in France
	'camp_site',  # string (basic|standard|serviced|deluxe) A classification scheme.
	'power_supply',  # (yes|no) Whether the site has electricity (e.g. for caravans)
	'addr',  # string Have they an address?
	'website',  # string Have they a website?
	'phone',  # string Have they a phone?
	'name'  # string Have they at least a name?
]

for place in places.values():
	print('---------------------------')
	# The camp that will be chosen for the current place (id)
	theBigWinner = None
	# If there are several winners for a criteria we must chose one of them
	itsATie = []
	usedCriteria = None

	for criteria in selectionCriteria:

		currentCriteriaWinners = []

		# Let's look at each camp
		for id, camp in place['camp_sites'].items():
			if criteria in camp:
				# Our criteria has been given, but it doesn't mean it is 'yes'
				if camp[criteria] is None or camp[criteria] == 'no':
					continue

				# We have a winner for the current criteria, so we add it to the winners list
				currentCriteriaWinners.append(id)

		# Now the list is filled with winners or not

		# If there is only one winner in the list, we have our camp, and we quit
		if len(currentCriteriaWinners) == 1:
			theBigWinner = currentCriteriaWinners[0]
			usedCriteria = criteria
			break

		# If there are several winners we need some more selection
		elif len(currentCriteriaWinners) > 1:
			# Stars go from (string) 1-5
			# Here all currentCriteriaWinners have stars
			# We select the one(s) which has/have the most
			if criteria == 'stars':
				bestRating = None
				campsWithMostStars = []  # IDs

				# Get the highest ranking
				for id in currentCriteriaWinners:
					stars = int(place['camp_sites'][id]['stars'])

					if bestRating is None:
						bestRating = stars
						continue

					if stars > bestRating:
						bestRating = stars

				# Get a list of camps with the highest ranking
				for id in currentCriteriaWinners:
					stars = int(place['camp_sites'][id]['stars'])

					if stars == bestRating:
						campsWithMostStars.append(id)

				# If there's only one with the highest ranking, we have our winner
				if len(campsWithMostStars) == 1:
					theBigWinner = campsWithMostStars[0]
					usedCriteria = criteria
					break
				# If there are multiple, we save them for next step
				else:
					itsATie = campsWithMostStars
					usedCriteria = criteria
					break

			# Camp rating are (string) basic|standard|serviced|deluxe
			# Here all currentCriteriaWinners have a rating
			# We select the one(s) which has/have the best
			elif criteria == 'camp_site':
				bestRating = None
				campsWithBestRating = []  # IDs

				# It's easier with ratings based on numbers, so we convert them
				campRatings = {
					'basic': 1,
					'standard': 2,
					'serviced': 3,
					'deluxe': 4
				}

				# Get the highest ranking
				for id in currentCriteriaWinners:
					rating = place['camp_sites'][id]['camp_site']  # Gives a string
					if rating in campRatings:
						rating = campRatings[rating] # Convert to value
					# If the rating isn't a valid rating it doesn't count
					# I saw 'reception' a few times, and it's not valid
					else:
						continue

					if bestRating is None:
						bestRating = rating
						continue

					if rating > bestRating:
						bestRating = rating

				# Get a list of camps with the highest rating
				for id in currentCriteriaWinners:
					rating = place['camp_sites'][id]['camp_site']
					if rating in campRatings:
						rating = campRatings[rating] # Convert to value
					else:
						continue

					if rating == bestRating:
						campsWithBestRating.append(id)

				# If there's only one with the highest ranking, we have our winner
				if len(campsWithBestRating) == 1:
					theBigWinner = campsWithBestRating[0]
					usedCriteria = criteria
					break
				# If there are multiple, we save them for next step
				else:
					itsATie = campsWithBestRating
					usedCriteria = criteria
					break

			# If there is still no clear winner, save them and quit for the next step
			else:
				itsATie = currentCriteriaWinners
				usedCriteria = criteria
				break

	if usedCriteria is not None:
		print('Criteria used: ' + usedCriteria)

	# If we found our winner we set it as selected camp site and go to the next place
	if theBigWinner is not None:
		place['selected_camp_site'] = theBigWinner
		continue

	# If not and none stand out, it's a tie
	if len(itsATie) == 0:
		itsATie = [i for i in place['camp_sites'].keys()]

	# All remaining camps will be decided by distance to the place visited

	# Will contain the id of the camp with the shortest distance
	shortestDistance = {'id': None, 'dist': None}

	for campID in itsATie:
		# Camp coordinates
		campCoord = place['camp_sites'][campID]['coordinates']
		# Place coordinates
		placeCoord = place['coordinates']

		dist = distanceBetween(campCoord, placeCoord)

		print(place['camp_sites'][campID]['name'] + ' - ' + str(dist) + 'm')

		# First iteration
		if shortestDistance['id'] is None:
			shortestDistance['id'] = campID
			shortestDistance['dist'] = dist
			continue

		# We found a shorter one
		if dist < shortestDistance['dist']:
			shortestDistance['id'] = campID
			shortestDistance['dist'] = dist

	# Now shortestDistance contains the one with the shortest distance to the place
	print('Shortest: ' + place['camp_sites'][shortestDistance['id']]['name'] + ' - ' + str(shortestDistance['dist']) + 'm')

	place['selected_camp_site'] = shortestDistance['id']


# Now we have a camp site for each place.

# In[78]:


for placeName, place in places.items():
	print('--- ' + placeName + ' ---')
	print('* ' + place['camp_sites'][place['selected_camp_site']]['name'])

	dist = distanceBetween(place['coordinates'],
	                       place['camp_sites'][place['selected_camp_site']]['coordinates'])
	dist = int(round(dist))

	print('* dist: ' + str(dist) + 'm')

	for key, val in place['camp_sites'][place['selected_camp_site']].items():
		if key == 'name':
			continue

		if val is None:
			val = '???'

		print('* ' + key + ': ' + str(val))


# Now we have the details we need for each place, we can generate some graphs
# 
# /!\ It takes about 20s per graph (16s for generating the graph, 4s for rendering) /!\

# In[79]:


for place in places.keys():
	print('---------------------------')
	print('Plotting the graph for {}...'.format(place))

	start = time.time()
	G = ox.graph_from_point(places[place]['coordinates'],
	                        distance=K,
	                        network_type='drive',
	                        name=place)
	end = time.time()

	print('Generation took: {}s'.format(end-start))

	start = time.time()
	fig, ax = ox.plot_graph(G,
	                        fig_height=15,
	                        axis_off=False,
	                        node_size=0,
	                        show=False,
	                        close=False)
	end = time.time()

	print('Plotting took: {}s'.format(end-start))

	coord = list(places[place]['coordinates'])

	# Now that we have the map, we can add other things on top

	# Place administrative boundaries
	gdf = ox.gdf_from_place(place)

	for geometry in gdf['geometry'].tolist():
		if isinstance(geometry, (Polygon, MultiPolygon)):
			if isinstance(geometry, Polygon):
				geometry = MultiPolygon([geometry])
			for polygon in geometry:
				patch = PolygonPatch(polygon,
				                     fc='#ffc0cb',
				                     ec='#ffc0cb',
				                     linewidth=3,
				                     alpha=0.3,
				                     zorder=-1)
				ax.add_patch(patch)

	# If no selected camp site, there are no camp sites at all
	if 'selected_camp_site' not in places[place]:
		print("It seems there isn't any camp site near " + place)
		# Saving the render, showing it & calling it a day
		plt.savefig('render/{}'.format(place))
		plt.show()
		continue

	# Nearby camp sites as red dots
	for camp in places[place]['camp_sites'].values():
		coord = list(camp['coordinates'])
		ax.scatter(coord[1], coord[0], c='red', s=100, zorder=5, alpha=0.6)

	# Selected camp site as green dot
	coord = places[place]['camp_sites'][places[place]['selected_camp_site']]['coordinates']
	selectedCampSiteCoord = coord # For later use, we save the tuple
	coord = list(coord)
	ax.scatter(coord[1], coord[0], c='green', s=300, zorder=10, alpha=1)

	# Adding the shortest route
	# From place
	source = ox.get_nearest_node(G, places[place]['coordinates'])
	# To camp site
	target = ox.get_nearest_node(G, selectedCampSiteCoord)

	route = nx.shortest_path(G, source=source, target=target)

	##########################################
	# Adapted from osmnx.plot_graph_routes() #
	##########################################

	# plot the route lines
	edge_nodes = list(zip(route[:-1], route[1:]))
	lines = []
	for u, v in edge_nodes:
		# if there are parallel edges, select the shortest in length
		data = min(G.get_edge_data(u, v).values(), key=lambda x: x['length'])

		# if it has a geometry attribute (ie, a list of line segments)
		if 'geometry' in data:
			# add them to the list of lines to plot
			xs, ys = data['geometry'].xy
			lines.append(list(zip(xs, ys)))
		else:
			# if it doesn't have a geometry attribute, the edge is a straight
			# line from node to node
			x1 = G.nodes[u]['x']
			y1 = G.nodes[u]['y']
			x2 = G.nodes[v]['x']
			y2 = G.nodes[v]['y']
			line = [(x1, y1), (x2, y2)]
			lines.append(line)

	# add the lines to the axis as a linecollection
	lc = LineCollection(lines, colors='blue', linewidths=4, alpha=0.5, zorder=3)
	ax.add_collection(lc)

	###################

	# Add route origin point (place)

	coord = list(list(places[place]['coordinates']))
	ax.scatter(coord[1], coord[0], c='blue', s=100, alpha=0.5, zorder=2)

	# Saving the render & showing it
	plt.savefig('render/{}'.format(place))
	plt.show()

