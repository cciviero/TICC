import math as m

latlonr = raw_input('Enter lat, lon, r?')

lat = float(latlonr.split()[0])
lon = float(latlonr.split()[1])
r = float(latlonr.split()[2])

lat_r = lat*m.pi/180.
lon_r = lon*m.pi/180.

print 'x, y, z'
print r*m.cos(lat_r)*m.cos(lon_r), r*m.cos(lat_r)*m.sin(lon_r), r*m.sin(lat_r)
