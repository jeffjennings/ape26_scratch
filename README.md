Step 1:
```
# Unchanging API
# Normal type hierarchy
ra, dec = ..., ...
rep = SphericalRepresentation(ra, dec)
c = ICRS(rep)  # the same. not yet deprecated.
sc = SkyCoord(c)
# Flexible inputs
c = ICRS(ra_arr, dec_arr)  # the same. not yet deprecated.
c = ICRS(rep)  # the same. not yet deprecated.
sc = SkyCoord(ra_arr, dec_arr, frame=ICRS())  # the same. not yet deprecated.
sc = SkyCoord(rep, frame=ICRS())  # the same. not yet deprecated.

# New API
frame = ICRSFrame()
sc = SkyCoord(rep, frame=frame)
sc = SkyCoord(ra_arr, dec_arr, frame=frame)
```

Step 2:
```
# Unchanging API
# Normal type hierarchy
ra, dec = ..., ...
rep = SphericalRepresentation(ra, dec)
c = ICRS(rep)  # the same. not yet deprecated.
sc = SkyCoord(c)
# Flexible inputs
c = ICRS(ra_arr, dec_arr)  # the same. not yet deprecated.
c = ICRS(rep)  # the same. not yet deprecated.
sc = SkyCoord(ra_arr, dec_arr, frame=ICRS())  # the same. not yet deprecated.
sc = SkyCoord(rep, frame=ICRS())  # the same. not yet deprecated.

# New API
frame = ICRSFrame()
c = Coordinate(rep, frame)
sc = SkyCoord(rep, frame=frame)
sc = SkyCoord(c.data, c.frame)  # (just showing the strict constructor)
sc = SkyCoord(c)  # flexible
sc = SkyCoord(ra_arr, dec_arr, frame=frame)  # flexible
```

Step 3:
```
# Unchanging API
# Normal type hierarchy
ra, dec = ..., ...
rep = SphericalRepresentation(ra, dec)
c = ICRS(rep)  # the same. not yet deprecated.
sc = SkyCoord(c)
# Flexible inputs
c = ICRS(ra_arr, dec_arr)  # the same. not yet deprecated.
c = ICRS(rep)  # the same. not yet deprecated.
sc = SkyCoord(ra_arr, dec_arr, frame=ICRS())  # redirects to SkyCoord(ra_arr, dec_arr, frame=ICRSFrame()) 
sc = SkyCoord(rep, frame=ICRS())  # redirects to SkyCoord(ra_arr, dec_arr, frame=ICRSFrame()) 

# New API
frame = ICRSFrame()
c = Coordinate(rep, frame)
sc = SkyCoord(rep, frame=frame)
sc = SkyCoord(c)  # flexible
sc = SkyCoord(ra_arr, dec_arr, frame=frame)  # flexible
```

Step 4:
```
# Unchanging API
# Normal type hierarchy
ra, dec = ..., ...
rep = SphericalRepresentation(ra, dec)

# Deprecated API
c = ICRS(rep)
c = ICRS(ra_arr, dec_arr)
sc = SkyCoord(c)
sc = SkyCoord(ra_arr, dec_arr, frame=ICRS())  # ICRS -> ICRSFrame 
sc = SkyCoord(rep, frame=ICRS())  # ICRS -> ICRSFrame

# New API
frame = ICRSFrame()
c = Coordinate(rep, frame)
sc = SkyCoord(rep, frame=frame)
sc = SkyCoord(c)  # flexible
sc = SkyCoord(ra_arr, dec_arr, frame=frame)  # flexible
```

End result:
```
# Unchanging API
ra, dec = ..., ...
rep = SphericalRepresentation(ra, dec)

frame = ICRSFrame()
c = Coordinate(rep, frame)
sc = SkyCoord(rep, frame=frame)
sc = SkyCoord(c)  # flexible
sc = SkyCoord(ra_arr, dec_arr, frame=frame) # flexible
```
