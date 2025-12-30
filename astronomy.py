''' astronomy.py -- hold constants and define some utility functions '''

from math import sqrt

class astronomy:
   ''' Class astronomy -- just holds some constants for now '''
# Universal constants:
   G              = 6.6743e-11
   mean_solar_day = 86400.

# kg
   m_sun     = 1.98847e30
   m_earth   = 5.9722e24
   m_jupiter = m_earth * 317.8
   m_moon    = m_earth *   0.0123

# meters
   radius_moon  = 1737.4e3
   radius_earth = 6371.1e3

# meters:
   rmoon   = 3.84399e8
   au      = 1.49597870700e11
   ly      = 9.46e15
   parsec  = 3.26*ly

   def __init___(self):
     return

#-----------------------------------------------------------------
# Utilities

def dist(x,y):
  ''' dist(x,y) -- compute the distance between two bodies x,y '''
#  dx = y.x - x.x
#  r = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2])
# Almost 3x faster this way
  r = sqrt((y.x[0]-x.x[0])**2 + (y.x[1]-x.x[1])**2 + (y.x[2]-x.x[2])**2 )
  return r

#keplerian speed for distance around body y:
def kepler(x, y):
  ''' kepler(x,y) -- compute the keplerian speed for body y around body x '''
  #debug: print("in kepler",flush=True)
  #debug: x.show()
  #debug: y.show()
  dx1 = y.x[0] - x.x[0]
  dx2 = y.x[1] - x.x[1]
  dx3 = y.x[2] - x.x[2]
  r = sqrt(dx1*dx1 + dx2*dx2 + dx3*dx3)
  return sqrt(astronomy.G*y.m/r)

# Roche limit for two real bodies:
def roche(x, y, rad = astronomy.radius_moon):
  ''' roche(x,y,rad) -- compute the roche limit for bodies x, y with radius rad '''
  mmax = max(x.m, y.m)
  mmin = min(x.m, y.m)
  return rad*(2.*mmax/mmin)**(1./3.)

#-----------------------------------------------------------------
#x = 'this' body, y = remote
def gravity(x, y, accel):
  ''' gravity(x,y,accel) -- compute gravitational acceleration of body x due to y ''' 
  dx1 = y.x1 - x.x1
  dx2 = y.x2 - x.x2
  dx3 = y.x3 - x.x3
  r  = sqrt(dx1*dx1 + dx2*dx2 + dx3*dx3)
  a0 = astronomy.G * y.m / r / r / r
  accel[0] = dx1 * a0
  accel[1] = dx2 * a0
  accel[2] = dx3 * a0
