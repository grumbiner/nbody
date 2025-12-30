''' earth plus some random moons. Number is an argument, followed by seed of rng '''

import sys

from math import ceil
import numpy as np

import astronomy
from point_mass import *

################################################################################
nbodies = int(sys.argv[1])
try:
  n = int(sys.argv[2])
except:
  n = 0
np.random.seed(n)

fout    = open("diagnostics."+sys.argv[2],"w", encoding="utf-8")
fdiag   = open("diag."+sys.argv[2],"w", encoding="utf-8")
fenergy = open("energy."+sys.argv[2], "w", encoding="utf-8")
fcom    = open("com."+sys.argv[2], "w", encoding="utf-8")


################################################################################
# Establish some objects and their initial conditions:
solar_system = point_system(nbody = nbodies)
#-------------- run control ---------------------
collision_tolerance = 2.5*astronomy.radius_moon

# Time step et al.:
dtref = 10.0
ratio = (24.*3600.)/dtref

freq  = int(3600./2./dtref + 0.5)
#diagnostic:
print("dt = ",dtref, "freq = ",freq, flush=True, file = fdiag)

nyears = 4
days = 365
#------------ end run control -------------------


bodies  = []
earth   = point_mass()
moon    = point_mass()

bodies.append(earth)
bodies.append(moon)

for i in range(2, nbodies):
  tmp = point_mass()
  bodies.append(tmp)
  del tmp

bodies[0].m   = astronomy.m_earth
for i in range(1, nbodies):
  bodies[i].m   = astronomy.m_moon

for i in range(0, nbodies):
  bodies[i].k   = bodies[i].m * astronomy.G

bodies[0].init_loc(0., 0., 0)
bodies[1].init_loc(astronomy.rmoon, 0., 0.)
bodies[1].u[1] = -kepler(bodies[1], bodies[0])

# start alternate moons at some multiples of real moon's distance and speed
for i in range(2, nbodies):
  tmp = np.random.uniform(low = 0.3, high = 3.0)
  bodies[i].init_loc(tmp*astronomy.rmoon, 0., 0.)

  tmp = np.random.uniform(low = 0.2, high = sqrt(2.) )
  bodies[i].u[1] = -tmp*kepler(bodies[i], bodies[0])
  print("body ",i,"start ",bodies[i].x[0]/astronomy.rmoon,
        "kepler multiplier: ",tmp, file = fout, flush=True)

# --- Make a system of the bodies ---------------------
for i in range(0, nbodies):
  solar_system.body[i] = copy.deepcopy(bodies[i])

# Center of Mass:
com = point_mass()
com = copy.deepcopy(solar_system.center_of_mass() )

# adopt comoving with center of mass, with center the center of mass at 0,0,0
solar_system.comoving()

# system energy
solar_system.conservation()
initial_totke = solar_system.ke()
initial_totpe = solar_system.PE()
initial_energy = initial_totke + initial_totpe

r0 = np.zeros((solar_system.nbody))
for i in range (0, solar_system.nbody):
  r0[i] = dist(solar_system.body[i],solar_system.com)

#-------------- ---------------    ---------------------

for i in range (0, int(days*int(ratio)*nyears) + 1):
  # adaptive dt -- note it is inverse cube in separation
  tmp = solar_system.closest()
  if (tmp < collision_tolerance):
    solar_system.collision(tolerance = collision_tolerance, fout = fdiag)
    initial_totke = solar_system.tke
    initial_totpe = solar_system.tpe
    initial_energy = initial_totke + initial_totpe
    #debug: print("collision at step ",i,"time ",i*dt, i*dt/60., i*dt/3600.)
    #debug: exit(1)

  kappa = 10./(tmp/astronomy.rmoon)**3
  mul = max(1., kappa/500.)
  #diagnostic: print(i, tmp, dtref/ceil(mul), mul, kappa, file=fout)

  dt = dtref/mul
  for ifast in range (0, ceil(mul)):
    solar_system.step(dt)

  #output:
  if (i%freq == 0):
    print(float(i)*(dtref/astronomy.mean_solar_day), end=" ")
    for j in range(0, solar_system.nbody):
      print(  (solar_system.body[j].x[0])/astronomy.rmoon,
              (solar_system.body[j].x[1])/astronomy.rmoon,
              (dist(solar_system.body[j], solar_system.com ) - r0[j])/astronomy.rmoon , end=" ")
    print(flush=True )

    totpe = solar_system.PE()
    totke = solar_system.ke()
    print(i, (totke+totpe-initial_energy)/initial_energy,
                solar_system.closest()/astronomy.rmoon, flush=True, file = fenergy)


#addl diagnostics:
# equilibrium tide magnitude
# angular momentum and conservation
