''' solar_system.py -- roughly the solar system with the sun, earth, jupiter, venus, and Saturn '''

from math import ceil
import astronomy
from point_mass import *

################################################################################
fout    = open("diagnostics7","w",encoding="utf-8")
fdiag   = open("diag7","w",encoding="utf-8")
fenergy = open("energy7", "w",encoding="utf-8")
#fcom    = open("com7", "w",encoding="utf-8")

np.random.seed(0)

#-------------- begin run control ---------------------
collision_tolerance = 2.5*astronomy.radius_moon

# Time step et al.:
dtref = 3600.0*6
ratio = (24.*3600.)/dtref

#Every half hour
#freq  = int(3600./2./dtref + 0.5)
#Every 6 hours
freq  = int(86400./dtref + 0.5)
#diagnostic:
print("dt = ",dtref, "freq = ",freq, "ratio = ",ratio, flush=True, file = fdiag)

# Run length ~= 1 hour on 1 core 58,056 years
nyears = int(59.*984)
days = int(365)

#--------------- specify solar system ----------------------------------
nbodies  = 5
solar_system = point_system(nbody = nbodies)

solar_system.body[0].m = astronomy.m_sun
solar_system.body[0].init_loc(0.*astronomy.au, 0., 0.)

solar_system.body[1].m = astronomy.m_earth
solar_system.body[1].init_loc(1.*astronomy.au, 0., 0.)

solar_system.body[2].m = astronomy.m_jupiter
solar_system.body[2].init_loc(5.2038*astronomy.au, 0., 0.)

solar_system.body[3].m = astronomy.m_earth * 0.815
solar_system.body[3].init_loc(0.723332*astronomy.au, 0., 0.)

solar_system.body[4].m = astronomy.m_earth * 95.159
solar_system.body[4].init_loc(9.5826*astronomy.au, 0., 0.)


# Center of Mass:
com = point_mass()
com = copy.deepcopy(solar_system.center_of_mass() )

for i in range(0, solar_system.nbody):
  solar_system.body[i].k = solar_system.body[i].m*astronomy.G
  if (i > 0):
    solar_system.body[i].u[1] = kepler(solar_system.body[i], com)
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

#-------------- end run control ---------------------

for i in range (0, int(days*int(ratio)*nyears) + 1):
  # adaptive dt -- note it is inverse cube in separation
  tmp   = solar_system.closest()
  kappa = (tmp/astronomy.au)**3
  mul   = max(1., 1.5e-5/kappa)

  dt = dtref/mul
  for ifast in range (0, ceil(mul)):
    solar_system.step(dt)

  #output:
  if (i%freq == 0):
    #time step diagnostic:
    print(i, tmp, dtref/ceil(mul), mul, kappa, file=fout)

    print(float(i)*(dtref/astronomy.mean_solar_day), end=" ")
    for j in range(0, solar_system.nbody):
      print(  (solar_system.body[j].x[0])/astronomy.au,
              (solar_system.body[j].x[1])/astronomy.au,
              (dist(solar_system.body[j], solar_system.com ) - r0[j])/astronomy.au , end=" ")
    print(flush=True )

    totpe = solar_system.PE()
    totke = solar_system.ke()
    print(i, (totke+totpe-initial_energy)/initial_energy,
                solar_system.closest()/astronomy.au, flush=True, file = fenergy)
