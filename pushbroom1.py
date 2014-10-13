import numpy as np

'''
Cartesian orbital elements calculated here
start with stupid circular orbit with no perturbations
c is camera, p is central body
'''
#oribtal elements
omega = np.random.random(1)*2*pi #argument of periapsis
I =  np.arccos(np.random.random(1)*2 - 1) #inclination
Omega = np.random.random(1)*2*pi #longitude of ascending node

#edist = lambda e: (e/0.0281)*np.exp((-e**2.)/2/0.0281) #rayleigh distributed eccentricity
#e = simpSample(edist,1,0,1)
e = 0

a = 1 #semi-major axis 1 distance unit

#generate anomaly (sampling orbit at some number of points)
nsamp = 1000
M = np.linspace(0,2*pi,nsamp)
#E = solveEccAnom(M, e) #eccentric 
E = M
r = a*(1-e*cos(E)) #orbital radius (DU)

# gravitational parameter:
#G = 6.67428e-11 #m^3/kg/s^2
#mAU = 1.495978707e11 #AU in meters
#mSun = 1.98892e30 #kg
#G = G/mAU**3.*86400**2.*mSun #AU^3/mSun/day^2
G = 1 #DU^3/MU/TU^2


#generate planetocentric orbital positions and velocities
A = np.vstack((a*(cos(Omega)*cos(omega) - sin(Omega)*cos(I)*sin(omega)),
     a*(sin(Omega)*cos(omega) + cos(Omega)*cos(I)*sin(omega)),
     a*sin(I)*sin(omega))).T

B = np.vstack((-a*sqrt(1-e**2)*(cos(Omega)*sin(omega) + sin(Omega)*cos(I)*cos(omega)),
      a*sqrt(1-e**2)*(-sin(Omega)*sin(omega) + cos(Omega)*cos(I)*cos(omega)),
      a*sqrt(1-e**2)*sin(I)*cos(omega))).T
 
mu = 1 #combined gravitational parameter
r_cp = A*np.tile(cos(E) - e,(3,1)).T+ B*np.tile(sin(E),(3,1)).T
v_cp = np.tile(sqrt(mu/a**3)/(1 - e*cos(E)),(3,1)).T*(-A*np.tile(sin(E),(3,1)).T+ B*np.tile(cos(E),(3,1)).T)

'''
Linear pushbroom model take 1
based on Gupta & Hartley 1997
'''

f = 1 #focal length/magnification
p_v = 0.1 #principal point offset in v dir

tmat = np.matrix([[1, 0, 0], [0, f, p_v], [0, 0, 1]])

#for now assume camera frame is aligned with world frame
camMat = np.zeros(r_cp.shape)   #[u,wvw]^T

for j,v in enumerate(v_cp):
    vmat = np.matrix([[1/v[0], 0, 0], [-v[1]/v[0], 1, 0], [-v[2]/v[0], 0, 1]])
    camMat[j] = (tmat*vmat*np.matrix(r_cp[j]).T).squeeze()


#(u,v) coords of centerline point over course of orbit
imcoords = np.vstack((camMat[:,0],camMat[:,1]/camMat[:,2]))
