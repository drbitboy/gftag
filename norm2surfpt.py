import os
import sys
import math
import numpy
import pprint
import random
import spiceypy as sp

doLog = "DOLOG" in os.environ

r = random.Random()

########################################################################
### Function to take two input vectors, v0 and v1, and return the
### vector [v0x*v1x,v1, v0y*v1y, v0z*v1z]
vXv = lambda v0,v1: sp.vpack(v0[0]*v1[0],v0[1]*v1[1],v0[2]*v1[2])


########################################################################
### Solve for ellipsoid surface point that has a given normal
def norm2surfpt(semi_axes,inputNormal):

  ### Ellipsoid semi-axes' lengths, per the ellipsoid formula:
  ###
  ###        2           2           2
  ###   / x \       / y \       / z \
  ###  ( --- )  +  ( --- )  +  ( --- ) =  1
  ###   \ a /       \ b /       \ c /
 
  ### Ensure abc components are positive
  abc = sp.vequ([abs(semi_axis) for semi_axis in semi_axes[:3]])

  ### Normalize input normal
  n = sp.vhat(inputNormal)

  ### Calculate two vectors:
  ### - [a*a, b*b, c*c]
  ### - [nx*nx, ny*ny, nz*nz]

  abcXabc = vXv(abc,abc)
  nXn = vXv(n,n)

  ### Direction of normal, N, at [x,y,z] is [ x/(a*a), y/(b*b), z/(c*c)]
  ### - That vector is not necessarily a unit vector
  ### - assume length of N is 1/k
  ### - So scaling unit normal by 1/k will yield N:
  ###     n/k = N = [x/(a*a), y/(b*b), z/(c*c)]
  ###
  ### so
  ###
  ###     nx*k = x/(a*a)
  ###     ny*k = y/(b*b)
  ###     nz*k = z/(c*c)
  ###
  ### and
  ###
  ###     x = nx*a*a/k
  ###     y = ny*b*b/k
  ###     z = nz*c*c/k
  ###
  ### Substituting back into ellipsoid formula:
  ###
  ###     (nx*a*a/k)^2/(a*a) + (ny*b*b/k)^2/(b*b) + (nz*c*c/k)/(c*c) = 1
  ###
  ### and Solving for k
  ###
  ###     (nx*a)^2 + (ny*b)^2 + (nz*c)^2 = k^2
  ###

  k2 = sp.vdot(abcXabc,nXn)
  k = math.sqrt(k2)

  ### Use k, abc and n to calculate surface point as vector xyz
  xyz=sp.vscl(1./k,vXv(abcXabc,n))

  if doLog:
    ### Calculate (x/a)^2 + (y/b)^2 + (z/c)^2; it should be = 1
    one=sp.vdot(vXv(xyz,xyz),1/abcXabc)
    pprint.pprint(dict(n=n,abc=abc,abcXabc=abcXabc,nXn=nXn,nMag=sp.vnorm(n),k=k,k2=k2,xyz=xyz,one=one))

  ### return surface point 
  return xyz


########################################################################
def dooneRandom():
  ### Random normal unit vector
  unitNormal = sp.vhat([r.uniform(-1,1) for i in xrange(3)])

  ### Random positive semi-axes' lengths
  semi_axes = [abs(r.uniform(.5,100)) for i in xrange(3)]

  ### Solve for surface point with that surface normal
  xyz = norm2surfpt(semi_axes,unitNormal)

  ### Use spiceypy.surfnm and xyz to calculate surface unit normal vector,
  ### at solved-for surface point vector xyz, to compare to input unit
  ### normal vector, and calculate error magnitude, errMag

  surfnm=sp.surfnm(semi_axes[0],semi_axes[1],semi_axes[2],xyz)
  errVec=sp.vsub(surfnm,unitNormal)
  errMag = sp.vnorm(errVec)

  if doLog:
    pprint.pprint(dict(errVec=errVec,errMag=errMag))

  ###return error magnitude
  return errMag


########################################################################
########################################################################
### Test code:  output only those error magnitudes above a very small
### level; in practice, this is typically one or two per 100,000 samples

if "__main__"==__name__ and sys.argv[1:]:
  successes = 0
  iterations =int(sys.argv[1])
  for i in xrange(iterations):
    errMag = dooneRandom()
    if abs(errMag>4.7e-16): print((i,errMag))
    if abs(errMag<1e-15): successes += 1
  print('Done:  %d successes and %d failure%s in %d iterations at the 1E-15 level'
        % (successes,iterations-successes
          , ((iterations-successes)!=1) and 's' or '',iterations
          ,))
