# Created to analyse the BOUT++ hermes-1 model plasma edge turbulence simulation output data

from boutdata import collect
from boututils.showdata import showdata
from boututils.plotdata import plotdata
import matplotlib.pyplot as plt
import numpy as np
import os

if not os.path.exists("results"):
	os.makedirs("results")
# Collect data from BOUT++
Ne = collect("Ne") # Collects electron density
Sn = collect("Sn") # Neutral density
Pe = collect("Pe") # Collects electron pressure
Phi = collect("phi") #Collects Electrostatic potential
Te = Pe/Ne # Collects electron temperature
Nnorm = collect("Nnorm") # density normalisation factor units will be number density m^-3
Tnorm = collect("Tnorm") # Temperature normalisation factor units will be in eV
Pnorm = Nnorm * Tnorm * 1.602e-19 # Pressure normalisation factor, units will be in pascals

# Average temperature over Z, so it's just a function of time and x:
tetx = np.mean(Te[:,:,0,:], axis=2)
# Average tetx over time, omitting the first 100 time points
tex = np.mean(tetx[100:,:], axis=0)
# Plot temperature, multiplying by tnorm to convert to eV
texTnorm = tex*Tnorm

# Similar for electrostatic potential
phitx = np.mean(Phi[:,:,0,:], axis=2)
phix = np.mean(phitx[100:,:], axis=0)
phixTnorm = phix * Tnorm

# Similar for density
ntx = np.mean(Ne[:,:,0,:], axis=2)
nx = np.mean(ntx[100:,:], axis=0)
nxNnorm = nx * Nnorm

# Similar for pressure
petx =np.mean(Pe[:,:,0,:], axis=2)
pex = np.mean(petx[100:,:], axis=0)
pexNorm = pex * Pnorm

# Coordinates
dx = collect("dx")
g11 = collect("g_11")  # Scale factors (e_x dot e_x)
dz = collect("dz")
g33 = collect("g_33")  # e_z dot e_z
rho_s = collect("rho_s0")  # Distance normalisation

# Grid spacing in meters
dxNorm = dx[0,0] * np.sqrt(g11[0,0]) * rho_s
dzNorm = dz * np.sqrt(g33[0,0]) * rho_s

xcoord1d = np.arange(Ne.shape[1])*dxNorm  # Number of points in x
zcoord1d = np.arange(Ne.shape[-1])*dzNorm

# 2D coordinates for plotting contours
x2d, z2d = np.meshgrid(xcoord1d, zcoord1d, indexing='ij')

# time
t_array = collect("t_array") # 1d array of normalised time
wci = collect("Omega_ci")
time = t_array / wci # Time in seconds
#print("Time duration: {} s".format(time[-1] - time[0]))

# Input power
Spe = collect("Spe")  # dimensionless
# Power per volume [W/m^3]
volume_power = (3./2) * Spe[:,0] * 1.602e-19 * Nnorm * Tnorm * wci
sink_invlpar = collect("sink_invlpar") 
Lpar = rho_s / sink_invlpar[:,0] # meters, connection length
# Power per area [W/m^2]
area_power = volume_power * Lpar

#Open file to output data.
f = open("results/Output.txt","w+")

# Input Power Flux
plt.figure()
plt.title("Input Power Flux")
plt.plot(xcoord1d, area_power * 1e-6) #why is power area multipled by 1e-6
plt.ylabel("Power Flux [MW/m$^{2}$]")
plt.xlabel("Radius [m]")
plt.savefig("results/InputPowerFlux.png")
#plt.show()

# Input Pressure
plt.figure()
plt.title("Input Pressure")
plt.plot(xcoord1d, Spe*Pnorm) #why is power area multipled by 1e-6
plt.ylabel("Pressure [Pa]")
plt.xlabel("Radius [m]")
plt.savefig("results/InputPe.png")
#plt.show()

# Input density
plt.figure()
plt.title("Input Density")
plt.plot(xcoord1d, Sn*Nnorm)
plt.ylabel("Density [m$^{-3}$]")
plt.xlabel("Radius [m]")
plt.savefig("results/InputNe.png")
#plt.show()

# Electrostatic potential contour plot
plt.figure()
plt.contourf(x2d, z2d, Phi[-1,:,0,:]*Tnorm, 50)
plt.colorbar()
plt.title(r"Electrostatic potential [eV]")
plt.xlabel("Radius [m]")
plt.ylabel("Distance [m]")
plt.savefig("results/ElectrostaticContour.png")
#plt.show()


# Electrostatic potential full width half maximum

plt.figure()
plt.plot(xcoord1d, phixTnorm) 
plt.title("Electrostatic potential spreading")
background = phixTnorm[0]  # Left boundary, use as estimate of background
maxind = np.argmax(phixTnorm) # x index where nx is maximum
half_max = background + 0.5 * (phixTnorm[maxind] - background)
plt.axhline(half_max, linestyle='--')
left = phixTnorm[:maxind] # All points to left of peak
right = phixTnorm[maxind:] # To the right
leftind = np.argmin( abs(half_max - left) )
plt.axvline(xcoord1d[leftind])
rightind = np.argmin( abs(half_max - right) ) + maxind
plt.axvline(xcoord1d[rightind])
FWHM = xcoord1d[rightind]-xcoord1d[leftind] # calculate full width at half maximum
innerwidth = xcoord1d[maxind]-xcoord1d[leftind] # calculate the inner width
outerwidth = xcoord1d[rightind]-xcoord1d[maxind] # calculate the outer width

plt.ylabel("Electrostatic potential [eV]")
plt.xlabel("Radius [m]")
plt.savefig("results/ElectrostaticSpreading.png")
f.write("Electrostatic potential spreading\n")
f.write("Peak value of " + "%s" % phixTnorm[maxind] + " eV" + " at x = " + "%s" % xcoord1d[maxind] + " m" + "\n")
f.write("Half max of " + "%s" % half_max + " eV" + " at x1 = " + "%s" % xcoord1d[leftind] + " m" + " and " "at x2 = " + "%s" % xcoord1d[rightind] +" m"+"\n")
f.write("FWHM of " + "%s" % FWHM +" m"+ "\n")
f.write("Inner width of " + "%s" % innerwidth +" m "+ "and "+ "Outer width of " + "%s" %outerwidth + " m"  "\n")
f.write("------------------------------------------------------\n")
#plt.show()


# Density contour plot
plt.figure()
plt.contourf(x2d, z2d, Ne[-1,:,0,:]*Nnorm, 50)
plt.colorbar()
plt.title(r"Electron density [m$^{-3}$]")
plt.xlabel("Radius [m]")
plt.ylabel("Distance [m]")
plt.savefig("results/NeContour.png")
#plt.show()


# Density full width half maximum

plt.figure()
plt.plot(xcoord1d, nxNnorm) 
plt.title("Density spreading")
background = nxNnorm[0]  # Left boundary, use as estimate of background
maxind = np.argmax(nxNnorm) # x index where nx is maximum
half_max = background + 0.5 * (nxNnorm[maxind] - background)
plt.axhline(half_max, linestyle='--')
left = nxNnorm[:maxind] # All points to left of peak
right = nxNnorm[maxind:] # To the right
leftind = np.argmin( abs(half_max - left) )
plt.axvline(xcoord1d[leftind])
rightind = np.argmin( abs(half_max - right) ) + maxind
plt.axvline(xcoord1d[rightind])
FWHM = xcoord1d[rightind]-xcoord1d[leftind] # calculate full width at half maximum
innerwidth = xcoord1d[maxind]-xcoord1d[leftind] # calculate the inner width
outerwidth = xcoord1d[rightind]-xcoord1d[maxind] # calculate the outer width

plt.ylabel("Density [m$^{-3}$]")
plt.xlabel("Radius [m]")
plt.savefig("results/NeSpreading.png")
f.write("Density spreading\n")
f.write("Peak value of " + "%s" % phixTnorm[maxind] + " m^-3" + " at x = " + "%s" % xcoord1d[maxind] + " m" + "\n")
f.write("Half max of " + "%s" % half_max + " m^-3" + " at x1 = " + "%s" % xcoord1d[leftind] + " m" + " and " "at x2 = " + "%s" % xcoord1d[rightind] +" m"+"\n")
f.write("FWHM of " + "%s" % FWHM +" m"+ "\n")
f.write("Inner width of " + "%s" % innerwidth +" m "+ "and "+ "Outer width of " + "%s" %outerwidth + " m"  "\n")
f.write("------------------------------------------------------\n")
#plt.show()


# Temperature contour plot
plt.figure()
plt.contourf(x2d, z2d, Te[-1,:,0,:]*Tnorm, 50)
plt.colorbar()
plt.title(r"Temperature [eV]")
plt.xlabel("Radius [m]")
plt.ylabel("Distance [m]")
plt.savefig("results/TeContour.png")
#plt.show()

# Temperature full width half maximum

plt.figure()
plt.plot(xcoord1d, texTnorm) 
plt.title("Temperature spreading")
background = texTnorm[0]  # Left boundary, use as estimate of background
maxind = np.argmax(texTnorm) # x index where nx is maximum
half_max = background + 0.5 * (texTnorm[maxind] - background)
plt.axhline(half_max, linestyle='--')
left = texTnorm[:maxind] # All points to left of peak
right = texTnorm[maxind:] # To the right
leftind = np.argmin( abs(half_max - left) )
plt.axvline(xcoord1d[leftind])
rightind = np.argmin( abs(half_max - right) ) + maxind
plt.axvline(xcoord1d[rightind])
FWHM = xcoord1d[rightind]-xcoord1d[leftind] # calculate full width at half maximum
innerwidth = xcoord1d[maxind]-xcoord1d[leftind] # calculate the inner width
outerwidth = xcoord1d[rightind]-xcoord1d[maxind] # calculate the outer width

plt.ylabel("Temperature [eV]")
plt.xlabel("Radius [m]")
plt.savefig("results/TeSpreading.png")
f.write("Temperature spreading\n")
f.write("Peak value of " + "%s" % phixTnorm[maxind] + " eV" + " at x = " + "%s" % xcoord1d[maxind] + " m" + "\n")
f.write("Half max of " + "%s" % half_max + " eV" + " at x1 = " + "%s" % xcoord1d[leftind] + " m" + " and " "at x2 = " + "%s" % xcoord1d[rightind] +" m"+"\n")
f.write("FWHM of " + "%s" % FWHM +" m"+ "\n")
f.write("Inner width of " + "%s" % innerwidth +" m "+ "and "+ "Outer width of " + "%s" %outerwidth + " m"  "\n")
f.write("------------------------------------------------------\n")
#plt.show()


# Pressure contour plot
plt.figure()
plt.contourf(x2d, z2d, Pe[-1,:,0,:]*Pnorm, 50)
plt.colorbar()
plt.title(r"Pressure [Pa]")
plt.xlabel("Radius [m]")
plt.ylabel("Distance [m]")
plt.savefig("results/PeContour.png")
#plt.show()

# Pressure full width half maximum

plt.figure()
plt.plot(xcoord1d, pexNorm) #changed from pex to pexNorm for normalised units
plt.title("Pressure spreading")
background = pexNorm[0]  # Left boundary, use as estimate of background
maxind = np.argmax(pexNorm) # x index where nx is maximum
half_max = background + 0.5 * (pexNorm[maxind] - background)
plt.axhline(half_max, linestyle='--')
left = pexNorm[:maxind] # All points to left of peak
right = pexNorm[maxind:] # To the right
leftind = np.argmin( abs(half_max - left) )
plt.axvline(xcoord1d[leftind])
rightind = np.argmin( abs(half_max - right) ) + maxind
plt.axvline(xcoord1d[rightind])
FWHM = xcoord1d[rightind]-xcoord1d[leftind] # calculate full width at half maximum
innerwidth = xcoord1d[maxind]-xcoord1d[leftind] # calculate the inner width
outerwidth = xcoord1d[rightind]-xcoord1d[maxind] # calculate the outer width

plt.ylabel("Pressure [Pa]")
plt.xlabel("Radius [m]")
plt.savefig("results/PeSpreading.png")
f.write("Pressure spreading\n")
f.write("Peak value of " + "%s" % phixTnorm[maxind] + " Pa" + " at x = " + "%s" % xcoord1d[maxind] + " m" + "\n")
f.write("Half max of " + "%s" % half_max + " Pa" + " at x1 = " + "%s" % xcoord1d[leftind] + " m" + " and " "at x2 = " + "%s" % xcoord1d[rightind] +" m"+"\n")
f.write("FWHM of " + "%s" % FWHM +" m"+ "\n")
f.write("Inner width of " + "%s" % innerwidth +" m "+ "and "+ "Outer width of " + "%s" %outerwidth + " m"  "\n")
f.write("------------------------------------------------------\n")
#plt.show()

# Elapsed simulation time
final_time = time[-1]
f.write("Final time " + "%s" % final_time +" seconds\n")
f.write("Final normalised time "+ "%s" % t_array[-1])


