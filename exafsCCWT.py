#!/usr/bin/env python3

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import griddata

"""
Author: Latif Ullah Khan
        Beamline Scientist
        BM08-XAFS/XRF Beamline, SESAME
Date:   December, 03, 2025
"""
# Header and Meta Information
print(" ")
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ")
print("%         CONTINUOUS CAUCHY WAVELET TRANSFORM     % ")
print("%                 OF EXAFS SIGNAL                 % ")
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ")
print(" ")
print("code freely downloaded from http://www.univ-mlv.fr/~farges/waw ")
print("(c) 2000, Univ. Marne la Vallee, France ")
print(" ")
print("Our reference paper for this code is : ")
print("  Munoz M., Argoul P. et Farges F. (2003) ")
print("  Continuous Cauchy wavelet transform analyses of EXAFS spectra: a qualitative approach.  ")
print("  American Mineralogist volume 88, pp. 694-700. ")
print(" ")

# User Preferences
n = 200  # Cauchy order
ri = 0.2  # minimum R-space distance
rf = 6.0  # maximum R-space distance
na = 200  # number of R-space intervals
skip = 0  # number of lines to skip in input file

# Input ascii file
file_in = 'znfoil_k.dat'
kold, xold = np.loadtxt(file_in, unpack=True, comments='#', skiprows=skip)

# Removing file extension
file_in_nodot = file_in.rsplit('.', 1)[0]

# EXAFS data interpolation
print('EXAFS data interpolation...')
nt = 256
#kfin = kold[-1]
kin = 3         #Initial K value
kfin = 11       #Final K value
pask = (kfin - kin) / nt
#pask = (kfin - kold[0]) / nt
#knew = np.arange(kold[0], kold[-1], pask)
knew = np.arange(kin, kfin, pask)
xnew = np.interp(knew, kold, xold)

# Wavelet transform analysis
z = 8
nk = len(knew)
ZF = z * nk
npt = ZF // 2
freq = 1 / pask * np.arange(0, npt) / ZF
omega = 2 * np.pi * freq

# Fourier Transform calculation
tff = np.fft.fft(xnew, ZF)
TF = tff[:npt] / np.max(np.abs(tff))

# Scale parameter
pasr = (rf - ri) / na
r = np.linspace(ri, rf, na)
a = n / (2 * r)

# Characteristic values of the Cauchy wavelet
derpha0 = n + 1
s = sum([math.log(y) for y in range(1, n)])
maxi = math.exp(math.log(2 * np.pi) + n * math.log(n) - n - s)

# Cauchy wavelet calculation
print('Cauchy wavelet calculation...')
filtre = np.zeros((na, npt))
for i in range(na):
    int_a = a[i] * omega
    filtre[i] = np.where(int_a == 0, 0, np.exp(math.log(2 * np.pi) - s + n * np.log(int_a) - int_a))

myttf = tff[:npt]

# Wavelet transform calculation
print('Wavelet transform calculation...')
to = [np.fft.ifft(np.conj(filtre[i]) * myttf, ZF) for i in range(na)]
to = np.array(to)

# Visualization
# Fourier Transform
plt.plot(freq * np.pi, np.abs(TF))
plt.xlabel('R / $\AA$')
plt.ylabel('FT Magnitude')
plt.title('Fourier Transform')
plt.show()

# EXAFS
plt.plot(knew, xnew)
plt.xlabel('k (Å$^{-1}$)')
plt.ylabel('χ(k)')
plt.title('Interpolated EXAFS Data')
plt.show()

# Continuous Cauchy Wavelet Transform
X, Y = np.meshgrid(knew, r)
Z = np.abs(to[:, :nk])  # Ensure to match the shape for interpolation

# Ensure the grid points and values are properly shaped
points = np.array([X.flatten(), Y.flatten()]).T
values = Z.flatten()

# Debugging prints
print(f"points shape: {points.shape}")
print(f"values shape: {values.shape}")

# Check if the number of points matches the number of values
if points.shape[0] != values.shape[0]:
    raise ValueError(f"The number of points ({points.shape[0]}) does not match the number of values ({values.shape[0]})")

# Perform interpolation
zi = griddata(points, values, (X, Y), method='linear')

plt.imshow(zi, extent=(knew.min(), knew.max(), r.max(), r.min()), interpolation='nearest', cmap=cm.jet, aspect='auto')
plt.colorbar()
plt.xlabel('k')
plt.ylabel('R')
plt.title('Continuous Cauchy Wavelet Transform')
plt.show()
