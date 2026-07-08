#!/usr/bin/env python3

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import griddata
import plotly.graph_objects as go

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
file_in = 'k2.dat'
kold, xold = np.loadtxt(file_in, unpack=True, comments='#', skiprows=skip)

# Removing file extension
file_in_nodot = file_in.rsplit('.', 1)[0]

# EXAFS data interpolation
print('EXAFS data interpolation...')
nt = 256
#kfin = kold[-1]
kin = 3         #Initial K value
kfin = 10.2       #Final K value
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
    int_a_safe = np.where(int_a == 0, 1, int_a)
    filtre[i] = np.where(int_a == 0, 0, np.exp(math.log(2 * np.pi) - s + n * np.log(int_a_safe) - int_a_safe))

myttf = tff[:npt]

# Wavelet transform calculation
print('Wavelet transform calculation...')
to = [np.fft.ifft(np.conj(filtre[i]) * myttf, ZF) for i in range(na)]
to = np.array(to)

# Visualization - Separate Figures

# Plot font sizes
label_fontsize = 14
tick_fontsize = 12
title_fontsize = 16

# Fourier Transform figure
fig_ft = plt.figure(figsize=(6, 4))
plt.plot(freq * np.pi, np.abs(TF))
plt.xlabel('R (Å)', fontsize=label_fontsize)
plt.ylabel('FT Magnitude', fontsize=label_fontsize)
plt.title('Fourier Transform', fontsize=title_fontsize)
plt.xlim(0, 6)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.tight_layout()
fig_ft.savefig(f'{file_in_nodot}_fourier_transform.png', dpi=300)
plt.show()
plt.close(fig_ft)

# EXAFS figure
fig_exafs = plt.figure(figsize=(6, 4))
plt.plot(knew, xnew)
plt.xlabel('k (Å$^{-1}$)', fontsize=label_fontsize)
plt.ylabel('χ(k)', fontsize=label_fontsize)
plt.title('Interpolated EXAFS Data', fontsize=title_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.tight_layout()
fig_exafs.savefig(f'{file_in_nodot}_exafs_interpolated.png', dpi=300)
plt.show()
plt.close(fig_exafs)

# Continuous Cauchy Wavelet Transform figure
X, Y = np.meshgrid(knew, r)
Z = np.abs(to[:, :nk])  # Ensure to match the shape for interpolation

points = np.array([X.flatten(), Y.flatten()]).T
values = Z.flatten()

print(f"points shape: {points.shape}")
print(f"values shape: {values.shape}")

if points.shape[0] != values.shape[0]:
    raise ValueError(f"The number of points ({points.shape[0]}) does not match the number of values ({values.shape[0]})")

zi = griddata(points, values, (X, Y), method='linear')

fig_wavelet = plt.figure(figsize=(7, 5))
im = plt.imshow(zi, extent=(knew.min(), knew.max(), r.min(), r.max()), origin='lower', interpolation='nearest', cmap=cm.jet, aspect='auto')
cbar = plt.colorbar(im, fraction=0.046, pad=0.02)
cbar.ax.tick_params(labelsize=tick_fontsize)
plt.xlabel('k (Å$^{-1}$)', fontsize=label_fontsize)
plt.ylabel('R (Å)', fontsize=label_fontsize)
plt.title('Continuous Cauchy Wavelet Transform', fontsize=title_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.tight_layout()
fig_wavelet.savefig(f'{file_in_nodot}_wavelet_transform.png', dpi=300)
plt.show()
plt.close(fig_wavelet)

# Interactive Plotly wavelet heatmap with dropdowns and colorbar slider
zmin = np.nanmin(zi)
zmax = np.nanmax(zi)
zmid = 0.5 * (zmin + zmax)
zspread = zmax - zmin
#colorscales = ['Jet', 'Viridis', 'Plasma', 'Cividis', 'HSV', 'YlGnBu', 'Portland', 'Greys', 'Hot', 'Rainbow', 'Earth', 'Electric', 'Picnic', 'Blackbody', 'Bluered', 'YlOrRd', 'YlGn']
colorscales = [
    # Rainbow-like
    'Jet', 'Rainbow', 'HSV', 'Portland',
    'Picnic', 'Earth',
    # Sequential
    'Viridis', 'Plasma', 'Inferno', 'Magma', 'Cividis',
    'Turbo', 'Electric', 'Hot', 'Blackbody',
    # Scientific / smooth
    'YlGnBu', 'YlOrRd', 'Blues', 'Greens',
    'Purples', 'Reds', 'Greys',
    # Diverging
    'RdBu', 'RdYlBu', 'RdYlGn',
    'Spectral', 'Balance',  
    # Cyclic
    'Phase', 'Twilight',
    # High contrast
    'Bluered'
]
contrast_levels = np.linspace(0.4, 1.0, 7)

fig_plotly = go.Figure(
    data=[go.Heatmap(
        z=zi,
        x=knew,
        y=r,
        colorscale='Jet',
        zmin=zmin,
        zmax=zmax,
        colorbar=dict(
            title=dict(text='Amplitude', font=dict(size=label_fontsize + 22)),
            tickfont=dict(size=tick_fontsize + 22),
            thickness=20,
            len=1.0,
            lenmode='fraction',
            y=1.0,
            yanchor='top',
            x=1.0
        )
    )]
)

color_buttons = [
    dict(
        method='restyle',
        label=scale,
        args=[{'colorscale': [scale]}]
    )
    for scale in colorscales
]

slider_steps = []
for level in contrast_levels:
    zmin_step = zmid - 0.5 * zspread * level
    zmax_step = zmid + 0.5 * zspread * level
    slider_steps.append(
        dict(
            label=f'{int(level * 100)}%',
            method='restyle',
            args=[{'zmin': [zmin_step], 'zmax': [zmax_step]}]
        )
    )

fig_plotly.update_layout(
    title=dict(text='Continuous Cauchy Wavelet Transform', x=0.5, xanchor='center', font=dict(size=title_fontsize + 4)),
    xaxis=dict(title=dict(text=r'k (Å<sup>-1</sup>)', font=dict(size=label_fontsize + 22)), tickfont=dict(size=tick_fontsize + 22), automargin=True),
    yaxis=dict(title=dict(text='R (Å)', font=dict(size=label_fontsize + 22)), tickfont=dict(size=tick_fontsize + 22), automargin=True),
    margin=dict(l=480, r=480, t=120, b=120),
    updatemenus=[
        dict(
            buttons=color_buttons,
            direction='down',
            pad={'r': 10, 't': 10},
            showactive=True,
            x=0.0,
            xanchor='left',
            y=1.15,
            yanchor='top',
            bgcolor='lightgrey',
            bordercolor='grey',
            borderwidth=1,
            active=0
        )
    ],
    sliders=[
        dict(
            active=len(slider_steps) - 1,
            currentvalue={'prefix': 'Contrast: '},
            pad={'t': 60},
            steps=slider_steps,
            bordercolor='lightgrey',
            borderwidth=1,
            x=0.0,
            xanchor='left',
            y=-0.10,
            yanchor='top'
        )
    ]
)

plotly_html = f'{file_in_nodot}_wavelet_transform_plotly.html'
fig_plotly.write_html(plotly_html, include_plotlyjs='cdn')
print(f'Interactive Plotly wavelet map saved to {plotly_html}')
fig_plotly.show()
