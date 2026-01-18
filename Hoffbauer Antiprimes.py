#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Hoffbauer Antiprimes - Berechnung mit SageMath
"""
# SageMath
from sage.all import *
# Sicherstellen, dass CC verfügbar ist (ComplexField für komplexe Zahlen)
if 'CC' not in dir():
    CC = ComplexField(53)  # Standard-Genauigkeit für komplexe Zahlen
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------- Familien (deine Definition) ----------
def fam_mod12(p):
    r = int(p % 12)
    if r == 1:  return "e"
    if r == 5:  return "a"
    if r == 7:  return "b"
    if r == 11: return "c"
    return "x"

def fam_phase(p):
    mp = {"e":0.0, "a":np.pi/2, "b":np.pi, "c":3*np.pi/2, "x":0.0}
    return mp[fam_mod12(p)]

# ---------- Hoffbauer-Hilfen ----------
def antiprime_from_gap(p_i, p_next):
    g = int(p_next - p_i)
    return int(p_i + g//2), g

def index_window(i, N, w):
    lo = max(0, i - w)
    hi = min(N-1, i + w)
    return range(lo, hi+1)

# ---------- A_theta aus Hoffbauer-Schnitt ----------
def A_theta_hoffbauer(theta, N=200, w=3, mu=0.35, kappa=0.55, eta=0.20,
                      use_family_phase=True, normalize=True):
    P = [int(nth_prime(i+1)) for i in range(N+1)]
    primes = P[:N]
    primes_next = P[1:N+1]

    gaps, antip = [], []
    for i in range(N):
        a_i, g_i = antiprime_from_gap(primes[i], primes_next[i])
        antip.append(a_i)
        gaps.append(g_i)

    lg = np.log(np.array(gaps, dtype=float))
    lg = (lg - lg.mean()) / (lg.std() + 1e-12)

    A = zero_matrix(CC, N, N)

    # Diagonal: Gap-Potential
    for i in range(N):
        A[i,i] = CC(mu * lg[i])

    # Kette: gerichteter Fluss (mit Holonomie theta)
    ph_edge = complex(np.cos(theta), np.sin(theta))
    for i in range(N-1):
        w_i = kappa / np.sqrt(gaps[i] + 1e-12)
        if use_family_phase:
            extra = np.exp(1j * (fam_phase(primes[i]) - fam_phase(primes[i+1])))
        else:
            extra = 1.0
        A[i, i+1] += CC(w_i * ph_edge * extra)
        A[i+1, i] += CC(w_i * np.conjugate(ph_edge) * np.conjugate(extra))

    # Hoffbauer-Fensterkopplungen über Antiprime-Nähe
    antip_arr = np.array(antip, dtype=float)
    scale = np.median(np.abs(np.diff(antip_arr))) + 1e-12

    for i in range(N):
        for j in index_window(i, N, w):
            if j == i:
                continue
            dist = abs(antip_arr[i] - antip_arr[j]) / scale
            wij = eta * np.exp(-dist) / np.sqrt((gaps[i]+1e-12)*(gaps[j]+1e-12))
            if use_family_phase:
                extra = np.exp(1j * (fam_phase(primes[i]) - fam_phase(primes[j])))
            else:
                extra = 1.0
            A[i,j] += CC(wij * ph_edge * extra)

    if normalize:
        nr = float(max(np.abs(np.array(A.list(), dtype=complex))))
        if nr > 0:
            A = (1.0/nr) * A
    return A

def make_chiral_dirac(A):
    N = A.nrows()
    Z = zero_matrix(CC, N, N)
    return block_matrix([[Z, A],
                         [A.conjugate_transpose(), Z]])

def det_phase(D, t=0.0, eps=1e-6):
    n = D.nrows()
    I = identity_matrix(CC, n)
    Z = D - t*I + (CC(0, eps))*I
    val = Z.det()
    return np.angle(complex(val))

def winding_number(N=200, thetas=361, t=0.0, eps=1e-6, w=3, params=None):
    th = np.linspace(0.0, 2.0*np.pi, thetas)
    phases = []
    for theta in th:
        A = A_theta_hoffbauer(theta, N=N, w=w, **(params or {}))
        D = make_chiral_dirac(A)
        phases.append(det_phase(D, t=t, eps=eps))
    phases = np.unwrap(np.array(phases))
    W = (phases[-1] - phases[0])/(2.0*np.pi)
    return float(W)

def spectral_flow(N=200, thetas=361, t=0.0, w=3, params=None):
    th = np.linspace(0.0, 2.0*np.pi, thetas)
    evs = []
    for theta in th:
        A = A_theta_hoffbauer(theta, N=N, w=w, **(params or {}))
        D = make_chiral_dirac(A) - t*identity_matrix(CC, 2*N)
        lam = np.array([float(rr) for rr in D.eigenvalues()])
        lam.sort()
        evs.append(lam)
    evs = np.array(evs)

    sf = 0
    for j in range(evs.shape[1]):
        s = np.sign(evs[:, j])
        s[s == 0] = 1
        flips = np.where(s[1:] != s[:-1])[0]
        for k in flips:
            if evs[k, j] < 0 and evs[k+1, j] > 0: sf += 1
            if evs[k, j] > 0 and evs[k+1, j] < 0: sf -= 1
    return int(sf)

# ---------- 1) Robustheits-Grid ----------
N0 = 200
thetas = 361

grid_w  = [2,3,4]
grid_eps = [1e-4, 1e-6, 1e-8]

rows = []
for w in grid_w:
    for eps in grid_eps:
        W = winding_number(N=N0, thetas=thetas, t=0.0, eps=eps, w=w)
        SF = spectral_flow(N=N0, thetas=thetas, t=0.0, w=w)
        rows.append({"N":N0, "w":w, "eps":eps, "W(t=0)":W, "SF(t=0)":SF})

df = pd.DataFrame(rows)
print(df)

# ---------- 2) "137-Umfeld"-Diagnose ----------
# Idee: wir betrachten zwei Cutoffs:
#   A) N=200 (enthält p=137 sicher)
#   B) N=200, aber wir "maskieren" die ersten m Knoten (setzen Kanten dort auf 0),
#      um zu sehen, ob die frühen Strukturen (inkl. 137-Region) topologisch dominieren.
def A_theta_masked(theta, N=200, mask_first=0, **kwargs):
    A = A_theta_hoffbauer(theta, N=N, **kwargs)
    if mask_first <= 0:
        return A
    # maskiere Zeilen/Spalten im A-Block -> entfernt frühe Knoten dynamisch
    for i in range(mask_first):
        for j in range(N):
            A[i,j] = 0
            A[j,i] = 0
    return A

def winding_number_masked(mask_first, N=200, w=3, eps=1e-6, thetas=361):
    th = np.linspace(0.0, 2.0*np.pi, thetas)
    phases = []
    for theta in th:
        A = A_theta_masked(theta, N=N, mask_first=mask_first, w=w)
        D = make_chiral_dirac(A)
        phases.append(det_phase(D, t=0.0, eps=eps))
    phases = np.unwrap(np.array(phases))
    return float((phases[-1]-phases[0])/(2*np.pi))

masks = [0, 20, 40, 60, 80]  # grob: "entferne" frühe Bereiche (137 liegt bei Index ~33)
Wmask = [winding_number_masked(m, N=N0, w=3, eps=1e-6, thetas=thetas) for m in masks]

plt.figure()
plt.plot(masks, Wmask, marker="o")
plt.xlabel("mask_first (Anzahl entfernter früher Knoten)")
plt.ylabel("Winding W bei t=0")
plt.title("Hoffbauer-Dirac: Sensitivität von W gegen frühe Regionen (inkl. 137)")
plt.grid(True)
plt.show()

# Plot Robustheits-Grid als Heatmap-ähnliche Darstellung (einfach)
pivot = df.pivot(index="w", columns="eps", values="W(t=0)")
plt.figure()
plt.imshow(pivot.values, aspect="auto")
plt.xticks(range(len(pivot.columns)), [str(e) for e in pivot.columns])
plt.yticks(range(len(pivot.index)), [str(w) for w in pivot.index])
plt.xlabel("eps")
plt.ylabel("w")
plt.title("W(t=0) über (w, eps)")
plt.colorbar()
plt.show()