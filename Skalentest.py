import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from sympy import next_prime, prime_pi

# --- Hilfsfunktionen (aus vorherigen Schritten) ---
def get_class(p):
    m = p % 12
    return {5:'A', 7:'B', 11:'C', 1:'D'}.get(m, '?')

def phi_a(t, a):
    if abs(t) < 1e-9: return 2.0 * a
    return (2.0/a) * (np.sin(a*t)/t)**2

# --- Der Kern: Optimierung für einen bestimmten Bereich ---
def find_optimal_epsilon(center_p, window=200, a=2.0):
    # 1. Daten holen
    p_start = next_prime(center_p - window)
    p_block = []
    curr = p_start
    for _ in range(window): # Hole ca 'window' Primzahlen
        p_block.append(curr)
        curr = next_prime(curr)
    
    # Simuliere Nullstellen (In Realität: zeros1.gz an die Stelle seeken)
    # t_n approx 2*pi*n / log(n)
    start_n = prime_pi(center_p)
    t_block = []
    for i, p in enumerate(p_block):
        n = start_n + i
        # Einfache Riemann-Siegel Näherung für die Simulation
        t_val = 14.1347 + (2 * np.pi * n) / np.log(n) 
        t_block.append(t_val)

    # 2. Zielfunktion (Trace-Fehler)
    rhs_target = 0 # Wir minimieren hier nur die LHS-Variation relativ zur Struktur
    # (Vereinfachung: Wir suchen das Eps, das die "interne" Resonanz maximiert)
    
    def target_func(eps):
        N = len(t_block)
        mid = N // 2
        target_class = get_class(p_block[mid])
        
        R = np.zeros((N, N))
        classes = [get_class(p) for p in p_block]
        
        for i in range(N):
            for j in range(i+1, N):
                if i == mid or j == mid: continue
                c1, c2 = classes[i], classes[j]
                res = False
                if {c1, c2} == {'B', 'C'} and target_class == 'A': res = True
                if {c1, c2} == {'A', 'C'} and target_class == 'B': res = True
                if {c1, c2} == {'A', 'B'} and target_class == 'C': res = True
                
                if res:
                    R[i, j] = 1.0; R[j, i] = -1.0
        
        # Singulärwerte
        A = np.diag(t_block)
        vals = np.linalg.svd(A + eps * R, compute_uv=False)
        
        # Wir messen hier die "Spektrale Schärfe" (Entropie-Minimierung)
        # Stabilized calculation to avoid overflow with large singular values
        beta = 0.01
        max_v = np.max(vals)
        # Z = sum(exp(beta*v) + exp(-beta*v))
        # LogZ = log(sum(exp(beta*v) + exp(-beta*v)))
        # Using log-sum-exp trick: log(sum(exp(x_i))) = M + log(sum(exp(x_i - M)))
        # Here x_i are beta*v and -beta*v
        combined_vals = np.concatenate([beta * vals, -beta * vals])
        M = np.max(combined_vals)
        logZ = M + np.log(np.sum(np.exp(combined_vals - M)))
        
        # dZ/Z = sum(v * (exp(beta*v) - exp(-beta*v))) / sum(exp(beta*v) + exp(-beta*v))
        # This is equivalent to E[v * tanh_like_weight]
        # weights = exp(beta*v - logZ) - exp(-beta*v - logZ)
        # Note: z_i = exp(beta*v_i) / Z and w_i = exp(-beta*v_i) / Z
        z_log_probs = beta * vals - logZ
        w_log_probs = -beta * vals - logZ
        exp_v = np.sum(vals * (np.exp(z_log_probs) - np.exp(w_log_probs)))
        
        S = logZ - beta * exp_v
        return S

    # Suche Minimum
    res = minimize_scalar(target_func, bounds=(0.0, 0.5), method='bounded')
    return res.x

# --- Der Scaling-Loop ---
scales = [1000, 5000, 10000, 50000, 100000, 500000, 1000000] # Bis 10^6
epsilons = []

print("Starte Renormierungsgruppen-Analyse...")
for s in scales:
    eps = find_optimal_epsilon(s)
    epsilons.append(eps)
    print(f"Höhe {s}: Optimales Epsilon = {eps:.5f}")

# Plot
plt.figure(figsize=(10, 6))
plt.semilogx(scales, epsilons, 'o-', color='crimson', linewidth=2)
plt.title(r'Renormierung der Bamberger Konstante $\epsilon$')
plt.xlabel('Größenordnung der Primzahlen (log scale)')
plt.ylabel(r'Optimale Kopplungsstärke $\epsilon_{opt}$')
plt.grid(True, which="both", ls="-", alpha=0.4)

# Trendlinie
z = np.polyfit(np.log(scales), epsilons, 1)
p = np.poly1d(z)
plt.plot(scales, p(np.log(scales)), "b--", label=f"Trend: {z[0]:.2e} * log(x)")

plt.legend()
plt.savefig('epsilon_scaling.png')
print("Analyse abgeschlossen.")