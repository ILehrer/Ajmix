import numpy as np
import pandas as pd
import filter
import matplotlib.pyplot as plt
import sys
from scipy.optimize import curve_fit

CHROM = int(sys.argv[1]) if len(sys.argv) > 1 else 18
USE_GLOBAL = len(sys.argv) > 2 and sys.argv[2] == 'global'
SMOOTH = True
EPS = 1e-10

# 10kb->100kb, 1kb increment
WINDOW_SIZES = np.arange(10_000, 101_000, 1000)

def load_recombo_map(chrom):
    rmap = pd.read_csv('sex-averaged_noncarrier.rmap.txt', sep='\t')
    rmap = rmap[rmap['chr'] == f'chr{chrom}'].copy()
    # Each bin is 0.01 Mb. 
    rmap['cM_delta'] = rmap['stdrate'] * 0.01 
    rmap['cM_cum'] = rmap['cM_delta'].cumsum()
    return rmap

def get_genetic_dist(pos_start, pos_end, rmap):
    """Interpolates cumulative cM to find Morgans between two physical points."""
    cM_start = np.interp(pos_start, rmap['pos'], rmap['cM_cum'])
    cM_end = np.interp(pos_end, rmap['pos'], rmap['cM_cum'])
    return (cM_end - cM_start) / 100.0

def run_hmm_iteration(df, rmap):
    n_obs = len(df)
    n_states = 2
    pi = np.array([0.5, 0.5])
    g = 30.0  # Initial guess
    
    pos_array = df['pos'].values
    dist_morgans = np.zeros(n_obs)
    for i in range(1, n_obs):
        if USE_GLOBAL:
            dist_morgans[i] = (pos_array[i] - pos_array[i-1]) * 1e-8
        else:
            dist_morgans[i] = get_genetic_dist(pos_array[i-1], pos_array[i], rmap)
    dist_morgans[0] = 0.0005 # Baseline jump for start of chromosome
    
    E = df[['E0_smooth', 'E1_smooth']].values + EPS
    for _ in range(10):
        # Forward Pass
        alpha = np.zeros((n_obs, n_states))
        c = np.zeros(n_obs)
        alpha[0] = pi * E[0]
        c[0] = 1.0 / (np.sum(alpha[0]) + EPS)
        alpha[0] *= c[0]
        for t in range(1, n_obs):
            p_switch = 0.5 * (1 - np.exp(-2 * g * dist_morgans[t]))
            A = np.array([[1-p_switch, p_switch], [p_switch, 1-p_switch]])
            alpha[t] = (alpha[t-1] @ A) * E[t]
            c[t] = 1.0 / (np.sum(alpha[t]) + EPS)
            alpha[t] *= c[t]
            
        # Backward Pass
        beta = np.zeros((n_obs, n_states))
        beta[-1] = np.ones(n_states) * c[-1]
        for t in range(n_obs - 2, -1, -1):
            p_switch = 0.5 * (1 - np.exp(-2 * g * dist_morgans[t+1]))
            A = np.array([[1-p_switch, p_switch], [p_switch, 1-p_switch]])
            beta[t] = (A @ (E[t+1] * beta[t+1])) * c[t]
            
        # Expectation & Maximization
        gamma = (alpha * beta)
        gamma /= np.sum(gamma, axis=1, keepdims=True)
        total_switches = 0
        sum_morgans = 0
        for t in range(n_obs - 1):
            p_switch = 0.5 * (1 - np.exp(-2 * g * dist_morgans[t+1]))
            A = np.array([[1-p_switch, p_switch], [p_switch, 1-p_switch]])
            xi_t = (alpha[t][:, None] * A * E[t+1] * beta[t+1])
            xi_t /= (np.sum(xi_t) + EPS)
            total_switches += (xi_t[0, 1] + xi_t[1, 0])
            sum_morgans += dist_morgans[t+1]
        g = total_switches / (2 * sum_morgans + EPS)
    return g

def exponential_decay(x, a, b, c):
    return a * np.exp(-b * x) + c

print(f"Sweeping Chr {CHROM}...")
raw_df = filter.get_df(CHROM)
rmap_data = None if USE_GLOBAL else load_recombo_map(CHROM)

# Precalc (log-)likelihoods for SNPs
af_eu, af_me = raw_df['af_eu'].values, raw_df['af_me'].values
obs = raw_df['dosage'].values
p_eu_all = np.vstack([(1-af_eu)**2, 2*af_eu*(1-af_eu), af_eu**2])
p_me_all = np.vstack([(1-af_me)**2, 2*af_me*(1-af_me), af_me**2])
raw_df['log_p_eu'] = np.log(p_eu_all[obs, np.arange(len(raw_df))] + EPS)
raw_df['log_p_me'] = np.log(p_me_all[obs, np.arange(len(raw_df))] + EPS)

results_g = []
mode_str = 'global' if USE_GLOBAL else 'local'

print(f"Starting {mode_str} sensitivity sweep (10kb to 100kb)...")
for w_size in WINDOW_SIZES:
    raw_df['window'] = raw_df['pos'] // w_size
    win_df = raw_df.groupby('window').agg({'pos': 'mean', 'log_p_eu': 'sum', 'log_p_me': 'sum'}).reset_index()
    max_log = win_df[['log_p_eu', 'log_p_me']].max(axis=1)
    
    # Apply 3-window smoothing kernel
    win_df['E0_smooth'] = np.exp(win_df['log_p_eu'] - max_log).rolling(3, center=True).mean().fillna(np.exp(win_df['log_p_eu'] - max_log))
    win_df['E1_smooth'] = np.exp(win_df['log_p_me'] - max_log).rolling(3, center=True).mean().fillna(np.exp(win_df['log_p_me'] - max_log))
    
    results_g.append(run_hmm_iteration(win_df, rmap_data))
    if w_size % 10_000 == 0:
        print(f"Window: {w_size/1000}kb | g: {results_g[-1]:.2f}")

# --- FITTING & PLOTTING ---
x_data = WINDOW_SIZES / 1000.0
y_data = np.array(results_g)

plt.figure(figsize=(10, 6))
plt.plot(x_data, y_data, 'bo', label='Data Points', markersize=4)

try:
    # Initial guess for [a, b, c]
    p0 = [y_data[0] - y_data[-1], 0.1, y_data[-1]]
    popt, _ = curve_fit(exponential_decay, x_data, y_data, p0=p0, maxfev=5000)
    
    x_fit = np.linspace(x_data.min(), x_data.max(), 500)
    y_fit = exponential_decay(x_fit, *popt)
    
    plt.plot(x_fit, y_fit, 'r-', label=f'Fit: {popt[0]:.1f}e^(-{popt[1]:.2f}w) + {popt[2]:.2f}')
    plt.axhline(popt[2], color='green', linestyle='--', label=f'Plateau (g={popt[2]:.2f})')
    print(f"\nConvergence identified at g = {popt[2]:.2f}")
except Exception as e:
    print(f"Curve fit failed: {e}")

plt.title(f"Sensitivity Analysis: {mode_str.capitalize()} Map - Chr {CHROM}")
plt.xlabel("Window Size (kb)")
plt.ylabel("Estimated Generations (g)")
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.5)

out_path = f'output/sensitivity_curve_{CHROM}_{mode_str}_fitted.png'
plt.savefig(out_path)
print(f"Saved to: {out_path}")
