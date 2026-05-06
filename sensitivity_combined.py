import numpy as np
import pandas as pd
import filter
import matplotlib.pyplot as plt
import sys
from scipy.optimize import curve_fit

# Configuration
CHROMS = [18, 19, 20, 21, 22]
USE_GLOBAL = len(sys.argv) > 1 and sys.argv[1] == 'global'
SMOOTH = True
EPS = 1e-10

# 10->100kb, 1kb inc
WINDOW_SIZES = np.arange(10_000, 101_000, 1000)

def load_recombo_map(chrom):
    rmap = pd.read_csv('sex-averaged_noncarrier.rmap.txt', sep='\t')
    rmap = rmap[rmap['chr'] == f'chr{chrom}'].copy()
    rmap['cM_delta'] = rmap['stdrate'] * 0.01
    rmap['cM_cum'] = rmap['cM_delta'].cumsum()
    return rmap

def get_genetic_dist(pos_start, pos_end, rmap):
    cM_start = np.interp(pos_start, rmap['pos'], rmap['cM_cum'])
    cM_end = np.interp(pos_end, rmap['pos'], rmap['cM_cum'])
    return (cM_end - cM_start) / 100.0

def run_combined_hmm(df, rmap_dict):
    n_obs = len(df)
    n_states = 2
    pi = np.array([0.5, 0.5])
    g = 30.0  # Initial guess

    pos_array = df['pos'].values
    chr_array = df['chr'].values
    dist_morgans = np.zeros(n_obs)

    # precalc distances with Chromosome Break logic
    for i in range(1, n_obs):
        if chr_array[i] != chr_array[i-1]:
            dist_morgans[i] = -1 # Flag for break
        elif USE_GLOBAL:
            dist_morgans[i] = (pos_array[i] - pos_array[i-1]) * 1e-8
        else:
            dist_morgans[i] = get_genetic_dist(pos_array[i-1], pos_array[i], rmap_dict[chr_array[i]])
    
    dist_morgans[0] = -1

    E = df[['E0_smooth', 'E1_smooth']].values + EPS

    for _ in range(10):
        # Forward Pass
        alpha = np.zeros((n_obs, n_states))
        c = np.zeros(n_obs)
        for t in range(n_obs):
            if dist_morgans[t] == -1:
                alpha[t] = pi * E[t]
            else:
                p_switch = 0.5 * (1 - np.exp(-2 * g * dist_morgans[t]))
                A = np.array([[1-p_switch, p_switch], [p_switch, 1-p_switch]])
                alpha[t] = (alpha[t-1] @ A) * E[t]
            c[t] = 1.0 / (np.sum(alpha[t]) + EPS)
            alpha[t] *= c[t]

        # Backward Pass
        beta = np.zeros((n_obs, n_states))
        beta[-1] = np.ones(n_states) * c[-1]
        for t in range(n_obs - 2, -1, -1):
            if dist_morgans[t+1] == -1:
                beta[t] = np.ones(n_states) * c[t]
            else:
                p_switch = 0.5 * (1 - np.exp(-2 * g * dist_morgans[t+1]))
                A = np.array([[1-p_switch, p_switch], [p_switch, 1-p_switch]])
                beta[t] = (A @ (E[t+1] * beta[t+1])) * c[t]

        # Expectation & Maximization
        total_switches = 0
        sum_morgans = 0
        for t in range(n_obs - 1):
            if dist_morgans[t+1] != -1:
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

# --- EXECUTION ---
print(f"Preparing Combined Data for Chromosomes {CHROMS}...")
all_raw = []
rmap_dict = {}

for chrom in CHROMS:
    raw = filter.get_df(chrom)
    if not USE_GLOBAL:
        rmap_dict[chrom] = load_recombo_map(chrom)
    
    # Likelihood pre-calc
    af_eu, af_me = raw['af_eu'].values, raw['af_me'].values
    obs = raw['dosage'].values
    p_eu_all = np.vstack([(1-af_eu)**2, 2*af_eu*(1-af_eu), af_eu**2])
    p_me_all = np.vstack([(1-af_me)**2, 2*af_me*(1-af_me), af_me**2])
    raw['log_p_eu'] = np.log(p_eu_all[obs, np.arange(len(raw))] + EPS)
    raw['log_p_me'] = np.log(p_me_all[obs, np.arange(len(raw))] + EPS)
    raw['chr'] = chrom
    all_raw.append(raw)

full_raw = pd.concat(all_raw)
results_g = []
mode_str = 'global' if USE_GLOBAL else 'local'

print(f"Starting combined sweep (10kb to 100kb)...")
for w_size in WINDOW_SIZES:
    full_raw['window'] = full_raw['pos'] // w_size
    # Group by both window and chromosome to keep them distinct
    win_df = full_raw.groupby(['chr', 'window']).agg({'pos': 'mean', 'log_p_eu': 'sum', 'log_p_me': 'sum'}).reset_index()
    
    max_log = win_df[['log_p_eu', 'log_p_me']].max(axis=1)
    win_df['E0'] = np.exp(win_df['log_p_eu'] - max_log)
    win_df['E1'] = np.exp(win_df['log_p_me'] - max_log)

    # Smooth within chromosome groups
    win_df['E0_smooth'] = win_df.groupby('chr')['E0'].transform(lambda x: x.rolling(3, center=True).mean().fillna(x))
    win_df['E1_smooth'] = win_df.groupby('chr')['E1'].transform(lambda x: x.rolling(3, center=True).mean().fillna(x))

    results_g.append(run_combined_hmm(win_df, rmap_dict))
    if w_size % 10_000 == 0:
        print(f"Window: {w_size/1000}kb | Combined g: {results_g[-1]:.2f}")

# --- FITTING & PLOTTING ---
x_data = WINDOW_SIZES / 1000.0
y_data = np.array(results_g)

plt.figure(figsize=(10, 6))
plt.plot(x_data, y_data, 'bo', label='Combined Data', markersize=5)

try:
    p0 = [y_data[0] - y_data[-1], 0.1, y_data[-1]]
    popt, _ = curve_fit(exponential_decay, x_data, y_data, p0=p0, maxfev=5000)
    x_fit = np.linspace(x_data.min(), x_data.max(), 500)
    y_fit = exponential_decay(x_fit, *popt)
    plt.plot(x_fit, y_fit, 'r-', label=f'Decay Fit (Plateau: {popt[2]:.2f})')
    plt.axhline(popt[2], color='green', linestyle='--', label=f'Asymptotic g={popt[2]:.2f}')
    print(f"\nConvergence identified at g = {popt[2]:.2f}")
    year = 2026 - int(popt[2] * 28)
    print(f"Historical Plateau Date: ~{year} AD")
except Exception as e:
    print(f"Curve fit failed: {e}")

plt.title(f"Combined Sensitivity Analysis (Chr 18-22) - {mode_str.capitalize()} Map")
plt.xlabel("Window Size (kb)")
plt.ylabel("Estimated Generations (g)")
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.savefig('output/combined_sensitivity_curve.png')
print('Saved to: output/combined_sensitivity_curve.png')

