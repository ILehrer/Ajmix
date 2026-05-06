import numpy as np
import pandas as pd
import filter
import matplotlib.pyplot as plt
import sys

CHROM = int(sys.argv[1]) if len(sys.argv) > 1 else 18
USE_GLOBAL = len(sys.argv) > 2 and sys.argv[2] == 'global'
WINDOW_INPUT = sys.argv[3] if len(sys.argv) > 3 and sys.argv[3].isdigit() else 0
SMOOTH = True

def load_recombo_map(chrom):
    """Integrates the 10kb rmap bins into a cumulative genetic map."""
    rmap = pd.read_csv('sex-averaged_noncarrier.rmap.txt', sep='\t')
    rmap = rmap[rmap['chr'] == f'chr{chrom}'].copy()
    # stdrate is cM/Mb. Each bin is 0.01 Mb. 
    rmap['cM_delta'] = rmap['stdrate'] * 0.01 
    rmap['cM_cum'] = rmap['cM_delta'].cumsum()
    return rmap

def get_genetic_dist(pos_start, pos_end, rmap):
    """Interpolates cumulative cM to find Morgans between two points."""
    cM_start = np.interp(pos_start, rmap['pos'], rmap['cM_cum'])
    cM_end = np.interp(pos_end, rmap['pos'], rmap['cM_cum'])
    return (cM_end - cM_start) / 100.0

def baum_welch_windowed(df, rmap, n_iter=10):
    n_obs = len(df)
    n_states = 2 # 0: EU, 1: ME
    eps = 1e-10
    pi = np.array([0.5, 0.5])
    g = 30.0   # Initial guess
    
    pos_array = df['pos'].values
    dist_morgans = np.zeros(n_obs)
    
    # Precalc distances
    for i in range(1, n_obs):
        if USE_GLOBAL:
            dist_morgans[i] = (pos_array[i] - pos_array[i-1]) * 1e-8
        else:
            dist_morgans[i] = get_genetic_dist(pos_array[i-1], pos_array[i], rmap)
    dist_morgans[0] = 0.0005 # Baseline jump

    E_cols = ['E0_smooth', 'E1_smooth'] if SMOOTH else ['E0', 'E1']
    E = df[E_cols].values + eps

    for iteration in range(n_iter):
        # Forward Pass
        alpha = np.zeros((n_obs, n_states))
        c = np.zeros(n_obs)
        alpha[0] = pi * E[0]
        c[0] = 1.0 / (np.sum(alpha[0]) + eps)
        alpha[0] *= c[0]
        
        for t in range(1, n_obs):
            p_switch = 0.5 * (1 - np.exp(-2 * g * dist_morgans[t]))
            A = np.array([[1-p_switch, p_switch], [p_switch, 1-p_switch]])
            alpha[t] = (alpha[t-1] @ A) * E[t]
            c[t] = 1.0 / (np.sum(alpha[t]) + eps)
            alpha[t] *= c[t]
            
        # Backward Pass
        beta = np.zeros((n_obs, n_states))
        beta[-1] = np.ones(n_states) * c[-1]
        for t in range(n_obs - 2, -1, -1):
            p_switch = 0.5 * (1 - np.exp(-2 * g * dist_morgans[t+1]))
            A = np.array([[1-p_switch, p_switch], [p_switch, 1-p_switch]])
            beta[t] = (A @ (E[t+1] * beta[t+1])) * c[t]
            
        # Expectation
        gamma = alpha * beta
        gamma /= np.sum(gamma, axis=1, keepdims=True)
        
        total_switches = 0
        sum_morgans = 0
        for t in range(n_obs - 1):
            p_switch = 0.5 * (1 - np.exp(-2 * g * dist_morgans[t+1]))
            A = np.array([[1-p_switch, p_switch], [p_switch, 1-p_switch]])
            xi_t = (alpha[t][:, None] * A * E[t+1] * beta[t+1])
            xi_t /= (np.sum(xi_t) + eps)
            total_switches += (xi_t[0, 1] + xi_t[1, 0])
            sum_morgans += dist_morgans[t+1]
            
        g = total_switches / (2 * sum_morgans + eps)
        print(f"Iteration {iteration+1} - g: {g:.2f}")

    return gamma, g

# Data Prep
print("Loading data...")
raw_df = filter.get_df(CHROM)
rmap_data = None if USE_GLOBAL else load_recombo_map(CHROM)

WINDOW_SIZE = 50_000 if not WINDOW_INPUT else int(WINDOW_INPUT) 
raw_df['window'] = raw_df['pos'] // WINDOW_SIZE
eps = 1e-10

# Likelihood math
af_eu, af_me = raw_df['af_eu'].values, raw_df['af_me'].values
obs = raw_df['dosage'].values
p_eu_all = np.vstack([(1-af_eu)**2, 2*af_eu*(1-af_eu), af_eu**2])
p_me_all = np.vstack([(1-af_me)**2, 2*af_me*(1-af_me), af_me**2])

raw_df['log_p_eu'] = np.log(p_eu_all[obs, np.arange(len(raw_df))] + eps)
raw_df['log_p_me'] = np.log(p_me_all[obs, np.arange(len(raw_df))] + eps)

windowed_df = raw_df.groupby('window').agg({'pos': 'mean', 'log_p_eu': 'sum', 'log_p_me': 'sum'}).reset_index()
max_log = windowed_df[['log_p_eu', 'log_p_me']].max(axis=1)
windowed_df['E0'] = np.exp(windowed_df['log_p_eu'] - max_log)
windowed_df['E1'] = np.exp(windowed_df['log_p_me'] - max_log)

# Optimal 3-window smoothing
windowed_df['E0_smooth'] = windowed_df['E0'].rolling(window=3, center=True).mean().fillna(windowed_df['E0'])
windowed_df['E1_smooth'] = windowed_df['E1'].rolling(window=3, center=True).mean().fillna(windowed_df['E1'])

mode_str = 'Global' if USE_GLOBAL else 'Map-based'
print(f"Running HMM ({mode_str})...")
gamma, final_g = baum_welch_windowed(windowed_df, rmap_data)
print(f'\nRESULT [{WINDOW_SIZE // 1000}kb]: {round(final_g, 2)} generations')
