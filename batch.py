import subprocess
import sys

CHROMS = [18, 19, 20, 21, 22]
# Range 45kb to 70kb in 2.5kb steps
WINDOWS = range(45000, 72500, 2500) 

print(f"{'Chr':<6} | {'50kb (Anch)':<12} | {'Range (45-70kb)':<20} | {'Max Delta':<10} | {'Year (50k)'}")
print("-" * 80)

for chrom in CHROMS:
    g_values = []
    g_50 = None
    
    for w in WINDOWS:
        cmd = [sys.executable, "windowed.py", str(chrom), "local", str(w)]
        process = subprocess.run(cmd, capture_output=True, text=True)

        for line in process.stdout.split('\n'):
            if line.startswith("RESULT"):
                parts = line.split()
                # Extract g-value
                g_val = float(parts[2])
                g_values.append(g_val)
                
                # Keep track of our 50kb anchor
                if w == 50000:
                    g_50 = g_val

    if g_values:
        g_min = min(g_values)
        g_max = max(g_values)
        delta = g_max - g_min
        year_50 = 2026 - int(g_50 * 28)

        # Formatting the output
        range_str = f"{g_min:.2f} - {g_max:.2f}"
        print(f"{chrom:<6} | {g_50:<12.2f} | {range_str:<20} | {delta:<10.2f} | ~{year_50} AD")
    else:
        print(f"{chrom:<6} | No data found.")
