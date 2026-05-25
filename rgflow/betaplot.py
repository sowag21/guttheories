import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Constants

# Initial values at MZ = 91.1876 GeV
mu0 = 91.1876
# We use the values you provided for g1, g2, g3
g1_0 =  0.461
g2_0 = 0.652
g3_0 = 1.22

# YOUR DEFINITION: alpha = g^2 / (4*pi)^2
alpha1_0 = (g1_0**2) / (4 * np.pi)**2
alpha2_0 = (g2_0**2) / (4 * np.pi)**2
alpha3_0 = (g3_0**2) / (4 * np.pi)**2

# Beta function coefficients (Standard Model 1-loop)
# Note: Since you removed the 1/(2*pi) from the inv_alpha function,
# we scale the b0 values by 1/(2*pi) to maintain physical consistency.
b1 = 2*(41/10) 
b2 = 2*(-19/6)  
b3 = 2*(-7) 

# YOUR DEFINITION: alpha^-1(mu) = alpha^-1(mu0) - b0 * ln(mu/mu0)
def inv_alpha(mu_ratio, b0, alpha0):
    return 1.0/alpha0 - b0 * np.log(mu_ratio)

# Define energy range mu/mu0 from 10^0 up to 10^17 GeV
mu_range = np.logspace(0, 16, 500)
energy_vals = mu_range * mu0

# Calculate running for each force
inv_alpha1 = inv_alpha(mu_range, b1, alpha1_0)
inv_alpha2 = inv_alpha(mu_range, b2, alpha2_0)
inv_alpha3 = inv_alpha(mu_range, b3, alpha3_0)

# Calculate Differences
# At unification, these should all equal zero
diff12 = inv_alpha1 - inv_alpha2
diff13 = inv_alpha1 - inv_alpha3
diff23 = inv_alpha2 - inv_alpha3


# We find where it crosses zero by interpolating the log of energy
log_e = np.log10(energy_vals)
gut_min = 10**np.interp(0, diff12[::-1], log_e[::-1]) # Red line crossing

# 2. Find the zero-crossing for the last line (diff23)
gut_max = 10**np.interp(0, diff23[::-1], log_e[::-1]) # Green line crossing

# Set up the plot
fig, ax = plt.subplots(figsize=(10, 8))

# Colors
colors = {'g1': "#d54143", 'g2': '#1f78b4', 'g3': '#33a02c'}

# Plot the lines
ax.plot(energy_vals, inv_alpha1, color=colors['g1'], linewidth=3, label=r'$\alpha_1^{-1}$ ($U(1)_Y$)')
ax.plot(energy_vals, inv_alpha2, color=colors['g2'], linewidth=3, label=r'$\alpha_2^{-1}$ ($SU(2)_L$)')
ax.plot(energy_vals, inv_alpha3, color=colors['g3'], linewidth=3, label=r'$\alpha_3^{-1}$ ($SU(3)_C$)')


# 2. Add the shaded band
ax.axvspan(gut_min, gut_max, color='gray', alpha=0.15, label='Potential GUT Scale')

# 3. Optional: Add a vertical text label for the band
ax.text(np.sqrt(gut_min * gut_max), 150, 'GUT SCALE', 
        fontsize=20, ha='center', fontweight='bold', color='gray', alpha=0.6)

# Formatting axes
ax.set_xscale('log')
ax.set_xlabel(r'Energy scale $\mu$ [GeV]', fontsize=32)
ax.set_ylabel(r'$\alpha_i^{-1}(\mu)$', fontsize=32)

# Set limits to see the "miss" near 10^15 GeV
ax.set_xlim(mu0, 1e18)
ax.set_ylim(0,1000)

# Tick formatting
ax.tick_params(axis='both', which='major', labelsize=26, width=1.5, length=8, direction='in', top=True, right=True)
ax.tick_params(axis='x', which='minor', width=1, length=4, direction='in', top=True)

# Legend
ax.legend(loc='upper right', fontsize=20, frameon=True, edgecolor='black', fancybox=False)

plt.tight_layout()
plt.savefig("sm_unification.pdf", format='pdf', bbox_inches='tight')
plt.show()


fig, ax = plt.subplots(figsize=(10, 8))

# Define colors consistent with your previous plots
color_12 = "#d54143" # Red for U(1) - SU(2)
color_23 = "#33a02c" # Green for SU(2) - SU(3)
color_13 = "#1f78b4" # Blue for U(1) - SU(3)

# Plotting the differences
ax.plot(energy_vals, diff12, color=color_12, linewidth=3, 
        label=r'$\alpha_1^{-1} - \alpha_2^{-1}$ ($U(1)_Y - SU(2)_L$)')
ax.plot(energy_vals, diff23, color=color_23, linewidth=3, 
        label=r'$\alpha_2^{-1} - \alpha_3^{-1}$ ($SU(2)_L - SU(3)_C$)')
ax.plot(energy_vals, diff13, color=color_13, linewidth=3, 
        label=r'$\alpha_1^{-1} - \alpha_3^{-1}$ ($U(1)_Y - SU(3)_C$)')


# Add a dashed line at 0 to represent perfect unification
ax.axhline(0, color='black', linewidth=2.5, linestyle='dotted', alpha=0.7)

# 2. Add the shaded band
ax.axvspan(gut_min, gut_max, color='gray', alpha=0.15, label='Potential GUT Scale')

# 3. Optional: Add a vertical text label for the band
ax.text(np.sqrt(gut_min * gut_max), 150, 'GUT SCALE', 
        fontsize=20, ha='center', fontweight='bold', color='gray', alpha=0.6)

# Formatting axes
ax.set_xscale('log')
ax.set_xlabel(r'Energy scale $\mu$ [GeV]', fontsize=32)
ax.set_ylabel(r'$\Delta \alpha_i^{-1}(\mu)$', fontsize=32)
ax.set_xlim(1e11, 1e18)
ax.set_ylim(-200, 200) # Centered around zero

# Hide the bottom-most y-tick label (-200) to prevent overlap
yticks = ax.get_yticks()
ax.set_yticks(yticks) # Lock in the tick positions
ylabels = [str(int(y)) for y in yticks] # Convert to strings
ylabels[0] = '' # Clear the first label
ax.set_yticklabels(ylabels)

# Tick formatting
ax.tick_params(axis='both', which='major', labelsize=26, width=1.5, length=8, 
               direction='in', top=True, right=True)
ax.tick_params(axis='x', which='minor', width=1, length=4, direction='in', top=True)

# Legend
ax.legend(loc='lower left', fontsize=18, frameon=True, edgecolor='black', fancybox=False)

plt.tight_layout()
plt.savefig("sm_unification_differences.pdf", format='pdf', bbox_inches='tight')
plt.show()










# Define the running coupling function derived from the user's equation
def alpha_g(mu_ratio, b0, alpha_0=0.1):
    return 1.0 / (1.0 / alpha_0 - b0 * np.log(mu_ratio))

# Define range for mu/mu_0
mu_ratio = np.logspace(0, 4, 5000)

# Evaluate for b0 = -1 and b0 = 1
alpha_b_neg1 = alpha_g(mu_ratio, b0=-1)
alpha_b_pos1 = alpha_g(mu_ratio, b0=1)

# Set up the plot
fig, ax = plt.subplots(figsize=(10, 8))

# Colors
color_af = '#1f78b4'  # Blue for b0 = -1
color_lp = '#33a02c'  # Green for b0 = 1

# Plot the lines
ax.plot(mu_ratio, alpha_b_neg1, color=color_af, linewidth=3)
ax.plot(mu_ratio, alpha_b_pos1, color=color_lp, linewidth=3, linestyle='--')

# Formatting axes
ax.set_xscale('log') 
ax.set_xlabel(r'$\mu / \mu_0$', fontsize=32)

# THE FIX: Added labelpad=25 to push the y-label to the left. 
# This adds invisible width to the bounding box to match the minus signs on the other plot!
ax.set_ylabel(r'$\alpha_g(\mu)$', fontsize=32, labelpad=25)

ax.set_xlim(1, 1.1*1e4)
ax.set_ylim(0, np.max(alpha_b_pos1) * 1.01)

# Tick formatting
ax.tick_params(axis='both', which='major', labelsize=26, width=1.5, length=8, direction='in', top=True, right=True)
ax.tick_params(axis='x', which='minor', width=1, length=4, direction='in', top=True)

ticks = [1, 10, 100, 1000, 10000]
labels = ['', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$'] 
ax.set_xticks(ticks)
ax.set_xticklabels(labels)

# Legend formatting (updated sizes to perfectly match the Beta function plot)
legend_elements = [
    Line2D([0], [0], color=color_af, linewidth=3, label=r'$b_0 = -1$'),
    Line2D([0], [0], color=color_lp, linewidth=3, label=r'$b_0 = 1$', linestyle='--')
]
ax.legend(handles=legend_elements, loc='upper left', fontsize=22, frameon=True, 
          edgecolor='black', fancybox=False, title='Running Coupling', title_fontsize=24)

plt.tight_layout()
plt.savefig("running_coupling_pretty.pdf", format='pdf', bbox_inches='tight')
plt.savefig("running_coupling_pretty.png", format='png', bbox_inches='tight', dpi=150)
plt.show()

fig, ax = plt.subplots(figsize=(10, 8))

# Colors
color_sm = '#1f78b4'  # Light blue for SM
color_bsm = '#33a02c'  # Red for nf=20

# Plot the lines
ax.plot(g_vals, beta_6, color=color_sm, linewidth=3)
ax.plot(g_vals, beta_27, color=color_bsm, linewidth=3, linestyle='--')

# Plot the beta = 0 boundary line
ax.axhline(0, color='black', linewidth=2.5, linestyle='dotted')

# Formatting axes
ax.set_xlabel('$g$', fontsize=32)
ax.set_ylabel(r'$\beta_{QCD}(g)$', fontsize=32)
ax.set_xlim(0, 3)

# Adjust y limits slightly for breathing room
ax.set_ylim(np.min(beta_6)*1.1, np.max(beta_27)*1.3)

ax.tick_params(axis='both', which='major', labelsize=26, width=1.5, length=8, direction='in', top=True, right=True)

# Legend formatting
legend_elements = [
    Line2D([0], [0], color=color_sm, linewidth=3, label=r'$n_f = 6$'),
    Line2D([0], [0], color=color_bsm, linewidth=3, label=r'$n_f = 27$'),
    Line2D([0], [0], color='black', linewidth=2.5, linestyle='dashed', label=r'$\beta=0$')
]
ax.legend(handles=legend_elements, loc='upper left', fontsize=22, frameon=True, 
          edgecolor='black', fancybox=False, title='SU(3) Beta function', title_fontsize=24)

plt.tight_layout()
plt.savefig("beta_function_pretty.pdf", format='pdf', bbox_inches='tight')
plt.savefig("beta_function_pretty.png", format='png', bbox_inches='tight', dpi=150) # Save png for preview if needed

# --- REPLACEMENT FOR THE FINAL PLOTTING BLOCK ---

# Define the beta function for alpha instead of g
# beta(alpha) = - (alpha^2 / 2*pi) * (11/3 N - 2/3 n_f)
def beta_alpha(alpha, N, n_f):
    b0 = (11/3 * N - 2/3 * n_f)
    return - (alpha**2 / (2 * np.pi)) * b0

# Set a sensible range for alpha (g=3 corresponds to alpha ~= 0.7)
alpha_vals = np.linspace(0, 0.8, 200)
N = 3 # SU(3)

beta_alpha_6 = beta_alpha(alpha_vals, N, 6)
beta_alpha_27 = beta_alpha(alpha_vals, N, 27)

# Set up the plot
fig, ax = plt.subplots(figsize=(10, 8))

# Colors
color_sm = '#1f78b4'  # Light blue for SM
color_bsm = '#33a02c' # Green for nf=27

# Plot the lines
ax.plot(alpha_vals, beta_alpha_6, color=color_sm, linewidth=3)
ax.plot(alpha_vals, beta_alpha_27, color=color_bsm, linewidth=3, linestyle='--')

# Plot the beta = 0 boundary line
ax.axhline(0, color='black', linewidth=2.5, linestyle='dotted')

# Formatting axes
ax.set_xlabel(r'$\alpha$', fontsize=32)
ax.set_ylabel(r'$\beta_{QCD}(\alpha)$', fontsize=32)
ax.set_xlim(0, 0.8)

# Adjust y limits slightly for breathing room
ax.set_ylim(np.min(beta_alpha_6)*1.1, np.max(beta_alpha_27)*1.3)

ax.tick_params(axis='both', which='major', labelsize=26, width=1.5, length=8, direction='in', top=True, right=True)

# Legend formatting
legend_elements = [
    Line2D([0], [0], color=color_sm, linewidth=3, label=r'$n_f = 6$'),
    Line2D([0], [0], color=color_bsm, linewidth=3, label=r'$n_f = 27$', linestyle='--'),
    Line2D([0], [0], color='black', linewidth=2.5, linestyle='dotted', label=r'$\beta=0$')
]
ax.legend(handles=legend_elements, loc='upper left', fontsize=22, frameon=True, 
          edgecolor='black', fancybox=False, title='SU(3) Beta function', title_fontsize=24)

plt.tight_layout()
plt.savefig("beta_function_alpha_pretty.pdf", format='pdf', bbox_inches='tight')
plt.savefig("beta_function_alpha_pretty.png", format='png', bbox_inches='tight', dpi=150)
plt.show()