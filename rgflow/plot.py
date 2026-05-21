import numpy as np
import matplotlib.pyplot as plt

# 1. Define the energy scale (Q^2)
# We start slightly above 1 to avoid a divide-by-zero error in our simplified log equations
Q2 = np.linspace(1.15, 10.5, 500)

# 2. Define the parameterized coupling equations (simplified for conceptual shape)
# QCD: Decreases logarithmically
alpha_qcd = 1.0 / np.log(Q2)

# QED: Increases, eventually diverging (conceptual Landau pole)
alpha_qed = 0.15 / (1 - 0.4 * np.log(1.1*Q2))

# 3. Setup the plot
fig, ax = plt.subplots(figsize=(8, 6))

# Plot the lines matching the image's style
ax.plot(Q2, alpha_qcd, 'k-', linewidth=2)   # Solid black line for QCD
ax.plot(Q2, alpha_qed, 'k:', linewidth=2)   # Dotted black line for QED

# 4. Add Text Annotations
ax.text(2.0, 4.5, r'$\alpha_{QCD}(q^2)$', fontsize=16)
ax.text(7.0, 4.5, r'$\alpha_{QED}(q^2)$', fontsize=16)

# 5. Format the axes (Remove numbers to make it a conceptual diagram)
ax.set_xticks([])
ax.set_yticks([])

# Add custom axis labels
ax.set_ylabel(r'$\alpha(q^2)$', fontsize=16, weight='bold')

ax.set_xlabel('$q^2$ ', fontsize=14)

# Add the top x-axis label (probing small distance scales)
ax_top = ax.twiny()
ax_top.set_xticks([])
ax_top.set_xlabel(r'large momentum transfer $\rightarrow$', fontsize=14, loc='right')

# Set axis limits to keep the view clean
ax.set_ylim(0, 5.5)
ax.set_xlim(0.8, 10.8)

# Emphasize the box border (spines) to match the image
for spine in ax.spines.values():
    spine.set_linewidth(1.5)

plt.tight_layout()
plt.savefig("qedqcdplot.pdf", format='pdf', bbox_inches='tight')
plt.show()