import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

# Define symbols and equations (same as before)
lambda1, lambda2, lam, A = sp.symbols('lambda1 lambda2 lam A')
N, x = sp.symbols('N x')

eq1 = (
    -lambda1
    - ((2 + 8/(x*N))*lambda1**2 + 2*lam**2
       - 6*(1 - 1/(x**2 * N**2))*A*lambda1
       + 3/sp.Integer(2) * A**2 * (1 + 1/(x*N) - 4/(x**2*N**2) + 2/(x**3*N**3)))
)

eq2 = (
    -lambda2
    - ((2 + 8/N)*lambda2**2 + 2*lam**2
       - 6*(1 - 1/(N**2))*A*lambda2
       + 3/sp.Integer(2) * A**2 * (1 + 1/N - 4/N**2 + 2/N**3))
)

eq3 = (
    -lam
    - lam * (
        -4/(N*sp.sqrt(x))*lam
        + 2*(1 + 1/(x*N))*lambda1
        + 2*(1 + 1/N)*lambda2
        - 3*A*(2 - 1/(x**2*N**2) - 1/N**2)
    )
)

constraint = 2*lambda2*(1 + 1/N) - sp.sqrt(x)*lam + 3*A*(1 - 1/N**2)

# Solver function (same as before)
def solve_full_system(N_val, x_val, guesses):
    solutions = []
    for g in guesses:
        try:
            sol = sp.nsolve(
                [eq1.subs({N: N_val, x: x_val}),
                 eq2.subs({N: N_val, x: x_val}),
                 eq3.subs({N: N_val, x: x_val}),
                 constraint.subs({N: N_val, x: x_val})],
                [lambda1, lambda2, lam, A],
                g,
                tol=1e-14,
                maxsteps=50
            )
            sol_floats = tuple(float(s) for s in sol)
            # avoid duplicates
            if all([abs(sol_floats[i]-s[i])>1e-6 for s in solutions for i in range(4)]):
                solutions.append(sol_floats)
        except:
            continue
    return solutions
# Scan ranges
N_values = np.linspace(3, 20, 50)  # adjust as needed
x_values = np.linspace(1, 5, 50)

# Store points where the system has a valid solution
plot_points = []

# Use reasonable initial guesses
guesses = [
    (0.1,0.1,0.1,1.0),
    (0.5,0.5,0.1,2.0),
    (1.0,1.0,0.0,3.0)
]

for N_val in N_values:
    for x_val in x_values:
        sols = solve_full_system(N_val, x_val, guesses)
        for sol in sols:
            lam1, lam2, lam_val, A_val = sol
            # record N, x, and A (or whatever quantity you want)
            plot_points.append((x_val, N_val, A_val))
plot_points = np.array(plot_points)
if plot_points.size > 0:
    x_plot = plot_points[:,0]
    N_plot = plot_points[:,1]
    A_plot = plot_points[:,2]  # optional, for color

    plt.figure(figsize=(6,5))
    plt.scatter(x_plot, N_plot, c=A_plot, cmap='viridis', s=30)
    plt.colorbar(label='A (from constraint)')
    plt.xlabel('x')
    plt.ylabel('N')
    plt.title('Valid points where m^2 ~ 0')
    plt.show()
else:
    print("No valid points found.")
