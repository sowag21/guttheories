import sympy as sp
"""
# -----------------------------
# Define symbols
# -----------------------------
lambda1, lambda2, lam, A = sp.symbols('lambda1 lambda2 lam A')
N, x = sp.symbols('N x')

# -----------------------------
# Define equations
# -----------------------------
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

# -----------------------------
# Solver function
# -----------------------------
def solve_full_system(N_val, x_val, guesses):
    
    #Solve the full system for given N, x using multiple initial guesses.
    #Returns a list of distinct solutions.
    
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

            # Convert to float
            sol_floats = tuple(float(s) for s in sol)

            # Check if solution is new (avoid duplicates)
            if all([abs(sol_floats[i] - s[i]) > 1e-6 for s in solutions for i in range(4)]):
                solutions.append(sol_floats)

        except:
            continue

    return solutions

# -----------------------------
# Example usage
# -----------------------------
N_val = 10.0   # Example N
x_val = 2.0    # Example x

# Multiple initial guesses to find all branches
guesses = [
    (0.1, 0.1, 0.1, 1.0),
    (1.0, 1.0, 0.1, 1.0),
    (-0.5, -0.5, 0.0, 1.0),
    (0.5, 0.5, 0.0, 4.0),
    (1.5, 1.5, 0.1, 3.0),
]

solutions = solve_full_system(N_val, x_val, guesses)

print(f"Found {len(solutions)} solution(s) for N={N_val}, x={x_val}:\n")
for sol in solutions:
    lam1, lam2, lam_val, A_val = sol
    # Compute constraint to verify
    con_val = float(constraint.subs({
        lambda1: lam1,
        lambda2: lam2,
        lam: lam_val,
        A: A_val,
        N: N_val,
        x: x_val
    }))
    print(f"λ1={lam1:.6f}, λ2={lam2:.6f}, λ={lam_val:.6e}, A={A_val:.6f}, constraint={con_val:.2e}")
"""
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

# -----------------------------
# Symbols and equations
# -----------------------------
lambda1, lambda2, lam, A_sym = sp.symbols('lambda1 lambda2 lam A')
N_sym, x_sym = sp.symbols('N x')

eq1 = (
    -lambda1
    - ((2 + 8/(x_sym*N_sym))*lambda1**2 + 2*lam**2
       - 6*(1 - 1/(x_sym**2 * N_sym**2))*A_sym*lambda1
       + 3/sp.Integer(2) * A_sym**2 * (1 + 1/(x_sym*N_sym) - 4/(x_sym**2*N_sym**2) + 2/(x_sym**3*N_sym**3)))
)

eq2 = (
    -lambda2
    - ((2 + 8/N_sym)*lambda2**2 + 2*lam**2
       - 6*(1 - 1/(N_sym**2))*A_sym*lambda2
       + 3/sp.Integer(2) * A_sym**2 * (1 + 1/N_sym - 4/N_sym**2 + 2/N_sym**3))
)

eq3 = (
    -lam
    - lam * (
        -4/(N_sym*sp.sqrt(x_sym))*lam
        + 2*(1 + 1/(x_sym*N_sym))*lambda1
        + 2*(1 + 1/N_sym)*lambda2
        - 3*A_sym*(2 - 1/(x_sym**2*N_sym**2) - 1/N_sym**2)
    )
)

constraint = 2*lambda2*(1 + 1/N_sym) - sp.sqrt(x_sym)*lam + 3*A_sym*(1 - 1/N_sym**2)

# -----------------------------
# Solver function
# -----------------------------
def solve_system_fixed_A(N_val, x_val, A_val):
    guesses = [
        (0.1,0.1,0.1),
        (0.5,0.5,0.1),
        (1.0,1.0,0.0),
        (-0.1,-0.1,0.0)
    ]
    for g in guesses:
        try:
            sol = sp.nsolve(
                [eq1.subs({N_sym: N_val, x_sym: x_val, A_sym: A_val}),
                 eq2.subs({N_sym: N_val, x_sym: x_val, A_sym: A_val}),
                 eq3.subs({N_sym: N_val, x_sym: x_val, A_sym: A_val})],
                [lambda1, lambda2, lam],
                g,
                tol=1e-10,
                maxsteps=50
            )
            lam1_val, lam2_val, lam_val = [float(s) for s in sol]
            # compute effective m^2 from constraint
            m2_val = float(constraint.subs({
                lambda1: lam1_val,
                lambda2: lam2_val,
                lam: lam_val,
                A_sym: A_val,
                N_sym: N_val,
                x_sym: x_val
            }))
            return m2_val
        except:
            continue
    return None  # no solution

# -----------------------------
# Scan grid
# -----------------------------
A_val = 0.1  # fixed value of A
N_values = np.linspace(3, 20, 50)
x_values = np.linspace(1, 5, 50)

M2_grid = np.zeros((len(N_values), len(x_values)))

for i, N_val in enumerate(N_values):
    for j, x_val in enumerate(x_values):
        m2 = solve_system_fixed_A(N_val, x_val, A_val)
        if m2 is not None:
            M2_grid[i,j] = m2
        else:
            M2_grid[i,j] = np.nan  # mark points with no solution

# -----------------------------
# Contour plot
# -----------------------------
plt.figure(figsize=(8,6))
X, Y = np.meshgrid(x_values, N_values)
cp = plt.contourf(X, Y, M2_grid, levels=[-np.inf, -1e-3, 1e-3, np.inf], colors=['red','yellow','green'], alpha=0.7)
plt.xlabel('x')
plt.ylabel('N')
plt.title(f'Regions of m^2 for A={A_val}')
plt.colorbar(cp, ticks=[-1,0,1], label='m^2 sign')
plt.show()
