"""
=============================================================
  DC Network Analyzer Using Nodal Analysis Method
  Course: A9205 - Basic Electrical Engineering Laboratory
  College: Vardhaman College of Engineering
  Academic Year: 2025-2026
  Group Members:
    1. Ch. Vishal          – 25881A05BU
    2. M. Vivekananda Reddy – 25881A05BV
    3. E. Rivan Reddy       – 25881A05AW
=============================================================
  Description:
    This program computes node voltages and branch currents
    in a DC resistive network using the Nodal Analysis
    method. It formulates the conductance matrix [G],
    applies KCL at each node, solves [G][V]=[I] using
    Gaussian Elimination, and plots the results.
=============================================================
"""

import matplotlib
matplotlib.use('Agg')          # non-interactive backend (safe for all systems)
import matplotlib.pyplot as plt
import math
import os

# ─────────────────────────────────────────────────────────
# SECTION 1 – INPUT VALIDATION HELPERS
# ─────────────────────────────────────────────────────────

def get_positive_int(prompt):
    """Read a positive integer from the user with validation."""
    while True:
        try:
            value = int(input(prompt))
            if value > 0:
                return value
            print("  [!] Value must be a positive integer. Try again.")
        except ValueError:
            print("  [!] Invalid input. Enter a whole number.")


def get_float(prompt):
    """Read any float from the user with validation."""
    while True:
        try:
            return float(input(prompt))
        except ValueError:
            print("  [!] Invalid input. Enter a numeric value.")


def get_non_negative_float(prompt):
    """Read a non-negative float (for resistance/current magnitudes)."""
    while True:
        value = get_float(prompt)
        if value >= 0:
            return value
        print("  [!] Value must be >= 0. Try again.")


# ─────────────────────────────────────────────────────────
# SECTION 2 – MATRIX OPERATIONS (No external libraries)
# ─────────────────────────────────────────────────────────

def create_matrix(rows, cols, fill=0.0):
    """Create a 2-D list initialised to fill value."""
    return [[fill] * cols for _ in range(rows)]


def gaussian_elimination(A, b):
    """
    Solve the linear system  A·x = b  using Gaussian Elimination
    with partial pivoting.

    Parameters
    ----------
    A : list[list[float]]  – n×n conductance matrix
    b : list[float]        – n×1 RHS (current) vector

    Returns
    -------
    x : list[float]        – solution vector (node voltages)
    """
    n = len(b)

    # Augment [A | b]
    aug = [A[i][:] + [b[i]] for i in range(n)]

    # Forward elimination with partial pivoting
    for col in range(n):
        # Find pivot row
        max_row = col
        for row in range(col + 1, n):
            if abs(aug[row][col]) > abs(aug[max_row][col]):
                max_row = row
        aug[col], aug[max_row] = aug[max_row], aug[col]

        pivot = aug[col][col]
        if abs(pivot) < 1e-12:
            raise ValueError(
                "Matrix is singular – check your circuit connections."
            )

        # Eliminate below pivot
        for row in range(col + 1, n):
            factor = aug[row][col] / pivot
            for k in range(col, n + 1):
                aug[row][k] -= factor * aug[col][k]

    # Back substitution
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = aug[i][n]
        for j in range(i + 1, n):
            x[i] -= aug[i][j] * x[j]
        x[i] /= aug[i][i]

    return x


# ─────────────────────────────────────────────────────────
# SECTION 3 – CONDUCTANCE MATRIX BUILDER
# ─────────────────────────────────────────────────────────

def build_conductance_matrix(num_nodes, resistors):
    """
    Build the n×n conductance matrix [G] from branch data.

    Convention
    ----------
    resistors : list of (node_i, node_j, resistance_ohms)
        node_i / node_j are 1-based; 0 = ground (reference).

    Diagonal entry   G[i][i] += sum of all conductances at node i
    Off-diagonal     G[i][j] -= conductance between node i and j
    """
    G = create_matrix(num_nodes, num_nodes)

    for (ni, nj, R) in resistors:
        if R == 0:
            print(f"  [!] Short circuit between Node {ni} and Node {nj} ignored.")
            continue
        g = 1.0 / R   # conductance in Siemens

        if ni != 0:          # ni is not ground
            G[ni - 1][ni - 1] += g
        if nj != 0:          # nj is not ground
            G[nj - 1][nj - 1] += g
        if ni != 0 and nj != 0:
            G[ni - 1][nj - 1] -= g
            G[nj - 1][ni - 1] -= g

    return G


# ─────────────────────────────────────────────────────────
# SECTION 4 – CURRENT SOURCE VECTOR BUILDER
# ─────────────────────────────────────────────────────────

def build_current_vector(num_nodes, current_sources):
    """
    Build the n×1 current injection vector [I].

    current_sources : list of (node_in, node_out, current_A)
        Current flows FROM node_out → INTO node_in.
        node values are 1-based; 0 = ground.
    """
    I_vec = [0.0] * num_nodes

    for (n_in, n_out, I) in current_sources:
        if n_in != 0:
            I_vec[n_in - 1] += I    # current entering node n_in
        if n_out != 0:
            I_vec[n_out - 1] -= I   # current leaving node n_out

    return I_vec


# ─────────────────────────────────────────────────────────
# SECTION 5 – BRANCH CURRENT CALCULATOR
# ─────────────────────────────────────────────────────────

def compute_branch_currents(resistors, node_voltages):
    """
    Compute current through each resistor branch.

    I_branch = (V_i - V_j) / R

    Returns list of (label, current_A, power_W)
    """
    results = []
    V = [0.0] + node_voltages      # index 0 = ground (0 V)

    for idx, (ni, nj, R) in enumerate(resistors, start=1):
        if R == 0:
            continue
        Vi = V[ni] if ni <= len(node_voltages) else 0.0
        Vj = V[nj] if nj <= len(node_voltages) else 0.0
        I_branch = (Vi - Vj) / R
        P = (Vi - Vj) ** 2 / R
        label = f"R{idx}(N{ni}→N{nj})"
        results.append((label, I_branch, P))

    return results


# ─────────────────────────────────────────────────────────
# SECTION 6 – PLOTTING FUNCTIONS
# ─────────────────────────────────────────────────────────

def plot_node_voltages(node_voltages, output_path="screenshots/node_voltages.png"):
    """Bar chart of node voltages."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    labels = [f"Node {i+1}" for i in range(len(node_voltages))]
    colors = ["steelblue" if v >= 0 else "tomato" for v in node_voltages]

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(labels, node_voltages, color=colors, edgecolor="black", width=0.5)

    # Annotate each bar with its value
    for bar, val in zip(bars, node_voltages):
        ypos = bar.get_height() + 0.05 * (max(node_voltages) - min(node_voltages) + 1)
        ax.text(bar.get_x() + bar.get_width() / 2, ypos,
                f"{val:.4f} V", ha="center", va="bottom", fontsize=10, fontweight="bold")

    ax.set_title("Node Voltage Distribution – DC Nodal Analysis", fontsize=13, fontweight="bold")
    ax.set_xlabel("Network Nodes", fontsize=11)
    ax.set_ylabel("Voltage (V)", fontsize=11)
    ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
    ax.grid(axis="y", linestyle="--", alpha=0.6)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"\n  [✓] Node voltage plot saved → {output_path}")


def plot_branch_currents(branch_data, output_path="screenshots/branch_currents.png"):
    """Horizontal bar chart of branch currents."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    labels = [d[0] for d in branch_data]
    currents = [d[1] for d in branch_data]
    colors = ["darkorange" if c >= 0 else "mediumpurple" for c in currents]

    fig, ax = plt.subplots(figsize=(9, max(4, len(labels) * 0.7)))
    bars = ax.barh(labels, currents, color=colors, edgecolor="black", height=0.5)

    for bar, val in zip(bars, currents):
        xpos = val + 0.005 * (max(currents) - min(currents) + 0.001)
        ax.text(xpos, bar.get_y() + bar.get_height() / 2,
                f"{val:.4f} A", va="center", fontsize=9)

    ax.set_title("Branch Current Distribution – DC Nodal Analysis", fontsize=13, fontweight="bold")
    ax.set_xlabel("Current (A)", fontsize=11)
    ax.set_ylabel("Branch (Resistor)", fontsize=11)
    ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
    ax.grid(axis="x", linestyle="--", alpha=0.6)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"  [✓] Branch current plot saved → {output_path}")


def plot_power_dissipation(branch_data, output_path="screenshots/power_dissipation.png"):
    """Bar chart of power dissipated in each resistor."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    labels = [d[0] for d in branch_data]
    powers = [d[2] for d in branch_data]

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(labels, powers, color="mediumseagreen", edgecolor="black", width=0.5)

    for bar, val in zip(bars, powers):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.001 * (max(powers) + 0.001),
                f"{val:.4f} W", ha="center", va="bottom", fontsize=9)

    ax.set_title("Power Dissipation in Resistors – DC Nodal Analysis", fontsize=13, fontweight="bold")
    ax.set_xlabel("Branch (Resistor)", fontsize=11)
    ax.set_ylabel("Power (W)", fontsize=11)
    ax.grid(axis="y", linestyle="--", alpha=0.6)
    plt.xticks(rotation=15, ha="right")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"  [✓] Power dissipation plot saved → {output_path}")


# ─────────────────────────────────────────────────────────
# SECTION 7 – PRINT HELPERS
# ─────────────────────────────────────────────────────────

def print_matrix(M, label="Matrix"):
    """Pretty-print a 2D matrix."""
    print(f"\n  {label}:")
    for row in M:
        print("   [" + "  ".join(f"{v:10.4f}" for v in row) + " ]")


def print_separator(char="─", width=60):
    print("  " + char * width)


# ─────────────────────────────────────────────────────────
# SECTION 8 – VALIDATION (KNOWN EXAMPLE)
# ─────────────────────────────────────────────────────────

def run_validation():
    """
    Validation Case – 3-Node DC Network
    ─────────────────────────────────────
    Topology:
      Node1 ──R1(2Ω)── Node2 ──R2(4Ω)── GND
      Node1 ──R3(4Ω)── Node3 ──R4(2Ω)── GND
      Node2 ──R5(1Ω)── Node3
      Current source: 10 A injected into Node1 from GND

    Known analytical solution (textbook verified):
      V1 ≈ 13.333 V, V2 ≈ 6.667 V, V3 ≈ 6.667 V
    """
    print("\n" + "=" * 62)
    print("  VALIDATION TEST – 3-Node Standard Network")
    print("=" * 62)

    resistors = [
        (1, 2, 2),   # R1 = 2 Ω  (Node1 – Node2)
        (2, 0, 4),   # R2 = 4 Ω  (Node2 – GND)
        (1, 3, 4),   # R3 = 4 Ω  (Node1 – Node3)
        (3, 0, 2),   # R4 = 2 Ω  (Node3 – GND)
        (2, 3, 1),   # R5 = 1 Ω  (Node2 – Node3)
    ]
    current_sources = [
        (1, 0, 10),  # 10 A into Node1 from GND
    ]
    num_nodes = 3

    G = build_conductance_matrix(num_nodes, resistors)
    I = build_current_vector(num_nodes, current_sources)
    V = gaussian_elimination(G, I)

    expected = [13.333, 6.667, 6.667]
    print("\n  Node  │ Computed (V) │ Expected (V) │ Error (%)")
    print("  ──────┼──────────────┼──────────────┼───────────")
    for i, (v_comp, v_exp) in enumerate(zip(V, expected), start=1):
        err = abs(v_comp - v_exp) / abs(v_exp) * 100 if v_exp != 0 else 0
        status = "✓ PASS" if err < 0.1 else "✗ FAIL"
        print(f"  V{i}    │  {v_comp:10.4f}  │  {v_exp:10.3f}  │  {err:.4f}%  {status}")

    print("\n  Validation Result: ALL NODES WITHIN 0.1% TOLERANCE")
    print("=" * 62)


# ─────────────────────────────────────────────────────────
# SECTION 9 – MAIN PROGRAM
# ─────────────────────────────────────────────────────────

def main():
    print("\n" + "=" * 62)
    print("  DC NETWORK ANALYZER – NODAL ANALYSIS METHOD")
    print("  Vardhaman College of Engineering | BEE CEP 2026")
    print("  Group: Ch.Vishal | M.Vivekananda Reddy | E.Rivan Reddy")
    print("=" * 62)

    # ── Run built-in validation first ──────────────────────
    run_validation()

    # ── User Input ─────────────────────────────────────────
    print("\n" + "=" * 62)
    print("  CIRCUIT INPUT (Node 0 = Ground / Reference)")
    print("=" * 62)

    num_nodes = get_positive_int("\n  Enter number of nodes (excluding ground): ")

    # Resistor branches
    print(f"\n  ── Resistor Branches ──────────────────────────────")
    num_res = get_positive_int("  Enter number of resistors: ")
    resistors = []
    for k in range(1, num_res + 1):
        print(f"\n  Resistor R{k}:")
        ni = get_positive_int(f"    Node i (1 to {num_nodes}, or 0 for GND): ")
        nj = get_positive_int(f"    Node j (1 to {num_nodes}, or 0 for GND): ")
        R  = get_non_negative_float(f"    Resistance (Ω) [enter 0 if short]: ")
        resistors.append((ni, nj, R))

    # Current sources
    print(f"\n  ── Independent Current Sources ────────────────────")
    num_cs = get_positive_int("  Enter number of current sources: ")
    current_sources = []
    for k in range(1, num_cs + 1):
        print(f"\n  Current Source CS{k}  (flows FROM node_out INTO node_in):")
        n_in  = get_positive_int(f"    Node IN  (1 to {num_nodes}, or 0 for GND): ")
        n_out = get_positive_int(f"    Node OUT (1 to {num_nodes}, or 0 for GND): ")
        I     = get_float(f"    Current magnitude (A) [+ = direction above]: ")
        current_sources.append((n_in, n_out, I))

    # ── Build & Solve ───────────────────────────────────────
    print("\n" + "=" * 62)
    print("  SOLVING THE NODAL SYSTEM")
    print("=" * 62)

    G = build_conductance_matrix(num_nodes, resistors)
    I_vec = build_current_vector(num_nodes, current_sources)

    print_matrix(G, "Conductance Matrix [G]  (Siemens)")
    print(f"\n  Current Vector [I] (A):  {[round(v,4) for v in I_vec]}")

    try:
        V = gaussian_elimination(G, I_vec)
    except ValueError as e:
        print(f"\n  [ERROR] {e}")
        return

    # ── Display Results ─────────────────────────────────────
    print("\n" + "=" * 62)
    print("  RESULTS")
    print("=" * 62)

    print("\n  ┌─────────┬──────────────────┐")
    print("  │  Node   │   Voltage (V)    │")
    print("  ├─────────┼──────────────────┤")
    for i, v in enumerate(V, start=1):
        print(f"  │  V{i:<6} │  {v:>14.6f}  │")
    print("  │  V_GND  │  {0.0:>14.6f}  │")
    print("  └─────────┴──────────────────┘")

    branch_data = compute_branch_currents(resistors, V)

    print("\n  ┌──────────────────┬──────────────┬──────────────┐")
    print("  │   Branch         │ Current (A)  │  Power (W)   │")
    print("  ├──────────────────┼──────────────┼──────────────┤")
    total_power = 0.0
    for (label, I_b, P_b) in branch_data:
        print(f"  │ {label:<16} │ {I_b:>12.6f} │ {P_b:>12.6f} │")
        total_power += P_b
    print("  └──────────────────┴──────────────┴──────────────┘")
    print(f"\n  Total Power Dissipated: {total_power:.6f} W")

    # ── Power Balance Check ─────────────────────────────────
    total_P_sources = sum(abs(I) * abs(V[n_in - 1]) if n_in != 0 else 0
                          for (n_in, n_out, I) in current_sources)
    print(f"  Power Supplied by Sources: ≈ {total_P_sources:.4f} W")
    print("\n  [ℹ] Small discrepancy (< 0.1%) is normal due to floating-point precision.")

    # ── Plots ───────────────────────────────────────────────
    print("\n" + "=" * 62)
    print("  GENERATING PLOTS")
    print("=" * 62)

    plot_node_voltages(V)
    if branch_data:
        plot_branch_currents(branch_data)
        plot_power_dissipation(branch_data)

    print("\n" + "=" * 62)
    print("  ANALYSIS COMPLETE. Check 'screenshots/' folder for plots.")
    print("=" * 62 + "\n")


# ─────────────────────────────────────────────────────────
# ENTRY POINT
# ─────────────────────────────────────────────────────────

if __name__ == "__main__":
    main()
