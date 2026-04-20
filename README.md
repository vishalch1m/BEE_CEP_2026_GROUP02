# BEE_CEP_2026_GROUP02
# DC Network Analyzer Using Nodal Analysis Method

**Course:** A9205 – Basic Electrical Engineering Laboratory (VCE-R25)  
**Class:** I B.Tech. II Semester | CSE – F  
**College:** Vardhaman College of Engineering  
**Academic Year:** 2025–2026  

---

## Group Members

| S.No | Name | Roll No |
|------|------|---------|
| 1    | Ch. Vishal | 25881A05BU |
| 2    | M. Vivekananda Reddy | 25881A05BV |
| 3    | E. Rivan Reddy | 25881A05AW |

---

## Problem Description

This project implements a **DC Network Analyzer** that computes unknown node voltages and branch currents in any DC resistive network using the **Nodal Analysis Method**.

Given:
- A DC circuit with resistors and independent current sources
- Number of nodes and their connectivity

The program finds:
- Voltage at each node (with respect to ground)
- Current through each resistor branch
- Power dissipated in each branch
- Graphical plots of all results

---

## Mathematical Formulation

Nodal Analysis is based on **Kirchhoff's Current Law (KCL)**:

> Sum of all currents leaving a node = 0

### System of Equations

For `n` nodes, the matrix equation is:

```
[G] × [V] = [I]
```

Where:

**Conductance Matrix [G] (n × n):**
- Diagonal:    `G[i][i] = Σ(1/R)` — sum of all conductances at node i
- Off-diagonal: `G[i][j] = −(1/R_ij)` — negative conductance between node i and j

**Voltage Vector [V]:**
- `V = [V1, V2, ..., Vn]` — unknown node voltages

**Current Vector [I]:**
- `I[i]` = algebraic sum of all current source injections at node i

### Solution Method

The system `[G][V] = [I]` is solved using **Gaussian Elimination with Partial Pivoting**:

```
Step 1: Build [G] from resistor branches
Step 2: Build [I] from current sources
Step 3: Solve [G][V]=[I] using Gaussian Elimination
Step 4: Compute branch currents: I_k = (V_i - V_j) / R_k
Step 5: Compute power: P_k = (V_i - V_j)² / R_k
```

---

## Input & Output Format

### Input (Console)
```
Number of nodes        : integer (excluding ground)
Number of resistors    : integer
For each resistor      : node_i, node_j, resistance (Ω)
Number of curr sources : integer
For each source        : node_in, node_out, current (A)
```

### Output (Console + Plots)
```
- Conductance matrix [G]
- Current vector [I]
- Node voltages (V)
- Branch currents (A)
- Branch power dissipation (W)
- Total power dissipated (W)
```

### Plots (saved to `screenshots/`)
- `node_voltages.png`    – Bar chart of node voltages
- `branch_currents.png`  – Horizontal bar chart of branch currents
- `power_dissipation.png`– Bar chart of power dissipated per resistor

---

## How to Run the Program

### Requirements
- Python 3.x
- matplotlib

### Install Dependency
```bash
pip install matplotlib
```

### Run
```bash
cd src
python dc_nodal_analyzer.py
```

---

## Sample Output

### Validation Test (Built-In)
```
Node  | Computed (V) | Expected (V) | Error (%)
V1    |    13.3333   |    13.333    |  0.0000% ✓ PASS
V2    |     6.6667   |     6.667    |  0.0000% ✓ PASS
V3    |     6.6667   |     6.667    |  0.0000% ✓ PASS
```

### User-Defined Circuit Example (2 nodes, 3 resistors, 1 source)
```
Enter number of nodes: 2
Resistor R1: Node 1 → Node 2, 4 Ω
Resistor R2: Node 2 → GND,    2 Ω
Resistor R3: Node 1 → GND,    8 Ω
Current Source: 3 A into Node 1

Result:
  V1 = 8.7273 V
  V2 = 5.8182 V
  Branch currents and power also displayed.
```

---

## Repository Structure

```
BEE_CEP_2026_GroupXX/
├── src/
│   └── dc_nodal_analyzer.py     ← Main program
├── report/
│   └── BEE_CEP_Report.pdf       ← Final report
├── screenshots/
│   ├── node_voltages.png
│   ├── branch_currents.png
│   └── power_dissipation.png
├── README.md
└── requirements.txt
```

---

## References

1. V. K. Mehta & Rohit Mehta, *Principles of Electrical Engineering*, S. Chand, 2014.
2. D. P. Kothari & I. J. Nagrath, *Basic Electrical Engineering*, Tata McGraw-Hill, 2010.
