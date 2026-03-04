#!/usr/bin/env python
# coding: utf-8

# In[2118]:


import re

# -----------------------------
# Input file
# -----------------------------
input_file = "enhancers_matrix_21_r_7.00000_g_0.40000_s1_0.990_s2_1.010_ma_5.00_mr_15.00_RC_25.00_AC5_t1.txt"

# -----------------------------
# Read file
# -----------------------------
with open(input_file, "r") as file:
    lines = file.readlines()

step_results = {}

i = 0
while i < len(lines):

    # -----------------------------
    # Detect step line
    # -----------------------------
    step_match = re.match(r"Step ([\d.]+):", lines[i])

    if not step_match:
        i += 1
        continue

    step_number = float(step_match.group(1))

    # -----------------------------
    # Read 3x23 matrix
    # -----------------------------
    matrix_lines = lines[i+1:i+4]
    matrix = [[float(x) for x in line.split()] for line in matrix_lines]

    num_cols = len(matrix[0])  

    # -----------------------------
    # extract boundaries
    # -----------------------------
    total_value = 0

    for n in range(1, num_cols-1):   # remove boundaries

        total_value += (
            n * matrix[0][n] +
            (n + 21) * matrix[1][n] +
            (n + 42) * matrix[2][n]
        )

    step_results[step_number] = total_value

    # move to next block
    i += 5


# -----------------------------
# Write output
# -----------------------------
output_file = "GE_distribution_lambda_5.txt"

with open(output_file, "w") as outfile:
    for step, value in step_results.items():
        outfile.write(f"{step} {value}\n")

print("Processing complete.")
print("Output written to:", output_file)


# In[2119]:


import pandas as pd
df=pd.read_csv('GE_distribution_lambda_5.txt',sep=' ',header=None)
df


# In[2120]:


ke = df[df.columns[1]].value_counts(normalize=True).sort_index()
ke


# In[2121]:


ke.to_csv('probability_values.txt',header=None,sep=' ')


# In[2122]:


import numpy as np


# ------------------------------------
# Generate rate arrays
# ------------------------------------
def generate_values(slopes, num_molecules, r, g, ma, mr, AC, RC, f):

    a = np.zeros(num_molecules)
    b = np.zeros(num_molecules)
    c = np.zeros(num_molecules)
    d = np.zeros(num_molecules)
    u = np.zeros(num_molecules)
    v = np.zeros(num_molecules)

    e_to_ea = np.zeros(num_molecules)
    ea_to_e = np.zeros(num_molecules)
    e_to_er = np.zeros(num_molecules)
    er_to_e = np.zeros(num_molecules)

    # Horizontal transitions (exclude boundaries)
    for i in range(1, num_molecules-1):

        a[i] = f * r * (slopes[0] ** (i - 1))
        b[i] = f * g * (slopes[1] ** (i - 1))

        c[i] = f * ma * r * (slopes[0] ** (i - 1))
        d[i] = f * g * (slopes[1] ** (i - 1))

        u[i] = f * r * (slopes[0] ** (i - 1))
        v[i] = f * mr * g * (slopes[1] ** (i - 1))

    # Vertical transitions
    for i in range(1, num_molecules-1):

        e_to_ea[i] = f * AC * r
        ea_to_e[i] = f * g

        e_to_er[i] = f * RC * r
        er_to_e[i] = f * g

    return a, b, c, d, u, v, e_to_ea, ea_to_e, e_to_er, er_to_e


# ------------------------------------
# Build W matrix
# ------------------------------------
def build_W(num_molecules, slopes, r, g, ma, mr, AC, RC, f):

    a, b, c, d, u, v, e_to_ea, ea_to_e, e_to_er, er_to_e = \
        generate_values(slopes, num_molecules, r, g, ma, mr, AC, RC, f)

    total_states = 3 * num_molecules
    W = np.zeros((total_states, total_states))

    def idx(row, col):
        return row * num_molecules + col

    # j = 2 : N-1 (MATLAB)
    for j in range(1, num_molecules-1):

        # -----------------------------
        # Activator row
        # -----------------------------
        current = idx(0, j)

        W[current, idx(1, j)] = ea_to_e[j]      # down
        W[current, idx(0, j+1)] = c[j]          # right
        W[current, idx(0, j-1)] = d[j-1]        # left

        # -----------------------------
        # Enhancer row
        # -----------------------------
        current = idx(1, j)

        W[current, idx(0, j)] = e_to_ea[j]      # up
        W[current, idx(2, j)] = e_to_er[j]      # down
        W[current, idx(1, j+1)] = a[j]          # right
        W[current, idx(1, j-1)] = b[j-1]        # left

        # -----------------------------
        # Repressor row
        # -----------------------------
        current = idx(2, j)

        W[current, idx(1, j)] = er_to_e[j]      # up
        W[current, idx(2, j+1)] = u[j]          # right
        W[current, idx(2, j-1)] = v[j-1]        # left

    return W


# ------------------------------------
# Parameters
# ------------------------------------
r = 7.0
g = 0.4
lam = 0.99
mu = 1.01

enh = 21
num_molecules = enh + 2   # 23 columns (with boundaries)

AC = 5
RC = 25.0
ma = 5.0
mr = 15.0

f = 1.0


# ------------------------------------
# Build matrix
# ------------------------------------
W = build_W(
    num_molecules=num_molecules,
    slopes=[lam, mu],
    r=r,
    g=g,
    ma=ma,
    mr=mr,
    AC=AC,
    RC=RC,
    f=f
)


print("W shape:", W.shape)
print(W)


# In[2123]:


import pandas as pd

# Load the data
df = pd.read_csv('probability_values.txt', header=None, sep=' ')


p = [0] * 3*enh


for index, row in df.iterrows():
    i = int(row[0])  
    if 1 <= i <= 3*enh:  
        p[i - 1] = row[1] 

for i in range(3*enh):
    print(f"p[{i + 1}] = {p[i]}")


# In[2124]:


import numpy as np

def extract_internal_block(W, enh):

    num_molecules = enh + 2
    total_states = 3 * num_molecules

    # Indices of internal columns (skip first and last column)
    internal_cols = range(1, num_molecules - 1)

    internal_indices = []

    for block in range(3):
        offset = block * num_molecules
        for col in internal_cols:
            internal_indices.append(offset + col)

    # Extract reduced matrix
    W_reduced = W[np.ix_(internal_indices, internal_indices)]

    return W_reduced


def calculate_entropy_production(W, p):

    num_states = W.shape[0]
    entropy_production = 0.0

    for i in range(num_states):
        for j in range(num_states-1):
            if W[i, j] > 0 and W[j, i] > 0 and p[i] > 0 and p[j] > 0:

                J_ij = W[i, j] * p[j] - W[j, i] * p[i]

                if J_ij != 0:
                    entropy_production += J_ij * np.log(
                        (W[i, j] * p[j]) / (W[j, i] * p[i])
                    )

    return entropy_production


# -------------------------
# W without boundary
# -------------------------
W_internal = extract_internal_block(W, enh)



# Compute entropy production
total_entropy_production = calculate_entropy_production(W_internal, p)

print("Total Entropy Production:", total_entropy_production)







