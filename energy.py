#!/usr/bin/env python
# coding: utf-8

# In[2118]:


import re

# Read data from file
with open("enhancers_matrix_21_r_0.40000_g_0.40000_s1_0.990_s2_1.010_ma_5.00_mr_15.00_RC_25.00_AC5_t1.txt", "r") as file:
    lines = file.readlines()

# Initialize variables
step_results = {}  

# Iterate through lines
i = 0
while i < len(lines):
    # Extract step number using regular expression
    step_number_match = re.match(r"Step ([\d.]+):", lines[i])
    if step_number_match:
        step_number = float(step_number_match.group(1))
    else:
        # Skip this line if it doesn't match the expected format
        i += 1
        continue

    # Extract matrix for this step
    matrix_lines = lines[i+1:i+4]
    matrix = [[float(x) for x in line.split()] for line in matrix_lines]

    # Calculate expression for this step
    results = [(n)* matrix[0][n] + (n+21) * matrix[1][n] + (n+42) * matrix[2][n] for n in range(len(matrix[0]))]

    step_results[step_number] = sum(results)

    # Move to next step
    i += 5  # Move to the next "Step X:" line (4 lines for the matrix + 1 blank line)

# Write results to file
with open("GE_distribution_lambda_5.txt", "w") as outfile:
    for step, result in step_results.items():
        outfile.write(f"{step} {result}\n")



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




#horizontal_transition(excluding bundaries)
    for i in range(2, num_molecules -1):
        a[i] = f * r * (slopes[0] ** (i - 2))
        b[i] = f * g * (slopes[1] ** (i - 2))
        c[i] = f * ma * r * (slopes[0] ** (i - 2))
        d[i] = f * g * (slopes[1] ** (i - 2))
        u[i] = f * r * (slopes[0] ** (i - 2))
        v[i] = f * mr * g * (slopes[1] ** (i - 2))

    # Vertical transitions
    for i in range(2, num_molecules):
        e_to_ea[i] = f * AC * r
        ea_to_e[i] = f * g
        e_to_er[i] = f * RC * r
        er_to_e[i] = f * g

    return a, b, c, d, u, v, e_to_ea, ea_to_e, e_to_er, er_to_e



def build_W(num_molecules, slopes, r, g, ma, mr, AC, RC, f):

    a, b, c, d, u, v, e_to_ea, ea_to_e, e_to_er, er_to_e = \
        generate_values(slopes, num_molecules, r, g, ma, mr, AC, RC, f)

    total_states = 3 * num_molecules
    W = np.zeros((total_states, total_states))

    def idx(row, col):
        return row * num_molecules + col

  
    for j in range(1, num_molecules-1):

        # -----------------------------
        # Activator row (row = 0)
        # -----------------------------
        current = idx(0, j)

        # down (to enhancer)
        W[current, idx(1, j)] = ea_to_e[j]

        # right
        W[current, idx(0, j+1)] = c[j]

        # left
        W[current, idx(0, j-1)] = d[j-1]

        # -----------------------------
        # Enhancer row (row = 1)
        # -----------------------------
        current = idx(1, j)

        # up (to activator)
        W[current, idx(0, j)] = e_to_ea[j]

        # down (to repressor)
        W[current, idx(2, j)] = e_to_er[j]

        # right
        W[current, idx(1, j+1)] = a[j]

        # left
        W[current, idx(1, j-1)] = b[j-1]

        # -----------------------------
        # Repressor row (row = 2)
        # -----------------------------
        current = idx(2, j)

        # up (to enhancer)
        W[current, idx(1, j)] = er_to_e[j]

        # right
        W[current, idx(2, j+1)] = u[j]

        # left
        W[current, idx(2, j-1)] = v[j-1]


    return W

r = 0.4
g = 0.4
lam = 0.99
mu = 1.01
enh = 21
num_molecules = enh + 2  # 23(including boundaries in gillespie)
AC = 5
RC = 25.0
ma = 5.0
mr = 15.0
f = 1.0
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


# print(W)


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







