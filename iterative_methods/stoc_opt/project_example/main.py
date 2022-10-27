# GENERALIZE BIN PROBLEM
# Imports
import time

import gurobipy as grb
import numpy as np

np.random.seed(1)


def solve():
    pass


# Defining initial variables 
n_bins = 10
n_compulsory_items = 3
n_non_compulsory_items = 3
n_items = n_compulsory_items + n_non_compulsory_items
# Generate cost
C = np.random.uniform(100, 120, n_bins)
# Generate item's profits
p = np.random.uniform(10, 12, n_items)
# Generate weights
w = np.random.uniform(1, 5, n_non_compulsory_items+n_compulsory_items)
W = np.random.uniform(20, 404 , n_bins)

model= grb.Model('gbpp')  # type: ignore

# Y = 1 if bin j is considered
Y = model.addVars(
    n_bins,
    vtype=grb.GRB.BINARY,
    name='Y'
)

X = model.addVars(
    n_items, n_bins,
    vtype=grb.GRB.BINARY,
    name='X'
)

bins = range(n_bins)
expr = sum(
    C[j]*Y[j] for j in bins
)

expr -= grb.quicksum(p[i]*X[i, j] for i in range(n_compulsory_items, n_items) for j in bins)

model.setObjective(expr, grb.GRB.MINIMIZE)

model.update()

model.addConstrs(
    (grb.quicksum(w[i]*X[i,j] for i in range(n_items))<= W[j]*Y[j] for j in bins),
    name="capacity_contraint"
)
model.addConstrs(
    (grb.quicksum(X[i,j] for j in bins) == 1 for i in range(n_compulsory_items)),
    name="compulsory_item"
)

# [0,1,...n compulsory items - 1,n_compulsory_items, ..., n_items]
model.addConstrs(
    (grb.quicksum(X[i,j] for j in bins) <= 1 for i in range(n_compulsory_items, n_items)),
    name="non_compulsory_item"
)
# model.setParam('MIPgap', gap)
# model.setParam(grb.GRB.Param.TimeLimit, time_limit)

model.setParam('OutputFlag', 1)

model.setParam('LogFile', 'gurobi.log')
model.write("model.lp")
start = time.time()
model.optimize()
end = time.time()
comp_time = end - start
print(f'Time: {comp_time}')

"""if model.status == grb.GRB.Status.OPTIMAL:
    sol = [Y[i].X for i in bins if Y[i].X > 0.5]
"""

for j in bins:
    if Y[j].X > 0.5:
        print(Y[j].X)
        print("Compulsory")
        for i in range(n_compulsory_items):
            print(X[i,j].X)
        print("No Compulsory")
        for i in range(n_compulsory_items, n_items):
            print(X[i, j].X)

