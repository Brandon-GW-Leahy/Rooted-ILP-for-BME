import random
import math
import numpy as np
import gurobipy as gp
from gurobipy import quicksum
from gurobipy import GRB
from openpyxl import Workbook




def sample_point_in_unit_ball():
    while True:
        x = random.uniform(-1, 1)
        y = random.uniform(-1, 1)
        z = random.uniform(-1, 1)
        if x * x + y * y + z * z <= 1:
            return (x, y, z)


def euclidean_distance(p, q):
    return math.sqrt(
        (p[0] - q[0]) ** 2 +
        (p[1] - q[1]) ** 2 +
        (p[2] - q[2]) ** 2
    )


def build_distance_matrix(n):
    # Set D = [[0 for _ in range(n)] for _ in range(n)]
    D = [[0.0 for _ in range(n)] for _ in range(n)]

    # Select n random points in the unit 3-ball
    points = [sample_point_in_unit_ball() for _ in range(n)]

    # Compute pairwise Euclidean distances
    for i in range(n):
        for j in range(i):
            d = euclidean_distance(points[i], points[j])
            D[i][j] = d
            D[j][i] = d

    return D





def generate_random_dissimilarity_matrix(n):
  # We generate a random binary tree
  def random_binary_tree(L):
    L = list(L)

    while len(L) > 1:
        # randomly pick two distinct elements
        x, y = random.sample(L, 2)

        # remove them
        L.remove(x)
        L.remove(y)

        # merge them into a binary node
        L.append([x, y])

    return L[0]

  L = list(range(n))
  B = random_binary_tree(L)

  # We assign bit strings
  def generate_codes(tree, prefix=None, codes=None):
    if codes is None:
        codes = {}
    if prefix is None:
        prefix = []

    if not isinstance(tree, list):
        codes[tree] = prefix
        return codes

    if len(tree) >= 1:
        generate_codes(tree[0], prefix + [0], codes)

    if len(tree) >= 2:
        generate_codes(tree[1], prefix + [1], codes)

    return codes

  B = generate_codes(B)

  # We compute the rooted path distance matrix $D_r$.
  D = [[0 for _ in range(n)] for _ in range(n)]

  for i in range(n):
    for j in range(i):
      l = 0
      while B[i][l] == B[j][l]:
        l = l + 1

      D[i][j] = len(B[i]) + len(B[j]) - 2*l
      D[j][i] = len(B[i]) + len(B[j]) - 2*l

  
  # We compute the unrooted path distance matrix Du.
  for i in range(n):
    for j in range(i):
      if B[i][0] != B[j][0]:
        D[i][j] = D[i][j] - 1
        D[j][i] = D[j][i] - 1
  

  # We add noise to the data
  for i in range(n):
    for j in range(i):
      noise = random.randint(0, 9)
      D[i][j] += noise
      D[j][i] = D[i][j]

  return D






# Rooted ILP 
def Rooted_model(D,L):
  n = len(D)

  model = gp.Model("Rooted_ILP")

  # Create the pairs of taxa (i,j) collated into a list E
  E = [(i, j) for j in range(n) for i in range(j)]

  # Create the binary numbers s_i
  s_il = {}
  for i in range(n):
    for l in range(L):
      s_il[(i, l)] = model.addVar(vtype=gp.GRB.BINARY, name=f"s_{i}_{l}")

  # Create the depth numbers d_i
  d_i = {}
  for i in range(n):
    d_i[i] = model.addVar(vtype=gp.GRB.INTEGER, lb=1, ub=L, name=f"d_{i}")

  # Create the variables z_el
  z_el = {}
  for e in E:
    for l in range(L):
      z_el[(e,l)] = model.addVar(vtype=gp.GRB.BINARY, name=f"z_{e}_{l}")

  # Create the variables b_il
  b_il = {}
  for i in range(n):
    for l in range(1,L+1):
      b_il[(i, l)] = model.addVar(vtype=gp.GRB.BINARY, name=f"b_{i}_{l}")

  # Create the variables x_el
  x_el = {}
  for e in E:
    for l in range(L):
      x_el[(e,l)] = model.addVar(vtype=gp.GRB.BINARY, name=f"x_{e}_{l}")

  # Create the variables y_el
  y_el = {}
  for e in E:
    for l in range(L):
      y_el[(e,l)] = model.addVar(vtype=gp.GRB.BINARY, name=f"y_{e}_{l}")

  # Create the variables c_el
  c_el = {}
  for e in E:
    for l in range(1,2*L+1):
      c_el[(e,l)] = model.addVar(vtype=gp.GRB.BINARY, name=f"c_{e}_{l}")

  # We add the objective function.
  model.setObjective(quicksum(D[e[0]][e[1]] * 2**(-l+1) * c_el[(e,l)] for e in E for l in range(1, 2*L + 1)), gp.GRB.MINIMIZE)

  # We add the constraints.
  # Constraint (1)
  for e in E:
    for l in range(L):
      model.addConstr(z_el[(e,l)] >= s_il[(e[0],l)] - s_il[(e[1],l)])
      model.addConstr(z_el[(e,l)] >= s_il[(e[1],l)] - s_il[(e[0],l)])
      model.addConstr(z_el[(e,l)] <= s_il[(e[0],l)] + s_il[(e[1],l)])
      model.addConstr(z_el[(e,l)] <= 2 - (s_il[(e[0],l)] + s_il[(e[1],l)]))

  # Constraint (2)
  for e in E:
    model.addConstr(quicksum(z_el[(e,l)] for l in range(0,L)) >= 1)

  # Constraint (3)
  for i in range(n):
    model.addConstr(quicksum(b_il[(i,l)] for l in range(1,L+1)) == 1)
    model.addConstr(quicksum(l*b_il[(i,l)] for l in range(1,L+1)) == d_i[i])

  # Constraint (4)
  model.addConstr(quicksum(b_il[(i,l)]*2**(-l) for i in range(n) for l in range(1,L+1)) == 1)

  # Constraint (5)
  for e in E:
    model.addConstr(x_el[(e,0)] == z_el[(e,0)])
    model.addConstr(y_el[(e,0)] == z_el[(e,0)])

  # Constraint (6)
  for e in E:
    for l in range(1,L):
      model.addConstr(x_el[(e,l)] >= x_el[(e,l-1)])
      model.addConstr(x_el[(e,l)] >= z_el[(e,l)])
      model.addConstr(x_el[(e,l)] <= x_el[(e,l-1)] + z_el[(e,l)])

  # Constraint (7)
  for e in E:
    for l in range(1,L):
      model.addConstr(y_el[(e,l)] <= z_el[(e,l)])
      model.addConstr(y_el[(e,l)] <= 1 - x_el[(e,l-1)])
      model.addConstr(y_el[(e,l)] >= z_el[(e,l)] - x_el[(e,l-1)])

  # Constraint (8)
  for e in E:
    model.addConstr(quicksum(l*y_el[(e,l)] for l in range(L)) <= d_i[e[0]]-1)
    model.addConstr(quicksum(l*y_el[(e,l)] for l in range(L)) <= d_i[e[1]]-1)

  # Constraint (9)
  for e in E:
    model.addConstr(quicksum(c_el[(e,l)] for l in range(1,2*L+1)) == 1)
    model.addConstr(quicksum(l*c_el[(e,l)] for l in range(1,2*L+1)) == d_i[e[0]] + d_i[e[1]] - 2*quicksum(l*y_el[(e,l)] for l in range(L)) - y_el[(e,0)])

  model.update()
  return model






# Path-Length ILP
# file: phylo_ilp.py
def Path_Length_model(D,L):
    n = len(D)
    Gamma = range(1, n + 1)
    L = range(2, L+1)

    model = gp.Model("Path_Length_ILP")

    # x_{ij}^l
    x = model.addVars([(i, j, l) for i in Gamma for j in Gamma if i != j for l in L], vtype=GRB.BINARY, name="x")

    # tau_{ij}
    tau = model.addVars([(i, j) for i in Gamma for j in Gamma if i != j], vtype=GRB.INTEGER, lb=2, ub=n - 1, name="tau")

    # y_{pq}^j
    nodes_no1 = [i for i in Gamma if i != 1]
    y = model.addVars([(j, p, q) for j in nodes_no1 for p in nodes_no1 for q in nodes_no1 if len({j, p, q}) == 3], vtype=GRB.BINARY, name="y")

    # Objective (11a)
    model.setObjective(gp.quicksum(D[i-1][j-1] * gp.quicksum((2 ** (-l)) * x[i, j, l] for l in L) for i in Gamma for j in Gamma if i != j), GRB.MINIMIZE)

    # (11b)
    for i in Gamma:
        for j in Gamma:
            if i != j:
                model.addConstr(gp.quicksum(x[i, j, l] for l in L) == 1)

    # (11c)
    for i in Gamma:
        for j in Gamma:
            if i < j:
                for l in L:
                    model.addConstr(x[i, j, l] == x[j, i, l])

    # (11d)
    for i in Gamma:
        for j in Gamma:
            if i != j:
                model.addConstr(tau[i, j] == gp.quicksum(l * x[i, j, l] for l in L))

    # (11e)
    for i in Gamma:
        model.addConstr(gp.quicksum((2 ** (-l)) * x[i, j, l] for j in Gamma if j != i for l in L) == 0.5)

    # (11f)
    for i in Gamma:
        for j in Gamma:
            if i != j and i != 1 and j != 1:
                model.addConstr(tau[1, i] + tau[1, j] - tau[i, j] >= 2)

    # (11g)
    for i in Gamma:
        for j in Gamma:
            if i != j and i != 1 and j != 1:
                model.addConstr(tau[1, i] + tau[i, j] - tau[1, j] >= 2)

    # (11h)
    for j in nodes_no1:
        for p in nodes_no1:
            for q in nodes_no1:
                if len({j, p, q}) == 3:
                    model.addConstr(y[j,p,q] + y[q,j,p] + y[p,j,q] == 1)

    # (11i)
    M = 2 * n - 2
    for j in nodes_no1:
        for p in nodes_no1:
            for q in nodes_no1:
                if len({j, p, q}) == 3:
                    model.addConstr(tau[1, j] + tau[p, q] >= 2 * (1 - y[q, j, p]) + tau[1, p] + tau[j, q] - M * y[j, p, q])

    # (11j)
    model.addConstr(gp.quicksum((l * 2**(- l)) * x[i, j, l] for i in Gamma for j in Gamma if i != j for l in L) == 2 * n - 3)

    model.update()
    return model










n = int(input("Number of taxa = "))
L = int(input("Maximum length of rooted tree = "))
r = 2*L - 1
x = int(input("Number of trials = "))
y = int(input("Time limit (hrs) = "))

wb = Workbook()
ws = wb.active
ws.title = "Results"
ws.append(["Best Rooted ILP Objective", "Rooted ILP Run Time", "Best Path-Length ILP Objective", "Path-Length ILP Run Time"])

for s in range(x):
  D = build_distance_matrix(n)

  time_limit_seconds = y * 3600



  model = Rooted_model(D,L)
  model.Params.TimeLimit = time_limit_seconds
  model.optimize()
  
  if model.SolCount > 0:
    our_obj = model.ObjVal
  else:
    our_obj = None
  our_time = model.Runtime




  model = Path_Length_model(D,r)
  model.Params.TimeLimit = time_limit_seconds
  model.optimize()
  
  if model.SolCount > 0:
    their_obj = model.ObjVal
  else:
    their_obj = None
  their_time = model.Runtime




  ws.append([our_obj, our_time, their_obj, their_time])



wb.save("ILP_comparison_results.xlsx")