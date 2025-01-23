#!/usr/bin/python3

import subprocess
import re
import numpy as np
from scipy.optimize import minimize
import sys
import time

num_params = 6
max_exec_time = 2 * 3600

prog = f"./{sys.argv[1]}"
pert_stddev = sys.argv[2]
rep = int(sys.argv[3])

np.random.seed(rep)

n = np.random.normal
n_stddev = 0.1

# assumed prior near the target mixture
def get_init_params():
  return np.array([1.5708 + n(0, n_stddev), 0.645 + n(0, n_stddev), -0.1113 + n(0, n_stddev), -0.285 + n(0, n_stddev), -0.645 + n(0, n_stddev), -0.589 + n(0, n_stddev)])

def external_function(params, executable):
    args = [executable, pert_stddev] + list(map(str, params))
    result = subprocess.run(args, capture_output=True, text=True)
    print(result.stdout)
    gradients = re.findall(r"dy/dx_.+: (.+)", result.stdout)
    exec_time = re.findall(r"took (.+)s", result.stdout)[0]
    y = re.findall(r"mean y: (.+)", result.stdout)[0]
    y_crisp = re.findall(r"crisp y: (.+)", result.stdout)[0]
    return np.array(list(map(float, gradients))), float(exec_time), float(y_crisp), float(y)

def optimize(dim, lr=1e-4, beta1=0.9, beta2=0.999, epsilon=1e-8, iterations=100000000, warmup=10):
    params = get_init_params()
    m, v = np.zeros(dim), np.zeros(dim)
    total_exec_time = 0
    with open(f"out_{sys.argv[1]}_{pert_stddev}_{rep}.csv", "w") as f:
      f.write(f"step,exec_time,y_crisp,y_crisp_scalar,y")
      for i in range(num_params):
        f.write(f",p{i}")
      f.write("\n")
      for step in range(iterations):
        gradients, exec_time, y_crisp, y = external_function(params, prog)
        total_exec_time += exec_time
        y_crisp_scalar = y_crisp
        print(y_crisp, y_crisp_scalar, y)
        f.write(f"{step},{total_exec_time},{y_crisp},{y_crisp_scalar},{y}")
        for p in params:
          f.write(f",{p}")
        f.write("\n")
        f.flush()
        params -= gradients * lr

        if total_exec_time > max_exec_time:
          return params

    return params

optimized_params = optimize(num_params)
