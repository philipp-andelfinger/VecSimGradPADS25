#!/usr/bin/python3

import os
import subprocess
import numpy as np

clang_bin = "clang++"

num_iter = 10
print(
    "stddev,nsample,flt_t,vec_width,spawn_period,t_scalar,t_vec,speedup,num_ego_updates,num_nb_accesses,ego_density,nb_density"
)

params = "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
print_trace = False

for stddev in [0, 1e-3, 1e-2, 1e-1]:
  for nsample in [128, 64, 32, 16, 8]:
      for flt_t in ["float", "double"]:
        for vec_width in [512]:
          for spawn_period in [1.0]:
            min_nb_samples = 1
            shared_flags = (
                f"-DNDEBUG -O3 -march=native -ffast-math -funroll-loops "
                f"-std=c++20 -Wall -DNENT=128 -DNREP=1 -DMIN_NB_SAMPLES={min_nb_samples} -DSPAWN_PERIOD={spawn_period} "
                f"-DFLT_T={flt_t} -DPRINT_TRACE={str(print_trace).lower()} -isystem ."
            )

            os.system(f"rm -f scalar vec vec_print_density")

            cmp_scalar = f"{clang_bin} evac.cpp {shared_flags} -DNSAMPLE={nsample} -o scalar"
            cmp_vec = f"{clang_bin} evac_vec.cpp {shared_flags} -DNSAMPLE={nsample} -mprefer-vector-width={vec_width} -DPRINT_DENSITY=false -o vec"
            cmp_vec_density = f"{clang_bin} evac_vec.cpp {shared_flags} -DNSAMPLE={nsample} -mprefer-vector-width={vec_width} -DPRINT_DENSITY=true -o vec_print_density"

            subprocess.run(f"{cmp_scalar} & {cmp_vec} & {cmp_vec_density} & wait", shell=True)

            if nsample == 64 and flt_t == "float":
              os.system("cp scalar scalar_opt; cp vec vec_opt")

            t_scalar, t_vec = 0, 0
            for i in range(1, num_iter + 1):
              np.random.seed(i)
              params = list(map(float, np.random.random(6) * 2 * np.pi - np.pi))
              scalar_out = subprocess.run(f"./scalar {stddev} {params}", shell=True, capture_output=True, text=True).stdout

              vec_out = subprocess.run(f"./vec {stddev} {params}", shell=True, capture_output=True, text=True).stdout 

              t_scalar += float( next((line.split()[2][:-1] for line in scalar_out.splitlines() if "took" in line), 0))
              t_vec += float(next((line.split()[2][:-1] for line in vec_out.splitlines() if "took" in line), 0))

            speedup = t_scalar / t_vec if t_vec > 0 else 0
            vec_density_out = subprocess.run(f"./vec_print_density {stddev} {params}", shell=True, capture_output=True, text=True).stdout 

            updates = int(next((line.split()[2] for line in vec_density_out.splitlines() if line.startswith("agent updates")), 0))
            nb_acc = int(next((line.split()[2] for line in vec_density_out.splitlines() if line.startswith("neighbor accesses")), 0))
            ego_density = float(next((line.split()[3] for line in vec_density_out.splitlines() if line.startswith("avg. ego density")), 0))
            nb_density = float(next((line.split()[3] for line in vec_density_out.splitlines() if line.startswith("avg. neighbor density")), 0))

            print(
                f"{stddev},{nsample},{flt_t},{vec_width},{spawn_period},{t_scalar},{t_vec},{speedup},{updates},{nb_acc},{ego_density},{nb_density}"
            )
