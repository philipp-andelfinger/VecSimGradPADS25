#!/bin/bash

#for start_rep in 1 11 21 31 41 51 61 71 81 91; do
for start_rep in 1; do
  end_rep=$(expr $start_rep + 9)

  for variant in scalar_opt vec_opt; do
    for stddev in 1e-3 1e-2 1e-1; do
      for rep in $(seq $start_rep $end_rep); do
        ./opt.py $variant $stddev $rep &
      done
    done
  done
  
  wait

done
