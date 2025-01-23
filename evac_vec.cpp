#include <math.h>
#include <unordered_set>
#include <vector>
#include <array>
#include <assert.h>
#include <stdio.h>
#include <cfloat>
#include <random>
#include <chrono>
#include "params.hpp"
#include "vfloat.hpp"

gvfloat2 init_p[nrep];
gvfloat init_v_desired[nrep];
int init_wp_offset[nrep][nent][nsample];


flt_t t_sim = 0.0;

vfloat angle(vfloat2& a, vfloat2& b) {
  vfloat angle_a = atan2(a[1], a[0]);
  vfloat angle_b = atan2(b[1], b[0]);

  vfloat angle_diff = angle_b - angle_a;

  vec_if (angle_diff > M_PI) {
    angle_diff = angle_diff - 2 * M_PI;
  } vec_fi

  vec_if (angle_diff <= -M_PI) {
    angle_diff = angle_diff + 2 * M_PI;
  } vec_fi

  return angle_diff;
}

vfloat2 compute_closest_point(flt_t o[2][2], vfloat2& p) {
  vfloat2 b_a_dist = { o[1][0] - o[0][0], o[1][1]- o[0][1] };
  vfloat2 p_a_dist = p - vfloat2(o[0][0], o[0][1]);

  vfloat lambda = dot(p_a_dist, b_a_dist) / (b_a_dist[0] * b_a_dist[0] + b_a_dist[1] * b_a_dist[1]);

  vfloat2 r;

  r = vfloat2(o[0][0], o[0][1]) + b_a_dist * lambda;
  
  vec_if (lambda <= 0.0) {
    r = { o[0][0], o[0][1] };
  } vec_fi

  vec_if (lambda >= 1.0) {
    r = { o[1][0], o[1][1] };
  } vec_fi

  return r;
}

vfloat2 compute_obstacle_force(vfloat2& p) {
  vfloat2 r = 0.0;
  for (size_t oid = 0; oid < num_obs; oid++) {
    vfloat2 closest_point = compute_closest_point(obstacles[oid], p);
    vfloat2 dist = p - closest_point;
    vfloat dist_norm = norm(dist);

    vfloat circle_dist = dist_norm - agent_radius;
    r += (dist / dist_norm) * exp(-1.0 * circle_dist / sigma);
  }

  return r;
}

vfloat waypoint_dist(vfloat2 p, vfloat2 next_waypoint) { return norm(p - next_waypoint); } 

vfloat2 left_normal(vfloat2& vec) { return {-1.0 * vec[1], vec[0]}; }

pair<int, int> getCellIndex(const flt_t x, const flt_t y) {
  return {static_cast<int>(x / cellSize),
          static_cast<int>(y / cellSize)};
}

vfloat sgn(vfloat& x) {
  vfloat r = -1.0;
  vec_if (x > 0) {
    r += 2;
  } vec_fi
  return r;
}

gvfloat all_evac_time[nrep];

flt_t log_gaussian_pdf(flt_t x, flt_t mean, flt_t stddev) {
    flt_t variance = stddev * stddev;
    return -0.5 * log(2 * M_PI * variance) - 0.5 * (x - mean) * (x - mean) / variance;
}
flt_t gaussian_overlap(flt_t mixture[num_in_comps][3]) {
    flt_t mu1 = mixture[0][1], sigma1 = mixture[0][2], mu2 = mixture[1][1], sigma2 = mixture[1][2];
    flt_t combined_variance = sigma1 * sigma1 + sigma2 * sigma2;
    flt_t diff_means = mu1 - mu2;
    flt_t exponent = - (diff_means * diff_means) / (2 * combined_variance);
    return (1.0 / sqrt(2 * M_PI * combined_variance)) * exp(exponent);
}

flt_t draw_from_mixture(flt_t mixture[num_in_comps][3], default_random_engine& rng) {
  uniform_real_distribution<flt_t> u_dist;
  flt_t u = uniform_real_distribution<flt_t>()(rng);

  for (int comp = 0; comp < num_in_comps; comp++) {
    if (mixture[comp][0] >= u) {
      normal_distribution<flt_t> norm_dist(mixture[comp][1], mixture[comp][2]);
      return norm_dist(rng);
    }
  }

  assert(false);
  return FLT_MAX;
}


void calc_likelihoods(flt_t ys[nsample])
{
  for (int sample = 0; sample < nsample; sample++) {
    ys[sample] = 0.0;
    for (int rep = 0; rep < nrep; rep++) {
      for (int aid = 0; aid < nent; aid++) {
        flt_t log_probs[num_out_comps];
        flt_t t = all_evac_time[rep][aid][sample];

        if (t == max_steps / delta_t)
          continue;

        for (int comp = 0; comp < num_out_comps; comp++) {
          log_probs[comp] = log(evac_time_ref[comp][0]) + log_gaussian_pdf(t, evac_time_ref[comp][1], evac_time_ref[comp][2]);
          ys[sample] += log_probs[comp];
        }
      }
    }
    ys[sample] /= -nrep;
  }
}

#if PRINT_DENSITY == true
  vector<double> ego_densities;
  vector<double> neighbor_densities;
  vector<double> neighbor_active_densities;
#endif

void evac(int rep)
{
  auto& evac_time = all_evac_time[rep];

  all_evac_time[rep] = max_steps / delta_t;

#if PRINT_TRACE == true
  FILE *trace_files[nsample];
  char filename[128];
  for (int sample = 0; sample < nsample; sample++) {
    snprintf(filename, sizeof(filename), "trace_vec_%d", sample);
    trace_files[sample] = fopen(filename, "w");
    if (trace_files[sample] == NULL) {
      perror("Error opening file");
      return;
    }
  }

  for (int sample = 0; sample < nsample; sample++) {
    fprintf(trace_files[sample], "width %.4f\n", scenario_width);
    for (size_t oid = 0; oid < num_obs; oid++) {
      auto& o = obstacles[oid];
      fprintf(trace_files[sample], "obstacle %.4f, %.4f; %4f, %.4f\n", o[0][0], o[0][1], o[1][0], o[1][1]);
    }
    fprintf(trace_files[sample], "waypoint tol %.4f\n", waypoint_tol);
    for (auto& wp : wp_alternatives) {
      fprintf(trace_files[sample], "waypoint %.4f, %.4f\n", wp[0][0], wp[0][1]);
    }
    fprintf(trace_files[sample], "t");
    for (int aid = 0; aid < nent; aid++)
      fprintf(trace_files[sample], ",a%d.active,a%d.y,a%d.x", aid, aid, aid);
    fprintf(trace_files[sample], "\n");
  }
#endif

  gvfloat2 p;
  gvfloat2 v;
  gvfloat2 a_old;
  
  gvfloat active = false;
  gvfloat2 next_waypoint;
  
  gvfloat t_spawn;
  gvfloat v_desired;



  int num_spawned_agents = 0;

  vfloat num_agents_near_wp[num_wp_alternatives];

  unordered_set<int> cand_nbs;
  for (int step = 0; step < max_steps; step++) {
    t_sim = step * delta_t;

    for (int wp_offset = 0; wp_offset < num_wp_alternatives; wp_offset++)
      num_agents_near_wp[wp_offset] = 0;

    for (int ego_aid = 0; ego_aid < nent; ego_aid++) {
      vec_if (active[ego_aid] > 0) {
        for (int wp_offset = 0; wp_offset < num_wp_alternatives; wp_offset++) {
          vfloat2 dist = { p[ego_aid][0] - wp_alternatives[wp_offset][1][0], p[ego_aid][1] - wp_alternatives[wp_offset][1][1] };
          vfloat dist_norm = norm(dist);
          vec_if (dist_norm < congestion_radius) {
            num_agents_near_wp[wp_offset] += 1;
          } vec_fi
        }
      } vec_fi
    }


    if (step % spawn_period == 0 && num_spawned_agents < nent) {
      int aid = num_spawned_agents;

      t_spawn[aid] = t_sim;
      v_desired[aid] = init_v_desired[rep][aid];
      p[aid] = init_p[rep][aid];

      for (int sample = 0; sample < nsample; sample++) {
        int wp_offset = init_wp_offset[rep][aid][sample];
        next_waypoint[aid][0][sample] = wp_alternatives[wp_offset][0][0];
        next_waypoint[aid][1][sample] = wp_alternatives[wp_offset][0][1];
      }


      v[aid] = {0.0, 0.0};
      a_old[aid] = {0.0, 0.0};
      active[aid] = true;
    
      num_spawned_agents++;
    }

    unordered_map<int, int> count_grid[numCells][numCells];
    vector<int> grid[numCells][numCells];
      
    for (int ego_aid = 0; ego_aid < nent; ego_aid++) {
      for (int sample = 0; sample < nsample; sample++) {
        if (!active[ego_aid][sample])
          continue;

        auto cell = getCellIndex(p[ego_aid][0][sample], p[ego_aid][1][sample]);
        if (cell.first < 0 || cell.first >= numCells ||
            cell.second < 0 || cell.second >= numCells) {
              active[ego_aid][sample] = false;
              continue;
        }
        count_grid[cell.first][cell.second][ego_aid]++;
      }
    }

    for (int x = 0; x < numCells; x++)
      for (int y = 0; y < numCells; y++)
        for (auto it : count_grid[x][y])
          if (it.second >= min_nb_samples)
            grid[x][y].push_back(it.first);
            

    for (int ego_aid = 0; ego_aid < nent; ego_aid++) {

      vec_if (active[ego_aid] > 0) {

#if PRINT_DENSITY == true
        ego_densities.push_back(cond_density());
#endif

        vfloat2 wp = next_waypoint[ego_aid];

        vfloat2 t_dist = wp - p[ego_aid];
        vfloat t_dist_norm = norm(wp - p[ego_aid]);
        vfloat2 e = t_dist / t_dist_norm;
        vfloat2 f_internal = e * v_desired[ego_aid] - v[ego_aid];
  
        bool visited_cell[numCells][numCells] = {};

        cand_nbs.clear();
        for (int sample = 0; sample < nsample; sample++) {
          auto egoCell = getCellIndex(p[ego_aid][0][sample], p[ego_aid][1][sample]);
          for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
              auto cx = egoCell.first + dx, cy = egoCell.second + dy;
              if (cx < 0 || cx >= numCells ||
                  cy < 0 || cy >= numCells ||
                  visited_cell[cx][cy])
                continue;

              visited_cell[cx][cy] = true;

              for (auto& other : grid[cx][cy])
                cand_nbs.insert(other);
            }
          }
        }


        vfloat2 f_interaction = {(flt_t)0.0, (flt_t)0.0};

        for (int other_aid : cand_nbs) {
          if (other_aid == ego_aid) continue;

          vec_if (active[other_aid] > (flt_t)0.0) {
#if PRINT_DENSITY == true
            neighbor_active_densities.push_back(cond_density());
#endif
            vfloat2 o_dist = p[other_aid] - p[ego_aid];
  
            vfloat o_dist_norm = norm(o_dist);

            vec_if (o_dist_norm < interaction_radius) {

#if PRINT_DENSITY == true
              neighbor_densities.push_back(cond_density());
#endif

              vfloat2 o_dir = o_dist / o_dist_norm;
  
              vfloat2 v_diff = v[ego_aid] - v[other_aid];
  
              vfloat2 int_v = lambda * v_diff + o_dir;
  
              vfloat int_norm = norm(int_v);
              vfloat2 int_dir = int_v / int_norm;
  
              vfloat theta = angle(int_dir, o_dir);
              
              vfloat theta_sign = (flt_t)0.0;
              vec_if (theta < (flt_t)0.0) {
                theta_sign = theta / abs(theta);
              } vec_fi

              vec_if (theta > (flt_t)0.0) {
                theta_sign = theta / abs(theta);
              } vec_fi

              vfloat B = gamma_ * int_norm;
  
              vfloat npBt = n_prime * B * theta;
              vfloat nBt = n * B * theta;

              vfloat2 force_v = (flt_t)-1.0 * int_dir * exp((flt_t)-1.0 * o_dist_norm / B - npBt * npBt);
              vfloat2 force_angle = (flt_t)-1.0 * left_normal(int_dir) * exp((flt_t)-1.0 * o_dist_norm / B - nBt * nBt) * theta_sign;
  
              f_interaction += force_v + force_angle;
            } vec_fi
          } vec_fi
        }
  
        vfloat2 f_obstacles = compute_obstacle_force(p[ego_aid]);
  
        vfloat2 a = f_internal * w_internal +
                 f_interaction * w_interaction +
                 f_obstacles * w_obstacles;

        vfloat2 p_new = p[ego_aid] + v[ego_aid] * delta_t + 0.5 * a_old[ego_aid] * delta_t * delta_t;
        vfloat2 v_new = v[ego_aid] + 0.5 * (a_old[ego_aid] + a) * delta_t;

        p[ego_aid] = p_new;
        v[ego_aid] = v_new;

        for (int wp_offset = 0; wp_offset < num_wp_alternatives; wp_offset++) {
          vfloat2 dist = { p[ego_aid][0] - wp_alternatives[wp_offset][0][0], p[ego_aid][1] - wp_alternatives[wp_offset][0][1] };
          vfloat dist_norm = norm(dist);
          vec_if (dist_norm < waypoint_tol) {
            active[ego_aid] = false;
            evac_time[ego_aid] = t_sim - t_spawn[ego_aid];
          } vec_fi
        }
      } vec_fi
      assert(flt_conds_vec[0][0] > 0.0);
    }

#if PRINT_TRACE == true
    for (int sample = 0; sample < nsample; sample++) {
      fprintf(trace_files[sample], "%.6f", t_sim);
      for (int ego_aid = 0; ego_aid < nent; ego_aid++)
        fprintf(trace_files[sample], ",%.0f,%.6f,%.6f", active[ego_aid][sample], p[ego_aid][0][sample], p[ego_aid][1][sample]);
      fprintf(trace_files[sample], "\n");
    }
#endif
  }

#if PRINT_TRACE == true
  for (int sample = 0; sample < nsample; sample++)
    fclose(trace_files[sample]);
#endif

}


int main(int argc, char **argv)
{
  flt_conds_vec.push_back(true);

  flt_t pert_stddev = atof(argv[1]);

  flt_t xs[num_in_comps][3];
  int i = 2;
  for (int comp = 0; comp < num_in_comps; comp++)
    for (int j = 0; j < 3; j++)
      xs[comp][j] = atof(argv[i++]);

  uniform_int_distribution<int> wp_dist(0, num_wp_alternatives - 1);
  normal_distribution<flt_t> v_desired_dist(v_desired_mean, v_desired_stddev);
  uniform_real_distribution<flt_t> spawn_p_dist(wall_margin_x * scenario_width + agent_radius * 2, (1 - wall_margin_x) * scenario_width - agent_radius * 2);

  normal_distribution<flt_t> pert_dist(0, 1);
  flt_t perturbations[nsample][num_in_comps][3];

  default_random_engine pert_rng;
  pert_rng.seed(random_device()());
  for (int sample = 0; sample < nsample; sample++)
    for (int comp = 0; comp < num_in_comps; comp++)
      for (int i = 0; i < 3; i++)
        perturbations[sample][comp][i] = (sample == 0 ? 0.0 : pert_dist(pert_rng));

  int sim_seed = 1;
  for (int sample = 0; sample < nsample; sample++) {
    flt_t x_perturbed[num_in_comps][3];
    for (int comp = 0; comp < num_in_comps; comp++) {
      project_param(xs[comp], perturbations[sample][comp], x_perturbed[comp], pert_stddev);
      x_perturbed[comp][0] += (comp > 0 ? x_perturbed[comp - 1][0] : 0.0);
    }

    for (int comp = 0; comp < num_in_comps; comp++) {
      x_perturbed[comp][0] /= x_perturbed[num_in_comps - 1][0];
      if (sample == 0) {
        printf("( ");
        for (int i = 0; i < 3; i++)
          printf("%.8f ", x_perturbed[comp][i]);
        printf(")\n");
      }
    }

    default_random_engine sim_rng(sim_seed);
    for (int rep = 0; rep < nrep; rep++) {

      for (int aid = 0; aid < nent; aid++) {
        init_v_desired[rep][aid][sample] = max(min_v_desired, min(max_v_desired, draw_from_mixture(x_perturbed, sim_rng)));
        init_p[rep][aid][0][sample] = 0.0;
        init_p[rep][aid][1][sample] = spawn_p_dist(sim_rng);
        init_wp_offset[rep][aid][sample] = wp_dist(sim_rng);
      }
    }

  }

  auto t0 = chrono::high_resolution_clock::now();
  for (int rep = 0; rep < nrep; rep++)
    evac(rep);
  auto t1 = chrono::high_resolution_clock::now();

  flt_t scalar_us = 1e-6 * chrono::duration_cast<chrono::microseconds>(t1 - t0).count();
  printf("vec took %.4fs\n", scalar_us);

#if PRINT_DENSITY == true
  printf("agent updates: %lu\n", ego_densities.size());
  printf("neighbor accesses: %lu\n", neighbor_densities.size());

  printf("avg. ego density: %.6f\n", std::accumulate(ego_densities.begin(), ego_densities.end(), 0.0) / ego_densities.size());
  printf("avg. neighbor active density: %.6f\n", std::accumulate(neighbor_active_densities.begin(), neighbor_active_densities.end(), 0.0) / neighbor_active_densities.size());
  printf("avg. neighbor density: %.6f\n", std::accumulate(neighbor_densities.begin(), neighbor_densities.end(), 0.0) / neighbor_densities.size());
#endif

  flt_t ys[nsample];

  calc_likelihoods(ys);

  flt_t mean = 0.0;
  for (int sample = 0; sample < nsample; sample++) {
    printf("y_%d: %.4f\n", sample, ys[sample]);
    if (sample > 0)
      mean += ys[sample] / (nsample - 1);
  }
  printf("mean y: %.4f\n", mean);
  printf("crisp y: %.4f\n", ys[0]);

  for (int comp = 0; comp < num_in_comps; comp++) {
    for (int i = 0; i < 3; i++) {
      flt_t deriv = 0.0;
      for (int sample = 1; sample < nsample; sample++) {
        const flt_t crisp_y = ys[0];
        deriv += ((ys[sample] - crisp_y) / pert_stddev * perturbations[sample][comp][i]) / nsample;
      }
      printf("dy/dx_%d: %.6f\n", comp * 3 + i, deriv);
    }
  }

  return 0;
}
