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
#include "vec2.hpp"

#define SGN(X) ((X > 0) ? 1 : -1)

#if PRINT_TRACE == true
const int trace_sample = 1;
#endif

typedef vec2<flt_t> flt2;
typedef vec2<int> int2;

struct InitAgentState {
  flt2 p;
  flt_t v_desired;
  int wp_offset;
  int waypoint_timer;
};

InitAgentState init_agent_states[nsample][nrep][nent];

flt_t sum_evac_time;
int num_evac;

flt_t t_sim = 0.0;

uniform_real_distribution<flt_t> u_dist;

pair<int, int> getCellIndex(const flt2 &pos) {
  return {static_cast<int>(pos[0] / cellSize),
          static_cast<int>(pos[1] / cellSize)};
}

flt_t norm(const flt2& vec) {
  flt2 s = vec * vec;

  return sqrt(s[0] + s[1]);
}

flt_t angle(flt2& a, flt2& b) {
  flt_t angle_a = atan2(a[1], a[0]);
  flt_t angle_b = atan2(b[1], b[0]);

  flt_t angle_diff = angle_b - angle_a;

  if (angle_diff > M_PI) angle_diff -= 2 * M_PI;
  else if (angle_diff <= -M_PI) angle_diff += 2 * M_PI;

  return angle_diff;
}

flt_t dot(flt2& a, flt2&b) {
  return a[0] * b[0] + a[1] * b[1];
}

flt2 compute_closest_point(flt_t o[2][2], flt2& p) {
  flt2 b_a_dist = { o[1][0] - o[0][0], o[1][1]- o[0][1] };
  flt2 p_a_dist = p - flt2({o[0][0], o[0][1]});

  flt_t lambda = dot(p_a_dist, b_a_dist) / (b_a_dist[0] * b_a_dist[0] + b_a_dist[1] * b_a_dist[1]);

  if (lambda <= 0)
    return { o[0][0], o[0][1] };
  else if (lambda >= 1)
    return { o[1][0], o[1][1] };
  return flt2({o[0][0], o[0][1]}) + b_a_dist * lambda;

}

flt2 compute_obstacle_force(flt2& p) {
  flt2 r = { 0.0, 0.0 };
  for (size_t oid = 0; oid < num_obs; oid++) {
    flt2 closest_point = compute_closest_point(obstacles[oid], p);
    flt2 dist = p - closest_point;
    flt_t dist_norm = norm(dist);

    flt_t circle_dist = dist_norm - agent_radius;
    r += (dist / dist_norm) * exp(-1.0 * circle_dist / sigma);
  }

  return r;
}

struct Agent {
  flt2 p;
  flt2 v;

  flt2 a_old;

  bool active = false;
  flt2 waypoint;
  flt2 prev_wp = { -1, -1 };

  flt_t t_spawn;

  flt_t v_desired;

  int waypoint_timer;

  Agent() { }

  Agent(int sample, int rep, int aid) {
    t_spawn = t_sim;
    v_desired = init_agent_states[sample][rep][aid].v_desired;
    p = init_agent_states[sample][rep][aid].p;
  }

  flt_t waypoint_dist() {
    flt2& wp = waypoint;
    if (wp[0] == -1)
      return abs(p[1] - wp[1]);
    if (wp[1] == -1)
      return abs(p[0] - wp[0]);
    
    flt2 dist = p - waypoint;
    return norm(dist);
  }
};

Agent agents[nent];

void spawn_agent(int sample, int aid, int rep) {
  if (aid >= nent)
    return;
  auto& ego = agents[aid];

  ego = Agent(sample, rep, aid);

  int wp_offset = init_agent_states[sample][rep][aid].wp_offset;
  ego.waypoint = {wp_alternatives[wp_offset][0][0], wp_alternatives[wp_offset][0][1]};

  ego.waypoint_timer = init_agent_states[sample][rep][aid].waypoint_timer;

  ego.v = {0.0, 0.0};
  ego.a_old = {0.0, 0.0};
  ego.active = true;
}

flt2 left_normal(flt2& vec) { return {-vec[1], vec[0]}; }

flt_t all_evac_time[nrep][nent][nsample];

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
      //printf("drawing from N(%.2f, %.2f)\n", mixture[comp][1], mixture[comp][2]);
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

        //printf("got one at %.2f\n", t);

        for (int comp = 0; comp < num_out_comps; comp++) {
          log_probs[comp] = log(evac_time_ref[comp][0]) + log_gaussian_pdf(t, evac_time_ref[comp][1], evac_time_ref[comp][2]);
          ys[sample] += log_probs[comp];
        }
      }
    }
    ys[sample] /= -nrep;
  }
}

void evac(int rep)
{
  auto& evac_time = all_evac_time[rep];

#if PRINT_TRACE == true
  fprintf(stderr, "width %.4f\n", scenario_width);
  for (size_t oid = 0; oid < num_obs; oid++) {
    auto& o = obstacles[oid];
    fprintf(stderr, "obstacle %.4f, %.4f; %4f, %.4f\n", o[0][0], o[0][1], o[1][0], o[1][1]);
  }
  fprintf(stderr, "waypoint tol %.4f\n", waypoint_tol);
  for (auto& wp : wp_alternatives) {
    fprintf(stderr, "waypoint %.4f, %.4f\n", wp[0][0], wp[0][1]);
  }
  fprintf(stderr, "t");
#endif

  for (int aid = 0; aid < nent; aid++)
    for (int sample = 0; sample < nsample; sample++)
      evac_time[aid][sample] = max_steps / delta_t;

  int num_active_agents = 0;
  for (int sample = 0; sample < nsample; sample++) {
    num_active_agents = 0;
    num_evac = 0;


    for (int aid = 0; aid < nent; aid++)
      agents[aid].active = false;

#if PRINT_TRACE == true
    if (sample == trace_sample) {
      for (int aid = 0; aid < nent; aid++)
        fprintf(stderr, ",a%d.active,a%d.y,a%d.x", aid, aid, aid);
      fprintf(stderr, "\n");
    }
#endif


    int num_agents_near_wp[num_wp_alternatives];
    vector<int> cand_nbs;

    bool active_agents_left = true;
    for (int step = 0; active_agents_left && step < max_steps; step++) {

      t_sim = step * delta_t;

      for (int wp_offset = 0; wp_offset < num_wp_alternatives; wp_offset++)
        num_agents_near_wp[wp_offset] = 0;

      for (int ego_aid = 0; ego_aid < nent; ego_aid++) {
        auto &ego = agents[ego_aid];
        if (!ego.active)
          continue;
        for (int wp_offset = 0; wp_offset < num_wp_alternatives; wp_offset++) {
          flt2 dist = ego.p - flt2({ wp_alternatives[wp_offset][1][0], wp_alternatives[wp_offset][1][1]});
          if (norm(dist) < congestion_radius)
            num_agents_near_wp[wp_offset]++;
        }
      }

      active_agents_left = false;
      vector<int> grid[numCells][numCells];

      for (int ego = 0; ego < nent; ego++) {
        if (!agents[ego].active)
          continue;
        auto cell = getCellIndex(agents[ego].p);

        if (cell.first < 0 || cell.first >= numCells ||
            cell.second < 0 || cell.second >= numCells) {
              agents[ego].active = false;
              continue;
        }

        grid[cell.first][cell.second].push_back(ego);
      }


      if (step % spawn_period == 0) {
        spawn_agent(sample, num_active_agents, rep);
        num_active_agents++;
      }

      for (int ego_aid = 0; ego_aid < nent; ego_aid++) {
        auto &ego = agents[ego_aid];
        if (!ego.active)
          continue;

        active_agents_left = true;

        flt2 f_internal = {0, 0};

        flt2 wp = ego.waypoint;

        for (int d = 0; d < 2; d++)
          if (wp[d] == -1)
            wp[d] = ego.p[d];

        flt2 t_dist = wp - ego.p;
        flt_t t_dist_norm = norm(t_dist);
        flt2 e = t_dist / t_dist_norm;
        f_internal = ego.v_desired * e - ego.v;

        auto egoCell = getCellIndex(ego.p);

        cand_nbs.clear();

        for (int dx = -1; dx <= 1; dx++) {
          for (int dy = -1; dy <= 1; dy++) {
            if (egoCell.first + dx < 0 || egoCell.first + dx >= numCells ||
                egoCell.second + dy < 0 || egoCell.second + dy >= numCells)
              continue;

            auto& cellList = grid[egoCell.first + dx][egoCell.second + dy];
            for (auto other_aid : cellList) {
               if (ego_aid != other_aid)
                 cand_nbs.push_back(other_aid);
             }
           }
        }

        flt2 f_interaction = {0, 0};
        for (auto other_aid : cand_nbs) {
          auto &other = agents[other_aid];
          if (other_aid == ego_aid || !other.active) continue;

          flt2 o_dist = other.p - ego.p;

          flt_t o_dist_norm = norm(o_dist);

          if (o_dist_norm >= interaction_radius)
            continue;

          flt2 o_dir = o_dist / o_dist_norm;

          flt2 v_diff = ego.v - other.v;

          flt2 int_v = lambda * v_diff + o_dir;

          flt_t int_norm = norm(int_v);
          flt2 int_dir = int_v / int_norm;

          flt_t theta = angle(int_dir, o_dir);
          int theta_sign = theta ? theta / abs(theta) : 0;

          flt_t B = gamma_ * int_norm;

          flt_t npBt = n_prime * B * theta;
          flt_t nBt = n * B * theta;
          flt2 force_v = -exp(-o_dist_norm / B - npBt * npBt) * int_dir;
          flt2 force_angle = -theta_sign * exp(-o_dist_norm / B - nBt * nBt) * left_normal(int_dir);

          f_interaction += force_v + force_angle;
        }

        flt2 f_obstacles = compute_obstacle_force(ego.p);

        flt2 a = f_internal * w_internal +
                 f_interaction * w_interaction +
                 f_obstacles * w_obstacles;

        flt2 p_new = ego.p + ego.v * delta_t + 0.5 * ego.a_old * delta_t * delta_t;
        flt2 v_new = ego.v + 0.5 * (ego.a_old + a) * delta_t;

        ego.p = p_new;
        ego.v = v_new;

        for (int wp_offset = 0; wp_offset < num_wp_alternatives; wp_offset++) {
          flt2 dist = ego.p - flt2({ wp_alternatives[wp_offset][0][0], wp_alternatives[wp_offset][0][1]});
          if (norm(dist) < waypoint_tol) {
            ego.active = false;
            evac_time[ego_aid][sample] = t_sim - ego.t_spawn;
          }
        }
      }


#if PRINT_TRACE == true
      if (sample == trace_sample) {
        fprintf(stderr, "%.6f", t_sim);
        for (auto &ego : agents) {
          fprintf(stderr, ",%d,%.6f,%.6f", ego.active, ego.p[0], ego.p[1]);
        }
        fprintf(stderr, "\n");
      }
#endif
    }
  }

}

int main(int argc, char **argv)
{
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
    if (pert_stddev == -1.0)
      exit(1);

    default_random_engine sim_rng(sim_seed);
    for (int rep = 0; rep < nrep; rep++) {

      for (int aid = 0; aid < nent; aid++) {
        init_agent_states[sample][rep][aid].v_desired = max(min_v_desired, min(max_v_desired, draw_from_mixture(x_perturbed, sim_rng)));
        init_agent_states[sample][rep][aid].p = {0, spawn_p_dist(sim_rng)};
        init_agent_states[sample][rep][aid].wp_offset = wp_dist(sim_rng);
      }
    }

  }

  auto t0 = chrono::high_resolution_clock::now();
  for (int rep = 0; rep < nrep; rep++)
    evac(rep);
  auto t1 = chrono::high_resolution_clock::now();

  flt_t scalar_us = 1e-6 * chrono::duration_cast<chrono::microseconds>(t1 - t0).count();
  printf("scalar took %.4fs\n", scalar_us);

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
