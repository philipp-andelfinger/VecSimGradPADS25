#pragma once

using namespace std;

typedef FLT_T flt_t;

constexpr int nent = NENT;
constexpr int nrep = NREP;
constexpr int nsample = NSAMPLE;

const int min_nb_samples = MIN_NB_SAMPLES;

const int num_in_comps = 2;
const int num_out_comps = 5;

const int num_params = num_in_comps * 3; // prob, mean, stddev

constexpr int32_t cceil(flt_t num) { return (int)num == num ?  num : (int)num + 1; }

constexpr flt_t scenario_width = 30;
constexpr flt_t interaction_radius = 5.0;
const flt_t delta_t = 0.1;

constexpr flt_t cellSize = interaction_radius;
constexpr int numCells = cceil((flt_t)scenario_width / cellSize);

const flt_t v_desired_mean = 1.29;
const flt_t v_desired_stddev = 0.5;
const flt_t min_v_desired = 0.5;
const flt_t max_v_desired = 2.5;

const flt_t min_v_desired_mean = 1.0;
const flt_t max_v_desired_mean = 2.0;

const flt_t lambda = 2.0;
const flt_t gamma_ = 0.35;
const flt_t n = 2;
const flt_t n_prime = 3;

const flt_t sigma = 0.8;

const flt_t w_internal = 1.0;
const flt_t w_interaction = 45;
const flt_t w_obstacles = 5;

const flt_t waypoint_tol = 2.0; // in this radius the wp is considered to be reached
const flt_t congestion_radius = 5; // agents in this radius are counted as part of the congestion around the waypoint
const flt_t agent_radius = 0.4;

const int max_steps = 120 / delta_t;

const int spawn_period = SPAWN_PERIOD / delta_t;

const flt_t door_width_m = 5;
const flt_t door_width = door_width_m / scenario_width;
const flt_t door_offset_front = 0.1;
const flt_t door_offset_sides = 0.6;
const flt_t wall_margin_y = 0.15;
const flt_t wall_margin_x = 0.15;
const flt_t door_wp_offset_m = 4;
const flt_t door_wp_offset = door_wp_offset_m / scenario_width;


constexpr size_t num_obs = 7;
flt_t obstacles[num_obs][2][2] = {
  {{(1 - wall_margin_y) * scenario_width, wall_margin_x * scenario_width}, {(1 - wall_margin_y) * scenario_width, (wall_margin_x + door_offset_front) * scenario_width}},
  {{(1 - wall_margin_y) * scenario_width, (wall_margin_x + door_offset_front + door_width) * scenario_width}, {(1 - wall_margin_y) * scenario_width, (1 - wall_margin_x - door_offset_front - door_width) * scenario_width}},
  {{(1 - wall_margin_y) * scenario_width, (1 - wall_margin_x - door_offset_front) * scenario_width}, {(1 - wall_margin_y) * scenario_width, (1 - wall_margin_x) * scenario_width}},

  {{0, wall_margin_x * scenario_width}, {(1 - wall_margin_y - door_offset_sides) * scenario_width, wall_margin_x * scenario_width}},
  {{(1 - wall_margin_y - door_offset_sides + door_width) * scenario_width, wall_margin_x * scenario_width}, {(1 - wall_margin_y) * scenario_width, wall_margin_x * scenario_width}},

  {{0, (1 - wall_margin_x) * scenario_width}, {(1 - wall_margin_y - door_offset_sides) * scenario_width, (1 - wall_margin_x) * scenario_width}},
  {{(1 - wall_margin_y - door_offset_sides + door_width) * scenario_width, (1 - wall_margin_x) * scenario_width}, {(1 - wall_margin_y) * scenario_width, (1 - wall_margin_x) * scenario_width}}
};

// y, x coordinates of waypoints and the expected center of congestion
constexpr size_t num_wp_alternatives = 4;
flt_t wp_alternatives[num_wp_alternatives][2][2] = {
  {{(1 - wall_margin_y + door_wp_offset) * scenario_width, (wall_margin_x + door_offset_front + door_width / 2) * scenario_width},
   {(1 - wall_margin_y - door_wp_offset) * scenario_width, (wall_margin_x + door_offset_front + door_width / 2) * scenario_width}},

  {{(1 - wall_margin_y + door_wp_offset) * scenario_width, (1 - wall_margin_x - door_offset_front - door_width / 2) * scenario_width},
   {(1 - wall_margin_y - door_wp_offset) * scenario_width, (1 - wall_margin_x - door_offset_front - door_width / 2) * scenario_width}},

  {{(1 - wall_margin_y - door_offset_sides + door_width / 2) * scenario_width, (wall_margin_x - door_wp_offset) * scenario_width},
   {(1 - wall_margin_y - door_offset_sides + door_width / 2) * scenario_width, (wall_margin_x + door_wp_offset) * scenario_width}},

  {{(1 - wall_margin_y - door_offset_sides + door_width / 2) * scenario_width, (1 - wall_margin_x + door_wp_offset) * scenario_width},
   {(1 - wall_margin_y - door_offset_sides + door_width / 2) * scenario_width, (1 - wall_margin_x - door_wp_offset) * scenario_width}}
};

void project_param(const flt_t x[3], const flt_t perts[3], flt_t x_perturbed[3], flt_t pert_stddev) { 
  x_perturbed[0] = 0.1 + 0.4 * (sin(x[0] + perts[0] * pert_stddev) + 1);
  x_perturbed[1] = min_v_desired_mean + (max_v_desired_mean - min_v_desired_mean) * 0.5 * (sin(x[1] + perts[1] * pert_stddev) + 1);
  x_perturbed[2] = 0.05 + 0.45 * 0.5 * (sin(x[2] + perts[2] * pert_stddev) + 1);
}

flt_t evac_time_ref[num_out_comps][3] = {
{ 0.3037560814867312, 13.292637935543487, 1.5973560091651458 },
{ 0.19951440954057273, 21.362474519814157, 1.842736985564384 },
{ 0.3480584155556005, 17.39034137814955, 1.5614024025011202 },
{ 0.11065641049868727, 26.393372726579265, 2.5319691700070313 },
{ 0.03801468291840818, 32.90404384510026, 5.204967488731871 }
};
