#ifndef NOPLAN_PLANNING_H
#define NOPLAN_PLANNING_H

#include "noplan_core.h"
#include "noplan_geometry.h"

/* ============================================================================
 * Nonholonomic Motion Planning — Planning Algorithms
 *
 * L5: Algorithms for nonholonomic planning
 *
 *   Core principle: Nonholonomic constraints make the connection between
 *   two configurations non-trivial. Standard holonomic planners (straight-line
 *   interpolation) produce inadmissible paths. Special algorithms are needed
 *   that respect the constraints.
 *
 *   Algorithm classes:
 *     - Sampling-based: RRT, PRM (with nonholonomic steering)
 *     - Geometric: Lie bracket steering, sinusoidal inputs, chained form
 *     - Optimal: Shooting methods, direct transcription
 *     - Search-based: Lattice planners, state-space discretization
 *
 * Key references:
 *   LaValle & Kuffner — "RRT-Connect" (2000)
 *   LaValle — "Planning Algorithms" (2006)
 *   Murray & Sastry — "Nonholonomic Motion Planning: Steering using Sinusoids" (1993)
 *   Lafferriere & Sussmann — "Motion Planning for Controllable Systems Without Drift" (1993)
 *   Choset et al. — "Principles of Robot Motion" (2005)
 * ============================================================================ */

/* ---------------------------------------------------------------------------
 * L5: Rapidly-exploring Random Tree (RRT) for Nonholonomic Systems
 * --------------------------------------------------------------------------- */

/**
 * RRT node: a configuration in the search tree with parent linkage
 * and the control used to reach it from the parent.
 */
typedef struct RRTNode {
    Config q;                         /* Configuration at this node */
    int parent_idx;                   /* Index of parent node (-1 for root) */
    double u_from_parent[NOPLAN_MAX_INPUTS]; /* Control from parent to this */
    double cost_from_root;            /* Accumulated cost */
    int depth;                        /* Tree depth */
} RRTNode;

/**
 * RRT planner state.
 * Maintains a tree rooted at q_start exploring the configuration space
 * using randomized control inputs applied over random durations.
 *
 * The EXTEND operation:
 *   1. Sample random configuration q_rand
 *   2. Find nearest node in tree q_near
 *   3. Apply random controls from q_near toward q_rand
 *   4. Add new node if the resulting q_new is collision-free
 */
typedef struct {
    RRTNode* nodes;                   /* Tree nodes */
    int n_nodes;                      /* Number of nodes */
    int capacity;                     /* Allocated node capacity */
    Config q_start;                   /* Root configuration */
    Config q_goal;                    /* Target configuration */
    NonholonomicSystem* system;       /* The control system */
    PlanningProblem* problem;         /* Full planning problem */
    double step_duration;             /* Control application duration per step */
    double goal_bias;                 /* Probability of sampling q_goal */
    double goal_tolerance;            /* Distance threshold for "reached goal" */
    int max_iterations;               /* Maximum extension attempts */
    bool goal_reached;                /* Whether goal has been reached */
    int goal_node_idx;                /* Index of node reaching goal */
} RRTPlanner;

/* RRT API */
RRTPlanner* noplan_rrt_create(NonholonomicSystem* sys,
                               PlanningProblem* prob,
                               double step_duration);
void    noplan_rrt_free(RRTPlanner* rrt);
int     noplan_rrt_extend(RRTPlanner* rrt);
bool    noplan_rrt_plan(RRTPlanner* rrt, int max_iters, Trajectory* traj);
int     noplan_rrt_nearest(const RRTPlanner* rrt, const Config* q);
void    noplan_rrt_extract_path(const RRTPlanner* rrt, int goal_idx,
                                 Trajectory* traj);
void    noplan_rrt_set_goal_bias(RRTPlanner* rrt, double bias);
bool    noplan_rrt_connect(RRTPlanner* rrt_a, RRTPlanner* rrt_b,
                            Trajectory* traj);

/* ---------------------------------------------------------------------------
 * L5: Probabilistic Roadmap (PRM) for Nonholonomic Systems
 * --------------------------------------------------------------------------- */

/**
 * PRM node in the roadmap graph.
 */
typedef struct PRMNode {
    Config q;                         /* Configuration */
    int* neighbors;                   /* Indices of connected neighbors */
    double* edge_costs;               /* Cost to each neighbor */
    int n_neighbors;                  /* Number of connections */
    int capacity;                     /* Allocated neighbor capacity */
} PRMNode;

/**
 * PRM planner.
 * Builds a graph in configuration space where edges are
 * nonholonomically feasible trajectories (computed by a steering method).
 *
 * Two phases:
 *   Construction:  Sample N configurations, connect k-nearest with steering.
 *   Query:         Find shortest path through roadmap from start to goal.
 */
typedef struct {
    PRMNode* nodes;                   /* Roadmap nodes */
    int n_nodes;                      /* Number of nodes */
    int capacity;                     /* Allocated capacity */
    NonholonomicSystem* system;       /* Control system */
    PlanningProblem* problem;         /* Planning problem */
    int k_nearest;                    /* Number of nearest neighbors to connect */
    double connection_radius;         /* Max distance for connection attempt */
    int max_steering_iters;           /* Max attempts to connect two nodes */
    bool is_built;                    /* Whether construction phase is done */
} PRMPlanner;

/* PRM API */
PRMPlanner* noplan_prm_create(NonholonomicSystem* sys,
                               PlanningProblem* prob,
                               int k_nearest, double connection_radius);
void    noplan_prm_free(PRMPlanner* prm);
void    noplan_prm_add_node(PRMPlanner* prm, const Config* q);
void    noplan_prm_connect_node(PRMPlanner* prm, int node_idx);
void    noplan_prm_build(PRMPlanner* prm, int n_samples);
bool    noplan_prm_query(const PRMPlanner* prm,
                          const Config* q_start,
                          const Config* q_goal,
                          Trajectory* traj);
void    noplan_prm_dijkstra(const PRMPlanner* prm,
                             int start_idx, int goal_idx,
                             int* path, int* path_len);

/* ---------------------------------------------------------------------------
 * L5: Sinusoidal Steering (Murray & Sastry, 1993)
 * --------------------------------------------------------------------------- */

/**
 * For systems in chained form, sinusoidal inputs at integrally
 * related frequencies produce net motion in the Lie bracket directions.
 *
 * Chained form (n-dimensional, 2-input):
 *   ż₁ = u₁
 *   ż₂ = u₂
 *   ż₃ = z₂ · u₁
 *   ż₄ = z₃ · u₁
 *   ...
 *   ż_n = z_{n-1} · u₁
 *
 * Steering strategy:
 *   - Use piecewise-constant u₁ to set the "phase"
 *   - Use sinusoids in u₂ to generate motion in the chained variables
 *   - Frequencies: u₂(t) = Σ a_i sin(2π·i·u₁·t + φ_i)
 *
 * This algorithm implements the Murray-Sastry steering method that
 * exploits the Lie bracket structure to achieve net displacement
 * in all n directions despite having only 2 inputs.
 */

/**
 * Sinusoidal steering result.
 */
typedef struct {
    double* u1_seq;                   /* Sequence of u₁ values */
    double* u2_seq;                   /* Sequence of u₂ values */
    double* times;                    /* Time points */
    int n_steps;                      /* Number of control steps */
    double duration;                  /* Total steering time */
    bool converged;                   /* Whether steering succeeded */
    double final_error;               /* ||q_final - q_goal|| */
} SinusoidalSteering;

SinusoidalSteering* noplan_steer_sinusoidal(
    const NonholonomicSystem* sys,
    const Config* q_start, const Config* q_goal,
    int n_frequencies, double amplitude, double duration);

void    noplan_steer_sinusoidal_free(SinusoidalSteering* ss);

/**
 * Apply sinusoidal steering and return the resulting trajectory.
 */
void    noplan_steer_sinusoidal_to_trajectory(
    const SinusoidalSteering* ss,
    const NonholonomicSystem* sys,
    Trajectory* traj);

/* ---------------------------------------------------------------------------
 * L5: Nilpotent Steering (Lafferriere & Sussmann, 1993)
 * --------------------------------------------------------------------------- */

/**
 * Using the nilpotent approximation of the system, solve the
 * motion planning problem by concatenating primitives that
 * produce motion along each bracket direction.
 *
 * The method:
 *   1. Compute nilpotent approximation at q_start
 *   2. Solve for the "Philip Hall" coordinates that produce
 *      the desired displacement
 *   3. Translate Philip Hall coordinates to piecewise-constant controls
 *      using the Campbell-Baker-Hausdorff formula
 *   4. Refine using Newton iteration on the original system
 */

/**
 * Philip Hall basis element for nilpotent Lie algebra.
 */
typedef struct {
    char label[256];                  /* Symbolic label, e.g., "[[g1,g2],g1]" */
    int degree;                       /* Homogeneous degree (weight) */
    TangentVector direction;          /* Vector field direction */
} PhilipHallElement;

/**
 * Nilpotent steering plan.
 */
typedef struct {
    int n_segments;                   /* Number of control segments */
    double** controls;                /* controls[seg][input] */
    double* durations;                /* Duration of each segment */
    int n_inputs;                     /* Number of control inputs */
} NilpotentSteering;

NilpotentSteering* noplan_steer_nilpotent(
    const NonholonomicSystem* sys,
    const Config* q_start, const Config* q_goal,
    int nilpotent_order, double eps);

void    noplan_steer_nilpotent_free(NilpotentSteering* ns);

void    noplan_steer_nilpotent_to_trajectory(
    const NilpotentSteering* ns,
    const NonholonomicSystem* sys,
    Trajectory* traj);

/* ---------------------------------------------------------------------------
 * L5: General Steering Functions
 * --------------------------------------------------------------------------- */

/**
 * Attempt to steer from q_from to q_to using random control sampling.
 * This is the "local planner" used inside RRT and PRM.
 *
 * Returns a short trajectory if successful, NULL if steering fails.
 */
Trajectory* noplan_steer_random(const NonholonomicSystem* sys,
                                 const Config* q_from,
                                 const Config* q_to,
                                 int max_attempts,
                                 double step_duration,
                                 double* best_distance);

/**
 * General steering dispatcher: choose the best steering method
 * based on the system type (chained form, nilpotent, or generic).
 */
Trajectory* noplan_steer_general(const NonholonomicSystem* sys,
                                  const Config* q_from,
                                  const Config* q_to,
                                  double eps);

/* ---------------------------------------------------------------------------
 * L5: Path Optimization
 * --------------------------------------------------------------------------- */

/**
 * Shortcut smoothing: attempt to replace sub-paths with shorter
 * direct connections. This is the classic "shortcut" heuristic
 * that improves RRT/PRM paths.
 */
void    noplan_path_shortcut(Trajectory* traj,
                              const NonholonomicSystem* sys,
                              const PlanningProblem* prob,
                              int n_iterations);

/**
 * Gradient-based trajectory optimization for nonholonomic systems.
 * Minimize a cost functional (e.g., path length + control effort)
 * subject to the nonholonomic constraints using a shooting method.
 *
 * The cost function: J = Σ (||q_k − q_{k-1}||² + λ·||u_k||²)
 * with constraints q_{k+1} = q_k + Σ g_i(q_k)·u_{k,i}·dt.
 */
void    noplan_trajectory_optimize(Trajectory* traj,
                                    const NonholonomicSystem* sys,
                                    int max_iters, double step_size);

/**
 * Reeds-Shepp optimal path length between two SE(2) configurations.
 * The Reeds-Shepp car can move forward and backward with minimum
 * turning radius.
 */
double  noplan_reeds_shepp_length(const Config* q_start,
                                   const Config* q_goal,
                                   double turning_radius);

/**
 * Compute Dubins path length (forward-only car).
 */
double  noplan_dubins_length(const Config* q_start,
                              const Config* q_goal,
                              double turning_radius);

/* ---------------------------------------------------------------------------
 * L5: Lattice-based Planning
 * --------------------------------------------------------------------------- */

/**
 * State lattice: a regular discretization of the configuration space
 * with pre-computed motion primitives connecting lattice points.
 *
 * Motion primitives are short, kinematically feasible trajectories
 * that connect one lattice state to another.
 */
typedef struct {
    Config* states;                   /* Lattice states */
    int n_states;                     /* Number of states */
    int* dimensions;                  /* Grid dimensions per axis */
    int n_dims;                       /* Number of discretized dimensions */
    Trajectory** primitives;          /* Motion primitives from each state */
    int* n_primitives;                /* Number of primitives per state */
    NonholonomicSystem* system;       /* The control system */
} StateLattice;

StateLattice* noplan_lattice_create(const NonholonomicSystem* sys,
                                     const double* bounds_min,
                                     const double* bounds_max,
                                     const int* grid_sizes,
                                     int n_dims);

void    noplan_lattice_free(StateLattice* lattice);

void    noplan_lattice_generate_primitives(StateLattice* lattice,
                                            double primitive_duration,
                                            int n_primitives_per_state);

bool    noplan_lattice_plan(const StateLattice* lattice,
                             const Config* q_start,
                             const Config* q_goal,
                             Trajectory* traj);

/* ---------------------------------------------------------------------------
 * Utility: Nearest Neighbor Search (k-d tree)
 * --------------------------------------------------------------------------- */

/**
 * Simple k-d tree for efficient nearest neighbor queries in
 * configuration space. Used by RRT and PRM.
 */
typedef struct KDNode {
    Config q;                         /* Stored configuration */
    int data_idx;                     /* Index in the original array */
    int split_dim;                    /* Dimension used for split */
    struct KDNode* left;
    struct KDNode* right;
} KDNode;

typedef struct {
    KDNode* root;
    int n_points;
    int dim;
    double (*distance_func)(const Config*, const Config*);
} KDTree;

KDTree* noplan_kdtree_create(int dim,
                              double (*dist_func)(const Config*,
                                                  const Config*));
void    noplan_kdtree_free(KDTree* tree);
void    noplan_kdtree_insert(KDTree* tree, const Config* q, int idx);
int     noplan_kdtree_nearest(const KDTree* tree, const Config* q,
                               double* best_dist);
void    noplan_kdtree_k_nearest(const KDTree* tree, const Config* q,
                                  int k, int* indices, double* dists);

/* ---------------------------------------------------------------------------
 * L5: Deterministic sampling sequences (Halton, Sobol)
 * --------------------------------------------------------------------------- */

/**
 * Generate a Halton sequence point in [0, 1]^dim.
 * Halton sequences provide low-discrepancy sampling for PRM construction.
 */
void    noplan_halton_sample(int index, int dim, double* point);

/**
 * Generate uniformly distributed random configuration within bounds.
 */
void    noplan_random_config(const Config* bound_min,
                              const Config* bound_max,
                              Config* q_rand);

#endif /* NOPLAN_PLANNING_H */
