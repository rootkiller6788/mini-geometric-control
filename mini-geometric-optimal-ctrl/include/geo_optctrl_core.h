/**
 * geo_optctrl_core.h ? Geometric Optimal Control: Core Definitions
 *
 * References:
 *   Jurdjevic (1997) "Geometric Control Theory"
 *   Agrachev & Sachkov (2004) "Control Theory from the Geometric Viewpoint"
 *   Bloch (2003) "Nonholonomic Mechanics and Control"
 *
 * This header defines the foundational types for geometric optimal control:
 *   - Smooth manifolds via atlas/chart representation
 *   - Tangent and cotangent bundles
 *   - Vector fields with Lie bracket operations
 *   - Control-affine systems on manifolds
 *   - Optimal control problem formulation
 *
 * Knowledge coverage: L1 Definitions, L3 Mathematical Structures
 */

#ifndef GEO_OPTCTRL_CORE_H
#define GEO_OPTCTRL_CORE_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>

/* -----------------------------------------------------------------------------
 * L1: Manifold dimension and chart definitions
 * ----------------------------------------------------------------------------- */

#define GEO_MAX_DIM 16
#define GEO_MAX_CHARTS 32
#define GEO_MAX_CONTROLS 8

typedef struct {
    int    dim;
    double center[GEO_MAX_DIM];
    double radius;
    int    chart_id;
    int    overlap_mask[GEO_MAX_CHARTS];
} GeoChart;

typedef struct {
    int      num_charts;
    int      manifold_dim;
    GeoChart charts[GEO_MAX_CHARTS];
} GeoAtlas;

/* -----------------------------------------------------------------------------
 * L1: Tangent and cotangent vector types
 * ----------------------------------------------------------------------------- */

typedef struct {
    double comp[GEO_MAX_DIM];
    int    dim;
    int    chart_id;
} GeoTangentVector;

typedef struct {
    double comp[GEO_MAX_DIM];
    int    dim;
    int    chart_id;
} GeoCotangentVector;

/* -----------------------------------------------------------------------------
 * L2: Vector field and Lie bracket
 * ----------------------------------------------------------------------------- */

typedef struct {
    int    dim;
    int    chart_id;
    void (*eval)(const double *x, double *result, void *ctx);
    void *ctx;
} GeoVectorField;

void geo_lie_bracket(const GeoVectorField *X, const GeoVectorField *Y,
                     const double *at_point, double *result,
                     int dim, double h);

void geo_iterated_lie_bracket(const GeoVectorField *X,
                              const GeoVectorField *Y,
                              int k, const double *at_point,
                              double *result, int dim, double h);

/* -----------------------------------------------------------------------------
 * L1: Control-affine system on a manifold
 * ----------------------------------------------------------------------------- */

typedef struct {
    int             state_dim;
    int             control_dim;
    GeoVectorField  drift;
    GeoVectorField  controls[GEO_MAX_CONTROLS];
    GeoAtlas       *atlas;
} GeoControlAffineSystem;

/* -----------------------------------------------------------------------------
 * L1: Optimal control problem (Bolza form)
 * ----------------------------------------------------------------------------- */

typedef double (*GeoRunningCost)(const double *x, const double *u,
                                  int state_dim, int control_dim, void *ctx);

typedef double (*GeoTerminalCost)(const double *x, int state_dim, void *ctx);

typedef struct {
    GeoControlAffineSystem *system;
    GeoRunningCost          running_cost;
    GeoTerminalCost         terminal_cost;
    void                   *cost_ctx;
    double                 *x0;
    double                 *xT;
    double                  T;
    int                     fixed_final_time;
    double                  control_box[GEO_MAX_CONTROLS][2];
} GeoOptimalControlProblem;

/* -----------------------------------------------------------------------------
 * L3: Product manifold
 * ----------------------------------------------------------------------------- */

typedef struct {
    int dim1;
    int dim2;
    GeoAtlas *atlas1;
    GeoAtlas *atlas2;
} GeoProductManifold;

void geo_product_embed(const double *x1, int dim1,
                       const double *x2, int dim2,
                       double *result);

void geo_product_project(const double *x, int dim_total,
                         double *x1, int dim1,
                         double *x2, int dim2);

/* -----------------------------------------------------------------------------
 * L2: Flow of a vector field (numerical)
 * ----------------------------------------------------------------------------- */

void geo_flow(const GeoVectorField *X, const double *x0, double t,
              int n_steps, double *result, int dim);

void geo_flow_jacobian(const GeoVectorField *X, const double *x0,
                       double t, int n_steps, double *jac,
                       int dim, double h);

/* -----------------------------------------------------------------------------
 * L3: Differential of a smooth function
 * ----------------------------------------------------------------------------- */

typedef double (*GeoScalarFunction)(const double *x, int dim, void *ctx);

void geo_differential(GeoScalarFunction f, void *ctx,
                      const double *at_p, double *result,
                      int dim, double h);

#endif /* GEO_OPTCTRL_CORE_H */
