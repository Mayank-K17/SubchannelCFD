/*============================================================================
 * SubChCFD (C++ port for code_saturne 9.x): code_saturne/code_saturne volume
 * coupling source terms.
 * Interface: subchcfd_coupling.h.
 *
 * Replaces the case's csc2ts_coarse.f90 together with the modified
 * csccel.f90, predvv.f90 and covofi.f90 (code_saturne 8.0). Remove those
 * from SRC/.
 *
 * Behavioural changes relative to the Fortran are marked "FIX(port)".
 *============================================================================*/
 
#include "subchcfd_coupling.h"
 
#include "subchcfd.h"
 
#include <ple_locator.h>
 
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
 
namespace subchcfd {
 
/*============================================================================
 * Access to the coupling locators
 *
 * cs_sat_coupling_t is opaque in 9.x. This mirrors the leading members of
 * struct _cs_sat_coupling_t in code_saturne 9.2.0 (src/base/cs_sat_coupling.cpp),
 * up to the two locators; a common initial sequence of identical members
 * has the same layout. Re-check against the source of any other release.
 *============================================================================*/
 
namespace {
 
struct sat_coupling_head_t {
  char                   *sat_name;
  cs_sat_coupling_tag_t  *tag_func;
  void                   *tag_context;
  int                     reverse;
  char                   *face_cpl_sel;
  char                   *cell_cpl_sel;
  char                   *face_loc_sel;
  char                   *cell_loc_sel;
  ple_locator_t          *localis_cel;
  ple_locator_t          *localis_fbr;
};
 
inline const sat_coupling_head_t *
cpl_head(int cpl_id)
{
  cs_sat_coupling_t *cpl = cs_sat_coupling_by_id(cpl_id);
  if (cpl == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "SubChCFD: code_saturne coupling %d is not defined.", cpl_id);
  return reinterpret_cast<const sat_coupling_head_t *>(cpl);
}
 
/* Fortran single-precision literal promoted to double */
constexpr double big_value = std::numeric_limits<double>::max();
 
/* Distant pressure gradient received at the last exchange, per coupling
 * (kept between exchange steps, as the Fortran rdpcpce) */
std::map<int, std::vector<cs_real_t>> &
last_dist_dp()
{
  static std::map<int, std::vector<cs_real_t>> m;
  return m;
}
 
/* Bounding boxes for mom_crr == 2, computed on the first time step of the
 * run and kept (Fortran "save" variables). FIX(port): stored per coupling;
 * the Fortran kept one set, overwritten by each coupling. */
struct boxes_t {
  double cen_min[3], cen_max[3], fb_min[3], fb_max[3];
};
 
std::map<int, boxes_t> &
saved_boxes()
{
  static std::map<int, boxes_t> m;
  return m;
}
 
/*----------------------------------------------------------------------------
 * Copy of code_saturne 9.2.0 _sat_coupling_compute_data_at_cells
 * (static in cs_sat_coupling.cpp)
 *----------------------------------------------------------------------------*/
 
void
compute_data_at_cells(cs_field_t        *f,
                      const int          reverse,
                      const cs_lnum_t    n_elts,
                      const cs_lnum_t   *elt_ids,
                      const cs_coord_t  *coords,
                      cs_real_t         *warray1,
                      cs_real_t         *warray2)
{
  if (reverse == 0) {
    const cs_real_3_t *cell_cen = cs_glob_mesh_quantities->cell_cen;
 
    if (f->dim == 1) {
      cs_real_3_t *grad = nullptr;
      CS_MALLOC(grad, cs_glob_mesh->n_cells_with_ghosts, cs_real_3_t);
      cs_field_gradient_scalar(f, true, 1, grad);
 
      for (cs_lnum_t e_id = 0; e_id < n_elts; e_id++) {
        cs_lnum_t c_id = elt_ids[e_id];
        cs_real_t dxyz[3];
        for (int i = 0; i < 3; i++)
          dxyz[i] = coords[3*e_id + i] - cell_cen[c_id][i];
        warray1[e_id] = f->val_pre[c_id] + cs_math_3_dot_product(dxyz, grad[c_id]);
      }
      CS_FREE(grad);
    }
    else if (f->dim == 3) {
      cs_real_33_t *grad = nullptr;
      CS_MALLOC(grad, cs_glob_mesh->n_cells_with_ghosts, cs_real_33_t);
      cs_field_gradient_vector(f, true, 1, grad);
 
      for (cs_lnum_t e_id = 0; e_id < n_elts; e_id++) {
        cs_lnum_t c_id = elt_ids[e_id];
        cs_real_t dxyz[3];
        for (int i = 0; i < 3; i++)
          dxyz[i] = coords[3*e_id + i] - cell_cen[c_id][i];
        for (int j = 0; j < 3; j++)
          warray1[3*e_id + j] =   f->val_pre[3*c_id + j]
                                + cs_math_3_dot_product(dxyz, grad[c_id][j]);
      }
      CS_FREE(grad);
    }
    else
      bft_error(__FILE__, __LINE__, 0,
                "%s: Exchange is possible only for scalars and vectors.\n",
                __func__);
  }
  else {
    const cs_real_t *cvar_rho = CS_F_(rho)->val;
    const cs_real_t *fluid_vol = cs_glob_mesh_quantities->cell_vol;
 
    if (f->dim == 1) {
      for (cs_lnum_t e_id = 0; e_id < n_elts; e_id++) {
        cs_lnum_t c_id = elt_ids[e_id];
        warray1[e_id] = fluid_vol[c_id]*cvar_rho[c_id]*f->val_pre[c_id];
        warray2[e_id] = fluid_vol[c_id]*cvar_rho[c_id];
      }
    }
    else if (f->dim == 3) {
      for (cs_lnum_t e_id = 0; e_id < n_elts; e_id++) {
        cs_lnum_t c_id = elt_ids[e_id];
        for (int i = 0; i < 3; i++) {
          warray1[3*e_id + i] = fluid_vol[c_id]*cvar_rho[c_id]
                                *f->val_pre[3*c_id + i];
          for (int j = 0; j < 3; j++)
            warray2[9*e_id + 3*i + j]
              = (i == j) ? fluid_vol[c_id]*cvar_rho[c_id] : 0.;
        }
      }
    }
    else
      bft_error(__FILE__, __LINE__, 0,
                "%s: Exchange is possible only for scalars and vectors.\n",
                __func__);
  }
}
 
/*----------------------------------------------------------------------------
 * Stock code_saturne 9.2.0 source terms (used when cfg::cpl is false)
 *----------------------------------------------------------------------------*/
 
void
stock_source_terms(const cs_field_t  *f,
                   int                reverse,
                   cs_lnum_t          n_cells_loc,
                   cs_lnum_t          n_g_cells_loc,
                   const cs_lnum_t   *cpl_cells_ids,
                   cs_lnum_t          n_cells_dist,
                   cs_lnum_t          n_g_cells_dist,
                   const cs_lnum_t   *dist_elt_ids,
                   const cs_real_t   *cw1_loc,
                   const cs_real_t   *cw1_dis,
                   const cs_real_t   *cw2_dis,
                   cs_real_t         *rhs,
                   cs_real_t         *fimp)
{
  const int dim = f->dim;
  const int dim2 = dim*dim;
  const cs_real_t *cvar_rho = CS_F_(rho)->val;
  const cs_real_t one_ov_xtau = 1./(100.0*cs_glob_time_step->dt_ref);
  const cs_real_t *cell_f_vol = cs_glob_mesh_quantities->cell_vol;
 
  if (reverse == 0 && n_g_cells_loc > 0) {
    for (cs_lnum_t e_id = 0; e_id < n_cells_loc; e_id++) {
      cs_lnum_t c_id = cpl_cells_ids[e_id];
      cs_real_t rovtau = cell_f_vol[c_id]*cvar_rho[c_id]*one_ov_xtau;
      for (int i = 0; i < dim; i++) {
        rhs[dim*c_id + i] += rovtau*cw1_loc[dim*e_id + i];
        fimp[dim2*c_id + dim*i + i] -= rovtau;
      }
    }
  }
 
  if (reverse == 1 && n_g_cells_dist > 0) {
    for (cs_lnum_t e_id = 0; e_id < n_cells_dist; e_id++) {
      cs_lnum_t c_id = dist_elt_ids[e_id];
      for (int i = 0; i < dim; i++) {
        rhs[dim*c_id + i] += cw1_dis[dim*e_id + i]*one_ov_xtau;
        for (int j = 0; j < dim; j++)
          fimp[dim2*c_id + dim*i + j] -= cw2_dis[dim2*e_id + dim*i + j]*one_ov_xtau;
      }
    }
  }
}
 
} // anonymous namespace
 
/*============================================================================
 * Coupling region bounding boxes
 *============================================================================*/
 
void
coupling_xyz(int     cpl_id,
             double  xyzcenmin[3],
             double  xyzcenmax[3],
             double  xyzfbmin[3],
             double  xyzfbmax[3])
{
  const sat_coupling_head_t *cpl = cpl_head(cpl_id);
  const cs_real_3_t *cell_cen = cs_glob_mesh_quantities_g->cell_cen;
 
  for (int k = 0; k < 3; k++) {
    xyzcenmin[k] =  big_value;  xyzcenmax[k] = -big_value;
    xyzfbmin[k]  =  big_value;  xyzfbmax[k]  = -big_value;
  }
 
  /* Local coupled cells (Fortran lelcpl) */
  if (cpl->localis_cel != nullptr) {
    const cs_lnum_t n = ple_locator_get_n_interior(cpl->localis_cel);
    const cs_lnum_t *ids = ple_locator_get_interior_list(cpl->localis_cel);
    for (cs_lnum_t i = 0; i < n; i++)
      for (int k = 0; k < 3; k++) {
        xyzcenmin[k] = std::min(xyzcenmin[k], cell_cen[ids[i]][k]);
        xyzcenmax[k] = std::max(xyzcenmax[k], cell_cen[ids[i]][k]);
      }
  }
 
  /* Distant coupled boundary faces located in the local mesh
   * (Fortran coocpl with ityvar = 2) */
  if (cpl->localis_fbr != nullptr) {
    const cs_lnum_t n = ple_locator_get_n_dist_points(cpl->localis_fbr);
    const ple_coord_t *xyz = ple_locator_get_dist_coords(cpl->localis_fbr);
    for (cs_lnum_t i = 0; i < n; i++)
      for (int k = 0; k < 3; k++) {
        xyzfbmin[k] = std::min(xyzfbmin[k], static_cast<double>(xyz[3*i + k]));
        xyzfbmax[k] = std::max(xyzfbmax[k], static_cast<double>(xyz[3*i + k]));
      }
  }
 
  cs_parall_min(3, CS_DOUBLE, xyzcenmin);
  cs_parall_max(3, CS_DOUBLE, xyzcenmax);
  cs_parall_min(3, CS_DOUBLE, xyzfbmin);
  cs_parall_max(3, CS_DOUBLE, xyzfbmax);
}
 
/*============================================================================
 * Pressure gradient exchange
 *============================================================================*/
 
bool
pressure_gradient_exchange_step()
{
  return    cfg::receive_dp
         && cs_glob_time_step->nt_cur % cfg::freceive_dp == 0;
}
 
void
exchange_pressure_gradient(int                     cpl_id,
                           std::vector<cs_real_t> &received)
{
  const sat_coupling_head_t *cpl = cpl_head(cpl_id);
  ple_locator_t *loc = cpl->localis_cel;
  if (loc == nullptr)
    return;
 
  const cs_lnum_t n_loc  = ple_locator_get_n_interior(loc);
  const cs_lnum_t n_dist = ple_locator_get_n_dist_points(loc);
  const cs_lnum_t *dist_ids = ple_locator_get_dist_locations(loc);
 
  cs_lnum_t n_g[2] = {n_loc, n_dist};
  cs_parall_sum(2, CS_LNUM_TYPE, n_g);
  if (n_g[0] == 0 && n_g[1] == 0)
    return;
 
  /* Local grad(p) at the cells containing the distant points */
  std::vector<cs_real_t> send(3*static_cast<size_t>(n_dist), 0.0);
  {
    cs_real_3_t *grad = nullptr;
    CS_MALLOC(grad, cs_glob_mesh->n_cells_with_ghosts, cs_real_3_t);
    cs_field_gradient_scalar(CS_F_(p), false, 1, grad);
    for (cs_lnum_t e = 0; e < n_dist; e++)
      for (int k = 0; k < 3; k++)
        send[3*e + k] = grad[dist_ids[e]][k];
    CS_FREE(grad);
  }
 
  received.assign(3*static_cast<size_t>(n_loc), 0.0);
 
  ple_locator_exchange_point_var(loc,
                                 (n_dist > 0) ? send.data() : nullptr,
                                 (n_loc > 0) ? received.data() : nullptr,
                                 nullptr,
                                 sizeof(cs_real_t),
                                 3,
                                 0);
}
 
/*============================================================================
 * SubChCFD coupling source terms (port of csc2ts_coarse)
 *============================================================================*/
 
void
coupling_source_terms(const cs_field_t  *f,
                      int                cpl_id,
                      cs_lnum_t          n_cpl_cells,
                      const cs_lnum_t    cpl_cells[],
                      const cs_real_t    dist_val[],
                      const cs_real_t    dist_dp[],
                      cs_real_t          rhs[],
                      cs_real_t          fimp[])
{
  if (!cfg::cpl)
    return;
 
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_real_t *cell_f_vol = cs_glob_mesh_quantities->cell_vol;   /* fluid */
  const cs_real_t *volume     = cs_glob_mesh_quantities_g->cell_vol; /* geom. */
  const cs_real_3_t *xyzcen   = cs_glob_mesh_quantities_g->cell_cen;
  const cs_time_step_t *ts = cs_glob_time_step;
  const double dtref = ts->dt_ref;
  const state_t &s = state();
 
  const cs_real_t *crom = CS_F_(rho)->val;
 
  /*--------------------------------------------------------------------------
   * Scalars: energy feedback
   *--------------------------------------------------------------------------*/
 
  if (f->dim == 1) {
 
    if (!cfg::feedback_energy)
      return;
 
    /* FIX(port): the Fortran applied this to every coupled scalar (all
       writing the single EnergySource field); restricted to the thermal
       variable. Other scalars receive no coupling term, as before when
       feedback_energy was off. */
    if (f != cs_thermal_model_field())
      return;
 
    const cs_real_t *cvara_s = f->val_pre;
    cs_real_t *mom_t         = detail::field_val("EnergySource");
    cs_real_t *received_temp = detail::field_val("received_T");
 
    const bool is_t = detail::thermal_is_temperature();
    const cs_real_t *cpro_cp = (CS_F_(cp) != nullptr) ? CS_F_(cp)->val : nullptr;
    const double cp0 = cs_glob_fluid_properties->cp0;
 
    const double xtau = 1.0/s.settings.lambda_energy*dtref;
 
    for (cs_lnum_t ipt = 0; ipt < n_cpl_cells; ipt++) {
      cs_lnum_t c = cpl_cells[ipt];
      double xcp = (!is_t) ? 1.0 : ((cpro_cp != nullptr) ? cpro_cp[c] : cp0);
      double xdis = dist_val[ipt];
      double xloc = cvara_s[c];
      received_temp[c] = xdis;
      mom_t[c] = (xdis - xloc)*crom[c]*xcp/xtau;
      rhs[c] += mom_t[c]*cell_f_vol[c];
    }
 
    /* Remove the volume-averaged source so the net energy input is zero */
    double sums[2] = {0.0, 0.0};   /* total source, total volume */
    for (cs_lnum_t ipt = 0; ipt < n_cpl_cells; ipt++) {
      cs_lnum_t c = cpl_cells[ipt];
      sums[0] += mom_t[c]*cell_f_vol[c];
      sums[1] += cell_f_vol[c];
    }
    cs_parall_sum(2, CS_DOUBLE, sums);
 
    if (cs_glob_rank_id <= 0)
      bft_printf("SubChCFD: coupling energy source total = %24.16e\n", sums[0]);
 
    const double mom_t_total = sums[0]/sums[1];
    for (cs_lnum_t ipt = 0; ipt < n_cpl_cells; ipt++) {
      cs_lnum_t c = cpl_cells[ipt];
      rhs[c] -= mom_t_total*cell_f_vol[c];
    }
 
    return;
  }
 
  /*--------------------------------------------------------------------------
   * Vectors: momentum feedback
   *--------------------------------------------------------------------------*/
 
  const cs_real_t *cvara_v = f->val_pre;
 
  cs_real_t *mom_s   = detail::field_val("MomentumSource");
  cs_real_t *mom_cx  = detail::field_val("MomSourceCX");
  cs_real_t *mom_hat = detail::field_val("MomSourceSHat");
  cs_real_t *mom_pt  = detail::field_val("MomSourcePotential");
  cs_real_t *feedback_zone = detail::field_val("FeedbackFlag");
  const cs_real_t *volume_porosi
    = cfg::prs ? detail::field_val("porosity") : nullptr;
 
  {
    cs_field_t *f_fb = cs_field_by_name("FeedbackFlag");
    const cs_lnum_t n = cs_mesh_location_get_n_elts(f_fb->location_id)[0];
    for (cs_lnum_t i = 0; i < n; i++)
      feedback_zone[i] = 1.0;
  }
  if (cfg::feedback_user)
    feedback_zone_user();
 
  const double xtau = 1.0/s.settings.lambda*dtref;
 
  /* mom_imp: only defined on coupled cells, stored per coupled cell */
  std::vector<double> mom_imp(3*static_cast<size_t>(n_cpl_cells), 0.0);
 
  for (cs_lnum_t ipt = 0; ipt < n_cpl_cells; ipt++) {
    cs_lnum_t c = cpl_cells[ipt];
    bool active = cfg::feedback && static_cast<int>(feedback_zone[c]) == 1;
    for (int k = 0; k < 3; k++) {
      double xdis = dist_val[3*ipt + k];
      double xloc = cvara_v[3*c + k];
      if (active) {
        mom_s[3*c + k] = (xdis - xloc)*crom[c]/xtau;
        mom_imp[3*ipt + k] = crom[c]/xtau;
        if (cfg::prs) {
          mom_s[3*c + k]     *= volume_porosi[c];
          mom_imp[3*ipt + k] *= volume_porosi[c];
        }
      }
      else {
        mom_s[3*c + k] = 0.0;
        mom_imp[3*ipt + k] = 0.0;
      }
    }
  }
 
  /* Source term correction */
  switch (cfg::mom_crr) {
 
  case 0:
    for (cs_lnum_t c = 0; c < n_cells; c++)
      for (int k = 0; k < 3; k++)
        mom_hat[3*c + k] = mom_s[3*c + k];
    break;
 
  case 1:   /* constant CX over the domain */
    {
      double sums[6] = {0, 0, 0, 0, 0, 0};   /* mom_c[3], itgtd_vol[3] */
      for (cs_lnum_t c = 0; c < n_cells; c++)
        for (int k = 0; k < 3; k++) {
          sums[k]   += mom_s[3*c + k]*volume[c];
          sums[3+k] += volume[c];
        }
      cs_parall_sum(6, CS_DOUBLE, sums);
      for (int k = 0; k < 3; k++) {
        double mom_c = sums[k]/sums[3+k];
        for (cs_lnum_t c = 0; c < n_cells; c++) {
          mom_cx[3*c + k]  = mom_c;
          mom_hat[3*c + k] = mom_s[3*c + k] - mom_c;
        }
      }
    }
    break;
 
  case 2:   /* constant CX over the coupling region */
    {
      for (cs_lnum_t c = 0; c < n_cells; c++)
        for (int k = 0; k < 3; k++)
          mom_hat[3*c + k] = 0.0;
 
      const int n_cpl = cs_sat_coupling_n_couplings();
      for (int icpl = 0; icpl < n_cpl; icpl++) {
 
        if (ts->nt_cur == ts->nt_prev + 1 || saved_boxes().count(icpl) == 0) {
          boxes_t b;
          coupling_xyz(icpl, b.cen_min, b.cen_max, b.fb_min, b.fb_max);
          saved_boxes()[icpl] = b;
        }
        const boxes_t &b = saved_boxes()[icpl];
 
        for (cs_lnum_t c = 0; c < n_cells; c++)
          for (int k = 0; k < 3; k++)
            mom_cx[3*c + k] = 0.0;
 
        double sums[6] = {0, 0, 0, 0, 0, 0};
        for (int k = 0; k < 3; k++)
          for (cs_lnum_t c = 0; c < n_cells; c++)
            if (xyzcen[c][k] > b.fb_min[k] && xyzcen[c][k] < b.fb_max[k]) {
              sums[k]   += mom_s[3*c + k]*volume[c];
              sums[3+k] += volume[c];
            }
        cs_parall_sum(6, CS_DOUBLE, sums);
 
        for (int k = 0; k < 3; k++) {
          double mom_c = sums[k]/sums[3+k];
          for (cs_lnum_t c = 0; c < n_cells; c++)
            if (xyzcen[c][k] > b.fb_min[k] && xyzcen[c][k] < b.fb_max[k]) {
              mom_cx[3*c + k]   = mom_c;
              mom_hat[3*c + k] += mom_s[3*c + k] - mom_c;
            }
        }
      }
    }
    break;
 
  default:  /* 3: to be added (mom_hat left unchanged, as in the Fortran) */
    break;
  }
 
  /* Pressure-gradient difference for the pressure reconstruction */
  if (cfg::ipre_rec == 2 && ts->nt_cur % cfg::fpre_rec == 0) {
 
    /* Fortran: inc = 0 (homogeneous BCs) unless iphydr == 1 or a free
       inlet is present; grdpor is not needed in 9.x (porosity is handled
       inside the gradient routines). */
    const cs_velocity_pressure_param_t *vp = cs_glob_velocity_pressure_param;
    int inc = (vp->iphydr == 1 || vp->iifren == 1) ? 1 : 0;
 
    cs_real_3_t *dp = nullptr;
    CS_MALLOC(dp, m->n_cells_with_ghosts, cs_real_3_t);
    cs_field_gradient_scalar(CS_F_(p), false, inc, dp);
 
    for (cs_lnum_t ipt = 0; ipt < n_cpl_cells; ipt++) {
      cs_lnum_t c = cpl_cells[ipt];
      bool active = cfg::feedback && static_cast<int>(feedback_zone[c]) == 1;
      for (int k = 0; k < 3; k++) {
        /* FIX(port): the Fortran else-branch wrote mom_pt(isou, iel) with
           isou out of range; all three components are zeroed here. */
        if (active && dist_dp != nullptr)
          mom_pt[3*c + k] = dp[c][k] - dist_dp[3*ipt + k];
        else
          mom_pt[3*c + k] = 0.0;
      }
    }
 
    CS_FREE(dp);
  }
 
  /* Add the source terms */
  for (cs_lnum_t ipt = 0; ipt < n_cpl_cells; ipt++) {
    cs_lnum_t c = cpl_cells[ipt];
    for (int k = 0; k < 3; k++) {
      rhs[3*c + k] += mom_hat[3*c + k]*cell_f_vol[c];
      if (cfg::imp_source) {
        double xloc = cvara_v[3*c + k];
        rhs[3*c + k] += mom_imp[3*ipt + k]*xloc*cell_f_vol[c];
        /* FIX(port): the Fortran overwrote the diagonal (crvimp = ...),
           discarding any implicit user source term (e.g. friction) on
           coupled cells, while keeping its explicit part. Added instead. */
        fimp[9*c + 3*k + k] -= mom_imp[3*ipt + k]*cell_f_vol[c];
      }
    }
  }
 
  /* Received pressure gradient, for checking the data passing */
  if (dist_dp != nullptr) {
    cs_real_t *delta_p = detail::field_val("delta_p");
    for (cs_lnum_t ipt = 0; ipt < n_cpl_cells; ipt++) {
      cs_lnum_t c = cpl_cells[ipt];
      for (int k = 0; k < 3; k++)
        delta_p[3*c + k] = dist_dp[3*ipt + k];
    }
  }
 
  CS_UNUSED(cpl_id);
}
 
} // namespace subchcfd
 
/*============================================================================
 * Replacement for cs_sat_coupling_exchange_at_cells (code_saturne 9.2.0).
 *
 * Same exchanges as the library version (so an unmodified distant instance
 * stays in step), plus the optional pressure-gradient exchange, then either
 * the SubChCFD source terms (cfg::cpl) or the stock ones.
 * Verbosity output of the stock version is omitted.
 *============================================================================*/
 
void
cs_sat_coupling_exchange_at_cells(cs_field_t  *f,
                                  cs_real_t   *rhs,
                                  cs_real_t   *fimp)
{
  using namespace subchcfd;
 
  if (!(f->type & CS_FIELD_VARIABLE))
    bft_error(__FILE__, __LINE__, 0,
              "'%s' Cannot be called for non variable fields.\n", __func__);
 
  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);
  if (eqp->verbosity > 0)
    bft_printf("\n  Code-code coupling (SubChCFD), add volume source terms "
               "for field %s\n\n", f->name);
 
  const int n_cpl = cs_sat_coupling_n_couplings();
 
  for (int cpl_id = 0; cpl_id < n_cpl; cpl_id++) {
 
    const sat_coupling_head_t *cpl = cpl_head(cpl_id);
    if (cpl->localis_cel == nullptr)
      continue;
 
    const int dim  = f->dim;
    const int dim2 = dim*dim;
    ple_locator_t *cells_locator = cpl->localis_cel;
 
    const cs_lnum_t n_cells_loc  = ple_locator_get_n_interior(cells_locator);
    const cs_lnum_t *cpl_cells_ids
      = ple_locator_get_interior_list(cells_locator);
    const cs_lnum_t n_cells_dist = ple_locator_get_n_dist_points(cells_locator);
 
    cs_lnum_t n_g[2] = {n_cells_loc, n_cells_dist};
    cs_parall_sum(2, CS_LNUM_TYPE, n_g);
    const cs_lnum_t n_g_cells_loc  = n_g[0];
    const cs_lnum_t n_g_cells_dist = n_g[1];
 
    const cs_lnum_t *dist_elt_ids = ple_locator_get_dist_locations(cells_locator);
    const cs_coord_t *dist_coords = ple_locator_get_dist_coords(cells_locator);
 
    std::vector<cs_real_t> cw1_dis(static_cast<size_t>(n_cells_dist)*dim, 0.0);
    std::vector<cs_real_t> cw2_dis(static_cast<size_t>(n_cells_dist)*dim2, 0.0);
    std::vector<cs_real_t> cw1_loc(static_cast<size_t>(n_cells_loc)*dim, 0.0);
    std::vector<cs_real_t> cw2_loc(static_cast<size_t>(n_cells_loc)*dim2, 0.0);
 
    if (n_g_cells_dist > 0 && cpl->reverse == 0)
      compute_data_at_cells(f, cpl->reverse, n_cells_dist, dist_elt_ids,
                            dist_coords, cw1_dis.data(), cw2_dis.data());
 
    if (n_g_cells_loc > 0 && cpl->reverse == 1)
      compute_data_at_cells(f, cpl->reverse, n_cells_loc, cpl_cells_ids,
                            dist_coords, cw1_loc.data(), cw2_loc.data());
 
    if (n_g_cells_loc > 0 || n_g_cells_dist > 0) {
      ple_locator_exchange_point_var(cells_locator,
                                     (n_cells_dist > 0) ? cw1_dis.data() : nullptr,
                                     (n_cells_loc > 0) ? cw1_loc.data() : nullptr,
                                     nullptr, sizeof(cs_real_t), dim,
                                     cpl->reverse);
 
      if (cpl->reverse == 1)
        ple_locator_exchange_point_var(cells_locator,
                                       (n_cells_dist > 0) ? cw2_dis.data() : nullptr,
                                       (n_cells_loc > 0) ? cw2_loc.data() : nullptr,
                                       nullptr, sizeof(cs_real_t), dim2,
                                       cpl->reverse);
    }
 
    /* Pressure gradient (velocity only, both instances) */
    const cs_real_t *dist_dp = nullptr;
    if (dim == 3 && pressure_gradient_exchange_step())
      exchange_pressure_gradient(cpl_id, last_dist_dp()[cpl_id]);
    if (dim == 3 && last_dist_dp().count(cpl_id) > 0
        && !last_dist_dp()[cpl_id].empty())
      dist_dp = last_dist_dp()[cpl_id].data();
 
    /* Source terms */
    if (cfg::cpl) {
      if (cpl->reverse == 0)
        coupling_source_terms(f, cpl_id, n_cells_loc, cpl_cells_ids,
                              cw1_loc.data(), dist_dp, rhs, fimp);
      else
        bft_error(__FILE__, __LINE__, 0,
                  "SubChCFD: reverse-mode couplings are not supported "
                  "with cfg::cpl.");
    }
    else
      stock_source_terms(f, cpl->reverse,
                         n_cells_loc, n_g_cells_loc, cpl_cells_ids,
                         n_cells_dist, n_g_cells_dist, dist_elt_ids,
                         cw1_loc.data(), cw1_dis.data(), cw2_dis.data(),
                         rhs, fimp);
  }
}
