#ifndef CS_SUBCHCFD_COUPLING_H
#define CS_SUBCHCFD_COUPLING_H
 
/*============================================================================
 * SubChCFD code_saturne/code_saturne volume coupling (subchcfd_coupling.cpp),
 * code_saturne 9.x.
 *
 * C++ port of the case's csc2ts_coarse.f90 (code_saturne 8.0): the source
 * terms the coarse SubChCFD instance adds to the momentum and energy
 * equations from the values received from the embedded fine-mesh instance.
 *
 * In 8.0 these were reached through modified predvv.f90 / covofi.f90 /
 * csccel.f90. In 9.x the core already calls
 *   cs_sat_coupling_exchange_at_cells(f, rhs, fimp)
 * from the velocity prediction and scalar solves. subchcfd_coupling.cpp
 * provides its own definition of that function, which takes precedence
 * over the library one (code_saturne must be built with shared libraries,
 * the default). It performs the same data exchange as the stock 9.2
 * function, so it stays compatible with an unmodified distant instance,
 * and then:
 *   - cfg::cpl == true : applies the SubChCFD source terms below
 *                        (csc2ts_coarse) instead of the stock relaxation;
 *   - cfg::cpl == false: applies the stock 9.2 relaxation.
 * When cfg::receive_dp is set, it also exchanges the pressure gradient with
 * the distant instance (see exchange_pressure_gradient()).
 *
 * The definition mirrors code_saturne 9.2.0's cs_sat_coupling.cpp and must
 * be re-checked when moving to another 9.x release.
 *============================================================================*/
 
#include "cs_headers.h"
 
#include <vector>
 
namespace subchcfd {
 
/*----------------------------------------------------------------------------
 * Bounding boxes of the coupling region for coupling cpl_id (0-based):
 *  - xyzcen*: local cells coupled with the distant instance;
 *  - xyzfb* : distant coupled boundary faces located in the local mesh.
 * Same quantities as the Fortran (lelcpl / coocpl); computed locally from
 * the coupling locators, so it does not need the distant instance to call
 * it (unlike the earlier C++ version, which exchanged boxes).
 *----------------------------------------------------------------------------*/
 
void coupling_xyz(int     cpl_id,
                  double  xyzcenmin[3],
                  double  xyzcenmax[3],
                  double  xyzfbmin[3],
                  double  xyzfbmax[3]);
 
/*----------------------------------------------------------------------------
 * SubChCFD coupling source terms (port of csc2ts_coarse) for one coupling.
 *
 * f            coupled variable (velocity, or a scalar)
 * cpl_id       coupling id (0-based)
 * n_cpl_cells  number of local coupled cells
 * cpl_cells    local coupled cell ids (0-based)
 * dist_val     distant values at the coupled cells
 *              (f->dim values per cell)
 * dist_dp      distant pressure gradient at the coupled cells
 *              (3 per cell), or nullptr if not exchanged
 * rhs, fimp    explicit / implicit source terms of the equation
 *              (layouts as in cs_user_source_terms)
 *----------------------------------------------------------------------------*/
 
void coupling_source_terms(const cs_field_t  *f,
                           int                cpl_id,
                           cs_lnum_t          n_cpl_cells,
                           const cs_lnum_t    cpl_cells[],
                           const cs_real_t    dist_val[],
                           const cs_real_t    dist_dp[],
                           cs_real_t          rhs[],
                           cs_real_t          fimp[]);
 
/*----------------------------------------------------------------------------
 * Symmetric pressure-gradient exchange on the cell locator of a coupling.
 *
 * Each instance sends grad(p) of the cell containing each distant point
 * (0th-order, no reconstruction) and receives the distant grad(p) at its
 * own coupled cells into received (3*n_cpl_cells values).
 *
 * Both coupled instances must call this at the same point: when the other
 * instance runs SubChCFD's subchcfd_coupling.cpp with the same
 * cfg::receive_dp / cfg::freceive_dp, this is automatic.
 *----------------------------------------------------------------------------*/
 
void exchange_pressure_gradient(int                     cpl_id,
                                std::vector<cs_real_t> &received);
 
/* True when the pressure gradient is exchanged at the current time step */
bool pressure_gradient_exchange_step();
 
} // namespace subchcfd
 
#endif /* SUBCHCFD_COUPLING_H */
 
