/*============================================================================
 * User functions for input of calculation parameters.
 *============================================================================*/

/* code_saturne version 8.0 */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_parameters.c
 *
 * \brief User functions for input of calculation parameters.
 *
 * See \ref parameters for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Select physical model options, including user fields.
 *
 * This function is called at the earliest stages of the data setup,
 * so field ids are not available yet.
 */
/*----------------------------------------------------------------------------*/
extern bool coupling();
extern bool porosity();
extern int iporosity_cal();

void
cs_user_model(void)
{

// Declarations and definitions of custom fields
 {

  // ID of the subchannel
    cs_parameters_add_property("SubchannelID",
                               1,
                               CS_MESH_LOCATION_CELLS);

  // Subchannel type (interior:1, edge:2, corner: 3, user ...)
    cs_parameters_add_property("SubchannelType",
                               1, CS_MESH_LOCATION_CELLS);

  // Pitch to diameter ratio
    cs_parameters_add_property("PitchtoDiameterRatio",
                               1, CS_MESH_LOCATION_CELLS);

  // Hydraulic diameter of the subchannel
    cs_parameters_add_property("HydraulicDiameter",
                               1,
                               CS_MESH_LOCATION_CELLS);

  // Total Surface area of the subchannel
    cs_parameters_add_property("TotalSurfaceArea",
                               1,
                               CS_MESH_LOCATION_CELLS);

  // Volume of the subchannel
    cs_parameters_add_property("SubchannelVolume",
                               1,
                               CS_MESH_LOCATION_CELLS);

  // Rod surface area
    cs_parameters_add_property("RodSurfaceArea",
                               1, 
                               CS_MESH_LOCATION_CELLS);

  // Housing wall surface area
    cs_parameters_add_property("HousingWallSurfaceArea",
                               1, 
                               CS_MESH_LOCATION_CELLS);

  // Average subchannel density of the coolant
    cs_parameters_add_property("AverageDensity",
                               1, 
                               CS_MESH_LOCATION_CELLS);

  // Average subchannel dynamic viscosity of the coolant
    cs_parameters_add_property("AverageDynamicViscosity",
                               1, 
                               CS_MESH_LOCATION_CELLS);

  // Average subchannel thermal conductivity of the coolant
    cs_parameters_add_property("AverageThermalConductivity",
                               1, 
                               CS_MESH_LOCATION_CELLS);

  // Average subchannel heat capacity of the coolant
    cs_parameters_add_property("AverageHeatCapacity",
                               1, 
                               CS_MESH_LOCATION_CELLS);

  // Average subchannel temperature
    cs_parameters_add_property("AverageTemperature",
                               1, 
                               CS_MESH_LOCATION_CELLS);

  // Subchannel bulk temperature
    cs_parameters_add_property("SubchannelBulkTemperature",
                               1, 
                               CS_MESH_LOCATION_CELLS);

  // Subchannel wall temperature
    cs_parameters_add_property("SubchannelWallTemperature",
                               1, 
                               CS_MESH_LOCATION_CELLS);

  // Average subchannel velocity
    cs_parameters_add_property("AverageVelocity",
                               3,
                               CS_MESH_LOCATION_CELLS);

  // Subchannel bulk velocity 
    cs_parameters_add_property("BulkVelocity",
                               3,
                               CS_MESH_LOCATION_CELLS);

  // Axial Reynolds number of sub channel
    cs_parameters_add_property("AxialReynoldsNumber",
                               1,
                               CS_MESH_LOCATION_CELLS);

  // Axial frictional factor of subchannel
    cs_parameters_add_property("AxialFrictionFactor",
                               1,
                               CS_MESH_LOCATION_CELLS);

  // Nusselt number of subchannel
    cs_parameters_add_property("NusseltNumber",
                               1,
                               CS_MESH_LOCATION_CELLS);

  // Heat transfer coefficient of the subchannel
    cs_parameters_add_property("HTC",
                               1,
                               CS_MESH_LOCATION_CELLS);

  // Wall temperatures
    cs_parameters_add_property("WallTemperature0",
                               1,
                               CS_MESH_LOCATION_BOUNDARY_FACES);

    cs_parameters_add_property("WallTemperature1",
                               1,
                               CS_MESH_LOCATION_BOUNDARY_FACES);

  // Wall thermal flux
    cs_parameters_add_property("WallThermalFlux",
                               1, 
                               CS_MESH_LOCATION_BOUNDARY_FACES);

  // temperature field in case of solving enthalpy
    if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_ENTHALPY) {
    cs_parameters_add_property("temperature",
                               1,
                               CS_MESH_LOCATION_CELLS);
    }

    cs_parameters_add_property("CHTInternalCoupling",
                               1,
                               CS_MESH_LOCATION_BOUNDARY_FACES);

  // Wall shear stress
    cs_parameters_add_property("WallShearStress",
                               1, 
                               CS_MESH_LOCATION_CELLS);

  // Properties defined for the coupling module 

    bool cpl;
    cpl = coupling();
    if (cpl) {
   // Additional momentum Source term for the coarse mesh model
    cs_parameters_add_property("MomentumSource",
                               3,
                               CS_MESH_LOCATION_CELLS);

    cs_parameters_add_property("MomSourceCX",
                               3,
                               CS_MESH_LOCATION_CELLS);

    cs_parameters_add_property("MomSourceSHat",
                               3,
                               CS_MESH_LOCATION_CELLS);

    cs_parameters_add_property("MomSourcePotential",
                               3,
                               CS_MESH_LOCATION_CELLS);

    cs_parameters_add_property("PressureOutput",
                               1,
                               CS_MESH_LOCATION_CELLS);

    cs_parameters_add_property("PressureGradient",
                               3,
                               CS_MESH_LOCATION_CELLS);

    cs_parameters_add_property("EnergySource",
                               1,
                               CS_MESH_LOCATION_CELLS);

    cs_parameters_add_property("FeedbackFlag",
                               1,
                               CS_MESH_LOCATION_CELLS);

    // temp check dp (received from the embedded model)
    cs_parameters_add_property("delta_p",
                               3,
                               CS_MESH_LOCATION_CELLS);

    cs_parameters_add_property("received_T",
                               1,
                               CS_MESH_LOCATION_CELLS);

    }

  // Properties defined for the porous module 

    bool prs;
    prs = porosity();
    if (prs) {
   // Field related to the porous model
    cs_parameters_add_property("PorousZone",
                               1,
                               CS_MESH_LOCATION_CELLS);

    cs_parameters_add_property("VolumePorosity",
                               1,
                               CS_MESH_LOCATION_CELLS);

    cs_parameters_add_property("AxialResistanceFriction",
                               1,
                               CS_MESH_LOCATION_CELLS);

    cs_parameters_add_property("AxialResistanceDrag",
                               1,
                               CS_MESH_LOCATION_CELLS);
    }

    int iporo_cal;
    iporo_cal = iporosity_cal();
    // Field related to the porous model (method II)
    if(prs && iporo_cal==3){
    // filter ID
    cs_parameters_add_property("FilterID",
                               1,
                               CS_MESH_LOCATION_CELLS);

    // Filtered density of the porous media
    cs_parameters_add_property("FilteredDensity",
                               1, 
                               CS_MESH_LOCATION_CELLS);

    // Filtered dynamic viscosity of the porous media
    cs_parameters_add_property("FilteredDynamicViscosity",
                               1, 
                               CS_MESH_LOCATION_CELLS);

    // Filtered thermal conductivity of the porous media
    cs_parameters_add_property("FilteredThermalConductivity",
                               1, 
                               CS_MESH_LOCATION_CELLS);
                               
    // Filtered heat capacity of the porous media
    cs_parameters_add_property("FilteredHeatCapacity",
                               1, 
                               CS_MESH_LOCATION_CELLS);

    // Filtered velocity
    cs_parameters_add_property("FilteredVelocity",
                               3,
                               CS_MESH_LOCATION_CELLS);

    }

     cs_parameters_create_added_properties();
  }

// Selection of the turbulence models
  {

  /*  Candidate turbulence models:
   *   CS_TURB_NONE: no turbulence model (laminar flow)
   *   CS_TURB_MIXING_LENGTH: mixing length model
   *   CS_TURB_K_EPSILON: standard k-epsilon model
   *   CS_TURB_K_EPSILON_LIN_PROD: k-epsilon model with
   *     Linear Production (LP) correction
   *   CS_TURB_K_EPSILON_LS: Launder-Sharma low Re k-epsilon model
   *   CS_TURB_RIJ_EPSILON_LRR: Rij-epsilon (LRR)
   *   CS_TURB_RIJ_EPSILON_SSG: Rij-epsilon (SSG)
   *   CS_TURB_RIJ_EPSILON_EBRSM: Rij-epsilon (EBRSM)
   *   CS_TURB_LES_SMAGO_CONST: LES (constant Smagorinsky model)
   *   CS_TURB_LES_SMAGO_DYN: LES ("classical" dynamic Smagorisky model)
   *   CS_TURB_LES_WALE: LES (WALE)
   *   CS_TURB_V2F_PHI: v2f phi-model
   *   CS_TURB_V2F_BL_V2K: v2f BL-v2-k
   *   CS_TURB_K_OMEGA: k-omega SST
   *   CS_TURB_SPALART_ALLMARAS: Spalart-Allmaras model */
  
  cs_turb_model_t *turb_model = cs_get_glob_turb_model();
  turb_model->iturb = CS_TURB_K_OMEGA;

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify general numerical and physical user parameters.
 *
 * At the calling point of this function, most model-related most variables
 * and other fields have been defined, so specific settings related to those
 * fields may be set here.
 *
 * At this stage, the mesh is not built or read yet, so associated data
 * such as field values are not accessible yet, though pending mesh
 * operations and some fields may have been defined.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/
extern int iporosity();

void
cs_user_parameters(cs_domain_t   *domain)
{

/*--------------------------------------------*/
// - 0 No porosity taken into account (Default)
// - 1 Porosity taken into account (isotropic)
// - 2 Porosity taken into account (anisotropic)

  {
    cs_glob_porous_model = 0;
    bool prs;
    prs = porosity();
    int iprs;
    iprs = iporosity();
    if(prs){
      cs_glob_porous_model = iprs;
      if(cs_glob_rank_id==0){
        printf(" Porous media module is used \n");
        if(cs_glob_porous_model==1)
          printf(" Porous model: %d (isotropic)\n", cs_glob_porous_model);
        if(cs_glob_porous_model==2)
          printf(" Porous model: %d (anisotropic)\n", cs_glob_porous_model);
      }
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define linear solver options.
 *
 * This function is called at the setup stage, once user and most model-based
 * fields are defined.
 *
 * Available native iterative linear solvers include conjugate gradient,
 * Jacobi, BiCGStab, BiCGStab2, and GMRES. For symmetric linear systems,
 * an algebraic multigrid solver is available (and recommended).
 *
 * External solvers may also be setup using this function, the cs_sles_t
 * mechanism allowing such through user-define functions.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_linear_solvers(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define time moments.
 *
 * This function is called at the setup stage, once user and most model-based
 * fields are defined, and before fine control of field output options
 * is defined.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_time_moments(void)
{

int start_time_step = 500;

  {
    int moment_f_id[] = {CS_F_(vel)->id};
    int moment_c_id[] = {-1};
    int n_fields = 1;
    cs_time_moment_define_by_field_ids("TimeAverageVelocity",
                                       n_fields,
                                       moment_f_id,
                                       moment_c_id,
                                       CS_TIME_MOMENT_MEAN,
                                       start_time_step, /* nt_start */
                                       -1,   /* t_start */
                                       CS_TIME_MOMENT_RESTART_AUTO,
                                       NULL);
  }

  {
    int moment_f_id[] = {CS_F_(p)->id};
    int moment_c_id[] = {-1};
    int n_fields = 1;
    cs_time_moment_define_by_field_ids("TimeAveragePressure",
                                       n_fields,
                                       moment_f_id,
                                       moment_c_id,
                                       CS_TIME_MOMENT_MEAN,
                                       start_time_step, /* nt_start */
                                       -1,   /* t_start */
                                       CS_TIME_MOMENT_RESTART_AUTO,
                                       NULL);
  }

  {
    cs_field_t *scal = cs_field_by_name("AverageVelocity");
    int moment_f_id[] = {scal->id};
    int moment_c_id[] = {-1};
    int n_fields = 1;
    cs_time_moment_define_by_field_ids("TimeAverageSubchannelVelocity",
                                       n_fields,
                                       moment_f_id,
                                       moment_c_id,
                                       CS_TIME_MOMENT_MEAN,
                                       start_time_step, /* nt_start */
                                       -1,   /* t_start */
                                       CS_TIME_MOMENT_RESTART_AUTO,
                                       NULL);
  }

  {
    bool cpl;
    cpl = coupling();
    if (cpl) {
      cs_field_t *scal = cs_field_by_name("PressureOutput");
      int moment_f_id[] = {scal->id};
      int moment_c_id[] = {-1};
      int n_fields = 1;
      cs_time_moment_define_by_field_ids("TimeAveragePressureOutput",
                                         n_fields,
                                         moment_f_id,
                                         moment_c_id,
                                         CS_TIME_MOMENT_MEAN,
                                         start_time_step, /* nt_start */
                                         -1,   /* t_start */
                                         CS_TIME_MOMENT_RESTART_AUTO,
                                         NULL);
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define internal coupling options.
 *
 * Options are usually defined using cs_internal_coupling_add_entity.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_internal_coupling(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify output user parameters.
 *
 * For CDO schemes, this function concludes the setup of properties,
 * equations, source terms...
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_finalize_setup(cs_domain_t   *domain)
{
  CS_UNUSED(domain);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculate the thermal flux at boundaries.
 *
 * A user field to store the calculated data needs to be defined
 * boundary temperature is can be forced to exist, however, boundary thermal
 * fluxes, for some reason, cannot be directly accessed in the code ...
 */
/*----------------------------------------------------------------------------*/

void
cal_boundary_thermal_flux(void)
{

cs_field_t *f = cs_thermal_model_field();
const char *name = f->name;

const int location_id = CS_MESH_LOCATION_BOUNDARY_FACES;
const cs_lnum_t n_elts = cs_mesh_location_get_n_elts(location_id)[0];
/*  cs_real_t *boundary_flux = NULL;*/
/*  BFT_MALLOC(boundary_flux, n_elts, cs_real_t);*/
cs_real_t *boundary_flux = cs_field_by_name("WallThermalFlux")->val;

cs_post_boundary_flux(f->name, n_elts, NULL, boundary_flux);
/*printf("%s \n", name);*/

}

/*----------------------------------------------------------------------------*/
END_C_DECLS
