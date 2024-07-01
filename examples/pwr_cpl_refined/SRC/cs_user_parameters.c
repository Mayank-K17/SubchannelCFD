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
extern set_interface_user();

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
                               
  // Axial Reynolds number of the sub channel
    cs_parameters_add_property("AxialReynoldsNumber",
                               1,
                               CS_MESH_LOCATION_CELLS);

  // Local Nusselt number of the wall surfacs
    cs_parameters_add_property("LocalNusseltNumber",
                               1,
                               CS_MESH_LOCATION_BOUNDARY_FACES);

  // temp check dp (sent to the coarse mesh baseline model)
    cs_parameters_add_property("delta_p",
                               3,
                               CS_MESH_LOCATION_CELLS);

    int interface_user;
    interface_user = set_interface_user();
    // Field related to the user defined interface boudnary conditions
    if(interface_user == 1){

    // interface velocity
      cs_parameters_add_property("InterfaceVelocity",
                                 3,
                                 CS_MESH_LOCATION_BOUNDARY_FACES);

    // interface bulk velocity
      cs_parameters_add_property("InterfaceBulkVelocity",
                                 3,
                                 CS_MESH_LOCATION_BOUNDARY_FACES);

    // interface k
      cs_parameters_add_property("InterfaceKineticEnergy",
                                 1,
                                 CS_MESH_LOCATION_BOUNDARY_FACES);

    // interface epilon
      cs_parameters_add_property("InterfaceEpsilon",
                                 1,
                                 CS_MESH_LOCATION_BOUNDARY_FACES);

    // interface omega
      cs_parameters_add_property("InterfaceOmega",
                                 1,
                                 CS_MESH_LOCATION_BOUNDARY_FACES); 

    // interface temperature
      cs_parameters_add_property("InterfaceTemperature",
                                 1,
                                 CS_MESH_LOCATION_BOUNDARY_FACES); 

    // interface bulk temperature
      cs_parameters_add_property("InterfaceBulkTemperature",
                                 1,
                                 CS_MESH_LOCATION_BOUNDARY_FACES); 

    // Blending factor for the remote and the local mesh
      cs_parameters_add_property("BlendingFactor",
                                 1,
                                 CS_MESH_LOCATION_BOUNDARY_FACES); 

    // check for wall distance calculation
      cs_parameters_add_property("WallDistanceCheck",
                                 1,
                                 CS_MESH_LOCATION_BOUNDARY_FACES);

    // check for wall heat flux calculation
      cs_parameters_add_property("WallHeatFluxCheck",
                                 1,
                                 CS_MESH_LOCATION_BOUNDARY_FACES); 

      cs_parameters_add_property("WallIndexCheck",
                                 1,
                                 CS_MESH_LOCATION_BOUNDARY_FACES); 

      cs_parameters_add_property("TplusCheck",
                                 1,
                                 CS_MESH_LOCATION_BOUNDARY_FACES);

      cs_parameters_add_property("TstarCheck",
                                 1,
                                 CS_MESH_LOCATION_BOUNDARY_FACES);

      cs_parameters_add_property("PrdtlCheck",
                                 1,
                                 CS_MESH_LOCATION_BOUNDARY_FACES);

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
  turb_model->iturb = CS_TURB_K_EPSILON_LS;
  
  }

  {
    cs_walldistance_activate();
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

void
cs_user_parameters(cs_domain_t   *domain)
{
  CS_NO_WARN_IF_UNUSED(domain);
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

int start_time_step = 1000;

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

END_C_DECLS
