/*============================================================================
 * User definition of porous media.
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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_porosity.c
 *
 * \brief User definitions of porous media.
 *
 * See \ref cs_porosity for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the porosity (volume factor \f$ \epsilon \f$
 *        when the porosity model is activated
 *        (\ref cs_glob_porous_model > 0).
 *
 * This function is called at the begin of the simulation only.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

extern bool porosity();
extern cal_porosity();

void
cs_user_porosity(cs_domain_t   *domain)
{
  /*!< [init_poro_mq] */
  
  bool prs;
  prs = porosity();
  if(!prs) return;
  
  cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  const cs_lnum_2_t *i_face_cells
    = (const cs_lnum_2_t *)m->i_face_cells;

  const cs_lnum_t *b_face_cells
    = (const cs_lnum_t *)m->b_face_cells;

  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)mq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)mq->b_face_cog;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)mq->cell_cen;
  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)mq->i_face_normal;
  cs_real_3_t *restrict i_f_face_normal
    = (cs_real_3_t *restrict)mq->i_f_face_normal;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)mq->b_face_normal;
  cs_real_3_t *restrict b_f_face_normal
    = (cs_real_3_t *restrict)mq->b_f_face_normal;

  const cs_real_t *i_f_face_surf = mq->i_f_face_surf;
  const cs_real_t *i_face_surf = mq->i_face_surf;
  /*!< [init_poro_mq] */

  /* Get the cell porosity field value */
  /*!< [init_poro_pro] */
  cal_porosity();
  cs_real_t *cpro_porosi = cs_field_by_name("porosity")->val;
  cs_real_t *volume_porosi = cs_field_by_name("VolumePorosity")->val;
  /*!< [init_poro_pro] */

  /* First, set cell porosity value; the fluid cell volume will be
   * automatically deduced */

  /*!< [set_poro_cells_1] */
  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++)
     cpro_porosi[cell_id] = volume_porosi[cell_id];

  /*!< [set_poro_cells_1] */

  /* Set interior face values */

  /*!< [set_poro_i_faces_1] */
  for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {
    
    cs_lnum_t ii = i_face_cells[face_id][0];
    cs_lnum_t jj = i_face_cells[face_id][1];
    
    cs_real_t face_porosity = 1.;
/*    face_porosity = 0.5*(cpro_porosi[ii]+cpro_porosi[jj]);*/

    for (int i = 0; i < 3; i++)
      i_f_face_normal[face_id][i] = face_porosity * i_face_normal[face_id][i];

    mq->i_f_face_surf[face_id] = cs_math_3_norm(i_f_face_normal[face_id]);

  }
  /*!< [set_poro_i_faces_1] */

  /* Set boundary face values */

  /*!< [set_poro_b_faces_1] */
  for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {

    cs_lnum_t ii = b_face_cells[face_id];
    
    cs_real_t face_porosity = 1.;
/*    face_porosity = cpro_porosi[ii];*/

    for (int i = 0; i < 3; i++)
      b_f_face_normal[face_id][i] = face_porosity * b_face_normal[face_id][i];

    mq->b_f_face_surf[face_id] = cs_math_3_norm(b_f_face_normal[face_id]);

  }
  /*!< [set_poro_b_faces_1] */

  /* Four set face factor */
  if (mq->i_f_face_factor != NULL) {
    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {
      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t face_porosity =
        CS_MAX(
            i_f_face_surf[face_id] / i_face_surf[face_id],
            cs_math_epzero);

      mq->i_f_face_factor[face_id][0] = cpro_porosi[ii] / face_porosity;
      mq->i_f_face_factor[face_id][1] = cpro_porosi[jj] / face_porosity;
    }
  }
}

//the following routines are used to derive mesh connectivity information for 
//the fortran routines
void n_face_vtx(cs_lnum_t *f_id, cs_lnum_t *n_face_vertices)
{
  cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  
  cs_lnum_t  n_i_faces = m->n_i_faces;
  cs_lnum_t  n_b_faces = m->n_b_faces;

  const cs_real_3_t *vtx_coord = (const cs_real_3_t *)m->vtx_coord;
  const cs_lnum_t *face_vtx_idx = (const cs_lnum_t *)m->i_face_vtx_idx;
  const cs_lnum_t *face_vtx = (const cs_lnum_t *)m->i_face_vtx_lst;
  const cs_real_3_t *face_cog = (const cs_real_3_t *)mq->i_face_cog;
  const cs_real_3_t *face_normal = (const cs_real_3_t *)mq->i_face_normal;
  
  cs_lnum_t fid = *f_id;
  
  cs_lnum_t s_id = face_vtx_idx[fid];
  cs_lnum_t e_id = face_vtx_idx[fid + 1];
  
  *n_face_vertices = e_id - s_id;

/*  printf("%d, %d\n",fid, s_id);  */
/*  printf("%d\n",fid);*/
  printf("%d\n",*n_face_vertices);

}


//the following routines are used to derive mesh quantities for 
//the fortran routines
void cord_face_vtx(cs_lnum_t *f_id, cs_lnum_t *n_face_vertices,
                   cs_real_t *cord)
{
  cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  
  cs_lnum_t  n_i_faces = m->n_i_faces;
  cs_lnum_t  n_b_faces = m->n_b_faces;

  const cs_real_3_t *vtx_coord = (const cs_real_3_t *)m->vtx_coord;
  const cs_lnum_t *face_vtx_idx = (const cs_lnum_t *)m->i_face_vtx_idx;
  const cs_lnum_t *face_vtx = (const cs_lnum_t *)m->i_face_vtx_lst;
  const cs_real_3_t *face_cog = (const cs_real_3_t *)mq->i_face_cog;
  const cs_real_3_t *face_normal = (const cs_real_3_t *)mq->i_face_normal;
  
  cs_lnum_t fid = *f_id;
  cs_lnum_t nfvtx = *n_face_vertices;
  
  cs_lnum_t s_id = face_vtx_idx[fid]; 
  cs_lnum_t idx[nfvtx];
  
  cs_real_3_t *cord_vtx = (cs_real_3_t *) cord;
  cs_lnum_t v0;
  
  for (cs_lnum_t i = 0; i < nfvtx; i++) {
    idx[i] = face_vtx[s_id+i];
    v0 = idx[i];
    for (cs_lnum_t j = 0; j < 3; j++)
      cord_vtx[i][j] = vtx_coord[v0][j];
  }

  printf("%d, %d\n",fid, idx[1]);  
//The first coordinate of the 4th vertex corresponds to cord_vtx(4,1) in the 
//Fortran routine
//  printf("%g\n",cord_vtx[0][3]);
  
}


//The following is a general testing routine 
void modifyarray(cs_real_t * insv)
{   
    cs_real_3_t *insv_3 = (cs_real_3_t *) insv;
    insv_3[0][1] = 1234.00;
}


//The following routine is to derive global cell id for the Fortran routine
void global_cell_id(cs_lnum_t *c_id, cs_gnum_t *c_g_id)
{
  cs_mesh_t *m = cs_glob_mesh;
  
  const cs_gnum_t *global_cell_num = (const cs_gnum_t *)m->global_cell_num;
  
  cs_lnum_t cid = *c_id;
  
  cs_gnum_t cgid = global_cell_num[cid-1];

  *c_g_id = cgid;
  

/*  printf("%d, %d, %d, %d\n",cs_glob_rank_id, cid, cgid, *c_g_id);*/
/*  printf("%d, %d, %d, %d\n",cs_glob_rank_id, cid, global_cell_num[cid-1], */
/*  m->global_cell_num[cid-1]);*/
  
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
