/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2009 Sandia National Laboratories.                          *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include <math.h>
#include "zz_const.h"
#include "zz_util_const.h"
#include "zoltan_dd.h"
#include "phg.h"
#include "matrix.h"

void
Zoltan_Matrix_Free(Zoltan_matrix *m)
{
  if (m->yend != m->ystart + 1)
    ZOLTAN_FREE(&m->yend);
  ZOLTAN_FREE(&m->ystart);
  ZOLTAN_FREE(&m->yGNO);
  ZOLTAN_FREE(&m->pinGNO);
  ZOLTAN_FREE(&m->pinwgt);
  ZOLTAN_FREE(&m->yGID);  
  ZOLTAN_FREE(&m->ybipart);
  ZOLTAN_FREE(&m->ypid);

  if (m->ddY != m->ddX && m->ddY != NULL)
    Zoltan_DD_Destroy(&m->ddY);
  if (m->ddX != NULL)
    Zoltan_DD_Destroy(&m->ddX);

  memset (m, 0, sizeof(Zoltan_matrix));
}

void
Zoltan_Matrix2d_Free(Zoltan_matrix_2d *m)
{
  Zoltan_Matrix_Free (&m->mtx);

  ZOLTAN_FREE(&m->dist_x);
  ZOLTAN_FREE(&m->dist_y);

  Zoltan_PHGComm_Destroy(m->comm);
}

void
Zoltan_Matrix_Reset(Zoltan_matrix* m)
{
  m->yGNO = NULL;
  m->ystart = NULL;
  m->yend = NULL;
  m->pinGNO = NULL;
  m->pinwgt = NULL;
  m->yGID = NULL;
}


void
Zoltan_Matrix2d_Init(Zoltan_matrix_2d *m)
{
  memset(m, 0, sizeof(Zoltan_matrix_2d));

  Zoltan_Distribute_Set(m, &Zoltan_Distribute_Origin, (void*)m);
}

int
Zoltan_Distribute_Set(Zoltan_matrix_2d* mat,
		      distFnct *hashDistFct, void * hashDistData)
{
  mat->hashDistData = hashDistData;
  mat->hashDistFct = hashDistFct;
  return (ZOLTAN_OK);
}


int
Zoltan_Matrix_Complete(ZZ* zz,Zoltan_matrix* m)
{
  static char *yo = "Zoltan_Matrix_Complete";
  int ierr = ZOLTAN_OK;

  if(m->completed)
    return (ZOLTAN_OK);

  ZOLTAN_TRACE_ENTER(zz, yo);

  if (m->yend != m->ystart +1 ) {/* Not in compact mode yet */
    int y;
    /* I have to also rewrites all the pinarrays */

    for (y = 1 ; y <= m->nY ; ++y) {
      int length;
      if (m->ystart[y] == m->yend[y-1]) /* No hole */
	continue;
      length = m->yend[y]-m->ystart[y];
      memcpy(m->pinGNO+m->yend[y-1], m->pinGNO+m->ystart[y], length*sizeof(int));
      memcpy(m->pinwgt+m->yend[y-1]*m->pinwgtdim,
	     m->pinGNO+m->ystart[y]*m->pinwgtdim, length*sizeof(float)*m->pinwgtdim);
      m->ystart[y] = m->yend[y-1];
      m->yend[y] = m->ystart[y] + length;
    }

    ZOLTAN_FREE(&m->yend);
    m->yend = m->ystart + 1;
  }

  /* Update data directories */
  m->yGID = ZOLTAN_MALLOC_GID_ARRAY(zz, m->nY);
  m->ypid = (int*) ZOLTAN_MALLOC(m->nY * sizeof(int));
  if (m->bipartite)
    m->ybipart = (int*) ZOLTAN_MALLOC(m->nY * sizeof(int));
  if (m->nY && ((m->yGID == NULL) || (m->ypid == NULL) || (m->bipartite && m->ybipart == NULL)))
    MEMORY_ERROR;

  /* Get Informations about Y */
  Zoltan_DD_Find (m->ddY, (ZOLTAN_ID_PTR)m->yGNO, m->yGID, (ZOLTAN_ID_PTR)m->ypid, m->ybipart,
		  m->nY, NULL);

  if (m->ddY != m->ddX) {
    Zoltan_DD_Destroy(&m->ddY);
    m->ddY = NULL;
  }

  m->completed = 1;
 End:
  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ierr);
}

/* Return an array of locally owned GID */
ZOLTAN_ID_PTR Zoltan_Matrix_Get_GID(ZZ* zz, Zoltan_matrix* m)
{
  ZOLTAN_ID_PTR yGID;

  yGID = ZOLTAN_MALLOC_GID_ARRAY(zz, m->nY);
  if (m->nY && yGID == NULL)
    return (NULL);

  /* Get Informations about Y */
  Zoltan_DD_Find (m->ddY, (ZOLTAN_ID_PTR)m->yGNO, yGID, NULL, NULL,
		  m->nY, NULL);

  return (yGID);
}

int
Zoltan_Matrix2d_adjproc (ZZ* zz, const Zoltan_matrix_2d * const mat, int **adjproc)
{
  static char *yo = "Zoltan_Matrix2d_adjproc";
  int ierr = ZOLTAN_OK;
  int i;
  ZOLTAN_TRACE_ENTER(zz, yo);

  *adjproc = (int*) ZOLTAN_MALLOC(mat->mtx.nPins*sizeof(int));
  if (mat->mtx.nPins && (*adjproc == NULL))
    MEMORY_ERROR;

  for (i = 0 ; i < mat->mtx.nPins ; ++i ) {
    (*adjproc)[i] = EDGE_TO_PROC_Y(mat, mat->mtx.pinGNO[i]);
  }

 End:
  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ierr);
}



#ifdef __cplusplus
}
#endif