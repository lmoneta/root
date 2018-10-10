
/*******************************************************************************
*
* file ranlxs.h
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef RANLXS_H
#define RANLXS_H

#ifdef __cplusplus
extern "C" {
#endif
   
extern void ranlxs(float r[],int n);
extern void rlxs_init(int level,int seed);
extern int rlxs_size(void);
extern void rlxs_get(int state[]);
extern void rlxs_reset(int state[]);

#ifdef __cplusplus
}
#endif

#endif
