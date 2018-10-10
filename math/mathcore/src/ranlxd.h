
/*******************************************************************************
*
* file ranlxd.h
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef RANLXD_H
#define RANLXD_H

#ifdef __cplusplus
extern "C" {
#endif

extern void ranlxd(double r[],int n);
extern void rlxd_init(int level,int seed);
extern int rlxd_size(void);
extern void rlxd_get(int state[]);
extern void rlxd_reset(int state[]);

#ifdef __cplusplus
}
#endif

#endif

