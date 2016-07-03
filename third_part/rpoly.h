/*
 * rpoly.h
 *
 *  Created on: 2012-1-2
 *      Author: hp
 */

#ifndef RPOLY_H_
#define RPOLY_H_

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>

extern "C" {

void quad(double a, double b1, double c, double *sr, double *si, double *lr,
double *li);
void fxshfr(int l2, int *nz);
void quadit(double *uu, double *vv, int *nz);
void realit(double sss, int *nz, int *iflag);
void calcsc(int *type);
void nextk(int *type);
void newest(int type, double *uu, double *vv);
void quadsd(int n, double *u, double *v, double *p, double *q, double *a,
double *b);

int rpoly_impl(double *op, int degree, double *zeror, double *zeroi, int info[]);

}

#endif /* RPOLY_H_ */
