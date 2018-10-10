/*
 *  Native routines registration, as per "Writing R extensions" and
 *  definition of native interfaces to two routines exported by
 *  package expint.
 *
 *  This is derived from code in packages zoo and xts.
 *
 *  Copyright (C) 2016 Vincent Goulet
 *  Copyright (C) 2010  Jeffrey A. Ryan jeff.a.ryan @ gmail.com
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License as
 *  published by the Free Software Foundation; either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 *  02110-1301, USA.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <expintAPI.h>		/* this is actually optional */
#include "expint.h"

/* Routine registration and native interfaces definitions. We prefix
 * names with pkg_ to avoid name clashes with expintAPI.h. */
void R_init_discreteHawkes(DllInfo *dll)
{
    /* native interfaces to routines from package expint */
    pkg_expint_E1 = (double(*)(double,int))    R_GetCCallable("expint", "expint_E1");
}
