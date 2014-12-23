//
//  profile.h
//  CalcIsoStruct
//
//  Created by Christian Gerber on 1/4/13.
//  Copyright (c) 2013 EAWAG. All rights reserved.
//

#ifndef CalcIsoStruct_profile_h
#define CalcIsoStruct_profile_h

#include "peak.h"

int calc_profile_with_trace(int n,
                            double* m,
                            double* a,
                            unsigned int tr_num,
                            double* trace,
                            double* profile_mass,
                            double* profile_a,
                            int *profile_n,
                            int res,
                            int profile_type,
                            double t);

int calc_profile(double* m,
                 double* a,
                 double* profile_mass,
                 double* profile_a,
                 unsigned int* num,
                 int res,
                 int profile_type,
                 double thres_profile);

#endif
