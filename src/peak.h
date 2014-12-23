//
//  peak.h
//  CalcIsoStruct
//
//  Created by Christian Gerber on 12/4/12.
//  Copyright (c) 2012 EAWAG. All rights reserved.
//

#ifndef CalcIsoStruct_peak_h
#define CalcIsoStruct_peak_h

#include "isotope.h"
#include "element.h"

typedef struct Peak Peak;

struct Peak{
    double abundance;
    double mass;
};

int peak_sort_by_mass(const void *a, const void *b);
int trace_sort_by_mass(const void *a, const void *b);
#endif
