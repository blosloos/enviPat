//
//  peak.c
//  CalcIsoStruct
//
//  Created by Christian Gerber on 12/4/12.
//  Copyright (c) 2012 EAWAG. All rights reserved.
//

#include <stdio.h>
#include "peak.h"


int peak_sort_by_mass(const void *a, const void *b)
{
    double y1 = ((const struct Peak*)a)->mass;
    double y2 = ((const struct Peak*)b)->mass;
    
    if (y1 > y2) {
        return 1;
    } else if (y1 < y2) {
        return -1;
    }
    return 0;
}

int trace_sort_by_mass(const void *a, const void *b)
{
    double y1 = *(double *)a;
    double y2 = *(double *)b;
    if (y1 > y2) {
        return 1;
    } else if (y1 < y2) {
        return -1;
    }
    return 0;
}
