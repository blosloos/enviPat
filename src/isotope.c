//
//  isotope.c
//  CalcIsoStruct
//
//  Created by Christian Gerber on 11/28/12.
//  Copyright (c) 2012 EAWAG. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "isotope.h"


int set_isotope(Isotope* isotope, char* symbol, char* iso_name, double mass, double abundance)
{
    strcpy(isotope->symbol, symbol);
    strcpy(isotope->isotope, iso_name);
    isotope->mass = mass;
    isotope->abundance = abundance;
	return 0;
}


int isotope_sort_by_abundance(const void *a, const void *b)
{
    long double y1 = ((const struct Isotope*)a)->abundance;
    long double y2 = ((const struct Isotope*)b)->abundance;
    
    if (y1 < y2) {
        return 1;
    } else if (y1 > y2) {
        return -1;
    }
    return 0;
}


int isotope2_sort_by_n_abundance_inc(const void *a, const void *b)
{
    size_t x1 = ((const struct Isotope2*)a)->amount;
    size_t x2 = ((const struct Isotope2*)b)->amount;
    
    long double y1 = ((const struct Isotope2*)a)->abundance;
    long double y2 = ((const struct Isotope2*)b)->abundance;
    
    if (x1*y1 > x2*y2) {
        return 1;
    } else if (x1*y1 < x2*y2) {
        return -1;
    }
    return 0;
    
}


int isotope2_sort_by_n_abundance_dec(const void *a, const void *b)
{
    size_t x1 = ((const struct Isotope2*)a)->amount;
    size_t x2 = ((const struct Isotope2*)b)->amount;
    
    long double y1 = ((const struct Isotope2*)a)->abundance;
    long double y2 = ((const struct Isotope2*)b)->abundance;
    
    if (x1*y1 < x2*y2) {
        return 1;
    } else if (x1*y1 > x2*y2) {
        return -1;
    }
    return 0;
}
