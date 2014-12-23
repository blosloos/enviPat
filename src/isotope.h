//
//  isotope.h
//  CalcIsoStruct
//
//  Created by Christian Gerber on 11/28/12.
//  Copyright (c) 2012 EAWAG. All rights reserved.
//

#ifndef isotope_h
#define isotope_h

#include "preferences.h"

typedef struct Isotope Isotope;
typedef struct Isotope2 Isotope2;

struct Isotope{
    char symbol[MAX_NAME_SIZE];
    char isotope[MAX_NAME_SIZE];
    double mass;
    double abundance;
};

struct Isotope2{
    char symbol[MAX_NAME_SIZE];
    char isotope[MAX_NAME_SIZE];
    double mass;
    double abundance;
    int element_nr;
    int iso_nr;
    int iso_e_nr;
    int amount;
};

int set_isotope(Isotope* isotope, char* symbol, char* iso_name, double mass, double abundance);
int isotope_sort_by_abundance(const void *a, const void *b);
int isotope2_sort_by_n_abundance_dec(const void *a, const void *b);
int isotope2_sort_by_n_abundance_inc(const void *a, const void *b);

#endif
