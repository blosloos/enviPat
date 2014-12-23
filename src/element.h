//
//  element.h
//  CalcIsoStruct
//
//  Created by Christian Gerber on 11/29/12.
//  Copyright (c) 2012 EAWAG. All rights reserved.
//

#ifndef CalcIsoStruct_element_h
#define CalcIsoStruct_element_h

#include "isotope.h"
#include "preferences.h"

typedef struct Element Element;

struct Element{
    Isotope isotopes[MAX_ISO_ELEM];
    char name[MAX_NAME_SIZE];
    int amount;
    int iso_amount;
};

int set_element(Element* element, Isotope* isotopes, char* name, int amount, int iso_amount);
void print_element(Element* element);
int elements_sort_by_isoamount_inc(const void *a, const void *b);
int elements_sort_by_isoamount_dec(const void *a, const void *b);

#endif
