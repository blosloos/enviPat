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
    size_t amount;
    int all_iso_calc_amount;
    unsigned short iso_amount;
};

int set_element(Element* element, Isotope* isotopes, char* name, size_t amount, unsigned short iso_amount);
void print_element(Element* element);
int elements_sort_by_isoamount_inc(const void *a, const void *b);
int elements_sort_by_isoamount_dec(const void *a, const void *b);
int combine_elements(Element* a, Element*b);


#endif
