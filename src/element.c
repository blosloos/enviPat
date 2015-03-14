//
//  element.c
//  CalcIsoStruct
//
//  Created by Christian Gerber on 11/29/12.
//  Copyright (c) 2012 EAWAG. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "element.h"
#include "isotope.h"


int set_element(Element* element, Isotope* isotopes, char* name, size_t amount, unsigned short iso_amount)
{
    element->amount = amount;
    element->iso_amount = iso_amount;
    strcpy(element->name, name);
    memmove(element->isotopes, isotopes, iso_amount*sizeof(Isotope));
    
    qsort(element->isotopes, iso_amount, sizeof(Isotope), isotope_sort_by_abundance);
	return 0;
}

int elements_sort_by_isoamount_inc(const void *a, const void *b)
{
    unsigned short y1 = ((const struct Element*)a)->iso_amount;
    unsigned short y2 = ((const struct Element*)b)->iso_amount;
    
    if (y1 > y2) {
        return 1;
    } else if (y1 < y2) {
        return -1;
    }
    return 0;
}

int elements_sort_by_isoamount_dec(const void *a, const void *b)
{
    unsigned short y1 = ((const struct Element*)a)->iso_amount;
    unsigned short y2 = ((const struct Element*)b)->iso_amount;
    
    if (y1 < y2) {
        return 1;
    } else if (y1 > y2) {
        return -1;
    }
    return 0;
}
