//
//  parse.h
//  CalcIsoStruct
//
//  Created by Christian Gerber on 11/28/12.
//  Copyright (c) 2012 EAWAG. All rights reserved.
//

#ifndef parse_h
#define parse_h

#include "isotope.h"
#include "element.h"

int parse_element(Element* elements,
                  char* e,
                  int amount,
                  char* iso_list);
int parse_sum_formula(Element* elements,
                      char* sum_formula,
                      unsigned short* element_amount,
                      unsigned short* mass_amount,
                      unsigned short* iso_amount,
                      char* iso_list);


#endif
