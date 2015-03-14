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
                  size_t amount,
                  char* iso_list);

int parse_sum_formula(Element* elements,
                      char* sum_formula,
                      unsigned short * element_amount,
                      unsigned short * iso_amount,
                      char* iso_list);

int parse_element_vector(  Element* element,
                         char* s,
                         size_t amount,
                         unsigned short iso_amount_global,
                         char* element_list[],
                         char* isotope_list[],
                         double isotope_mass_list[],
                         double isotope_abundance_list[]
                         );

int parse_sum_formula_vector(
                             Element* elements,
                             char* sum_formula,
                             unsigned short * element_amount,
                             unsigned short * iso_amount,
                             unsigned short iso_amount_global,
                             char* element_list[],
                             char* isotope_list[],
                             double isotope_mass_list[],
                             double isotope_abundance_list[]
                             );

int parse_element_vector_R(  Element* element,
                           char* s,
                           size_t amount,
                           unsigned short iso_amount_global,
                           char* element_list,
                           char* isotope_list,
                           double isotope_mass_list[],
                           double isotope_abundance_list[]
                           );
int parse_sum_formula_vector_R(
                               Element* elements,
                               char* sum_formula,
                               unsigned short * element_amount,
                               unsigned short * iso_amount,
                               unsigned short iso_amount_global,
                               char* element_list,
                               char* isotope_list,
                               double isotope_mass_list[],
                               double isotope_abundance_list[]
                               );



#endif
