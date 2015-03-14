//
//  parse.c
//  CalcIsoStruct
//
//  Created by Christian Gerber on 11/28/12.
//  Copyright (c) 2012 EAWAG. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stddef.h> /* For int */   

#include "parse.h"
#include "isotope.h"
#include "element.h"


int parse_element(  Element* element,
                    char* s,
                    size_t amount,
                    char* iso_list)
{
    
    int number = 0;
    char* symbol = (char*)calloc(MAX_NAME_SIZE, sizeof(char));
    char* iso = (char*)calloc(MAX_NAME_SIZE, sizeof(char));

    Isotope* isotopes = (Isotope*)calloc(MAX_ISO_ELEM, sizeof(Isotope));
    
    double mass = 0.0;
    double abundance = 0.0;
    char* tmp = (char*)calloc(128,sizeof(char));
    size_t i = 0;
    size_t j = 0;
    unsigned short k = 0;
    int found = 0;
    
    while(strcmp((iso_list + i), "@") != 0)
    {
        tmp[j] = iso_list[i];
        if (!strncmp((tmp + j), "$", sizeof(char)) || !strcmp((iso_list + i + 1), "@")) {
            sscanf(tmp, " %d %s %s %lf %lf ",&number, symbol, iso, &mass, &abundance);
            if (strcmp(s, symbol) == 0 && abundance != 0.0){
                set_isotope((isotopes + k), symbol, iso, mass, abundance);
                k++;
                
                found = 1;
            }else if(found && strcmp(s, symbol) != 0){
                free(symbol);
                free(iso);
                free(tmp);
                
                if (k > 0) {
                    set_element(element, isotopes, s, amount, k);
                    free(isotopes);
                    return 0;
                }
                free(isotopes);
                return 1;
            }
            j = 0;
        }
        j++;
        i++;
        
    }
    free(symbol);
    free(iso);
    free(tmp);
    
    if (k > 0) {
        set_element(element, isotopes, s, amount, k);
        free(isotopes);
        return 0;
    }
    free(isotopes);
    return 1;
}


int parse_sum_formula(  Element* elements,
                        char* sum_formula,
                        unsigned short * element_amount,
                        unsigned short * iso_amount,
                        char* iso_list)
{
    if(!sum_formula) return 1;
    
    char symbol[MAX_NAME_SIZE];
    char num[MAX_NAME_SIZE];
    size_t m = 0;
    unsigned long number = 0;
    unsigned short i = 0;

    while(*sum_formula)
    {
        m = 0;
        if (ispunct(*sum_formula)){
            if (isdigit(*(sum_formula + 1)) ) {
                symbol[m] = *sum_formula;
                symbol[m + 1]='\0';
                m++;
                sum_formula++;
                while (isdigit(*sum_formula)) {
                    symbol[m] = *sum_formula;
                    symbol[m+1]='\0';
                    m++;
                    sum_formula++;
                }
                symbol[m] = *sum_formula;
                symbol[m+1]='\0';
                m++;
                sum_formula++;
                if (!isalpha(*(sum_formula))) {
                    return 1;
                }
            }
            else{
                sum_formula++;
                continue;
            }
        }
        while( isalpha(*sum_formula) )
        {
            symbol[m] = *sum_formula;
            symbol[m+1] = '\0';
            m++;
            sum_formula++;
            if (isupper(*sum_formula)) {
                if (!parse_element((elements + i), symbol, 1, iso_list)) {
                    m=0;
                    i++;
                }
            }

        }

        m=0;
        if (isdigit(*sum_formula)) {
            while( isdigit( *sum_formula ) )
            {
                num[m]=*sum_formula;
                m++;
                sum_formula++;
            }
            num[m]='\0';
            number = (unsigned long)atol(num);
        }else{
            number = 1;
        }

        if( number <= 0){
            return 1;
        }
        
        if(parse_element((elements + i), symbol, number, iso_list)){
            return 1;
        }
        i++;
    }

    if (i == 0) {
        return 1;
    }
    *(element_amount) = i;
    for (ptrdiff_t  j = 0; j < *element_amount; j++) {
        *iso_amount += (elements + j)->iso_amount;
    }
    return 0;
}

int parse_element_vector(  Element* element,
                         char* s,
                         size_t amount,
                         unsigned short iso_amount_global,
                         char* element_list[],
                         char* isotope_list[],
                         double isotope_mass_list[],
                         double isotope_abundance_list[]
                         )
{
    unsigned short elem_iso_count = 0;
    strcpy(element->name, s);
    element->amount = amount;
    //size_t size = sizeof(isotope_abundance_list) / sizeof(isotope_abundance_list[0]);
    
    
    int found_element = 0;
    for (int i = 0; i < iso_amount_global; i++) {
        if(strcmp(element_list[i], s) == 0){
            if(isotope_abundance_list[i] > 0.0){
                strcpy(element->isotopes[elem_iso_count].isotope, isotope_list[i]);
                element->isotopes[elem_iso_count].mass = isotope_mass_list[i];
                element->isotopes[elem_iso_count].abundance = isotope_abundance_list[i];
                elem_iso_count++;
            }
            found_element = 1;
        }else{
            if (found_element) {
                break;
            }
        }
    }
    element->iso_amount = elem_iso_count;
    qsort(element->isotopes, elem_iso_count, sizeof(Isotope), isotope_sort_by_abundance);
    return 0;
}

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
                           )
{
#if USE_IN_R == 0
    
    
    // additional calculation for profiling purposes
    char list[iso_amount_global * MAX_NAME_SIZE];
    size_t size_list = sizeof(list);
    
    char* i_el_ = (char*)malloc(size_list);
    char* i_il_ = (char*)malloc(size_list);
    char test[10] = "123456789\0";
    
    for (int i = 0; i < iso_amount_global; i++) {
        strcpy(&i_el_[i*MAX_NAME_SIZE], test);
        strcpy(&i_il_[i*MAX_NAME_SIZE], test);
    }
    
    free(i_el_);
    free(i_il_);
    
#endif
    
    if(!sum_formula) return 1;
    
    char symbol[MAX_NAME_SIZE];
    char num[MAX_NAME_SIZE];
    size_t m = 0;
    unsigned long number = 0;
    unsigned short i = 0;
    
    while(*sum_formula)
    {
        m = 0;
        if (ispunct(*sum_formula)){
            if (isdigit(*(sum_formula + 1)) ) {
                symbol[m] = *sum_formula;
                symbol[m + 1]='\0';
                m++;
                sum_formula++;
                while (isdigit(*sum_formula)) {
                    symbol[m] = *sum_formula;
                    symbol[m+1]='\0';
                    m++;
                    sum_formula++;
                }
                symbol[m] = *sum_formula;
                symbol[m+1]='\0';
                m++;
                sum_formula++;
                if (!isalpha(*(sum_formula))) {
                    return 1;
                }
            }
            else{
                sum_formula++;
                continue;
            }
        }
        while( isalpha(*sum_formula) )
        {
            symbol[m] = *sum_formula;
            symbol[m+1] = '\0';
            m++;
            sum_formula++;
            if (isupper(*sum_formula)) {
                if (!parse_element_vector((elements + i),
                                          symbol,
                                          1,
                                          iso_amount_global,
                                          element_list,
                                          isotope_list,
                                          isotope_mass_list,
                                          isotope_abundance_list
                                          )
                    ) {
                    m=0;
                    i++;
                }
            }
            
        }
        
        m=0;
        if (isdigit(*sum_formula)) {
            while( isdigit( *sum_formula ) )
            {
                num[m]=*sum_formula;
                m++;
                sum_formula++;
            }
            num[m]='\0';
            number = (unsigned long)atol(num);
        }else{
            number = 1;
        }
        
        if( number <= 0){
            return 1;
        }
        
        if(parse_element_vector((elements + i),
                                symbol,
                                number,
                                iso_amount_global,
                                element_list,
                                isotope_list,
                                isotope_mass_list,
                                isotope_abundance_list)
           ){
            return 1;
        }
        i++;
    }
    
    if (i == 0) {
        return 1;
    }
    *(element_amount) = i;
    for (ptrdiff_t  j = 0; j < *element_amount; j++) {
        *iso_amount += (elements + j)->iso_amount;
    }
    return 0;
}



int parse_element_vector_R(  Element* element,
                         char* s,
                         size_t amount,
                         unsigned short iso_amount_global,
                         char* element_list,
                         char* isotope_list,
                         double isotope_mass_list[],
                         double isotope_abundance_list[]
                         )
{
    unsigned short elem_iso_count = 0;
    strcpy(element->name, s);
    element->amount = amount;
    
    int found_element = 0;
    for (int i = 0; i < iso_amount_global; i++) {
        if(strcmp(&element_list[i*MAX_NAME_SIZE], s) == 0){
            if(isotope_abundance_list[i] > 0.0){
                strcpy(element->isotopes[elem_iso_count].isotope, &isotope_list[i*MAX_NAME_SIZE]);
                element->isotopes[elem_iso_count].mass = isotope_mass_list[i];
                element->isotopes[elem_iso_count].abundance = isotope_abundance_list[i];
                elem_iso_count++;
            }
            found_element = 1;
        }else{
            if (found_element) {
                break;
            }
        }
    }
    element->iso_amount = elem_iso_count;
    qsort(element->isotopes, elem_iso_count, sizeof(Isotope), isotope_sort_by_abundance);
    return 0;
}

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
                             )
{
    
    if(!sum_formula) return 1;
    
    char symbol[MAX_NAME_SIZE];
    char num[MAX_NAME_SIZE];
    size_t m = 0;
    unsigned long number = 0;
    unsigned short i = 0;
    
    while(*sum_formula)
    {
        m = 0;
        if (ispunct(*sum_formula)){
            if (isdigit(*(sum_formula + 1)) ) {
                symbol[m] = *sum_formula;
                symbol[m + 1]='\0';
                m++;
                sum_formula++;
                while (isdigit(*sum_formula)) {
                    symbol[m] = *sum_formula;
                    symbol[m+1]='\0';
                    m++;
                    sum_formula++;
                }
                symbol[m] = *sum_formula;
                symbol[m+1]='\0';
                m++;
                sum_formula++;
                if (!isalpha(*(sum_formula))) {
                    return 1;
                }
            }
            else{
                sum_formula++;
                continue;
            }
        }
        while( isalpha(*sum_formula) )
        {
            symbol[m] = *sum_formula;
            symbol[m+1] = '\0';
            m++;
            sum_formula++;
            if (isupper(*sum_formula)) {
                if (!parse_element_vector_R((elements + i),
                                          symbol,
                                          1,
                                          iso_amount_global,
                                          element_list,
                                          isotope_list,
                                          isotope_mass_list,
                                          isotope_abundance_list
                                          )
                    ) {
                    m=0;
                    i++;
                }
            }
            
        }
        
        m=0;
        if (isdigit(*sum_formula)) {
            while( isdigit( *sum_formula ) )
            {
                num[m]=*sum_formula;
                m++;
                sum_formula++;
            }
            num[m]='\0';
            number = (unsigned long)atol(num);
        }else{
            number = 1;
        }
        
        if( number <= 0){
            return 1;
        }
        
        if(parse_element_vector_R((elements + i),
                                symbol,
                                number,
                                iso_amount_global,
                                element_list,
                                isotope_list,
                                isotope_mass_list,
                                isotope_abundance_list)
           ){
            return 1;
        }
        i++;
    }
    
    if (i == 0) {
        return 1;
    }
    *(element_amount) = i;
    for (ptrdiff_t  j = 0; j < *element_amount; j++) {
        *iso_amount += (elements + j)->iso_amount;
    }
    return 0;
}
