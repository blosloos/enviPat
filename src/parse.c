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

#include "parse.h"
#include "isotope.h"
#include "element.h"


int parse_element(  Element* element,
                    char* s,
                    int amount,
                    char* iso_list)
{
    
    int number = 0;
    
    char* symbol = (char*)malloc(MAX_NAME_SIZE * sizeof(char));
    char* iso = (char*)malloc(MAX_NAME_SIZE * sizeof(char));

    Isotope* isotopes = (Isotope*)malloc(MAX_ISO_ELEM * sizeof(Isotope));
    
    double mass = 0.0;
    double abundance = 0.0;
    
    char* tmp = (char*)calloc(128,sizeof(char));
    
    int i = 0;
    int j = 0;
    int k = 0;
    while(strcmp((iso_list + i), "@") != 0)
    {
        tmp[j] = iso_list[i];
        if (!strncmp((tmp + j), "$", sizeof(char)) || !strcmp((iso_list + i + 1), "@")) {
            sscanf(tmp, " %d %s %s %lf %lf ",&number, symbol, iso, &mass, &abundance);
            if (strcmp(s, symbol) == 0 && abundance != 0.0){
                set_isotope((isotopes + k), symbol, iso, mass, abundance);
                k++;
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
                        unsigned short* element_amount,
                        unsigned short* mass_amount,
                        unsigned short* iso_amount,
                        char* iso_list)
{
    if(!sum_formula) return 1;
    
    char symbol[MAX_NAME_SIZE];
    char num[MAX_NAME_SIZE];
    int m = 0;
    int number = 0;
    int i = 0;

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
            number = atoi(num);
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
    for (unsigned short j = 0; j < *element_amount; j++) {
        *mass_amount += (elements + j)->amount;
        *iso_amount += (elements + j)->iso_amount;
    }
    return 0;
}
