//
//  combination.c
//  CalcIsoStruct
//
//  Created by Christian Gerber on 11/29/12.
//  Copyright (c) 2012 EAWAG. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "element.h"
#include "combination.h"
#include "parse.h"
#include "preferences.h"
#include "isotope.h"

#if USE_IN_R == 1
    #include <R.h>
    #include <Rdefines.h>
#endif


/*************************   algorithm 3  *************************************************************************
 
 calculates combinations of the isotopes of the given element
 
 combination:         	contains calculated combinations
 element:    			element struct with isotope information
 mass_amount: 		amount of overall elements within the chemical formula. It is used for calculating a appropriate threshold.
 threshold:      		threshold value to reject combinations with a lower abundance than the threshold value
 
 ***************************************************************************************************************/
int create_combination_algo_3(	Combination_3* combination,
                              Element *element,
                              double threshold)
{
    
#if SHOW_DETAILS == 1
    int update_steps_3 = 0;
#endif
#if USE_REALLOC == 1
    size_t allocated_mem_3 = ALLOC_START_3;
#endif
    
    long double abundance = 1.0;
    double mass = 0.0;
    size_t c = 0;
    size_t last_c = 1;
    unsigned short s_a = element->iso_amount;
    int sum_form[s_a];
    const size_t sum_size = sizeof(sum_form);
    element->all_iso_calc_amount = 0;
    
    // berechnung des monoisotopischen peaks
    for (ptrdiff_t i  = 0; i< element->amount; i++) {
        
        // isotope des elementes sind nach absteigender abundance sortiert
        // deshalb kann muss man keinen index fuer die pointer setzen
        mass += (element->isotopes)->mass;
        abundance *= (element->isotopes)->abundance;
    }
    
    // initialisieren des containers mit anzahl 0 und dem ersten element
    combination->amount = 0;
    
    // die maximal vorkommende abundance waehrend der berechnung der element kombinationen
    combination->max_abundance = abundance;
    combination->element = *element;
    combination->compounds->mass = mass;
    combination->compounds->abundance = abundance;
    
    // counter ist 0 weil nur isotope mit der haeufigsten abundance vorkommen
    combination->compounds->counter = 0;
    
    for (int i = 0; i<s_a; i++) {
        combination->compounds->sum[i] = 0;
    }
    
    combination->compounds->sum[0] = (unsigned int)element->amount;
    
    // zahlvariable fuer das argument rica um die anzahl berechnungen zu erhalten
    element->all_iso_calc_amount++;
    
    // zaehler fuer die anzahl element kombinationen
    c++;
    
    
    // erste for schlaufe geht durch alle isotope eines elementes,
    // faengt bei 1 an, weil der monositopische mit isotopen an der stelle 0 im
    // struct element berechnet wurde
    for (ptrdiff_t j = 1; j < element->iso_amount; j++) {
        
        
        // last_c ist die groesse des kontainers, jede bereits gespeicherte
        // kombination wird nun mit dem j-ten isotop erweitert
        for (ptrdiff_t l = 0; l < last_c; l++) {
            
            // ausgangsmasse und abundance von der vorhergehenden element kombination
            mass = (combination->compounds + l)->mass;
            abundance = (combination->compounds + l)->abundance;
            memmove(sum_form, (combination->compounds + l)->sum, sum_size);
            
            
            // for-schlaufe fuer die anzahl elemente, die ->counter struct variable steht fuer die anzahl isotope
            // die nicht fuer die berechnung des monoistopischen benutzt werden. also wieviel man noch ersetzen kann.
            // element->amount ist die gesamt anzahl elemente, die im molekuel vorkommen
            for (ptrdiff_t k = (combination->compounds + l)->counter; k < element->amount; k++) {
                
                // um die masse und abundance zu berechnen, ist immer nur ein isotop das sich aendert,
                // die summe des j-ten isotopes wir um ein inkrementiert wobei spaeter das monoisotopic isotop
                // um 1 dekrementiert
                sum_form[j]++;
                mass -= element->isotopes->mass;
                mass += (element->isotopes + j)->mass;
                abundance *= (element->isotopes + j)->abundance * (sum_form[0]);
                abundance /= (element->isotopes )->abundance * sum_form[j];
                sum_form[0]--;
                
                
                // die berechneten werte der element kombination werden nun im element kontainer abgespeichert
                (combination->compounds + c)->mass = mass;
                (combination->compounds + c)->abundance =  abundance;
                // counter erhoeht sich, da ein isotop, das zur berechnung des monoisotopischen weniger ist und
                // ein anderes mit einer geringeren abundance dessen platz eingenommen hat.
                (combination->compounds + c)->counter = (unsigned int)k + 1;
                memmove((combination->compounds + c)->sum, sum_form, sum_size);
                
                c++;
                element->all_iso_calc_amount++;
                
                if (c >= MAX_COMPOUNDS_3) {
                    return 3003;
                }
#if USE_REALLOC == 1
                if (c == allocated_mem_3) {
#if SHOW_DETAILS == 1
                    update_steps_3++;
#endif
                    allocated_mem_3 *= ALLOC_FACTOR_3;
                    size_t alloc_ = allocated_mem_3 * sizeof(Compound);
                    combination->compounds = (Compound*)realloc(combination->compounds, alloc_);
                    if (combination->compounds == NULL) {
                        return 3913;
                    }
                }
#endif
                if (combination->max_abundance < abundance) {
                    combination->max_abundance = abundance;
                }
            }
        }
        last_c = c - 1;
    }
    combination->amount = c;
#if SHOW_DETAILS == 1
#if USE_IN_R == 1
    Rprintf("\n    %s, comb_realloc_steps: %d, ", element->name, update_steps_3);
#else
    printf("\n    %s, comb_realloc_steps: %d, ", element->name, update_steps_3);
#endif
    
#endif
   	return 0;
}



/*************************   algorithm 3  *************************************************************************
 
 this function combines the combinations of each element with the maximum combination of the other elements.
 if the resulting abundance value is lower than the threshold value, this combination is deleted, because it is impossible to
 reach the threshold value.
 
 combinations:        contains calculated combinations of each element
 threshold:      		threshold value to reject combinations with a lower abundance than the threshold value
 element_amount:	number of different elements
 
 ***************************************************************************************************************/
int clean_combination_algo_3( Combination_3* combinations, long double threshold, size_t element_amount){
    
    
    long double clean_abundance = 1.0;
    long double clean_abundance_other = 1.0;
    for (ptrdiff_t b = 0; b < element_amount; b++) {
        
        //qsort((combinations + b)->compounds, (combinations + b)->amount, sizeof(Compound), compound_sort_by_abundance_dec);
    
        clean_abundance_other = 1.0;
        for (ptrdiff_t d = 0; d < element_amount; d++) {
            if (d != b) {
                clean_abundance_other *= (combinations + d)->max_abundance;
            }
        }
        
        for (ptrdiff_t c = (ptrdiff_t)(combinations + b)->amount - 1; c >= 0; c--) {
            clean_abundance = 1.0;
            clean_abundance *= ((combinations + b)->compounds + c)->abundance * clean_abundance_other;
            //printf("clean:  %lf\n",((combinations + b)->compounds + c)->abundance);
            if ( clean_abundance < threshold) {
                if (c == (combinations + b)->amount - 1) {
                    (combinations + b)->amount--;
                    
                }
                else{
                    *((combinations + b)->compounds + c) = *((combinations + b)->compounds + (combinations + b)->amount - 1);
                    (combinations + b)->amount--;
                }
            }
        }
    }
    
    return 0;
}


int combine_combinations_algo_3(Combination_3* combinations,
                               double threshold,
                               unsigned short  element_amount,
                               size_t* peak_amount,
                               int  peak_limit,
                               unsigned short iso_amount,
                               long double max_abundance,
                               int rtm, double **m_, double **a_, int **cc_){
    
#if USE_REALLOC == 1
    size_t allocated_mem_pl = ALLOC_START_PL;
#endif
    
    size_t tracking[element_amount];
    tracking[0] = 0;
    for (unsigned short i = 0; i < element_amount; i++) {
        tracking[i] = 0;
        qsort((combinations + i)->compounds, (combinations + i)->amount, sizeof(Compound), compound_sort_by_abundance_dec);
    }
    
    const unsigned short iso_a = iso_amount;
    int sum[iso_a];
    size_t iso_sum_size = sizeof(sum);
    
    double mass = 0.0;
    long double abundance = 1.0;
    size_t v = 0;
    
    while (tracking[0] < (combinations)->amount) {
        
        mass = 0.0;
        abundance = 1.0;
        unsigned int cc_count = 0;
        unsigned int cc_tmp[iso_a];
        
        for (unsigned short j = 0; j < element_amount; j++) {
            
            const unsigned short e_a = (combinations + j)->element.iso_amount;
            int sum_e_iso[e_a];
            size_t e_iso_size = sizeof(sum_e_iso);
            // the last element
            if (j == element_amount - 1) {
                double tmp_mass = 0.0;
                long double tmp_abundance = 1.0;
                //double test_abundance = 1.0;
                for (ptrdiff_t h = 0; h < (combinations + j)->amount; h++) {
                    tmp_mass = mass;
                    tmp_abundance = abundance;
                    tmp_mass += ((combinations + j)->compounds + tracking[j])->mass;
                    //test_abundance = ((combinations + j)->compounds + tracking[j])->abundance;
                    tmp_abundance *= ((combinations + j)->compounds + tracking[j])->abundance;
                    if ((100/max_abundance)* tmp_abundance >= threshold) {
                        *(*m_ + v) = tmp_mass;
                        if (rtm == 3 || rtm == 4) {
                            *(*a_ + v) = (double)tmp_abundance;
                        }else{
                            //printf("\n%zu  %lf  max_abundance: %lf", v, tmp_abundance, max_abundance);
                            *(*a_ + v) = (double)((100/max_abundance) * tmp_abundance);
                        }
                        memmove(*cc_ + v * iso_a, cc_tmp, iso_sum_size);
                        memmove(*cc_ + v * iso_a + cc_count, combinations[j].compounds[tracking[j]].sum, e_iso_size);
                        v++;
                        if (v >= peak_limit) {
                            *(peak_amount) = v;
                            return 3201;
                        }
                        
#if USE_REALLOC == 1
                        if (v == allocated_mem_pl) {
                            
                            allocated_mem_pl *= ALLOC_FACTOR_PL;
                            *m_ = (double*) realloc(*m_, allocated_mem_pl * sizeof(double));
                            if (*m_ == NULL) {
                                *(peak_amount) = v;
                                return 32901;
                            }
                            *a_ = (double*) realloc(*a_, allocated_mem_pl * sizeof(double));
                            if (*a_ == NULL) {
                                *(peak_amount) = v;
                                return 32902;
                            }
                            *cc_ = (int*)realloc(*cc_, allocated_mem_pl * iso_a * sizeof(int));
                            if (*cc_ == NULL) {
                                *(peak_amount) = v;
                                return 32903;
                            }
                        }
#endif
                        tracking[j]++;
                    }else{
                        tracking[j]++;
                        break;
                    }
                }
            }else{
                if (tracking[j] < (combinations + j)->amount) {
                    mass += ((combinations + j)->compounds + tracking[j])->mass;
                    abundance *= ((combinations + j)->compounds + tracking[j])->abundance;
                    
                    for (unsigned short u = 0; u < (combinations + j)->element.iso_amount; u++) {
                        cc_tmp[cc_count] = ((combinations + j)->compounds + tracking[j])->sum[u];
                        cc_count++;
                        if (cc_count >= iso_a) {
                            *(peak_amount) = v;
                            return 3202;
                        }
                    }
                }else{
                    break;
                }
                
            }
        }
        for (ptrdiff_t k = element_amount - 2; k >= 0; k--) {
            if (tracking[k] < combinations[k].amount) {
                tracking[k]++;
                for (ptrdiff_t l = k + 1; l < element_amount; l++) {
                    tracking[l] = 0;
                }
                break;
            }
        }
    }
    
    *(peak_amount) = v;
    return 0;
}


int calc_pattern_algo_3(Element *elements
                        , size_t  *peak_amount
                        , double t
                        , unsigned short iso_amount
                        , unsigned short element_amount
                        , long double a_monoisotopic
                        , int p_l
                        , char* l_n
                        , int rtm
                        , double **m_
                        , double **a_
                        , int **cc_
                        )
{
    
    long double max_abundance = 1.0;
    
    Combination_3* combinations = (Combination_3*)calloc(element_amount, sizeof(Combination_3));
    if(combinations == NULL){
        return 3901;
    }
#if USE_REALLOC == 1
    for (ptrdiff_t  i = 0; i < element_amount; i++) {
        Compound* comp;
        size_t alloc_;
        if (ALLOC_START_3 > MAX_COMPOUNDS_3) {
            alloc_ = MAX_COMPOUNDS_3 * sizeof(Compound);
            comp = (Compound*)malloc(alloc_);
        }else{
            alloc_ = ALLOC_START_3 * sizeof(Compound);
            comp = (Compound*)malloc(alloc_);
        }
        
        if(comp == NULL){
            for (ptrdiff_t  j = 0; j < i; j++) {
                free((combinations + j)->compounds);
            }
            free(combinations);
            return 3902;
        }
        combinations[i].compounds = comp;
    }
#endif
    
    
    for (ptrdiff_t  j = 0; j < element_amount; j++) {
        int msg = create_combination_algo_3((combinations + j), elements + j, t );
        if (msg) {
#if USE_REALLOC == 1
            for (ptrdiff_t  i = 0; i < element_amount; i++) {
                free((combinations + i)->compounds);
            }
#endif
            free(combinations);
            return msg;
        }

        max_abundance *= (combinations + j)->max_abundance;
        //printf("\nalgo3 %s  max_a %lf", (combinations + j)->element.name,(combinations + j)->max_abundance);
    }
    //printf("\nalgo3 max abundance  %lf", max_abundance);
    CompoundMulti* monoisotopic = (CompoundMulti*)calloc(1,sizeof(CompoundMulti));
    if(monoisotopic == NULL){
#if USE_REALLOC == 1
        for (ptrdiff_t  i = 0; i < element_amount; i++) {
            free((combinations + i)->compounds);
        }
#endif
        free(combinations);
        return 3903;
    }
    
    calc_monoisotopic(elements, element_amount, monoisotopic);
    if(rtm == 1 || rtm == 4){
        max_abundance = monoisotopic->abundance;
    }else if(rtm == 2){
        max_abundance = 100;
    }
    free(monoisotopic);
    
    if (element_amount == 0) {
#if USE_REALLOC == 1
        for (ptrdiff_t  i = 0; i < element_amount; i++) {
            free((combinations + i)->compounds);
        }
#endif
        free(combinations);
        return 3001;
    }
    
    if (rtm != 2) {
        clean_combination_algo_3(combinations, t * max_abundance / 100 , element_amount);
    }else{
        clean_combination_algo_3(combinations, t , element_amount);
    }
    
    
//    for (int j = 0; j < element_amount; j++) {
//        printf("\nalgo3 %s   amount %zu", combinations[j].element.name, combinations->amount);
//        for (int i = 0; i < (combinations+j)->amount; i++) {
//            printf("\n  %lf", combinations[j].compounds[i].abundance);
//        }
//    }
    

    int msg = combine_combinations_algo_3(combinations, t, element_amount, peak_amount, p_l, iso_amount, max_abundance, rtm, m_, a_, cc_);
    if (msg) {
#if USE_REALLOC == 1
        for (ptrdiff_t  i = 0; i < element_amount; i++) {
            free((combinations + i)->compounds);
        }
#endif
        free(combinations);
        return 3004;
    }
    
    
#if USE_REALLOC == 1
    for (ptrdiff_t  i = 0; i < element_amount; i++) {
        free((combinations + i)->compounds);
    }
#endif
    free(combinations);
    return 0;
}



/*************************   algorithm 1  *************************************************************************
 
 this functions calculates the maximum possible abundance
 
 ***************************************************************************************************************/
int calc_combination_max_abundance(
                                   Combination_1* combination,
                                   Element *element,
                                   Isotope2 *isotopes,
                                   double threshold,
                                   CombinationMulti_1_A* A,
                                   long double mono_abundance,
                                   int rtm)
{
    
#if SHOW_DETAILS == 1
    combination->maxA = 0;
    combination->maxA2 = 0;
    combination->realloc_steps = 0;
    combination->A_realloc_steps = 0;
    combination->A2_realloc_steps = 0;
#endif
    
#if USE_REALLOC == 1
    size_t allocated_mem_1 = ALLOC_START_1;
    size_t allocated_mem_1_A = ALLOC_START_1_A;
    size_t allocated_mem_1_A2 = ALLOC_START_1_A2;
#endif
    
    size_t e_iso_a = element->iso_amount;
    int e_iso_sum[e_iso_a];
    const size_t sum_size = sizeof(e_iso_sum);
    
    Compound* monoisotopic = (Compound*)calloc(1,sizeof(Compound));
	if(monoisotopic == NULL){
        return 10902;
    }
    
    Compound* current_highest = (Compound*)calloc(1,sizeof(Compound));
    if(current_highest == NULL){
        free(monoisotopic);
        return 10903;
    } 
    
    Compound* current = (Compound*)calloc(1,sizeof(Compound));
    if(current == NULL){
        free(monoisotopic);
        free(current_highest);
        return 10904;
    }    
    
    CombinationMulti_1_C* C = (CombinationMulti_1_C*)calloc(1, sizeof(CombinationMulti_1_C));
    if(C == NULL){
        free(monoisotopic);
        free(current_highest);
        free(current);
        return 10905;
    }
    
    calc_monoisotopic_single(element, monoisotopic);
    
//    long double t_abundance = 1.0;
//    for (unsigned int i = 0; i < element->amount; i++) {
//        t_abundance *= element->isotopes->abundance;
//    }
//    monoisotopic->abundance = t_abundance;
    
    A->amount = 0;
    A->max_abundance = 1.0;
    C->amount = 0;
    C->max_abundance = 1.0;
    
    *current = *monoisotopic;

    element->all_iso_calc_amount = 0;
    size_t c = 0;
    
    long double max_a = monoisotopic->abundance;// mono_abundance;
    long double comb_max_a = max_a;
    
    combination->compounds->mass = monoisotopic->mass;
    combination->compounds->abundance = monoisotopic->abundance;
    memmove(combination->compounds->sum, monoisotopic->sum, sum_size);
    combination->a2_amount = 0;
    element->all_iso_calc_amount++;
    c++;
    
    free(monoisotopic);
    size_t  iso_nr_max = 0;
    
    while (current->abundance != -1.0) {
        *current_highest = *current;
        C->amount = 0;
        iso_nr_max = 0;
        
        for (unsigned short j = current->indicator_iso; j < element->iso_amount - 1; j++) {
            if ( current->counter < element->amount ) {
                
                Isotope *isotope = element->isotopes;
                unsigned short  iso_e_nr = (isotopes + j)->iso_e_nr;
                C->compounds[C->amount] = *current;
                Compound *comp = &C->compounds[C->amount];
                
                comp->counter++;
                comp->indicator_iso = j;
                comp->sum[iso_e_nr]++;
                comp->mass -= isotope->mass;
                comp->mass += (isotope + iso_e_nr)->mass;
                comp->abundance *= (isotope + iso_e_nr)->abundance;
                comp->abundance *= comp->sum[0];
                comp->abundance /= isotope->abundance * comp->sum[0 + iso_e_nr];

//                long double aa = 1.0;
//
//                for (int m = 0 ; m < comp->sum[0 + iso_e_nr];  m++) {
//                    aa *= (isotope + iso_e_nr)->abundance;
//                    aa *= element->amount - m;
//                    aa /= m + 1;
//                }
//
//                for (int n = 0; n < comp->sum[0]; n++) {
//                    aa*=isotope->abundance;
//                }
//                
//                //comp->abundance = aa;// l_a;
               comp->sum[0]--;
//                
//                long double diff = comp->abundance - aa;

                //printf("%zu comp %4.10Lf aa %4.10Lf d %10.120Lf\n", c, comp->abundance, aa, diff);
                
                if (current_highest->abundance < comp->abundance) {
                    *current_highest = *comp;
                }
                
                if (comp->abundance >= current->abundance) {
                    iso_nr_max = j;
                }
                
                element->all_iso_calc_amount++;
                C->amount++;
            }
        }
        
        if (comb_max_a < current_highest->abundance ) {
            
            if (!(rtm == 1 || rtm == 4)) {
                max_a = current_highest->abundance;
            }
            comb_max_a = current_highest->abundance;
        }
        
        if(current_highest->abundance > current->abundance){
            for (ptrdiff_t v = (ptrdiff_t)C->amount - 1; v >= 0 ; v--) {
                if ( C->compounds[v].abundance != current_highest->abundance
                    ) {
                    if(C->compounds[v].indicator_iso <= iso_nr_max) {
#if USE_REALLOC == 1
                        if (A->amount == allocated_mem_1_A) {
#if SHOW_DETAILS == 1
                            combination->A_realloc_steps += 1;
#endif
                            allocated_mem_1_A *= ALLOC_FACTOR_1_A;
                            A->compounds = (Compound*) realloc(A->compounds, allocated_mem_1_A * sizeof(Compound));
                            if (A->compounds == NULL) {
                                free(current_highest);
                                free(current);
                                free(C);
                                return 10905;
                            }
                        }
#endif
                        A->compounds[A->amount] = C->compounds[v];
                        A->amount++;
#if SHOW_DETAILS == 1
                        if (combination->maxA < A->amount) {
                            combination->maxA = A->amount;
                        }
#endif
                        if (A->amount >= MAX_COMPOUNDS_1_A) {;
                            free(current_highest);
                            free(current);
                            free(C);
                            return 1001;
                        }
                    }else{
                        if ( (100/ max_a) * C->compounds[v].abundance >= threshold) {
#if USE_REALLOC == 1
                            if (combination->a2_amount == allocated_mem_1_A2) {
#if SHOW_DETAILS == 1
                                combination->A2_realloc_steps += 1;
#endif
                                allocated_mem_1_A2 *= ALLOC_FACTOR_1_A2;
                                combination->a2_list = (Compound*) realloc(combination->a2_list, allocated_mem_1_A2 * sizeof(Compound));
                                if (combination->a2_list == NULL) {
                                    free(C);
                                    free(current_highest);
                                    free(current);
                                    return 10906;
                                }
                            }
#endif
                            combination->a2_list[combination->a2_amount] = C->compounds[v];
                            combination->a2_amount++;
#if SHOW_DETAILS == 1
                            if (combination->maxA2 < combination->a2_amount) {
                                combination->maxA2 = combination->a2_amount;
                            }
#endif
                            if (combination->a2_amount >= MAX_COMPOUNDS_1_A2) {
                                combination->amount = c;
                                free(current_highest);
                                free(current);
                                free(C);
                                return 1002;
                            }
                        }
                    }
                }
            }
            
            *current = *current_highest;
            if ((100/ max_a) * current_highest->abundance >= threshold) {
#if USE_REALLOC == 1
                if (c == allocated_mem_1) {
#if SHOW_DETAILS == 1
					combination->realloc_steps += 1;
#endif
                    allocated_mem_1 *= ALLOC_FACTOR_1;
                    size_t alloc_ = allocated_mem_1 * sizeof(Compound);
                    combination->compounds = (Compound*) realloc(combination->compounds, alloc_);
                    if (combination->compounds == NULL) {
                        free(current_highest);
                        free(current);
                        free(C);
                        return 10907;
                    }
                }
#endif
                (combination->compounds + c)->mass = current_highest->mass;
                (combination->compounds + c)->abundance =  current_highest->abundance;
                memmove((combination->compounds + c)->sum, current_highest->sum, sum_size);
                c++;
                if (c >= MAX_COMPOUNDS_1) {
                    free(current_highest);
                    free(current);
                    free(C);
                    return 1003;
                }
            }
        }
        else{
            for (ptrdiff_t v = 0; v < C->amount; v++) {
                if ( (100/ max_a) * C->compounds[v].abundance >= threshold) {
#if USE_REALLOC == 1
                    if (combination->a2_amount == allocated_mem_1_A2) {
#if SHOW_DETAILS == 1
                        combination->A2_realloc_steps += 1;
#endif
                        allocated_mem_1_A2 *= ALLOC_FACTOR_1_A2;
                        size_t alloc_ = allocated_mem_1_A2 * sizeof(Compound);
                        combination->a2_list = (Compound*) realloc(combination->a2_list, alloc_);
                        if (combination->a2_list == NULL) {
                            free(current_highest);
                            free(current);
                            free(C);
                            return 10908;
                        }
                    }
#endif
                    combination->a2_list[combination->a2_amount] = C->compounds[v];
                    combination->a2_amount++;
                    
#if SHOW_DETAILS == 1
                    if (combination->maxA2 < combination->a2_amount) {
                        combination->maxA2 = combination->a2_amount;
                    }
#endif
                    if (combination->a2_amount >= MAX_COMPOUNDS_1_A2) {
                        combination->amount = c;
                        free(current_highest);
                        free(current);
                        free(C);
                        return 1002;
                    }
                }
            }
            
            if (A->amount > 0) {
                Compound *a_c = &A->compounds[A->amount - 1];
                *current = *a_c;
                if ((100/ max_a) * a_c->abundance >= threshold) {
#if USE_REALLOC == 1
                    if (c == allocated_mem_1) {
#if SHOW_DETAILS == 1
						combination->realloc_steps += 1;
#endif
                        allocated_mem_1 *= ALLOC_FACTOR_1;
                        size_t alloc_ = allocated_mem_1 * sizeof(Compound);
                        combination->compounds = (Compound*) realloc(combination->compounds, alloc_);
                        if (combination->compounds == NULL) {
                            free(current_highest);
                            free(current);
                            free(C);
                            return 10907;
                        }
                    }
#endif
                    (combination->compounds + c)->mass = a_c->mass;
                    (combination->compounds + c)->abundance =  a_c->abundance;
                    memmove((combination->compounds + c)->sum, a_c->sum, sum_size);
                    c++;
                    
                    if (c >= MAX_COMPOUNDS_1) {
                        free(current_highest);
                        free(current);
                        free(C);
                        return 1003;
                    }
                }
                A->amount--;
            }
            else if(A->amount == 0){
                break;
            }
            else {
                current->abundance = -1.0;
            }
        }
    }

    combination->max_abundance = comb_max_a;
    combination->amount = c;
    free(current_highest);
    free(current);
    free(C);
    
    return 0;
}

#if WITH_C_LIST == 0
int create_combination_algo_1(Combination_1 *combination,
                              Element *element,
                              Isotope2 *isotopes,
                              long double clean_abundance,
                              long double threshold,
                              int peak_limit,
                              CombinationMulti_1_A* A
                              )
{
#if USE_REALLOC == 1
    size_t allocated_mem_1 = ALLOC_START_1;
    if (combination->amount >= ALLOC_START_1) {
#if SHOW_DETAILS == 1
        combination->realloc_steps += 1;
#endif
        allocated_mem_1 = ALLOC_FACTOR_1 * combination->amount;
        size_t alloc_tmp = allocated_mem_1 * sizeof(Compound);
        combination->compounds = (Compound*)realloc(combination->compounds, alloc_tmp);
        if (combination->compounds == NULL) {
            return 11900;
        }
    }
    
    size_t allocated_mem_1_A = ALLOC_START_1_A;
    if (A->amount >= ALLOC_START_1_A) {
#if SHOW_DETAILS == 1
        combination->A_realloc_steps += 1;
#endif
        allocated_mem_1_A = ALLOC_FACTOR_1_A * A->amount;
        size_t alloc_tmp = allocated_mem_1_A * sizeof(Compound);
        A->compounds = (Compound*) realloc(A->compounds, alloc_tmp);
        if (A->compounds == NULL) {
            return 11900;
        }

    }
    
    size_t allocated_mem_1_A2 = ALLOC_START_1_A2;
    if (combination->a2_amount >= ALLOC_START_1_A2) {
#if SHOW_DETAILS == 1
        combination->A2_realloc_steps += 1;
#endif
        allocated_mem_1_A2 = ALLOC_FACTOR_1_A2 * combination->a2_amount;
        size_t alloc_tmp = allocated_mem_1_A2 * sizeof(Compound);
        combination->a2_list = (Compound*) realloc(combination->a2_list, alloc_tmp);
        if (combination->a2_list == NULL) {
            return 11900;
        }

    }
#endif
    
    size_t e_iso_a = element->iso_amount;
    int e_iso_sum[e_iso_a];
    const size_t sum_size = sizeof(e_iso_sum);
    
    Compound* current = (Compound*)calloc(1,sizeof(Compound));
    if(current == NULL){
        return 11903;
    }
    A->amount = 0;
    A->max_abundance = 1.0;
    combination->element = *element;
    size_t c = combination->amount;
    
    if(combination->a2_amount > 0){
        *current = combination->a2_list[combination->a2_amount - 1];
        combination->a2_amount--;
    }else{
        free(current);
        
        if (c == 0 || combination->max_abundance == 1.0) {
            //return 1100;
        }
#if SHOW_DETAILS == 1
#if USE_IN_R == 1
        Rprintf("\n    %s, maxA: %d, maxA2: %d, comb_realloc_steps: %d, A_realloc_steps: %d, A2_realloc_steps: %d, updated_steps: %d, ", element->name, combination->maxA, combination->maxA2, combination->realloc_steps, combination->A_realloc_steps, combination->A2_realloc_steps, element->all_iso_calc_amount);
#else
        printf("\n    %s, maxA: %zu, maxA2: %zu, comb_realloc_steps: %d, A_realloc_steps: %d, A2_realloc_steps: %d, updated_steps: %d, ", element->name, combination->maxA, combination->maxA2, combination->realloc_steps, combination->A_realloc_steps, combination->A2_realloc_steps, element->all_iso_calc_amount);
#endif
#endif
        return 0;
    }
    
    if (clean_abundance * current->abundance >= threshold && current->mass > 1.0) {
        (combination->compounds + c)->mass = current->mass;
        (combination->compounds + c)->abundance =  current->abundance;
        memmove((combination->compounds + c)->sum, current->sum, sum_size);
        c++;
        //element->all_iso_calc_amount++;
    } 
    
    Compound *comp = (Compound*)calloc(1,sizeof(Compound));
    
    while (current->abundance != -1.0) {
        for (unsigned short j = current->indicator_iso; j < element->iso_amount - 1; j++) {
            if ( current->counter < element->amount ) {
                
                Isotope *isotope = element->isotopes;
                unsigned short  iso_e_nr = (isotopes + j)->iso_e_nr;
                *comp = *current;
                
                comp->counter++;
                comp->indicator_iso = j;
                comp->sum[iso_e_nr]++;
                comp->mass -= isotope->mass;
                comp->mass += (isotope + iso_e_nr)->mass;
                comp->abundance *= ( isotope + iso_e_nr)->abundance * (comp->sum[0]);
                comp->abundance /= ( isotope )->abundance * comp->sum[0 + iso_e_nr];
                comp->sum[0]--;
                //element->all_iso_calc_amount++;
                
                if (comp->abundance > current->abundance) {
#if USE_REALLOC == 1
                        if (A->amount == allocated_mem_1_A) {
#if SHOW_DETAILS == 1
                            combination->A_realloc_steps += 1;
#endif
                            allocated_mem_1_A *= ALLOC_FACTOR_1_A;
                            size_t alloc_ = allocated_mem_1_A * sizeof(Compound);
                            A->compounds = (Compound*) realloc(A->compounds, alloc_);
                            if (A->compounds == NULL) {
                                free(current);
                                free(comp);
                                return 11904;
                            }
                        }
#endif
                        A->compounds[A->amount] = *comp;
                        A->amount++;
#if SHOW_DETAILS == 1
                        if (combination->maxA < A->amount) {
                            combination->maxA = A->amount;
                        }
#endif
                        if (A->amount >= MAX_COMPOUNDS_1_A) {
                            combination->amount = c;
                            free(current);
                            free(comp);
                            return 1101;
                        }
                }else{
					if ( clean_abundance * comp->abundance >= threshold ) {
#if USE_REALLOC == 1
	                    if (combination->a2_amount == allocated_mem_1_A2) {
#if SHOW_DETAILS == 1
                            combination->A2_realloc_steps += 1;
#endif
	                        allocated_mem_1_A2 *= ALLOC_FACTOR_1_A2;
	                        size_t alloc_ = allocated_mem_1_A2 * sizeof(Compound);
	                        combination->a2_list = (Compound*) realloc(combination->a2_list, alloc_ );
	                        if (combination->a2_list == NULL) {
	                            free(current);
                                free(comp);
	                            return 11907;
	                        }
	                    }
#endif
	                    combination->a2_list[combination->a2_amount] = *comp;
	                    combination->a2_amount++;
#if SHOW_DETAILS == 1
	                    if (combination->maxA2 < combination->a2_amount) {
	                        combination->maxA2 = combination->a2_amount;
	                    }
#endif
	                    if (combination->a2_amount >= MAX_COMPOUNDS_1_A2) {
	                        combination->amount = c;
	                        free(current);
                            free(comp);
	                        return 1102;
	                    }
					}
				}
                element->all_iso_calc_amount++;
                if (combination->max_abundance < comp->abundance) {
					combination->max_abundance = comp->abundance;
				}
            }
        }
        
        if (c >= MAX_COMPOUNDS_1 || combination->amount >= MAX_COMPOUNDS_1) {
            combination->amount = c;
            free(current);
            free(comp);
            return 1103;
        }

		if (A->amount > 0) {
			Compound *a_c = &A->compounds[A->amount - 1];
			*current = *a_c;
			
			if (clean_abundance * a_c->abundance >= threshold && c < MAX_COMPOUNDS_1) {
#if USE_REALLOC == 1
				if (c == allocated_mem_1) {
#if SHOW_DETAILS == 1
                    combination->realloc_steps += 1;
#endif
					allocated_mem_1 *= ALLOC_FACTOR_1;
					size_t alloc_ = allocated_mem_1 * sizeof(Compound);
					combination->compounds = (Compound*) realloc(combination->compounds, alloc_);
					if (combination->compounds == NULL) {
						free(current);
                        free(comp);
						return 11909;
					}
				}
#endif
				(combination->compounds + c)->mass = a_c->mass;
				(combination->compounds + c)->abundance =  a_c->abundance;
				memmove((combination->compounds + c)->sum, a_c->sum, sum_size);
				c++;
				
				if (c >= MAX_COMPOUNDS_1) {
					free(current);
                    free(comp);
					return 1103;
				}
			}
			A->amount--;
		}else if(A->amount == 0 && combination->a2_amount > 0){
			Compound *a2 = &combination->a2_list[combination->a2_amount - 1];
			*current = *a2;
			if ( clean_abundance * a2->abundance >= threshold && c < MAX_COMPOUNDS_1) {
#if USE_REALLOC == 1
				if (c == allocated_mem_1) {
#if SHOW_DETAILS == 1
                    combination->realloc_steps += 1;
#endif
					allocated_mem_1 *= ALLOC_FACTOR_1;
					size_t alloc_ = allocated_mem_1 * sizeof(Compound);
					combination->compounds = (Compound*) realloc(combination->compounds, alloc_);
					if (combination->compounds == NULL) {
						free(current);
                        free(comp);
						return 11909;
					}
				}
#endif
				(combination->compounds + c)->mass = a2->mass;
				(combination->compounds + c)->abundance =  a2->abundance;
				memmove((combination->compounds + c)->sum, a2->sum, sum_size);
				c++;
				
				if (c >= MAX_COMPOUNDS_1) {
					free(current);
                    free(comp);
					return 1103;
				}
			}
			combination->a2_amount--;
		}else {
			current->abundance = -1.0;
		}
    }
    
#if SHOW_DETAILS == 1
    #if USE_IN_R == 1
        Rprintf("\n    %s, maxA: %d, maxA2: %d, comb_realloc_steps: %d, A_realloc_steps: %d, A2_realloc_steps: %d, updated_steps: %d, ", element->name, combination->maxA, combination->maxA2, combination->realloc_steps, combination->A_realloc_steps, combination->A2_realloc_steps, element->all_iso_calc_amount);
    #else
        printf("\n    %s, maxA: %zu, maxA2: %zu, comb_realloc_steps: %d, A_realloc_steps: %d, A2_realloc_steps: %d, updated_steps: %d, ", element->name, combination->maxA, combination->maxA2, combination->realloc_steps, combination->A_realloc_steps, combination->A2_realloc_steps, element->all_iso_calc_amount);
    #endif
#endif

    combination->amount = c;
    free(current);
    free(comp);
   	return 0;
}

#else // WITH_C_LIST == 1

int create_combination_algo_1(Combination_1 *combination,
                              Element *element,
                              Isotope2 *isotopes,
                              long double clean_abundance,
                              long double threshold,
                              int peak_limit,
                              CombinationMulti_1_A* A
                              )
{
#if USE_REALLOC == 1
    size_t allocated_mem_1 = ALLOC_START_1;
    if (combination->amount >= ALLOC_START_1) {
#if SHOW_DETAILS == 1
        combination->realloc_steps += 1;
#endif
        allocated_mem_1 = ALLOC_FACTOR_1 * combination->amount;
        size_t alloc_tmp = allocated_mem_1 * sizeof(Compound);
        combination->compounds = (Compound*)realloc(combination->compounds, alloc_tmp);
        if (combination->compounds == NULL) {
            return 11900;
        }
    }
    
    size_t allocated_mem_1_A = ALLOC_START_1_A;
    if (A->amount >= ALLOC_START_1_A) {
#if SHOW_DETAILS == 1
        combination->A_realloc_steps += 1;
#endif
        allocated_mem_1_A = ALLOC_FACTOR_1_A * A->amount;
        size_t alloc_tmp = allocated_mem_1_A * sizeof(Compound);
        A->compounds = (Compound*) realloc(A->compounds, alloc_tmp);
        if (A->compounds == NULL) {
            return 11900;
        }
        
    }
    
    size_t allocated_mem_1_A2 = ALLOC_START_1_A2;
    if (combination->a2_amount >= ALLOC_START_1_A2) {
#if SHOW_DETAILS == 1
        combination->A2_realloc_steps += 1;
#endif
        allocated_mem_1_A2 = ALLOC_FACTOR_1_A2 * combination->a2_amount;
        size_t alloc_tmp = allocated_mem_1_A2 * sizeof(Compound);
        combination->a2_list = (Compound*) realloc(combination->a2_list, alloc_tmp);
        if (combination->a2_list == NULL) {
            return 11900;
        }
        
    }
#endif
    
    size_t e_iso_a = element->iso_amount;
    int e_iso_sum[e_iso_a];
    const size_t sum_size = sizeof(e_iso_sum);
    
    Compound* current_highest = (Compound*)calloc(1,sizeof(Compound));
    if(current_highest == NULL){
        return 11902;
    }
    
    Compound* current = (Compound*)calloc(1,sizeof(Compound));
    if(current == NULL){
        free(current_highest);
        return 11903;
    }
    CombinationMulti_1_C* C = (CombinationMulti_1_C*)calloc(1, sizeof(CombinationMulti_1_C));
    if(C == NULL){
        free(current_highest);
        free(current);
        return 10905;
    }
    
    create_isotope_list_single(*element, isotopes);
    A->amount = 0;
    A->max_abundance = 1.0;
    C->amount = 0;
    C->max_abundance = 1.0;
    combination->element = *element;
    
    size_t iso_nr_max = 0;
    size_t c = combination->amount;
    
    if(combination->a2_amount > 0){
        *current = combination->a2_list[combination->a2_amount - 1];
        combination->a2_amount--;
    }else{
        free(C);
        free(current);
        free(current_highest);
        
        if (c == 0 || combination->max_abundance == 1.0) {
            //return 1100;
        }
#if SHOW_DETAILS == 1
#if USE_IN_R == 1
        Rprintf("\n    %s, maxA: %d, maxA2: %d, comb_realloc_steps: %d, A_realloc_steps: %d, A2_realloc_steps: %d, updated_steps: %d, ", element->name, combination->maxA, combination->maxA2, combination->realloc_steps, combination->A_realloc_steps, combination->A2_realloc_steps, element->all_iso_calc_amount);
#else
        printf("\n    %s, maxA: %zu, maxA2: %zu, comb_realloc_steps: %d, A_realloc_steps: %d, A2_realloc_steps: %d, updated_steps: %d, ", element->name, combination->maxA, combination->maxA2, combination->realloc_steps, combination->A_realloc_steps, combination->A2_realloc_steps, element->all_iso_calc_amount);
#endif
#endif
        return 0;
    }
    
    if (clean_abundance * current->abundance >= threshold && current->mass > 1.0) {
        (combination->compounds + c)->mass = current->mass;
        (combination->compounds + c)->abundance =  current->abundance;
        memmove((combination->compounds + c)->sum, current->sum, sum_size);
        c++;
    }
    
    while (current->abundance != -1.0) {
        *current_highest = *current;
        C->amount = 0;
        iso_nr_max = 0;
        
        for (unsigned short j = current->indicator_iso; j < element->iso_amount - 1; j++) {
            if ( current->counter < element->amount ) {
                
                Isotope *isotope = element->isotopes;
                size_t  iso_e_nr = (isotopes + j)->iso_e_nr;
                C->compounds[C->amount] = *current;
                Compound *comp = &C->compounds[C->amount];
                
                comp->counter++;
                comp->indicator_iso = j;
                comp->sum[iso_e_nr]++;
                
                comp->mass -= isotope->mass;
                comp->mass += (isotope + iso_e_nr)->mass;
                comp->abundance *= ( isotope + iso_e_nr)->abundance * (comp->sum[0]);
                comp->abundance /= ( isotope )->abundance * comp->sum[0 + iso_e_nr];
                
                comp->sum[0]--;
                
                if (current_highest->abundance < comp->abundance) {
                    *current_highest = *comp;
                }
                
                if (comp->abundance >= current->abundance) {
                    iso_nr_max = j;
                }
                
                element->all_iso_calc_amount++;
                C->amount++;
            }
        }
        
        if (c >= MAX_COMPOUNDS_1 || combination->amount >= MAX_COMPOUNDS_1) {
            combination->amount = c;
            free(C);
            free(current);
            free(current_highest);
            return 1103;
        }
        
        if (combination->max_abundance < current_highest->abundance) {
            combination->max_abundance = current_highest->abundance;
        }
        
        if(current_highest->abundance > current->abundance){
            for (ptrdiff_t v = (ptrdiff_t)(C->amount - 1); v >= 0 ; v--) {
                if ( C->compounds[v].abundance != current_highest->abundance
                    ) {
                    if(C->compounds[v].indicator_iso <= iso_nr_max) {
#if USE_REALLOC == 1
                        if (A->amount == allocated_mem_1_A) {
#if SHOW_DETAILS == 1
                            combination->A_realloc_steps += 1;
#endif
                            allocated_mem_1_A *= ALLOC_FACTOR_1_A;
                            size_t alloc_ = allocated_mem_1_A * sizeof(Compound);
                            A->compounds = (Compound*) realloc(A->compounds, alloc_);
                            if (A->compounds == NULL) {
                                free(C);
                                free(current_highest);
                                free(current);
                                return 11904;
                            }
                        }
#endif
                        A->compounds[A->amount] = C->compounds[v];
                        A->amount++;
#if SHOW_DETAILS == 1
                        if (combination->maxA < A->amount) {
                            combination->maxA = A->amount;
                        }
#endif
                        if (A->amount >= MAX_COMPOUNDS_1_A) {
                            combination->amount = c;
                            free(C);
                            free(current);
                            free(current_highest);
                            return 1101;
                        }
                    }else{
                        if(clean_abundance * C->compounds[v].abundance >= threshold) {
#if USE_REALLOC == 1
                            if (combination->a2_amount == allocated_mem_1_A2) {
#if SHOW_DETAILS == 1
                                combination->A2_realloc_steps += 1;
#endif
                                allocated_mem_1_A2 *= ALLOC_FACTOR_1_A2;
                                size_t alloc_ = allocated_mem_1_A2 * sizeof(Compound);
                                combination->a2_list = (Compound*) realloc(combination->a2_list, alloc_);
                                if (combination->a2_list == NULL) {
                                    free(C);
                                    free(current_highest);
                                    free(current);
                                    return 11905;
                                }
                            }
#endif
                            combination->a2_list[combination->a2_amount] = C->compounds[v];
                            combination->a2_amount++;
#if SHOW_DETAILS == 1
                            if (combination->maxA2 < combination->a2_amount) {
                                combination->maxA2 = combination->a2_amount;
                            }
#endif
                            if (combination->a2_amount >= MAX_COMPOUNDS_1_A2) {
                                combination->amount = c;
                                free(C);
                                free(current);
                                free(current_highest);
                                return 1102;
                            }
                        }
                    }
                }
            }
            
            *current = *current_highest;
            if (clean_abundance * current_highest->abundance >= threshold && c < MAX_COMPOUNDS_1) {
                
#if USE_REALLOC == 1
                if (c == allocated_mem_1) {
#if SHOW_DETAILS == 1
                    combination->realloc_steps += 1;
#endif
                    allocated_mem_1 *= ALLOC_FACTOR_1;
                    size_t alloc_ = allocated_mem_1 * sizeof(Compound);
                    combination->compounds = (Compound*) realloc(combination->compounds, alloc_);
                    if (combination->compounds == NULL) {
                        free(current_highest);
                        free(current);
                        free(C);
                        return 11909;
                    }
                }
#endif
                (combination->compounds + c)->mass = current_highest->mass;
                (combination->compounds + c)->abundance =  current_highest->abundance;
                memmove((combination->compounds + c)->sum, current_highest->sum, sum_size);
                c++;
                
                if (c >= MAX_COMPOUNDS_1) {
                    free(current_highest);
                    free(current);
                    free(C);
                    return 1103;
                }
            }
        }
        else{
            for (ptrdiff_t v = 0; v < C->amount; v++) {
                if ( clean_abundance * C->compounds[v].abundance >= threshold ) {
#if USE_REALLOC == 1
                    if (combination->a2_amount == allocated_mem_1_A2) {
#if SHOW_DETAILS == 1
                        combination->A2_realloc_steps += 1;
#endif
                        allocated_mem_1_A2 *= ALLOC_FACTOR_1_A2;
                        size_t alloc_ = allocated_mem_1_A2 * sizeof(Compound);
                        combination->a2_list = (Compound*) realloc(combination->a2_list, alloc_ );
                        if (combination->a2_list == NULL) {
                            free(C);
                            free(current_highest);
                            free(current);
                            return 11907;
                        }
                    }
#endif
                    combination->a2_list[combination->a2_amount] = C->compounds[v];
                    combination->a2_amount++;
#if SHOW_DETAILS == 1
                    if (combination->maxA2 < combination->a2_amount) {
                        combination->maxA2 = combination->a2_amount;
                    }
#endif
                    if (combination->a2_amount >= MAX_COMPOUNDS_1_A2) {
                        combination->amount = c;
                        free(C);
                        free(current);
                        free(current_highest);
                        return 1102;
                    }
                }
            }
            
            if (A->amount > 0) {
                Compound *a_c = &A->compounds[A->amount - 1];
                *current = *a_c;
                
                if (clean_abundance * a_c->abundance >= threshold && c < MAX_COMPOUNDS_1) {
#if USE_REALLOC == 1
                    if (c == allocated_mem_1) {
#if SHOW_DETAILS == 1
                        combination->realloc_steps += 1;
#endif
                        allocated_mem_1 *= ALLOC_FACTOR_1;
                        size_t alloc_ = allocated_mem_1 * sizeof(Compound);
                        combination->compounds = (Compound*) realloc(combination->compounds, alloc_);
                        if (combination->compounds == NULL) {
                            free(C);
                            free(current_highest);
                            free(current);
                            return 11909;
                        }
                    }
#endif
                    (combination->compounds + c)->mass = a_c->mass;
                    (combination->compounds + c)->abundance =  a_c->abundance;
                    memmove((combination->compounds + c)->sum, a_c->sum, sum_size);
                    c++;
                    
                    if (c >= MAX_COMPOUNDS_1) {
                        free(C);
                        free(current_highest);
                        free(current);
                        return 1103;
                    }
                }
                A->amount--;
            }else if(A->amount == 0 && combination->a2_amount > 0){
                Compound *a2 = &combination->a2_list[combination->a2_amount - 1];
                *current = *a2;
                if ( clean_abundance * a2->abundance >= threshold && c < MAX_COMPOUNDS_1) {
                    
#if USE_REALLOC == 1
                    if (c == allocated_mem_1) {
#if SHOW_DETAILS == 1
                        combination->realloc_steps += 1;
#endif
                        allocated_mem_1 *= ALLOC_FACTOR_1;
                        size_t alloc_ = allocated_mem_1 * sizeof(Compound);
                        combination->compounds = (Compound*) realloc(combination->compounds, alloc_);
                        if (combination->compounds == NULL) {
                            free(C);
                            free(current_highest);
                            free(current);
                            return 11909;
                        }
                    }
#endif
                    (combination->compounds + c)->mass = a2->mass;
                    (combination->compounds + c)->abundance =  a2->abundance;
                    memmove((combination->compounds + c)->sum, a2->sum, sum_size);
                    c++;
                    
                    if (c >= MAX_COMPOUNDS_1) {
                        free(C);
                        free(current_highest);
                        free(current);
                        return 1103;
                    }
                }
                combination->a2_amount--;
            }else {
                current->abundance = -1.0;
            }
        }
    }
    
#if SHOW_DETAILS == 1
#if USE_IN_R == 1
    Rprintf("\n    %s, maxA: %d, maxA2: %d, comb_realloc_steps: %d, A_realloc_steps: %d, A2_realloc_steps: %d, updated_steps: %d, ", element->name, combination->maxA, combination->maxA2, combination->realloc_steps, combination->A_realloc_steps, combination->A2_realloc_steps, element->all_iso_calc_amount);
#else
    printf("\n    %s, maxA: %zu, maxA2: %zu, comb_realloc_steps: %d, A_realloc_steps: %d, A2_realloc_steps: %d, updated_steps: %d, ", element->name, combination->maxA, combination->maxA2, combination->realloc_steps, combination->A_realloc_steps, combination->A2_realloc_steps, element->all_iso_calc_amount);
#endif
#endif
    
    combination->amount = c;
    free(C);
    free(current);
    free(current_highest);
   	return 0;
}
#endif

int combine_combinations_algo_1(Combination_1* combinations,
                                double threshold,
                                unsigned short  element_amount,
                                size_t* peak_amount,
                                int  peak_limit,
                                unsigned short iso_amount,
                                long double max_abundance,
                                int rtm,
                                double **m_,
                                double **a_,
                                int **cc_){
    
    
#if USE_REALLOC == 1
    size_t allocated_mem_pl = ALLOC_START_PL;
#endif
    
    if (element_amount < 1) {
        return 1200;
    }
    size_t tracking[element_amount];
    tracking[0] = 0;
    for (ptrdiff_t i = 0; i < element_amount; i++) {
        tracking[i] = 0;
        qsort((combinations + i)->compounds, (combinations + i)->amount, sizeof(Compound), compound_sort_by_abundance_dec);
    }
    
    const size_t iso_a = iso_amount;
    int sum[iso_a];
    size_t iso_sum_size = sizeof(sum);
    
    double mass = 0.0;
    long double abundance = 1.0;
    size_t v = 0;
    
    while (tracking[0] < (combinations)->amount) {
        
        mass = 0.0;
        abundance = 1.0;
        size_t cc_count = 0;
        unsigned int cc_tmp[iso_a];
        
        for (ptrdiff_t j = 0; j < element_amount; j++) {
            
            const size_t e_a = (combinations + j)->element.iso_amount;
            unsigned int sum_e_iso[e_a];
            size_t e_iso_size = sizeof(sum_e_iso);
            // the last element
            if (j == element_amount - 1) {
                double tmp_mass = 0.0;
                long double tmp_abundance = 1.0;
                for (ptrdiff_t h = 0; h < (combinations + j)->amount; h++) {
                    tmp_mass = mass;
                    tmp_abundance = abundance;
                    tmp_mass += ((combinations + j)->compounds + tracking[j])->mass;
                    tmp_abundance *= ((combinations + j)->compounds + tracking[j])->abundance;
                    if ((100/max_abundance)* tmp_abundance >= threshold) {
                        *(*m_ + v) = tmp_mass;
                        if (rtm == 3 || rtm == 4) {
                            *(*a_ + v) = (double)tmp_abundance;
                        }else{
                            *(*a_ + v) = (double)((100/max_abundance) * tmp_abundance);
                        }
                        
                        memmove(*cc_ + v * iso_a, cc_tmp, iso_sum_size);
                        memmove(*cc_ + v * iso_a + cc_count, combinations[j].compounds[tracking[j]].sum, e_iso_size);
                        
                        v++;
                        if (v >= peak_limit) {
                            *(peak_amount) = v;
                            return 1201;
                        }
#if USE_REALLOC == 1
                        if (v == allocated_mem_pl) {
                            allocated_mem_pl *= ALLOC_FACTOR_PL;
                            *m_ = (double*) realloc(*m_, allocated_mem_pl * sizeof(double));
                            if (*m_ == NULL) {
                                *(peak_amount) = v;
                                return 12901;
                            }
                            *a_ = (double*) realloc(*a_, allocated_mem_pl * sizeof(double));
                            if (*a_ == NULL) {
                                *(peak_amount) = v;
                                return 12902;
                            }
                            *cc_ = (int*)realloc(*cc_, allocated_mem_pl * iso_a * sizeof(int));
                            if (*cc_ == NULL) {
                                *(peak_amount) = v;
                                return 12903;
                            }
                        }
#endif
                        tracking[j]++;
                    }else{
                        tracking[j]++;
                        break;
                    }
                }
            }else{
                if (tracking[j] < (combinations + j)->amount) {
                    mass += ((combinations + j)->compounds + tracking[j])->mass;
                    abundance *= ((combinations + j)->compounds + tracking[j])->abundance;
                    
                    for (ptrdiff_t u = 0; u < (combinations + j)->element.iso_amount; u++) {
                        cc_tmp[cc_count] = ((combinations + j)->compounds + tracking[j])->sum[u];
                        cc_count++;
                        if (cc_count >= iso_a) {
                            *(peak_amount) = v;
                            return 1202;
                        }
                    }
                }else{
                    break;
                }
                
            }
        }
        for (ptrdiff_t k = element_amount - 2; k >= 0; k--) {
            if (tracking[k] < combinations[k].amount) {
                tracking[k]++;
                for (ptrdiff_t l = k + 1; l < element_amount; l++) {
                    tracking[l] = 0;
                }
                break;
            }
        }
    }
    
    *(peak_amount) = v;
    return 0;
}


int calc_pattern_algo_1(Element *elements,
                        size_t  *p_a,
                        double threshold,
                        unsigned short i_a,
                        unsigned short e_a,
                        long double mono_abundance,
                        int  peak_limit,
                        int rtm,
                        double **m_,
                        double **a_,
                        int **cc_
                        )
{

    size_t element_amount = e_a;
    //size_t iso_amount = i_a;
    long double a_max = 1.0;
    
    CombinationMulti_1_A* A = (CombinationMulti_1_A*)calloc(1, sizeof(CombinationMulti_1_A));
    if (A == NULL) {
		return 1900;
	}

	Isotope2** isotopes = (Isotope2**)malloc(element_amount * sizeof(Isotope2*));
    if(isotopes == NULL){
        free(A);
        return 1901;
    }
    Combination_1* combinations = (Combination_1*)calloc(element_amount, sizeof(Combination_1));
    if(combinations == NULL){
        free(A);
        free(isotopes);
        return 1902;
    }
    
#if USE_REALLOC == 1
    Compound* compmult_a;
    if (ALLOC_START_1_A > MAX_COMPOUNDS_1_A) {
		size_t alloc_ = MAX_COMPOUNDS_1_A * sizeof(Compound);
        compmult_a = (Compound*)malloc(alloc_);
    }else{
        size_t alloc_ = ALLOC_START_1_A * sizeof(Compound);
        compmult_a = (Compound*)malloc(alloc_);
    }
    if(compmult_a == NULL){
        free(A);
        free(combinations);
        free(isotopes);
        *p_a = 0;
        return 1903;
    }
    A->compounds = compmult_a;
    
    for (ptrdiff_t  i = 0; i < element_amount; i++) {
        Compound* comp;
        if (ALLOC_START_1 > MAX_COMPOUNDS_1) {
			size_t alloc_ = MAX_COMPOUNDS_1 * sizeof(Compound);
            comp = (Compound*)malloc(alloc_);
        }else{
            size_t alloc_ = ALLOC_START_1 * sizeof(Compound);
            comp = (Compound*)malloc(alloc_);
        }
        if(comp == NULL){
            for (ptrdiff_t  j = 0; j < i; j++) {
                if ((combinations + j)->compounds != NULL) {
                    free((combinations + j)->compounds);
                }
            }
            for (ptrdiff_t  j = 0; j < i; j++) {
                if ((combinations + j)->a2_list != NULL) {
                    free((combinations + j)->a2_list);
                }
            }
            free(A->compounds);
            free(A);
            free(combinations);
            free(isotopes);
            return 1904;
        }
        (combinations + i)->compounds = comp;
        
        Compound* compmult = NULL;
        if (ALLOC_START_1_A2 > MAX_COMPOUNDS_1_A2) {
			size_t alloc_ = MAX_COMPOUNDS_1_A2 * sizeof(Compound);
            compmult = (Compound*)malloc(alloc_);
        }else{
            size_t alloc_ = ALLOC_START_1_A2 * sizeof(Compound);
            compmult = (Compound*)malloc(alloc_);
        }
        
        if(compmult == NULL){
            for (ptrdiff_t  j = 0; j <= i; j++) {
                if ((combinations + j)->compounds != NULL) {
                    free((combinations + j)->compounds);
                }
            }
            for (ptrdiff_t  j = 0; j < i; j++) {
                if ((combinations + j)->a2_list != NULL) {
                    free((combinations + j)->a2_list);
                }
            }

            free(A->compounds);
            free(A);
            free(combinations);
            free(isotopes);
            return 1905;
        }
        (combinations + i)->a2_list = compmult;
    }
#endif
    
    for (ptrdiff_t  i = 0; i < element_amount; i++) {
		isotopes[i] = (Isotope2*)malloc((elements + i)->iso_amount * sizeof(Isotope2));
        if (isotopes[i] == NULL) {
#if USE_REALLOC == 1
            for (ptrdiff_t  i = 0; i < element_amount; i++) {
                free((combinations + i)->compounds);
                free((combinations + i)->a2_list);
            }
            free(A->compounds);
#endif
            for (ptrdiff_t  j = 0; j < i; j++) {
                free(isotopes[j]);
            }
            free(combinations);
            free(A);
            free(isotopes);
            return 1906;
        }
		create_isotope_list_single(*(elements + i), isotopes[i]);
        int msg = calc_combination_max_abundance((combinations + i), (elements + i), isotopes[i], threshold, A, mono_abundance, rtm);
        if (msg) {
#if USE_REALLOC == 1
            for (ptrdiff_t  j = 0; j < element_amount; j++) {
                free((combinations + j)->compounds);
                free((combinations + j)->a2_list);
            }
            free(A->compounds);
#endif
            for (ptrdiff_t  j = 0; j <= i; j++) {
				free(isotopes[j]);
            }
            free(combinations);
            free(A);
            free(isotopes);
            return msg;
        }
        a_max *= (combinations + i)->max_abundance;
    }
    
    if(rtm == 1 || rtm == 4){
        long double clean_abundance = 1.0;
        for (ptrdiff_t  j = 0; j < element_amount; j++) {
            long double combination_mono_abundance = 1.0;
            for(ptrdiff_t k = 0; k < (elements+j)->amount; k++){
                combination_mono_abundance *= (elements + j)->isotopes[0].abundance;
            }
            clean_abundance = a_max/mono_abundance;// combination_mono_abundance;

            int result = create_combination_algo_1(combinations + j, elements + j, isotopes[j], clean_abundance, (threshold * mono_abundance)/100.0, peak_limit,A);
            if(result){
#if USE_REALLOC == 1
                for (ptrdiff_t  i = 0; i < element_amount; i++) {
                    free((combinations + i)->compounds);
                    free((combinations + i)->a2_list);
                }
                free(A->compounds);
#endif
                
                for (ptrdiff_t  i = 0; i < element_amount; i++) {
                    free(isotopes[i]);
                }
                free(combinations);
                free(A);
                free(isotopes);
                return result;
            }
        }
        a_max = mono_abundance;
        
    }else if(rtm == 0 || rtm == 3){
        for (ptrdiff_t  j = 0; j < element_amount; j++) {
            long double clean_abundance = a_max/(combinations + j)->max_abundance;
            int result = create_combination_algo_1(combinations + j, elements + j, isotopes[j], clean_abundance, (threshold * a_max)/100.0, peak_limit,A);
            if(result){
#if USE_REALLOC == 1
                for (ptrdiff_t  i = 0; i < element_amount; i++) {
                    free((combinations + i)->compounds);
                    free((combinations + i)->a2_list);
                }
                free(A->compounds);
#endif
                for (ptrdiff_t  i = 0; i < element_amount; i++) {
					free(isotopes[i]);
				}
                free(combinations);
                free(A);
                free(isotopes);
                return result;
            }
        }
    }else if(rtm == 2){
        for (ptrdiff_t  j = 0; j < element_amount; j++) {
            long double clean_abundance = a_max/(combinations + j)->max_abundance;
            int result = create_combination_algo_1(combinations + j, elements + j, isotopes[j], clean_abundance, threshold, peak_limit,A);
            if(result){
#if USE_REALLOC == 1
                for (ptrdiff_t  i = 0; i < element_amount; i++) {
                    free((combinations + i)->compounds);
                    free((combinations + i)->a2_list);
                }
                free(A->compounds);
#endif
				for (ptrdiff_t  i = 0; i < element_amount; i++) {
					free(isotopes[i]);
				}
				free(combinations);
                free(A);
                free(isotopes);
                return result;
            }
        }
        a_max = 100;
        
    }else{
#if USE_REALLOC == 1
        for (ptrdiff_t  i = 0; i < element_amount; i++) {
            free((combinations + i)->compounds);
            free((combinations + i)->a2_list);
        }
        free(A->compounds);
#endif
		for (ptrdiff_t  i = 0; i < element_amount; i++) {
			free(isotopes[i]);
		}
        free(combinations);
        free(A);
        free(isotopes);
        return 22;
    }
    
    
#if USE_REALLOC == 1
    for (ptrdiff_t  i = 0; i < element_amount; i++) {
        free((combinations + i)->a2_list);
    }
    free(A->compounds);
#endif
    free(A);
    for (ptrdiff_t  i = 0; i < element_amount; i++) {
        free(isotopes[i]);
    }
    free(isotopes);
    
    int err_combination = combine_combinations_algo_1(combinations, threshold, e_a, p_a, peak_limit, i_a, a_max,rtm, m_, a_, cc_);
    if(err_combination){
#if USE_REALLOC == 1
        for (ptrdiff_t  i = 0; i < element_amount; i++) {
            free((combinations + i)->compounds);
        }
#endif
        free(combinations);
        return err_combination;
    }
    
#if USE_REALLOC == 1
    for (ptrdiff_t  i = 0; i < element_amount; i++) {
        free((combinations + i)->compounds);
    }
#endif
    free(combinations);
    return 0;
}

// algo 2 ////////////////////////////////////////////////////////////////////////////////////////////////////////
int calc_pattern_algo_2(
                        long double* max_a,
                        Element *elements,
                        unsigned short element_amount,
                        double threshold,
                        size_t * peak_amount,
                        int peak_limit,
                        int *iso_amount_stats,
                        int rtm,
                        double **m_,
                        double **a_,
                        int **cc_
                        )
{
#if USE_REALLOC == 1
    size_t allocated_mem_pl = ALLOC_START_PL;
#endif
    
#if SHOW_DETAILS == 1
    size_t max_a2_list_2 = 0;
    size_t max_a_list_2 = 0;
    int A_realloc_steps = 0;
    int A2_realloc_steps = 0;
#endif
    
    Isotope2 *isotopes = (Isotope2*)calloc(MAX_ISO_SIZE, sizeof(Isotope2));
    if (isotopes == NULL) {
        return 2901;
    }
    CompoundMulti* monoisotopic = (CompoundMulti*)calloc(1,sizeof(CompoundMulti));
    if(monoisotopic == NULL){
        free(isotopes);
        return 2902;
    }
    CombinationMulti_2_A* A = (CombinationMulti_2_A*)calloc(1,sizeof(CombinationMulti_2_A));
    if(A == NULL){
        free(isotopes);
        free(monoisotopic);
        return 2903;
    }
    CombinationMulti* A2 = (CombinationMulti*)calloc(1,sizeof(CombinationMulti));
    if(A2 == NULL){
        free(isotopes);
        free(monoisotopic);
        free(A);
        return 2904;
    }
    CombinationMulti_2_C* C = (CombinationMulti_2_C*)calloc(1,sizeof(CombinationMulti_2_C));
    if(C == NULL){
        free(isotopes);
        free(monoisotopic);
        free(A);
        free(A2);
        return 2905;
    }
    CompoundMulti* current_highest = (CompoundMulti*)calloc(1,sizeof(CompoundMulti));
    if(current_highest == NULL){
        free(isotopes);
        free(monoisotopic);
        free(A);
        free(A2);
        free(C);
        return 2906;
    }
    CompoundMulti* current = (CompoundMulti*)calloc(1,sizeof(CompoundMulti));
    if(current == NULL){
        free(isotopes);
        free(monoisotopic);
        free(A);
        free(A2);
        free(C);
        free(current_highest);
        return 2907;
    }
#if USE_REALLOC == 1
    size_t allocated_mem_2_A = ALLOC_START_2_A;
    size_t allocated_mem_2_A2 = ALLOC_START_2_A2;
    
    CompoundMulti* compmult;
    if (ALLOC_START_2_A2 > MAX_COMPOUNDS_2_A2) {
        size_t alloc_ = MAX_COMPOUNDS_2_A2 * sizeof(CompoundMulti);
        compmult = (CompoundMulti*)malloc(alloc_);
    }else{
        size_t alloc_ = ALLOC_START_2_A2 * sizeof(CompoundMulti);
        compmult = (CompoundMulti*)malloc(alloc_);
    }
    if(compmult == NULL){
        free(isotopes);
        free(monoisotopic);
        free(A);
        free(A2);
        free(C);
        free(current_highest);
        free(current);
        return 2908;
    }
    A2->compounds = compmult;
    
    CompoundMulti* compmult_a;
    if (ALLOC_START_2_A > MAX_COMPOUNDS_2_A) {
        size_t alloc_ = MAX_COMPOUNDS_2_A * sizeof(CompoundMulti);
        compmult_a = (CompoundMulti*)malloc(alloc_);
    }else{
        size_t alloc_ = ALLOC_START_2_A * sizeof(CompoundMulti);
        compmult_a = (CompoundMulti*)malloc(alloc_);
    }
    if(compmult_a == NULL){
        free(A2->compounds);
        free(isotopes);
        free(monoisotopic);
        free(A);
        free(A2);
        free(C);
        free(current_highest);
        free(current);
        return 2909;
    }
    A->compounds = compmult_a;
#endif
    
    unsigned short iso_c = 0;
    
    create_isotope_list(elements, element_amount, isotopes, &iso_c);
    calc_monoisotopic(elements, element_amount, monoisotopic);
    
    A->amount = 0;
    A->max_abundance = 1.0;
    A2->amount = 0;
    A2->max_abundance = 1.0;
    C->amount = 0;
    C->max_abundance = 1.0;
    
    unsigned short iso_amount = iso_c + element_amount;
    if (iso_amount < 1) {
#if USE_REALLOC == 1
        free(A->compounds);
        free(A2->compounds);
#endif
        free(A);
        free(A2);
        free(C);
        free(isotopes);
        free(monoisotopic);
        free(current);
        free(current_highest);
        return 2004;
    }

    size_t  iso_nr_max = 0;
    size_t c = 0;
    int  c_stats = 0;
    size_t  h = 0;

    unsigned int sum_s[iso_amount];
    size_t size_sum = sizeof(sum_s);

    **m_ = monoisotopic->mass;
    **a_ = (double)monoisotopic->abundance;
    *max_a = monoisotopic->abundance;
    memmove(*cc_, monoisotopic->sum, size_sum);
    c++;
    
    *current = *monoisotopic;
    size_t  iso_pos[MAX_ELEMENTS];
    for (unsigned short  d = 0; d < element_amount; d++) {
        iso_pos[d] = 0;
        for (ptrdiff_t  b = 0; b < d; b++) {
            iso_pos[d] += (elements + b)->iso_amount;
        }
    }
    
    while (current->abundance != -1.0) {
        *current_highest = *current;
        C->amount = 0;
        iso_nr_max = 0;
        
        for (unsigned short j = current->indicator_iso; j < iso_c; j++) {
            h = (isotopes + j)->element_nr;
            
            if ( current->counter[h] < (elements + h)->amount ) {
                
                Isotope *isotope = (elements + h)->isotopes;
                size_t  iso_e_nr = (isotopes + j)->iso_e_nr;
                C->compounds[C->amount] = *current;
                CompoundMulti *comp = &C->compounds[C->amount];
                
                comp->counter[h]++;
                comp->indicator_iso = j;
                comp->sum[iso_pos[h] + iso_e_nr]++;
                comp->mass -= isotope->mass;
                comp->mass += (isotope + iso_e_nr)->mass;
                comp->abundance *= ( isotope + iso_e_nr)->abundance * (comp->sum[iso_pos[h]]);
                comp->abundance /= ( isotope )->abundance * comp->sum[iso_pos[h] + iso_e_nr];
                comp->sum[iso_pos[h]]--;
                
                
                //printf("\n%zu  %lf  max_abundance: %lf", c, comp->abundance, max_a);
                elements->all_iso_calc_amount++;
                if (current_highest->abundance < comp->abundance) {
                    *current_highest = *comp;
                }
                if (comp->abundance >= current->abundance) {
                    iso_nr_max = j;
                }
                c_stats++;
                if (c_stats > MAX_ITERATION_2) {
#if USE_REALLOC == 1
                    free(A->compounds);
                    free(A2->compounds);
#endif
                    free(A);
                    free(A2);
                    free(C);
                    free(isotopes);
                    free(monoisotopic);
                    free(current);
                    free(current_highest);
                    return 2006;
                }
                C->amount++;
            }
        }
        
        if(current_highest->abundance > current->abundance){
            if (*max_a < current_highest->abundance && (rtm == 0 || rtm == 3)) {
                *max_a = current_highest->abundance;
            }else if(rtm == 2){
                *max_a = 100;
            }
            
            for (ptrdiff_t v = (ptrdiff_t)C->amount - 1; v >= 0 ; v--) {
                if ( C->compounds[v].abundance != current_highest->abundance
                    ) {
                    if(C->compounds[v].indicator_iso <= iso_nr_max) {
                        A->compounds[A->amount] = C->compounds[v];
                        A->amount++;
#if SHOW_DETAILS == 1
                        if (max_a_list_2 < A->amount) {
                            max_a_list_2 = A->amount;
                        }
#endif
                        if (A->amount >= MAX_COMPOUNDS_2_A) {
                            *peak_amount = c;
#if USE_REALLOC == 1
                            free(A->compounds);
                            free(A2->compounds);
#endif
                            free(A);
                            free(A2);
                            free(C);
                            free(isotopes);
                            free(monoisotopic);
                            free(current);
                            free(current_highest);
                            return 2001;
                        }
#if USE_REALLOC == 1
                            size_t a_a = A->amount;
                            if (a_a >= allocated_mem_2_A) {
#if SHOW_DETAILS == 1
                                A_realloc_steps++;
#endif
                                allocated_mem_2_A *= ALLOC_FACTOR_2_A;
                                size_t alloc_ = allocated_mem_2_A * sizeof(CompoundMulti);
                                A->compounds = (CompoundMulti*) realloc(A->compounds, alloc_);
                                if (A->compounds == NULL) {
                                    free(A->compounds);
                                    free(A2->compounds);
                                    free(A);
                                    free(A2);
                                    free(C);
                                    free(isotopes);
                                    free(monoisotopic);
                                    free(current);
                                    free(current_highest);
                                    return 2910;
                                }
                            }
#endif
                    }else{
                        if((100/ *max_a) * C->compounds[v].abundance >= threshold) {
                            A2->compounds[A2->amount] = C->compounds[v];
                            A2->amount++;
                            
#if SHOW_DETAILS == 1
                            if (max_a2_list_2 < A2->amount) {
                                max_a2_list_2 = A2->amount;
                            }
#endif
                            
                            if (A2->amount >= MAX_COMPOUNDS_2_A2) {
                                *peak_amount = c;
#if USE_REALLOC == 1
                                free(A->compounds);
                                free(A2->compounds);
#endif
                                free(A);
                                free(A2);
                                free(C);
                                free(isotopes);
                                free(monoisotopic);
                                free(current);
                                free(current_highest);
                                return 2002;
                            }
#if USE_REALLOC == 1
                            size_t a2_a = A2->amount;
                            if (a2_a >= allocated_mem_2_A2) {
#if SHOW_DETAILS == 1
                                A2_realloc_steps++;
#endif
                                allocated_mem_2_A2 *= ALLOC_FACTOR_2_A2;
                                size_t alloc_ = allocated_mem_2_A2 * sizeof(CompoundMulti);
                                A2->compounds = (CompoundMulti*) realloc(A2->compounds, alloc_);
                                if (A2->compounds == NULL) {
                                    free(A->compounds);
                                    free(A2->compounds);
                                    free(A);
                                    free(A2);
                                    free(C);
                                    free(isotopes);
                                    free(monoisotopic);
                                    free(current);
                                    free(current_highest);
                                    return 2911;
                                }
                            }
#endif
                        }
                    }
                }else{
                    C->compounds[C->amount] = C->compounds[v];
                    *current_highest = C->compounds[C->amount];
                }
            }
            
            *current = *current_highest;
            
            if ((100/ *max_a) * current_highest->abundance >= threshold) {
                *(*m_ + c) = current_highest->mass;
                *(*a_ + c) = (double)current_highest->abundance;
                memmove((*cc_ + c * iso_amount), current_highest->sum, size_sum);
                c++;
                
                if (c >= peak_limit) {
                    *peak_amount = c;
#if USE_REALLOC == 1
                    free(A->compounds);
                    free(A2->compounds);
#endif
                    free(A);
                    free(A2);
                    free(C);
                    free(isotopes);
                    free(monoisotopic);
                    free(current);
                    free(current_highest);
                    return 2003;
                }
#if USE_REALLOC == 1
                if (c >= allocated_mem_pl) {
                    allocated_mem_pl *= ALLOC_FACTOR_PL;
                    size_t m_n =  allocated_mem_pl;
                    *m_ = (double*) realloc(*m_, m_n * sizeof(double));
                    if (*m_ == NULL) {
                        free(A->compounds);
                        free(A2->compounds);
                        free(A);
                        free(A2);
                        free(C);
                        free(isotopes);
                        free(monoisotopic);
                        free(current);
                        free(current_highest);
                        return 2912;
                    }
                    *a_ = (double*) realloc(*a_, m_n * sizeof(double));
                    if (*a_ == NULL ) {
                        free(A->compounds);
                        free(A2->compounds);
                        free(A);
                        free(A2);
                        free(C);
                        free(isotopes);
                        free(monoisotopic);
                        free(current);
                        free(current_highest);
                        return 2913;
                    }
                    *cc_ = (int*)realloc(*cc_, m_n * iso_amount * sizeof(int));
                    if (*cc_ == NULL ) {
                        free(A->compounds);
                        free(A2->compounds);
                        free(A);
                        free(A2);
                        free(C);
                        free(isotopes);
                        free(monoisotopic);
                        free(current);
                        free(current_highest);
                        return 2914;
                    }
                }
#endif
            }
        }
        else{
            for (ptrdiff_t v = 0; v < C->amount; v++) {
                if ( (100/ *max_a) * C->compounds[v].abundance >= threshold ) {
                    A2->compounds[A2->amount] = C->compounds[v];
                    A2->amount++;
                    
#if SHOW_DETAILS == 1
                    if (max_a2_list_2 < A2->amount) {
                        max_a2_list_2 = A2->amount;
                    }
#endif
                    
                    if (A2->amount >= MAX_COMPOUNDS_2_A2) {
                        *peak_amount = c;
#if USE_REALLOC == 1
                        free(A->compounds);
                        free(A2->compounds);
#endif
                        free(A);
                        free(A2);
                        free(C);
                        free(isotopes);
                        free(monoisotopic);
                        free(current);
                        free(current_highest);
                        return 2002;
                    }
#if USE_REALLOC == 1
                    size_t a2_a = A2->amount;
                    if (a2_a >= allocated_mem_2_A2) {
#if SHOW_DETAILS == 1
                        A2_realloc_steps++;
#endif
                        allocated_mem_2_A2 *= ALLOC_FACTOR_2_A2;
                        size_t alloc_ = allocated_mem_2_A2 * sizeof(CompoundMulti);
                        A2->compounds = (CompoundMulti*) realloc(A2->compounds, alloc_);
                        if (A2->compounds == NULL) {
                            free(A->compounds);
                            free(A2->compounds);
                            free(A);
                            free(A2);
                            free(C);
                            free(isotopes);
                            free(monoisotopic);
                            free(current);
                            free(current_highest);
                            return 2915;
                        }
                    }
#endif
                }
            }
            
            if (A->amount > 0) {
                CompoundMulti *a_c = &A->compounds[A->amount - 1];
                *current = *a_c;
                
                if ((100/ *max_a) * a_c->abundance >= threshold) {
                    *(*m_ + c) = a_c->mass;
                    *(*a_ + c) = (double)a_c->abundance;
                    memmove((*cc_ + c * iso_amount), a_c->sum, size_sum);
                    c++;
                    
                    if (c >= peak_limit) {
                        *peak_amount = c;
#if USE_REALLOC == 1
                        free(A->compounds);
                        free(A2->compounds);
#endif
                        free(A);
                        free(A2);
                        free(C);
                        free(isotopes);
                        free(monoisotopic);
                        free(current);
                        free(current_highest);
                        return 2003;
                    }
#if USE_REALLOC == 1
                    if (c >= allocated_mem_pl) {
                        allocated_mem_pl *= ALLOC_FACTOR_PL;
                        size_t m_n =  allocated_mem_pl;
                        *m_ = (double*) realloc(*m_, m_n * sizeof(double));
                        if (*m_ == NULL) {
                            free(A->compounds);
                            free(A2->compounds);
                            free(A);
                            free(A2);
                            free(C);
                            free(isotopes);
                            free(monoisotopic);
                            free(current);
                            free(current_highest);
                            return 2912;
                        }
                        *a_ = (double*) realloc(*a_, m_n * sizeof(double));
                        if (*a_ == NULL ) {
                            free(A->compounds);
                            free(A2->compounds);
                            free(A);
                            free(A2);
                            free(C);
                            free(isotopes);
                            free(monoisotopic);
                            free(current);
                            free(current_highest);
                            return 2913;
                        }
                        *cc_ = (int*)realloc(*cc_, m_n * iso_amount * sizeof(int));
                        if (*cc_ == NULL ) {
                            free(A->compounds);
                            free(A2->compounds);
                            free(A);
                            free(A2);
                            free(C);
                            free(isotopes);
                            free(monoisotopic);
                            free(current);
                            free(current_highest);
                            return 2914;
                        }
                    }
#endif
                }
                A->amount--;
            }else if(A->amount == 0 && A2->amount > 0){
                CompoundMulti *a2 = &A2->compounds[A2->amount - 1];
                *current = *a2;
                
                if ( (100/ *max_a) * a2->abundance >= threshold ) {
                    *(*m_ + c) = a2->mass;
                    *(*a_ + c) = (double)a2->abundance;
                    memmove((*cc_ + c * iso_amount), a2->sum, size_sum);
                    c++;
                    
                    if (c >= peak_limit) {
                        *peak_amount = c;
#if USE_REALLOC == 1
                        free(A->compounds);
                        free(A2->compounds);
#endif
                        free(A);
                        free(A2);
                        free(C);
                        free(isotopes);
                        free(monoisotopic);
                        free(current);
                        free(current_highest);
                        return 2003;
                    }
#if USE_REALLOC == 1
                    if (c >= allocated_mem_pl) {
                        allocated_mem_pl *= ALLOC_FACTOR_PL;
                        size_t m_n =  allocated_mem_pl;
                        *m_ = (double*) realloc(*m_, m_n * sizeof(double));
                        if (*m_ == NULL) {
                            free(A->compounds);
                            free(A2->compounds);
                            free(A);
                            free(A2);
                            free(C);
                            free(isotopes);
                            free(monoisotopic);
                            free(current);
                            free(current_highest);
                            return 2912;
                        }
                        *a_ = (double*) realloc(*a_, m_n * sizeof(double));
                        if (*a_ == NULL ) {
                            free(A->compounds);
                            free(A2->compounds);
                            free(A);
                            free(A2);
                            free(C);
                            free(isotopes);
                            free(monoisotopic);
                            free(current);
                            free(current_highest);
                            return 2913;
                        }
                        *cc_ = (int*)realloc(*cc_, m_n * iso_amount * sizeof(int));
                        if (*cc_ == NULL ) {
                            free(A->compounds);
                            free(A2->compounds);
                            free(A);
                            free(A2);
                            free(C);
                            free(isotopes);
                            free(monoisotopic);
                            free(current);
                            free(current_highest);
                            return 2914;
                        }
                    }
#endif
                }
                A2->amount--;
                
            }else {
                current->abundance = -1.0;
            }
        }
    }
    *peak_amount =  c;
    *iso_amount_stats = c_stats;
    
#if SHOW_DETAILS == 1
#if USE_IN_R == 1
    Rprintf("\n    maxA: %zu, maxA2: %zu, A_realloc_steps: %d, A2_realloc_steps: %d, ", max_a_list_2, max_a2_list_2, A_realloc_steps,A2_realloc_steps);
#else
    printf("\n    maxA: %zu, maxA2: %zu, A_realloc_steps: %d, A2_realloc_steps: %d, ", max_a_list_2, max_a2_list_2, A_realloc_steps,A2_realloc_steps);

#endif
#endif

#if USE_REALLOC == 1
    free(A->compounds);
    free(A2->compounds);
#endif
    free(A);
    free(A2);
    free(C);
    free(isotopes);
    free(monoisotopic);
    free(current);
    free(current_highest);
   	return 0;
}


// algo4 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int calc_pattern_algo_4(Element *elements,
                        size_t  *peak_amount,
                        double t,
                        unsigned short iso_amount,
                        unsigned short element_amount,
                        long double a_monoisotopic,
                        int p_l,
                        char* l_n,
                        int rtm,
                        double **m_,
                        double **a_,
                        int **cc_,
                        long double *max_a)
{
//    for (int i = 1; i < 2/*ceil(element_amount/ALGO4_BLOCK_SIZE)*/; i++) {
//        Element* tmp_element = (Element*)calloc(1, sizeof(Element));
//        *tmp_element = *(elements + i);
//        elements[i] = elements[element_amount - i];
//        elements[element_amount - i] = *tmp_element;
//    }
    
    long double max_abundance = 1.0;
    
    unsigned short container_factor = ALGO4_BLOCK_SIZE;
    if (ALGO4_BLOCK_SIZE <= 0) {
        container_factor = 1;
    }
    
    if (ALGO4_BLOCK_SIZE > element_amount) {
        container_factor = element_amount;
    }
    
    int sum_amt = 0;
    
    for (int i = 0; i < element_amount; i++) {
        sum_amt += elements[i].amount;
    }
    
    unsigned short container_border = 0;
    unsigned short container_amount = 0;
    unsigned short containter_amount_tmp = 0;
    
    //if (sum_amt < 1000 && element_amount <= 5) {
//    if (0) {
//        // puts all elements in a group of ALGO4_BLOCK_SIZE into a container
//        
//        double div = (double)element_amount/(double)container_factor;
//        unsigned short container_amount_orig = (unsigned short)ceil(div);
//        container_amount = 1;
//        
//        if (container_amount_orig > container_amount) {
//            container_amount = container_amount_orig;
//        }
//    }else{
    
        // take the first elements with 3 or less isotopes into a single element container, the other elements into ALGO4_BLOCK_SIZE containers
        container_border = 3;
        for (int i = 0; i < element_amount; i++) {
            if (elements[i].iso_amount <= container_border) {
                container_amount++;
            }else{
                containter_amount_tmp++;
            }
        }
        double div = (double)containter_amount_tmp/(double)container_factor;
        container_amount += (unsigned short)ceil(div);
//    }
        
    Combination_4* combinations = (Combination_4*)malloc(sizeof(Combination_4[container_amount]));
    
    if(combinations == NULL){
        return 4901;
    }
    Isotope2 *isotopes = (Isotope2*)calloc(MAX_ISO_SIZE, sizeof(Isotope2));
    if (isotopes == NULL) {
        free(combinations);
        return 4902;
    }
    CombinationMulti_2_A* A = (CombinationMulti_2_A*)calloc(1,sizeof(CombinationMulti_2_A));
    if(A == NULL){
        free(isotopes);
        free(combinations);
        return 4903;
    }
    CombinationMulti* A2 = (CombinationMulti*)calloc(1,sizeof(CombinationMulti));
    if(A2 == NULL){
        free(isotopes);
        free(A);
        free(combinations);
        return 4904;
    }
    CombinationMulti_2_C* C = (CombinationMulti_2_C*)calloc(1,sizeof(CombinationMulti_2_C));
    if(C == NULL){
        free(isotopes);
        free(A);
        free(A2);
        free(combinations);
        return 4905;
    }
    
#if USE_REALLOC == 1
    
    CompoundMulti* compmult;
    if (ALLOC_START_4_A2 > MAX_COMPOUNDS_4_A2) {
        size_t alloc_ = MAX_COMPOUNDS_4_A2 * sizeof(CompoundMulti);
        compmult = (CompoundMulti*)malloc(alloc_);
    }else{
        size_t alloc_ = ALLOC_START_4_A2 * sizeof(CompoundMulti);
        compmult = (CompoundMulti*)malloc(alloc_);
    }
    if(compmult == NULL){
        free(isotopes);
        free(A);
        free(A2);
        free(C);
        free(combinations);
        return 4906;
    }
    A2->compounds = compmult;
    
    CompoundMulti* compmult_a;
    if (ALLOC_START_4_A > MAX_COMPOUNDS_4_A) {
        size_t alloc_ = MAX_COMPOUNDS_4_A * sizeof(CompoundMulti);
        compmult_a = (CompoundMulti*)malloc(alloc_);
    }else{
        size_t alloc_ = ALLOC_START_4_A * sizeof(CompoundMulti);
        compmult_a = (CompoundMulti*)malloc(alloc_);
    }
    if(compmult_a == NULL){
        free(A2->compounds);
        free(isotopes);
        free(A);
        free(A2);
        free(C);
        free(combinations);
        return 4907;
    }
    A->compounds = compmult_a;
    
    for (ptrdiff_t  i = 0; i < container_amount; i++) {
        Compound* comp;
        if (ALLOC_START_4 > MAX_COMPOUNDS_4) {
            comp = (Compound*)calloc(MAX_COMPOUNDS_4, sizeof(Compound));
            
        }else{
            comp = (Compound*)calloc(ALLOC_START_4, sizeof(Compound));
        }
        
        if(comp == NULL){
            for (ptrdiff_t  j = 0; j < i; j++) {
                free((combinations + j)->compounds);
            }
            free(A->compounds);
            free(A2->compounds);
            free(isotopes);
            free(A);
            free(A2);
            free(combinations);
            free(C);
            return 4908;
        }
        (combinations + i)->compounds = comp;
    }
    
#endif
    CompoundMulti* monoisotopic = (CompoundMulti*)calloc(1,sizeof(CompoundMulti));
    calc_monoisotopic(elements, element_amount, monoisotopic);
    a_monoisotopic = monoisotopic->abundance;
    free(monoisotopic);
    unsigned short element_index = 0;
    for (ptrdiff_t  j = 0; j < container_amount; j++) {
        
        int msg = 0;
        
        if (elements[element_index].iso_amount > container_border) {
            container_factor = ALGO4_BLOCK_SIZE;
        }else{
            container_factor = 1;
        }
        
        if (j == container_amount - 1) {
            unsigned short element_range = element_amount - element_index;
            msg = create_combination_algo_4((combinations + j), (elements + element_index), element_range, isotopes, A, A2, C, t, rtm, a_monoisotopic);
            if (element_range == 0) {
                container_amount--;
                break;
            }

        }else{
            unsigned short element_range = container_factor;
            if (element_index + container_factor >= element_amount) {
                element_range = element_amount - element_index;
            }

            if (element_range == 0) {
                container_amount--;
                break;
            }

            msg = create_combination_algo_4((combinations + j), (elements + element_index), element_range, isotopes, A, A2, C, t, rtm, a_monoisotopic);
        }
        
        if (msg) {
#if USE_REALLOC == 1
            for (ptrdiff_t  i = 0; i < container_amount; i++) {
                free((combinations + i)->compounds);
            }
            free(A->compounds);
            free(A2->compounds);
#endif
            free(isotopes);
            free(A);
            free(A2);
            free(combinations);
            free(C);
            return msg;
        }
        max_abundance *= (combinations + j)->max_abundance;
        element_index += container_factor;
    }
    
    *max_a = max_abundance;
    
#if USE_REALLOC == 1
    free(A->compounds);
    free(A2->compounds);
#endif
    free(isotopes);
    free(A);
    free(A2);
    free(C);
    
    if(rtm == 1 || rtm == 4){
        max_abundance = a_monoisotopic;
    }else if(rtm == 2){
        max_abundance = 100;
    }
    
    if (container_amount == 0) {
#if USE_REALLOC == 1
        for (ptrdiff_t  i = 0; i < container_amount; i++) {
            free((combinations + i)->compounds);
        }
#endif
        free(combinations);
        return 4004;
    }
    
    if (rtm != 2) {
        clean_combination_algo_4(combinations, t * max_abundance / 100 , container_amount);
    }else{
        clean_combination_algo_4(combinations, t , container_amount);
    }
    
    int msg = combine_combinations_algo_4(combinations, t, container_amount, peak_amount, p_l, iso_amount, max_abundance, rtm, m_, a_, cc_);
    if (msg) {
#if USE_REALLOC == 1
        for (ptrdiff_t  i = 0; i < container_amount; i++) {
            free((combinations + i)->compounds);
        }
#endif
        free(combinations);
        return msg;
    }
    
    
#if USE_REALLOC == 1
    for (ptrdiff_t  i = 0; i < container_amount; i++) {
        free((combinations + i)->compounds);
    }
#endif
    free(combinations);
    return 0;
}

// algo 2 ////////////////////////////////////////////////////////////////////////////////////////////////////////
int create_combination_algo_4(
                        Combination_4* combination
                        ,Element *elements
                        ,unsigned short element_amount
                        ,Isotope2 *isotopes
                        ,CombinationMulti_2_A* A
                        ,CombinationMulti* A2
                        ,CombinationMulti_2_C* C
                        ,double threshold
                        ,int rtm
                        ,long double mono_a
                        )
{
#if SHOW_DETAILS == 1
    size_t max_a2_list_4 = 0;
    size_t max_a_list_4 = 0;
    int A_realloc_steps = 0;
    int A2_realloc_steps = 0;
    int realloc_steps_4 = 0;
    
#endif
    
    CompoundMulti* monoisotopic = (CompoundMulti*)calloc(1,sizeof(CompoundMulti));
    if(monoisotopic == NULL){
        return 4909;
    }
    CompoundMulti* current = (CompoundMulti*)calloc(1,sizeof(CompoundMulti));
    if(current == NULL){
        free(monoisotopic);
        return 4911;
    }
#if USE_REALLOC == 1
    size_t allocated_mem_4 = ALLOC_START_4;
    if (ALLOC_START_4 > MAX_COMPOUNDS_4) {
        allocated_mem_4 = MAX_COMPOUNDS_4;
    }
    size_t allocated_mem_4_A = ALLOC_START_4_A;
    size_t allocated_mem_4_A2 = ALLOC_START_4_A2;
    
#endif
    
    unsigned short iso_c = 0;
    create_isotope_list(elements, element_amount, isotopes, &iso_c);
    calc_monoisotopic(elements, element_amount, monoisotopic);
    
    *current = *monoisotopic;
    
    A->amount = 0;
    A->max_abundance = 1.0;
    A2->amount = 0;
    A2->max_abundance = 1.0;
    C->amount = 0;
    C->max_abundance = 1.0;
    
    unsigned short iso_amount = iso_c + element_amount;
    if (iso_amount < 1) {
        free(current);
        free(monoisotopic);
        return 4005;
    }
    
    size_t iso_nr_max = 0;
    size_t c = 0;
    size_t h = 0;
    
    unsigned int sum_s[iso_amount];
    size_t size_sum = sizeof(sum_s);
    
    combination->compounds->mass = current->mass;
    combination->compounds->abundance = current->abundance;
    memmove(combination->compounds->sum, current->sum, size_sum);
    c++;
    elements->all_iso_calc_amount++;
    
    long double max_a = mono_a;// monoisotopic->abundance;
    long double comb_max_a = max_a;
    
    size_t  iso_pos[MAX_ELEMENTS];
    for (unsigned short  d = 0; d < element_amount; d++) {
        iso_pos[d] = 0;
        for (ptrdiff_t  b = 0; b < d; b++) {
            iso_pos[d] += (elements + b)->iso_amount;
        }
    }
    
    CompoundMulti *current_highest = NULL;
    
    while (current->abundance != -1.0) {
        current_highest = current;
        C->amount = 0;
        iso_nr_max = 0;
        
        for (unsigned short j = current->indicator_iso; j < iso_c; j++) {
            h = (isotopes + j)->element_nr;

            
            if ( current->counter[h] < (elements + h)->amount ) {

                Isotope *isotope = (elements + h)->isotopes;
                size_t  iso_e_nr = (isotopes + j)->iso_e_nr;
 
                C->compounds[C->amount] = *current;
                CompoundMulti *comp = &C->compounds[C->amount];
                
                comp->counter[h]++;
                comp->indicator_iso = j;
                comp->sum[iso_pos[h] + iso_e_nr]++;
                comp->mass -= isotope->mass;
                comp->mass += (isotope + iso_e_nr)->mass;
                comp->abundance *= ( isotope + iso_e_nr)->abundance * (comp->sum[iso_pos[h]]);
                comp->abundance /= ( isotope )->abundance * comp->sum[iso_pos[h] + iso_e_nr];
                comp->sum[iso_pos[h]]--;
                
                if (current_highest->abundance < comp->abundance) {
                    current_highest = comp;
                }
                if (comp->abundance >= current->abundance) {
                    iso_nr_max = j;
                }
                elements->all_iso_calc_amount++;
                if (elements->all_iso_calc_amount > MAX_ITERATION_4) {
                    free(current);
                    free(monoisotopic);
                    return 4006;
                }
                C->amount++;
            }
        }
        
        if (comb_max_a < current_highest->abundance) {
            comb_max_a = current_highest->abundance;
            if (rtm == 0 || rtm == 3) {
                max_a = current_highest->abundance;
            }
        }else if(rtm == 2){
            max_a = 100;
        }
        
        if(current_highest->abundance > current->abundance){
            for (ptrdiff_t v = (ptrdiff_t)C->amount - 1; v >= 0 ; v--) {
                if ( C->compounds[v].abundance != current_highest->abundance
                    ) {
                    if(C->compounds[v].indicator_iso <= iso_nr_max) {
                        A->compounds[A->amount] = C->compounds[v];
                        A->amount++;
                        
#if USE_REALLOC == 1
                        if (A->amount == allocated_mem_4_A) {
#if SHOW_DETAILS == 1
                            realloc_steps_4++;
#endif
                            allocated_mem_4_A *= ALLOC_FACTOR_4_A;
                            A->compounds = (CompoundMulti*) realloc(A->compounds, allocated_mem_4_A * sizeof(CompoundMulti));
                            if (A->compounds == NULL) {
                                //free(current_highest);
                                free(current);
                                free(monoisotopic);
                                return 4913;
                            }
                        }
#endif
                        
#if SHOW_DETAILS == 1
                        if (max_a_list_4 < A->amount) {
                            max_a_list_4 = A->amount;
                        }
#endif
                        if (A->amount >= MAX_COMPOUNDS_4_A) {
                            free(current);
                            free(monoisotopic);
                            return 4001;
                        }
                    }else{
                        if((100/ max_a) * C->compounds[v].abundance >= threshold) {
                            A2->compounds[A2->amount] = C->compounds[v];
                            A2->amount++;
                            if (A2->amount >= MAX_COMPOUNDS_4_A2) {
                                free(current);
                                free(monoisotopic);
                                return 4002;
                            }

#if USE_REALLOC == 1
                            size_t a2_a = A2->amount;
                            if (a2_a >= allocated_mem_4_A2) {
#if SHOW_DETAILS == 1
                                A2_realloc_steps++;
#endif
                                allocated_mem_4_A2 *= ALLOC_FACTOR_4_A2;
                                size_t alloc_ = allocated_mem_4_A2 * sizeof(CompoundMulti);
                                A2->compounds = (CompoundMulti*) realloc(A2->compounds, alloc_);
                                if (A2->compounds == NULL) {
                                    free(current);
                                    free(monoisotopic);
                                    return 4914;
                                }
                            }
#endif
#if SHOW_DETAILS == 1
                            if (max_a2_list_4 < A2->amount) {
                                max_a2_list_4 = A2->amount;
                            }
#endif
                        }
                    }
                }else{
                    C->compounds[C->amount] = C->compounds[v];
                    *current_highest = C->compounds[C->amount];
                }
            }
            
            *current = *current_highest;
            if ((100/ max_a) * current_highest->abundance >= threshold) {
                (combination->compounds + c)->mass = current_highest->mass;
                (combination->compounds + c)->abundance = current_highest->abundance;
                memmove((combination->compounds + c)->sum, current_highest->sum, size_sum);
                c++;
                if (c >= MAX_COMPOUNDS_4) {
                    free(current);
                    free(monoisotopic);
                    return 4003;
                }
                
#if USE_REALLOC == 1
                if (c == allocated_mem_4) {
#if SHOW_DETAILS == 1
                    realloc_steps_4++;
#endif
                    allocated_mem_4 *= ALLOC_FACTOR_4;
                    size_t alloc_ = allocated_mem_4 * sizeof(Compound);
                    combination->compounds = (Compound*)realloc(combination->compounds, alloc_);
                    if (combination->compounds == NULL) {
                        free(current);
                        free(monoisotopic);
                        return 4918;
                    }
                }
#endif
            }
        }
        else{
            for (ptrdiff_t v = 0; v < C->amount; v++) {
                if ( (100/ max_a) * C->compounds[v].abundance >= threshold ) {
                    A2->compounds[A2->amount] = C->compounds[v];
                    A2->amount++;
                    if (A2->amount >= MAX_COMPOUNDS_4_A2) {
                        free(current);
                        free(monoisotopic);
                        return 4002;
                    }
                    
#if USE_REALLOC == 1
                    size_t a2_a = A2->amount;
                    if (a2_a >= allocated_mem_4_A2) {
#if SHOW_DETAILS == 1
                        A2_realloc_steps++;
#endif
                        allocated_mem_4_A2 *= ALLOC_FACTOR_4_A2;
                        size_t alloc_ = allocated_mem_4_A2 * sizeof(CompoundMulti);
                        A2->compounds = (CompoundMulti*) realloc(A2->compounds, alloc_);
                        if (A2->compounds == NULL) {
                            free(current);
                            free(monoisotopic);
                            return 4916;
                        }
                    }
#endif
                    
#if SHOW_DETAILS == 1
                    if (max_a2_list_4 < A2->amount) {
                        max_a2_list_4 = A2->amount;
                    }
#endif
                }
            }
            
            if (A->amount > 0) {
                CompoundMulti *a_c = &A->compounds[A->amount - 1];
                *current = *a_c;
                A->amount--;
                
                if ((100/ max_a) * a_c->abundance >= threshold) {
                    (combination->compounds + c)->mass = a_c->mass;
                    (combination->compounds + c)->abundance = a_c->abundance;
                    memmove((combination->compounds + c)->sum, a_c->sum, size_sum);
                    c++;
                    if (c >= MAX_COMPOUNDS_4) {
                        free(current);
                        free(monoisotopic);
                        return 4003;
                    }
#if USE_REALLOC == 1
                    if (c == allocated_mem_4) {
#if SHOW_DETAILS == 1
                        realloc_steps_4++;
#endif
                        allocated_mem_4 *= ALLOC_FACTOR_4;
                        size_t alloc_ = allocated_mem_4 * sizeof(Compound);
                        combination->compounds = (Compound*)realloc(combination->compounds, alloc_);
                        if (combination->compounds == NULL) {
                            free(current);
                            free(monoisotopic);
                            return 4917;
                        }
                    }
#endif
                }
            }else if(A->amount == 0 && A2->amount > 0){
                CompoundMulti *a2 = &A2->compounds[A2->amount - 1];
                *current = *a2;
                A2->amount--;
            
                if ( (100/ max_a) * a2->abundance >= threshold ) {
                    (combination->compounds + c)->mass = a2->mass;
                    (combination->compounds + c)->abundance = a2->abundance;
                    memmove((combination->compounds + c)->sum, a2->sum, size_sum);
                    c++;
                    if (c >= MAX_COMPOUNDS_4) {
                        free(current);
                        free(monoisotopic);
                        return 4003;
                    }
                    
#if USE_REALLOC == 1
                    if (c == allocated_mem_4) {
#if SHOW_DETAILS == 1
                        realloc_steps_4++;
#endif
                        allocated_mem_4 *= ALLOC_FACTOR_4;
                        size_t alloc_ = allocated_mem_4 * sizeof(Compound);
                        combination->compounds = (Compound*)realloc(combination->compounds, alloc_);
                        if (combination->compounds == NULL) {
                            free(current);
                            free(monoisotopic);
                            return 4918;
                        }
                    }
#endif
                }
            }else {
                current->abundance = -1.0;
            }
        }
    }
    
    combination->amount = c;
    combination->max_abundance = comb_max_a;
    combination->element.iso_amount = iso_amount;
    strcpy (combination->element.name,elements->name);
    for (int i = 1; i < element_amount; i++) {
        strcat(combination->element.name,(elements+i)->name);
    }
    
#if SHOW_DETAILS == 1
#if USE_IN_R == 1
    Rprintf("\n    container: %s, maxA: %zu, maxA2: %zu, realloc_steps: %d, A_realloc_steps: %d, A2_realloc_steps: %d, update_steps: %d, ", combination->element.name, max_a_list_4, max_a2_list_4, realloc_steps_4, A_realloc_steps,A2_realloc_steps, elements->all_iso_calc_amount);
#else
    printf("\n    container: %s, maxA: %zu, maxA2: %zu, realloc_steps: %d, A_realloc_steps: %d, A2_realloc_steps: %d, update_steps: %d, ",combination->element.name, max_a_list_4, max_a2_list_4, realloc_steps_4, A_realloc_steps,A2_realloc_steps,elements->all_iso_calc_amount);
    
#endif
#endif
    free(current);
    free(monoisotopic);
   	return 0;
}

// cleaning compounds for algo4 ///////////////////////////////////////////////////////////////////////////
int clean_combination_algo_4(Combination_4* combinations, long double threshold, size_t element_amount){
    
    
    long double clean_abundance = 1.0;
    long double clean_abundance_other = 1.0;
    for (ptrdiff_t b = 0; b < element_amount; b++) {
        
        //qsort((combinations + b)->compounds, (combinations + b)->amount, sizeof(Compound), compound_sort_by_abundance_dec);
        
        clean_abundance_other = 1.0;
        for (ptrdiff_t d = 0; d < element_amount; d++) {
            if (d != b) {
                clean_abundance_other *= (combinations + d)->max_abundance;
            }
        }
        
        for (ptrdiff_t c = (ptrdiff_t)(combinations + b)->amount - 1; c >= 0; c--) {
            clean_abundance = 1.0;
            clean_abundance *= ((combinations + b)->compounds + c)->abundance * clean_abundance_other;
            
            if ( clean_abundance < threshold) {
                if (c == (combinations + b)->amount - 1) {
                    (combinations + b)->amount--;
                    
                }
                else{
                    *((combinations + b)->compounds + c) = *((combinations + b)->compounds + (combinations + b)->amount - 1);
                    (combinations + b)->amount--;
                }
            }
        }
    }
    
    return 0;
}


int combine_combinations_algo_4(Combination_4* combinations,
                                double threshold,
                                unsigned short  element_amount,
                                size_t* peak_amount,
                                int  peak_limit,
                                unsigned short iso_amount,
                                long double max_abundance,
                                int rtm, double **m_, double **a_, int **cc_){
    
#if USE_REALLOC == 1
    size_t allocated_mem_pl = ALLOC_START_PL;
#endif
    
    size_t tracking[element_amount];
    tracking[0] = 0;
    for (ptrdiff_t i = 0; i < element_amount; i++) {
        tracking[i] = 0;
        qsort((combinations + i)->compounds, (combinations + i)->amount, sizeof(Compound), compound_sort_by_abundance_dec);
    }
    
    const unsigned short iso_a = iso_amount;
    int sum[iso_a];
    size_t iso_sum_size = sizeof(sum);
    
    double mass = 0.0;
    long double abundance = 1.0;
    size_t v = 0;
    
    while (tracking[0] < (combinations)->amount) {
        
        mass = 0.0;
        abundance = 1.0;
        unsigned int cc_count = 0;
        unsigned int cc_tmp[iso_a];
        
        for (unsigned short j = 0; j < element_amount; j++) {
            
            const unsigned short e_a = (combinations + j)->element.iso_amount;
            int sum_e_iso[e_a];
            size_t e_iso_size = sizeof(sum_e_iso);
            // the last element
            if (j == element_amount - 1) {
                double tmp_mass = 0.0;
                long double tmp_abundance = 1.0;
                for (ptrdiff_t h = 0; h < (combinations + j)->amount; h++) {
                    tmp_mass = mass;
                    tmp_abundance = abundance;
                    tmp_mass += ((combinations + j)->compounds + tracking[j])->mass;
                    tmp_abundance *= ((combinations + j)->compounds + tracking[j])->abundance;
                    if ((100/max_abundance)* tmp_abundance >= threshold) {
                        *(*m_ + v) = tmp_mass;
                        if (rtm == 3 || rtm == 4) {
                            *(*a_ + v) = (double)tmp_abundance;
                        }else{
                            *(*a_ + v) = (double)((100/max_abundance) * tmp_abundance);
                        }
                        memmove(*cc_ + v * iso_a, cc_tmp, iso_sum_size);
                        memmove(*cc_ + v * iso_a + cc_count, combinations[j].compounds[tracking[j]].sum, e_iso_size);
                        v++;
                        if (v >= peak_limit) {
                            *(peak_amount) = v;
                            return 4201;
                        }
                        
#if USE_REALLOC == 1
                        if (v == allocated_mem_pl) {
                            
                            allocated_mem_pl *= ALLOC_FACTOR_PL;
                            *m_ = (double*) realloc(*m_, allocated_mem_pl * sizeof(double));
                            if (*m_ == NULL) {
                                *(peak_amount) = v;
                                return 4919;
                            }
                            *a_ = (double*) realloc(*a_, allocated_mem_pl * sizeof(double));
                            if (*a_ == NULL) {
                                *(peak_amount) = v;
                                return 4920;
                            }
                            *cc_ = (int*)realloc(*cc_, allocated_mem_pl * iso_a * sizeof(int));
                            if (*cc_ == NULL) {
                                *(peak_amount) = v;
                                return 4921;
                            }
                        }
#endif
                        tracking[j]++;
                    }else{
                        tracking[j]++;
                        break;
                    }
                }
            }else{
                if (tracking[j] < (combinations + j)->amount) {
                    mass += ((combinations + j)->compounds + tracking[j])->mass;
                    abundance *= ((combinations + j)->compounds + tracking[j])->abundance;
                    
                    for (unsigned short u = 0; u < (combinations + j)->element.iso_amount; u++) {
                        cc_tmp[cc_count] = ((combinations + j)->compounds + tracking[j])->sum[u];
                        cc_count++;
                        if (cc_count >= iso_a) {
                            *(peak_amount) = v;
                            return 4202;
                        }
                    }
                }else{
                    break;
                }
                
            }
        }
        for (ptrdiff_t k = element_amount - 2; k >= 0; k--) {
            if (tracking[k] < combinations[k].amount) {
                tracking[k]++;
                for (ptrdiff_t l = k + 1; l < element_amount; l++) {
                    tracking[l] = 0;
                }
                break;
            }
        }
    }
    
    *(peak_amount) = v;
    return 0;
}

int compound_sort_by_abundance_dec(const void *a, const void *b)
{
    long double y1 = ((const struct Compound*)a)->abundance;
    long double y2 = ((const struct Compound*)b)->abundance;
    
    if (y1 < y2) {
        return 1;
    } else if (y1 > y2) {
        return -1;
    }
    return 0;
}

int compound_sort_by_abundance_inc(const void *a, const void *b)
{
    long double y1 = ((const struct Compound*)a)->abundance;
    long double y2 = ((const struct Compound*)b)->abundance;
    
    if (y1 > y2) {
        return 1;
    } else if (y1 < y2) {
        return -1;
    }
    return 0;
}

int compoundmulti_sort_by_abundance_dec(const void *a, const void *b)
{
    long double y1 = ((const struct CompoundMulti*)a)->abundance;
    long double y2 = ((const struct CompoundMulti*)b)->abundance;
    
    if (y1 < y2) {
        return 1;
    } else if (y1 > y2) {
        return -1;
    }
    return 0;
}

int compoundmulti_sort_by_abundance_inc(const void *a, const void *b)
{
    long double y1 = ((const struct CompoundMulti*)a)->abundance;
    long double y2 = ((const struct CompoundMulti*)b)->abundance;
    
    if (y1 > y2) {
        return 1;
    } else if (y1 < y2) {
        return -1;
    }
    return 0;
}

int combinations_sort_by_amount_inc(const void *a, const void *b)
{
    double y1 = ((const struct Combination_1*)a)->amount;
    double y2 = ((const struct Combination_1*)b)->amount;
    
    if (y1 > y2) {
        return 1;
    } else if (y1 < y2) {
        return -1;
    }
    return 0;
}

int combinations_sort_by_amount_dec(const void *a, const void *b)
{
    double y1 = ((const struct Combination_1*)a)->amount;
    double y2 = ((const struct Combination_1*)b)->amount;
    
    if (y1 < y2) {
        return 1;
    } else if (y1 > y2) {
        return -1;
    }
    return 0;
}

void calc_monoisotopic_single(Element* element, Compound *monoisotopic){
    
    long double abundance = 1.0;
    double mass = 0.0;
    int  pos = 0;
    
    for (ptrdiff_t g = 0; g < element->amount; g++) {
        mass += (element->isotopes)->mass;
        abundance *= (element->isotopes)->abundance;
        monoisotopic->sum[pos]++;
    }
    
    monoisotopic->counter = 0;
    monoisotopic->mass = mass;
    monoisotopic->abundance = abundance;
    monoisotopic->indicator_iso = 0;
}


void create_isotope_list_single(Element element, Isotope2 *isotopes){
    
    size_t i = 0;
    for (unsigned short j = 1; j < element.iso_amount; j++) {
        (isotopes + i)->element_nr = j;
        (isotopes + i)->iso_e_nr = j;
        (isotopes + i)->amount = element.amount;
        (isotopes + i)->abundance = element.isotopes[j].abundance;
        (isotopes + i)->mass = element.isotopes[j].mass;
        strcpy((isotopes + i)->symbol, element.isotopes[j].symbol);
        strcpy((isotopes + i)->isotope, element.isotopes[j].isotope);
        i++;
    }
    qsort(isotopes, i, sizeof(Isotope2), isotope2_sort_by_n_abundance_dec);
}

void create_isotope_list(Element *elements, size_t element_amount, Isotope2 *isotopes, unsigned short *iso_c){
    *iso_c = 0;
    for (unsigned short i = 0; i < element_amount; i++) {
        for (unsigned short j = 1; j < (elements + i)->iso_amount; j++) {
            (isotopes + *iso_c)->element_nr = i;
            (isotopes + *iso_c)->iso_e_nr = j;
            (isotopes + *iso_c)->amount = (elements + i)->amount;
            (isotopes + *iso_c)->abundance = (((elements + i)->isotopes + j)->abundance);
            (isotopes + *iso_c)->mass = (((elements + i)->isotopes + j)->mass);
            strcpy((isotopes + *iso_c)->symbol, ((elements + i)->isotopes + j)->symbol);
            strcpy((isotopes + *iso_c)->isotope, ((elements + i)->isotopes + j)->isotope);
            *(iso_c) += 1;
        }
    }
    qsort(isotopes, *iso_c, sizeof(Isotope2), isotope2_sort_by_n_abundance_dec);
}

void calc_monoisotopic(Element* elements, size_t element_amount, CompoundMulti *monoisotopic){
    long double abundance = 1.0;
    double mass = 0.0;
    size_t  pos = 0;
    for (unsigned short i  = 0; i < element_amount; i++) {
        for (ptrdiff_t g = 0; g < (elements + i)->amount; g++) {
            mass += ((elements + i)->isotopes)->mass;
            abundance *= ((elements + i)->isotopes)->abundance;
            monoisotopic->sum[pos]++;
        }
        pos+= (elements + i)->iso_amount;
        monoisotopic->counter[i] = 0;
    }
    monoisotopic->mass = mass;
    monoisotopic->abundance = abundance;
    monoisotopic->indicator_iso = 0;
}


