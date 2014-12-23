//
//  main.c
//
//  Created by Christian Gerber on 11/28/12.
//  Copyright (c) 2012 EAWAG. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>

#include "isotope.h"
#include "parse.h"
#include "element.h"
#include "combination.h"
#include "peak.h"
#include "profile.h"
#include "n-tuple.h"
#include "preferences.h"


int calc_algo_1(Element *elements, double *mass, double *a, int *cc, unsigned int *p_a, double threshold, int i_a,int e_a, double mono_abundance){

    unsigned int peak_limit = MAX_PEAKS;
    unsigned int peak_amount = 0;
    
    unsigned short element_amount = e_a;
    unsigned short iso_amount = i_a;
    
    double a_max = 1.0;
    peak_amount = 0;
    
    CombinationMulti_A* A = (CombinationMulti_A*)malloc(sizeof(CombinationMulti_A));
    CombinationMulti_C* C = (CombinationMulti_C*)malloc(sizeof(CombinationMulti_C));
    Combination2* combinations = (Combination2*)calloc(element_amount,sizeof(Combination2));

	if(mono_abundance > 0.0){
		
		for (unsigned short i = 0; i < element_amount; i++) {
			calc_combination_max_abundance_mono((combinations + i), (elements + i), threshold,A,C,mono_abundance);
			a_max *= (combinations + i)->max_abundance;	
		} 
		double clean_abundance = 1.0;
		for (unsigned short j = 0; j < element_amount; j++) {
			
			double combination_mono_abundance = 1.0;
			for(int k = 0; k < (elements+j)->amount; k++){
				combination_mono_abundance *= (elements + j)->isotopes[0].abundance;
			}
			
			clean_abundance = a_max/combination_mono_abundance;
			if(create_combination_algo_1(combinations + j, elements + j, clean_abundance, (threshold * mono_abundance)/100.0, peak_limit,A,C)){
				Rprintf("ERROR: could not create combinations\n");
				free(combinations);
				free(A);
				free(C);
				return 1;
			}
			//Rprintf("\ncombination %s    amount %d\n    threshold %f", (combinations +j)->element.name, (combinations +j)->amount,(threshold * mono_abundance)/100.0);
			
			//for(int l = 0; l < combinations[j].amount; l++){
				//Rprintf("\n\t%d,  %s   %f       %f",l,combinations[j].element.name,combinations[j].compounds[l].mass,combinations[j].compounds[l].abundance);
			//}
			
		}
		
		if(combine_combinations_algo_1(combinations, threshold, element_amount, mass, a, cc, &peak_amount, peak_limit, iso_amount, mono_abundance)){
			Rprintf("ERROR: could not combine combinations\n");
			free(combinations);
			free(A);
			free(C);
			return 1;
		}
		
	}else{
		
		for (unsigned short i = 0; i < element_amount; i++) {
			calc_combination_max_abundance((combinations + i), (elements + i), threshold,A,C);
			a_max *= (combinations + i)->max_abundance;	
		}  

	    for (unsigned short j = 0; j < element_amount; j++) {
			double clean_abundance = a_max/(combinations + j)->max_abundance;
			if(create_combination_algo_1(combinations + j, elements + j, clean_abundance, (threshold * a_max)/100.0, peak_limit,A,C)){
				Rprintf("ERROR: could not create combinations\n");
				free(combinations);
				free(A);
				free(C);
				return 1;
			}
			//Rprintf("\ncombination %s    amount %d\n   threshold %f", (combinations +j)->element.name, (combinations +j)->amount, clean_abundance,(threshold * a_max)/100.0);
			
			//for(int l = 0; l < combinations[j].amount; l++){
				//Rprintf("\n\t%d,  %s   %f       %f",l,combinations[j].element.name,combinations[j].compounds[l].mass,combinations[j].compounds[l].abundance);
			//}
		}
		if(combine_combinations_algo_1(combinations, threshold, element_amount, mass, a, cc, &peak_amount, peak_limit, iso_amount, a_max)){
			Rprintf("ERROR: could not combine combinations\n");
			free(combinations);
			free(A);
			free(C);
			return 1;
		}
	}
    
    *p_a = peak_amount;
    
    free(combinations);
    free(A);
    free(C);
    return 0;
}

SEXP iso_ppm_Call(SEXP start, SEXP end, SEXP ppm_R) {
    
    double s;
    double e;
    double ppm;
    unsigned int storage_limit = 2e7;
    
    double *tr;
    
    PROTECT(start = AS_NUMERIC(start));
    PROTECT(end = AS_NUMERIC(end));
    PROTECT(ppm_R = AS_NUMERIC(ppm_R));
    
    s = NUMERIC_VALUE(start);
    e = NUMERIC_VALUE(end);
    ppm = NUMERIC_VALUE(ppm_R);
    
    SEXP trace_R;
    PROTECT(trace_R = NEW_NUMERIC(storage_limit));
    tr = NUMERIC_POINTER(trace_R);
    
    *tr = s;
    
    int n = 0;
    while (*(tr + n) < e ) {
        double tmp = *(tr + n) + *(tr + n) * ppm/1e6;
        n++;
        *(tr + n) = tmp;
        if (n > storage_limit - 1) {
            Rprintf("ERROR: too many mass points for ppm trace\n");
            UNPROTECT(4);
            return R_NilValue;
        }
    }
    SETLENGTH(trace_R, n + 1);
    UNPROTECT(4);
    return trace_R;
}

/*************************
 
 calculates isotopic pattern with algorithm 2
 
 sum:           character array of sum formula
 peak_limit:    maximum allowed peak amount
 thres:         only peaks equal or above this threshold are returned
 iso_list:      character array of elements with their isotopes
 
 *************************/
SEXP iso_pattern_Call_2(SEXP sum, SEXP peak_limit, SEXP threshold, SEXP iso_list, SEXP rel_to_mono) {
    
    PROTECT(sum = AS_CHARACTER(sum));
    PROTECT(iso_list = AS_CHARACTER(iso_list));
    PROTECT(peak_limit = AS_INTEGER(peak_limit));
    PROTECT(threshold = AS_NUMERIC(threshold));
    PROTECT(rel_to_mono = AS_INTEGER(rel_to_mono));
    
    char *s = R_alloc(strlen(CHARACTER_VALUE(sum)), sizeof(char));
    char *i_l = R_alloc(strlen(CHARACTER_VALUE(iso_list)), sizeof(char));
    int p_l = INTEGER_VALUE(peak_limit);
    double t = NUMERIC_VALUE(threshold);
    int r_m = INTEGER_VALUE(rel_to_mono);
    
    UNPROTECT(5);

    if (p_l > MAX_PEAKS) {
        Rprintf("ERROR: peak limit of %d exeeds MAX peak limit of %d\n", p_l, MAX_PEAKS);
        return R_NilValue;
    }
    
    strcpy(s, CHARACTER_VALUE(sum));
    strcpy(i_l, CHARACTER_VALUE(iso_list));

    unsigned short element_amount = 0;
    unsigned short mass_amount = 0;
    unsigned short iso_amount = 0;
    
    Element* elements = (Element*)malloc(MAX_ELEMENTS * sizeof(Element));
    if ( parse_sum_formula(elements, s, &element_amount, &mass_amount, &iso_amount, i_l) ) {
        Rprintf("ERROR: cannot parse sum formula with the given isolist\n");
        free(elements);
        return R_NilValue;
    }

    double *mass = (double*)malloc(p_l * sizeof(double));
    double *a = (double*)malloc(p_l * sizeof(double));

    double max_a = 0.0;
    int *cc = (int*)malloc(p_l * MAX_ISO_SIZE * sizeof(int));
    unsigned int peak_amount = 0;

    CompoundMulti* monoisotopic = (CompoundMulti*)calloc(1,sizeof(CompoundMulti)); 
    calc_monoisotopic(elements, element_amount, monoisotopic);

	if(r_m){
	    if(calc_pattern_algo_2_mono(mass, a ,cc, &max_a, elements, element_amount, t, &peak_amount, p_l,monoisotopic->abundance)){
	        Rprintf("ERROR: could not create combinations\n");
	        free(elements);
	        free(mass);
	        free(a);
	        free(cc);
	        return R_NilValue;
	    }
	}else{
	    if(calc_pattern_algo_2(mass, a ,cc, &max_a, elements, element_amount, t, &peak_amount, p_l)){
	        Rprintf("ERROR: could not create combinations\n");
	        free(elements);
	        free(mass);
	        free(a);
	        free(cc);
	        return R_NilValue;
		}
	
	}
    
    SEXP iso_pattern = PROTECT(allocVector(VECSXP, 3 + iso_amount));
    SEXP mass_R;
    SEXP a_R;

    PROTECT(mass_R = NEW_NUMERIC(peak_amount));
    PROTECT(a_R = NEW_NUMERIC(peak_amount));
    double *p_m_R = NUMERIC_POINTER(mass_R);
    double *p_a_R = NUMERIC_POINTER(a_R);
    

    if(r_m){
		max_a = monoisotopic->abundance;
	}
    free(monoisotopic);
     
    int c = 0;
    for (int i = 0; i < peak_amount; i++) {
        double tmp = *(a + i);
        *(a + i) = (100/ max_a) * tmp;
        if (*(a + i) >= t) {
            p_m_R[c] = *(mass + i);
            p_a_R[c] = *(a + i);
            c++;
        }
    }
    
    int r = 0;
    for (int j = 0; j < iso_amount; j++) {
        SEXP compound_R;
        PROTECT(compound_R = NEW_INTEGER(c));
        int *compound = INTEGER_POINTER(compound_R);
        r = 0;
        for (int k = 0; k < peak_amount; k++) {
            if (*(a + k) >= t) {
                *(compound + r) = *(cc + k * MAX_ISO_SIZE + j);
                r++;
            }
        }
        /*SETLENGTH(compound_R, iso_amount);*/
        SET_VECTOR_ELT(iso_pattern, j + 2, compound_R);
        UNPROTECT(1);
    }
    
    SETLENGTH(mass_R, r);
    SETLENGTH(a_R, r);
    SET_VECTOR_ELT(iso_pattern, 0, mass_R);
    SET_VECTOR_ELT(iso_pattern, 1, a_R);
    
    SEXP list_names;
    int n_c = 2;
    char* names = (char*)malloc((iso_amount + 2) * MAX_NAME_SIZE * sizeof(char));
    strncpy(names, "mass\0", 5 * sizeof(char));
    strncpy(names + MAX_NAME_SIZE, "abundance\0", 10 * sizeof(char));
    
    for (unsigned short e = 0; e < element_amount; e++) {
        for (unsigned short ee = 0; ee < (elements + e)->iso_amount; ee++) {
            strncpy(names + n_c * MAX_NAME_SIZE, ((elements + e)->isotopes +ee)->isotope, MAX_NAME_SIZE * sizeof(char));
            n_c++;
        }
    }
    
    PROTECT(list_names = allocVector(STRSXP, iso_amount + 3));
    for(int o = 0; o < n_c; o++){
        char tmp[ MAX_NAME_SIZE ];
        strncpy(tmp, names + o * MAX_NAME_SIZE, MAX_NAME_SIZE * sizeof(char));
        SET_STRING_ELT(list_names, o,  mkChar(tmp));
    }
    SET_STRING_ELT(list_names, n_c,  mkChar("NAMES"));
    setAttrib(iso_pattern, R_NamesSymbol, list_names);
    SET_VECTOR_ELT(iso_pattern, n_c, list_names);
    
    free(elements);
    free(mass);
    free(a);
    free(cc);
    free(names);
    UNPROTECT(4);
    return iso_pattern;
}

/*************************
 
 calculates isotopic pattern with algorithm 3
 
 sum:           character array of sum formula
 peak_limit:    maximum allowed peak amount
 threshold:     only peaks equal or above this threshold are returned
 iso_list:      character array of elements with their isotopes
 
 *************************/
SEXP iso_pattern_Call_3(SEXP sum, SEXP peak_limit, SEXP threshold, SEXP iso_list, SEXP rel_to_mono) {

    char *s;
    char *i_l;
    int p_l;
    double t;
    int r_m;
    
    PROTECT(sum = AS_CHARACTER(sum));
    PROTECT(iso_list = AS_CHARACTER(iso_list));
    PROTECT(peak_limit = AS_INTEGER(peak_limit));
    PROTECT(threshold = AS_NUMERIC(threshold));
    PROTECT(rel_to_mono = AS_INTEGER(rel_to_mono));
    
    s = R_alloc(strlen(CHARACTER_VALUE(sum)), sizeof(char));
    i_l = R_alloc(strlen(CHARACTER_VALUE(iso_list)), sizeof(char));
    p_l = INTEGER_VALUE(peak_limit);
    t = NUMERIC_VALUE(threshold);
    r_m = INTEGER_VALUE(rel_to_mono);
    
    UNPROTECT(5);
    
    int max_p = MAX_PEAKS;
    if (p_l > max_p) {
        Rprintf("ERROR: peak limit of %d exeeds MAX peak limit of %d\n", p_l, max_p);
        return R_NilValue;
    }

    strcpy(s, CHARACTER_VALUE(sum));
    strcpy(i_l, CHARACTER_VALUE(iso_list));

    Element* elements = NULL;
    
    unsigned short element_amount = 0;
    unsigned short mass_amount = 0;
    unsigned short iso_amount = 0;
    
    elements = (Element*)malloc(MAX_ELEMENTS * sizeof(Element));
    if ( parse_sum_formula(elements, s, &element_amount, &mass_amount, &iso_amount, i_l) ) {
        Rprintf("ERROR: cannot parse sum formula with the given isolist\n");
        return R_NilValue;
    }

       
    double a_monoisotopic = 0.0;
    for (unsigned short i = 0; i < element_amount; i++) {
        a_monoisotopic += pow((elements + i)->isotopes[0].abundance, (elements + i)->amount);
    }
    if (iso_amount > MAX_ISO_SIZE || iso_amount == 0) {
        free(elements);
        return R_NilValue;
    }
        
    Combination* combinations = NULL;
    
    if (element_amount > 0) {
        combinations = (Combination*)malloc(element_amount * sizeof(Combination));
    }else{
        Rprintf("ERROR: cannot create combinations of the single elements\n");
        free(combinations);
        free(elements);
        return R_NilValue;
    }
    
    double* mass;
    double* a;
    int *cc = (int*)malloc(p_l * iso_amount * sizeof(int));
    unsigned int peak_amount = 0;
    mass = (double*)malloc(p_l * sizeof(double));
    a = (double*)malloc(p_l * sizeof(double));
    
    double max_abundance = 1.0;
    double max_mass = 0.0;
    
    for (unsigned short j = 0; j < element_amount; j++) {
        if (create_combinations_algo_3((combinations + j), elements + j, mass_amount, t * a_monoisotopic / 100)) {
            Rprintf("ERROR: cannot create combinations of the single elements\n");
            UNPROTECT(4);
            free(cc);
            free(mass);
            free(a);
            free(combinations);
            free(elements);
            return R_NilValue;
        }
        //Rprintf("\ncombination %s    amount %d\n", (combinations +j)->element.name, (combinations +j)->amount);
	
        max_abundance *= (combinations + j)->max_abundance;
        max_mass += (combinations + j)->max_mass;
    }
    
    CompoundMulti* monoisotopic = (CompoundMulti*)calloc(1,sizeof(CompoundMulti)); 
    calc_monoisotopic(elements, element_amount, monoisotopic);
    //Rprintf("\nAbundance %f", monoisotopic->abundance);
    if(r_m){
		max_abundance = monoisotopic->abundance;
	}
	free(monoisotopic);
    
    unsigned int ar[element_amount];
    const size_t len = element_amount;
    size_t sizes[element_amount];

    clean_combinations_algo_3(combinations, t * max_abundance / 100 , element_amount);

    //for (unsigned short j = 0; j < element_amount; j++) {
        //Rprintf("\nAfter Cleaning       combination %s    amount %d", (combinations +j)->element.name, (combinations +j)->amount);
        //for(int l = 0; l < combinations[j].amount; l++){
			//Rprintf("\n\t%d,  %s   %f       %f",l,combinations[j].element.name,combinations[j].compounds[l].mass,combinations[j].compounds[l].abundance);
		//}
    //}
    
    double mass_;
    double abundance_;
    for (int l = 0; l < element_amount; l++) {
        sizes[l] = (combinations + l)->amount;
        ar[l] = 0;
    }
    
    SEXP list_names;

    char* l_n = (char*)calloc((iso_amount + 2) * MAX_NAME_SIZE, sizeof(char));
    strncpy(l_n, "mass\0", 5 * sizeof(char));
    strncpy(l_n + MAX_NAME_SIZE, "abundance\0", 10* sizeof(char));
    
    int v = 0;
    int index = 0;
    do {
        mass_ = 0.0;
        abundance_ = 1.0;
        index = 0;
        int com_tmp[iso_amount];
        for (int b = 0; b < element_amount; b++) {
            if (ar[b] <=  (combinations + b)->amount) {
                mass_ += ((combinations + b)->compounds + ar[b])->mass;
                abundance_ *= ((combinations + b)->compounds + ar[b])->abundance;
                
                for (int bb = 0; bb < (combinations + b)->element.iso_amount; bb++) {
                    com_tmp[index] = ((combinations + b)->compounds + ar[b])->sum[bb];
                    if (v == 0) {
                        strncpy(l_n + (index + 2) * MAX_NAME_SIZE, (((combinations + b)->element.isotopes + bb))->isotope, MAX_NAME_SIZE*sizeof(char));
                    }
                    index++;
                }
                
            }else{
                Rprintf("ERROR: Threshold too high, cannot combine all elements\n");
                free(cc);
                free(l_n);
                free(mass);
                free(a);
                free(combinations);
                free(elements);
                return R_NilValue;
            }
        }
        if ( (100/max_abundance)* abundance_ > t) {
            *(mass + v) = mass_;
            *(a + v) = (100/max_abundance)* abundance_;
            
            for (int iu = 0; iu < index; iu++) {
                *(cc + v * iso_amount + iu) = com_tmp[iu];
            }
            v++;
        }
        
        if ( v > p_l ) {
            Rprintf("ERROR: reached peak limit of %d, calculation stopped\n", p_l);
            free(cc);
            free(l_n);
            free(mass);
            free(a);
            free(combinations);
            free(elements);
            return R_NilValue;
        }
    } while (MBnext_n_tuple(ar, len, sizes));

    free(combinations);
    free(elements);
    
    
    peak_amount = v;
    
    SEXP iso_pattern = PROTECT(allocVector(VECSXP, 3 + iso_amount));
    SEXP mass_R;
    SEXP a_R;
    double *p_m_R;
    double *p_a_R;
    PROTECT(mass_R = NEW_NUMERIC(peak_amount));
    PROTECT(a_R = NEW_NUMERIC(peak_amount));
    p_m_R = NUMERIC_POINTER(mass_R);
    p_a_R = NUMERIC_POINTER(a_R);
    
    for (int k = 0; k < peak_amount; k ++) {
        p_m_R[k] = *(mass + k);
        p_a_R[k] = *(a + k);
    }
    
    SET_VECTOR_ELT(iso_pattern, 0, mass_R);
    SET_VECTOR_ELT(iso_pattern, 1, a_R);

    for (int q = 0; q < iso_amount; q++) {
        SEXP compound_R;
        int *compound;
        PROTECT(compound_R = NEW_INTEGER(peak_amount));
        compound = INTEGER_POINTER(compound_R);
        
        for (int r = 0; r < peak_amount; r++) {
            compound[r] = *(cc + r * iso_amount + q);
        }
        
        SET_VECTOR_ELT(iso_pattern, q + 2, compound_R);
        UNPROTECT(1);
    }

    PROTECT(list_names = allocVector(STRSXP, iso_amount + 3));
    for(int o = 0; o < iso_amount + 2; o++){
        char tmp[ MAX_NAME_SIZE ];
        strncpy(tmp, l_n + o * MAX_NAME_SIZE, MAX_NAME_SIZE * sizeof(char));
        SET_STRING_ELT(list_names, o,  mkChar(tmp));
    }
    SET_STRING_ELT(list_names, iso_amount + 2,  mkChar("NAMES"));
    setAttrib(iso_pattern, R_NamesSymbol, list_names);
    SET_VECTOR_ELT(iso_pattern, iso_amount + 2, list_names);

    free(mass);
    free(a);
    free(l_n);
    free(cc);
    
    UNPROTECT(4);
    return iso_pattern;
}
 
/*************************
 
 calculates isotopic pattern with algorithm 1
 
 sum:           character array of sum formula
 peak_limit:    maximum allowed peak amount
 threshold:     only peaks equal or above this threshold are returned
 iso_list:      character array of elements with their isotopes
 
 *************************/
SEXP iso_pattern_Call(SEXP sum, SEXP peak_limit, SEXP threshold, SEXP iso_list, SEXP rel_to_mono) {
    
    char *s;
    char *i_l;
    int p_l;
    double t;
    int r_m;
    
    PROTECT(sum = AS_CHARACTER(sum));
    PROTECT(iso_list = AS_CHARACTER(iso_list));
    PROTECT(peak_limit = AS_INTEGER(peak_limit));
    PROTECT(threshold = AS_NUMERIC(threshold));
    PROTECT(rel_to_mono = AS_INTEGER(rel_to_mono));
    
    s = R_alloc(strlen(CHARACTER_VALUE(sum)), sizeof(char));
    i_l = R_alloc(strlen(CHARACTER_VALUE(iso_list)), sizeof(char));
    p_l = INTEGER_VALUE(peak_limit);
    t = NUMERIC_VALUE(threshold);
    r_m = INTEGER_VALUE(rel_to_mono);
    
    UNPROTECT(5);
    
    int max_p = MAX_PEAKS;
    if (p_l > max_p) {
        Rprintf("\n ERROR: peak limit of %d exeeds MAX peak limit of %d\n", p_l, max_p);
        return R_NilValue;
    }
    
    strcpy(s, CHARACTER_VALUE(sum));
    strcpy(i_l, CHARACTER_VALUE(iso_list));
    
    Element* elements = NULL;
    
    unsigned short element_amount = 0;
    unsigned short mass_amount = 0;
    unsigned short iso_amount = 0;

    elements = (Element*)malloc(MAX_ELEMENTS * sizeof(Element));
    
    if ( parse_sum_formula(elements, s, &element_amount, &mass_amount, &iso_amount, i_l) ) {
        Rprintf("\n ERROR: cannot parse sum formula with the given isolist\n");
        return R_NilValue;
    }
    
    qsort(elements, element_amount, sizeof(Element), elements_sort_by_isoamount_inc);
    
    double* mass;
    double* a;
    int *cc = (int*)malloc(p_l * iso_amount * sizeof(int));
    unsigned int peak_amount = 0;
    mass = (double*)malloc(p_l * sizeof(double));
    a = (double*)malloc(p_l * sizeof(double));
    
    SEXP list_names;
    char* l_n = (char*)calloc((iso_amount + 2) * MAX_NAME_SIZE, sizeof(char));
    strncpy(l_n, "mass\0", 5 * sizeof(char));
    strncpy(l_n + MAX_NAME_SIZE, "abundance\0", 10* sizeof(char));
    
    CompoundMulti* monoisotopic = (CompoundMulti*)calloc(1,sizeof(CompoundMulti)); 
    calc_monoisotopic(elements, element_amount, monoisotopic);
    //Rprintf("\nAbundance %f", monoisotopic->abundance);
    double mono_abundance = -1;
    if(r_m){
		mono_abundance = monoisotopic->abundance;
	}
    free(monoisotopic);
    
    if(calc_algo_1(elements, mass, a, cc, &peak_amount, t, iso_amount,element_amount, mono_abundance)){
        Rprintf("\n ERROR: Unable to calculate isotope pattern\n");
        return R_NilValue;
    }

    int index = 0;
    for (int b = 0; b < element_amount; b++) {
        for (int bb = 0; bb < (elements+b)->iso_amount; bb++) {
            strncpy(l_n + (index + 2) * MAX_NAME_SIZE, ((elements + b)->isotopes + bb)->isotope, MAX_NAME_SIZE*sizeof(char));
            index++;
        }
    }

    free(elements);
    
    SEXP iso_pattern = PROTECT(allocVector(VECSXP, 3 + iso_amount));
    SEXP mass_R;
    SEXP a_R;
    double *p_m_R;
    double *p_a_R; 
    PROTECT(mass_R = NEW_NUMERIC(peak_amount));
    PROTECT(a_R = NEW_NUMERIC(peak_amount));
    p_m_R = NUMERIC_POINTER(mass_R);
    p_a_R = NUMERIC_POINTER(a_R);
    
    for (int k = 0; k < peak_amount; k ++) {
        p_m_R[k] = *(mass + k);
        p_a_R[k] = *(a + k);
    }
    
    SET_VECTOR_ELT(iso_pattern, 0, mass_R);
    SET_VECTOR_ELT(iso_pattern, 1, a_R);
    
    for (int q = 0; q < iso_amount; q++) {
        SEXP compound_R;
        int *compound;
        PROTECT(compound_R = NEW_INTEGER(peak_amount));
        compound = INTEGER_POINTER(compound_R);
        
        for (int r = 0; r < peak_amount; r++) {
            compound[r] = *(cc + r * iso_amount + q);
        }
        
        SET_VECTOR_ELT(iso_pattern, q + 2, compound_R);
        UNPROTECT(1);
    }
    
    PROTECT(list_names = allocVector(STRSXP, iso_amount + 3));
    for(int o = 0; o < iso_amount + 2; o++){
        char tmp[ MAX_NAME_SIZE ];
        strncpy(tmp, l_n + o * MAX_NAME_SIZE, MAX_NAME_SIZE * sizeof(char));
        SET_STRING_ELT(list_names, o,  mkChar(tmp));
    }
    SET_STRING_ELT(list_names, iso_amount + 2,  mkChar("NAMES"));
    setAttrib(iso_pattern, R_NamesSymbol, list_names);
    SET_VECTOR_ELT(iso_pattern, iso_amount + 2, list_names);
    
    free(mass);
    free(a);
    free(l_n);
    free(cc);
    
    UNPROTECT(4);
    return iso_pattern;
}


/*************************
 
 calculates points of interest; centroid,intensoid,valley with a given profile
 ATTENTION sort profile according to the mass before executing this function (increasing order)
 
 profile_mass:          mass values array of the profile
 profile_abundance:     abundance values array of the profile
 type:                  return type of point of interest
    centroid:                   		type = 0
    local maxima(intensoid):    	type = 1
    valley:                     			type = 2
 
 *************************/
SEXP iso_centroid_Call(SEXP profile_mass, SEXP profile_abundance, SEXP type) {

    PROTECT(profile_mass = AS_NUMERIC(profile_mass));
    PROTECT(profile_abundance = AS_NUMERIC(profile_abundance));
    PROTECT(type = AS_INTEGER(type));
    
    double *p_m = NUMERIC_POINTER(profile_mass);
    double *p_a = NUMERIC_POINTER(profile_abundance);
    int t = INTEGER_VALUE(type);
    unsigned int n = LENGTH(profile_mass);

    double* c_m = (double*)malloc(n * sizeof(double));
    double* c_a = (double*)malloc(n * sizeof(double));
    
    double centroid = 0.0;
    double sum_sticks = 0.0;
    
    double upper_sum_a = 0.0;
    double lower_sum_a = 0.0;
    double max_centroid = 0.0;
    double max_intensoid = 0.0;
    int centroid_count = 0;
    
    int j = 0;
    unsigned short step = 1;
    
    for (unsigned int i = step; i < n - step; i++) {
        if (t == 1)
        {
            if (    (
                        *(p_a + i) > *(p_a + i + step) && *(p_a + i) > *(p_a + i - step)
                    )
                    &&
                    (
                        *(p_m + i - step) < *(p_m + i) < *(p_m + i + step)
                    )
                
                )
            {
                *(c_a + j) = *(p_a + i);
                *(c_m + j) = *(p_m + i);
                j++;
                
                if (*(c_a + j) > max_intensoid) {
                    max_intensoid = *(c_a + j);
                }
            }
        }
        else{
            
            centroid += *(p_a + i) * *(p_m + i);
            sum_sticks += *(p_a + i);
            
            double diff_mass = *(p_m + i) - *(p_m + i - 1);
            upper_sum_a += diff_mass * *(p_a + i);
            lower_sum_a += diff_mass * *(p_a + i - 1);
            
            if (( (*(p_a + i) <= *(p_a + i + step) && *(p_a + i) < *(p_a + i - step)) || i == n - step - 1)) {
                
                double centroid_temp = (upper_sum_a + lower_sum_a)/2.0;
                if (t == 2) {
                    *(c_a + j) = *(p_a + i);
                    *(c_m + j) = *(p_m + i);
                    j++;
                }
                
                if (max_centroid < centroid_temp) {
                    max_centroid = centroid_temp;
                }
                
                if (t == 0) {
                    if (sum_sticks > 0.0 ) {
						*(c_a + j) = centroid_temp;
						*(c_m + j) = centroid/sum_sticks;
						j++;
						centroid_count++;
                    }
                }
                centroid = 0.0;
                sum_sticks = 0.0;
                
                upper_sum_a = 0.0;
                lower_sum_a = 0.0;
            }
        }
    }
    
    if(t == 0){
        for (int l = 0; l < centroid_count; l++) {
            double temp = *(c_a +  l) * (100/max_centroid);
            *(c_a + l) = temp;
        }
    }
    
    UNPROTECT(3);

    SEXP POI = PROTECT(allocVector(VECSXP, 2));
    SEXP POI_m_R;
    SEXP POI_a_R;
    double *POI_m;
    double *POI_a;
    PROTECT(POI_m_R = NEW_NUMERIC(j));
    PROTECT(POI_a_R = NEW_NUMERIC(j));
    POI_m = NUMERIC_POINTER(POI_m_R);
    POI_a = NUMERIC_POINTER(POI_a_R);
    
    for (int k = 0; k < j; k ++) {
        POI_m[k] = *(c_m + k);
        POI_a[k] = *(c_a + k);
    }
    
    SET_VECTOR_ELT(POI, 0, POI_m_R);
    SET_VECTOR_ELT(POI, 1, POI_a_R);
    
    UNPROTECT(3);
    
    free(c_m);
    free(c_a);
    
    return POI;
}


/*************************
 
 calculates profile with a given isotope pattern and a given mass array.
 the mass array, called trace in this function, determines the position
 of the profile sticks.
 
 profile_type: 
    Gaussian        0
    Cauchy-Lorentz  1
 mass:          mass values array of the isotope pattern
 abundance:     abundance values array of the isotope pattern
 trace:         mass position of profile sticks
 resolution:    
 threshold:
    threshold == 0:     reduce calculation around current stick with mass value m to the intervall witin the maximum peak
                        [m - Gaussian/Cauchy-Lorentz(max_abundance), m + Gaussian/Cauchy-Lorentz(max_abundance)]
    threshold > 0:      reduce calculation around current stick with mass m to the intervall
                        [m - threshold, m + threshold]
 
 *************************/
SEXP iso_profile_with_trace_Call(
                            SEXP profile_type,
                            SEXP mass,
                            SEXP abundance,
                            SEXP trace,
                            SEXP resolution,
                            SEXP threshold
                            ) {
    int p_t;
    double *m;
    double *a;
    double *tr;
    int r;
    double t;
    
    PROTECT(profile_type = AS_INTEGER(profile_type));
    PROTECT(mass = AS_NUMERIC(mass));
    PROTECT(abundance = AS_NUMERIC(abundance));
    PROTECT(trace = AS_NUMERIC(trace));
    PROTECT(resolution = AS_INTEGER(resolution));
    PROTECT(threshold = AS_NUMERIC(threshold));
    
    p_t = NUMERIC_VALUE(profile_type);
    m = NUMERIC_POINTER(mass);
    a = NUMERIC_POINTER(abundance);
    tr = NUMERIC_POINTER(trace);
    r = INTEGER_VALUE(resolution);
    t = NUMERIC_VALUE(threshold);
    
    unsigned int tr_num = LENGTH(trace);
    unsigned int n = LENGTH(mass);

    double* p_m;
    double* p_a;
    p_m = (double*)malloc(LENGTH(trace) * sizeof(double));
    p_a = (double*)malloc(LENGTH(trace) * sizeof(double));
    
    int p_n = 0;
 
    calc_profile_with_trace(n, m, a, tr_num, tr, p_m, p_a, &p_n, r, p_t, t);
    
    UNPROTECT(6);

    SEXP profile = PROTECT(allocVector(VECSXP, 2));
    SEXP profile_mass_R;
    SEXP profile_a_R;
    double *p_m_R;
    double *p_a_R;
    
    PROTECT(profile_mass_R = NEW_NUMERIC(p_n));
    PROTECT(profile_a_R = NEW_NUMERIC(p_n));
    p_m_R = NUMERIC_POINTER(profile_mass_R);
    p_a_R = NUMERIC_POINTER(profile_a_R);
    
    for (int k = 0; k < p_n; k ++) {
        p_m_R[k] = *(p_m + k);
        p_a_R[k] = *(p_a + k);
    }
    
    SET_VECTOR_ELT(profile, 0, profile_mass_R);
    SET_VECTOR_ELT(profile, 1, profile_a_R);
    
    UNPROTECT(3);
    
    free(p_m);
    free(p_a);
    
    return profile;
}
