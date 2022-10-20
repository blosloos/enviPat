//
//  main.c
//  enviPat_test
//
//  Created by chrigugu on 9/23/14.
//  Copyright (c) 2014 cg. All rights reserved.
//


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stddef.h>
#include <time.h>

#include "preferences.h"
#include "isotope.h"
#include "parse.h"
#include "element.h"
#include "combination.h"
#include "peak.h"
#include "profile.h"

#if USE_IN_R == 1
    #include <R.h>
    #include <Rdefines.h>
#else
    #include "data.h"
#endif


int alloc_peaks(int p_l, size_t iso_amount, double **m_, double **a_, int **cc_){

#if USE_REALLOC == 1
    if (ALLOC_START_PL > p_l) {
        size_t alloc_ = sizeof(double[p_l]);
        *m_ = (double*)malloc(alloc_);
        if(*m_ == NULL){
            return 1;
        }
        *a_ = (double*)malloc(alloc_);
        if(*a_ == NULL){
            free(*m_);
            return 1;
        }
        alloc_ = sizeof(int[p_l * iso_amount]);
        *cc_ = (int*)malloc(alloc_);
        if(*cc_ == NULL){
            free(*m_);
            free(*a_);
            return 1;
        }
    }else{
        size_t alloc_ = ALLOC_START_PL * sizeof(double);
        *m_ = (double*)malloc(alloc_);
        if(*m_ == NULL){
            return 1;
        }
        alloc_ = ALLOC_START_PL * sizeof(double);
        *a_ = (double*)malloc(alloc_);
        if(*a_ == NULL){
            free(*m_);
            return 1;
        }
        alloc_ = ALLOC_START_PL * iso_amount * sizeof(int);
        *cc_ = (int*)malloc(alloc_);
        if(*cc_ == NULL){
            free(*m_);
            free(*a_);
            return 1;
        }
    }

#else
    size_t alloc_ = sizeof(double[p_l]);
    *m_ = (double*)malloc(alloc_);
    if(*m_ == NULL){
        return 1;
    }
    *a_ = (double*)malloc(alloc_);
    if(*a_ == NULL){
        free(*m_);
        return 1;
    }
    alloc_ = sizeof(int[p_l * iso_amount]);
    *cc_ = (int*)malloc(alloc_);
    if(*cc_ == NULL){
        free(*m_);
        free(*a_);
        return 1;
    }
#endif
    return 0;
}


#if USE_IN_R == 1

    SEXP iso_ppm_Call(SEXP start, SEXP end, SEXP ppm_R) {

        double s;
        double e;
        double ppm;
        unsigned int  storage_limit = 2e7;
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
                Rprintf("\nERROR: too many mass points for ppm trace");
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
    SEXP iso_pattern_Call_2(SEXP sum, SEXP peak_limit, SEXP threshold, SEXP iso_list, SEXP rel_to_mono, SEXP return_iso_calc_amount) {

        char *s;
        char *i_l;
        int p_l;
        double t;
        int rtm;
        int rica;

        double* m_;
        double* a_;
        int* cc_;

        PROTECT(sum = AS_CHARACTER(sum));
        PROTECT(iso_list = AS_CHARACTER(iso_list));
        PROTECT(peak_limit = AS_INTEGER(peak_limit));
        PROTECT(threshold = AS_NUMERIC(threshold));
        PROTECT(rel_to_mono = AS_INTEGER(rel_to_mono));
        PROTECT(return_iso_calc_amount = AS_INTEGER(return_iso_calc_amount));

        s = R_alloc(strlen(CHARACTER_VALUE(sum)), sizeof(char));
        i_l = R_alloc(strlen(CHARACTER_VALUE(iso_list)), sizeof(char));
        p_l = INTEGER_VALUE(peak_limit);
        t = NUMERIC_VALUE(threshold);
        rtm = INTEGER_VALUE(rel_to_mono);
        rica = INTEGER_VALUE(return_iso_calc_amount);

        if (p_l >= INT_LIMIT || p_l < 1) {

            Rprintf("\ninvalid peak limit");
			UNPROTECT(6);
            return R_NilValue;
        }

        strcpy(s, CHARACTER_VALUE(sum));
        strcpy(i_l, CHARACTER_VALUE(iso_list));

        UNPROTECT(6);

#if SHOW_DETAILS == 1
        Rprintf("\n%s, algo2, ",s);
        clock_t start, end;
        double cpu_time_used = 0.0;
        start = clock();
#endif

        unsigned short element_amount = 0;
        unsigned short iso_amount = 0;

        long double max_a = 0.0;
        size_t  peak_amount = 0;
        int iso_count_stats = 0;

        if ( rtm > 4 || rtm < 0) {
            Rprintf("\nERROR:  wrong value for rtm");
            return R_NilValue;
        }
        size_t alloc_ = MAX_ELEMENTS * sizeof(Element);
        Element* elements = (Element*)malloc(alloc_);
        if (elements == NULL) {
            Rprintf("\nERROR: cannot allocate memory for elements pointer");
            return R_NilValue;
        }

        // parse the given molecul s /////////////////////////////////////////////////////////////////////
        if ( parse_sum_formula(elements, s, &element_amount, &iso_amount, i_l) ) {
            Rprintf("\nERROR: cannot parse sum formula with the given isolist");
            free(elements);
            return R_NilValue;
        }

        if (iso_amount >= MAX_ISO_SIZE
            || iso_amount <= 0
            || element_amount <= 0
            || element_amount >= MAX_ELEMENTS
            ) {
            free(elements);
            return R_NilValue;
        }

        // allocate pointer for peaks ///////////////////////////////////////////////////////////////
        int msg_peaks = alloc_peaks(p_l, iso_amount, &m_, &a_, &cc_);
        if (msg_peaks) {
            Rprintf("\nERROR: pointer allocation, error code: %d", msg_peaks);
            free(elements);
            return R_NilValue;
        }

        // calc pattern algo 2 //////////////////////////////////////////////////////////////////////
        int msg = calc_pattern_algo_2(&max_a, elements, element_amount, t, &peak_amount, p_l, &iso_count_stats,rtm, &m_, &a_, &cc_);

        if(msg || peak_amount < 1){
            Rprintf("\nERROR: could not create combinations, error code: %d", msg);
            free(elements);
            free(m_);
            free(a_);
            free(cc_);
            return R_NilValue;
        }

        // return data to R ////////////////////////////////////////////////////////////////////////
        if (rica) {
            SEXP out = PROTECT(allocVector(INTSXP, 1));
            INTEGER(out)[0] = iso_count_stats;
            free(elements);
            free(m_);
            free(a_);
            free(cc_);
            UNPROTECT(1);
            return out;
        }else{

            SEXP iso_pattern = PROTECT(allocVector(VECSXP, 3 + iso_amount));
            SEXP mass_R;
            SEXP a_R;

            PROTECT(mass_R = NEW_NUMERIC((long)peak_amount));
            PROTECT(a_R = NEW_NUMERIC((long)peak_amount));
            double *p_m_R = NUMERIC_POINTER(mass_R);
            double *p_a_R = NUMERIC_POINTER(a_R);

            size_t c = 0;
            for (ptrdiff_t i = 0; i < peak_amount; i++) {
                double tmp = *(a_ + i);
                if (rtm != 2) {
                    *(a_ + i) = (100/ (double)max_a) * tmp;
                }

                if (*(a_ + i) > t) {
                    p_m_R[c] = *(m_ + i);

                    if (rtm == 3 || rtm == 4) {
                        p_a_R[c] = tmp;
                    }else{
                        p_a_R[c] = *(a_ + i);
                    }

                    c++;
                }
            }
            int r = 0;
            for (ptrdiff_t j = 0; j < iso_amount; j++) {
                SEXP compound_R;
                double *compound;
                PROTECT(compound_R = NEW_NUMERIC((long)c));
                compound = NUMERIC_POINTER(compound_R);
                r = 0;
                for (ptrdiff_t k = 0; k < peak_amount; k++) {
                    if (*(a_ + k) > t) {
                        *(compound + r) = *(cc_ + k * iso_amount + j);
                        r++;
                    }
                }
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
            if (names == NULL) {
                Rprintf("\nERROR: allocation names pointer failed");
                free(elements);
                free(m_);
                free(a_);
                free(cc_);
				UNPROTECT(3);
                return R_NilValue;
            }

            strncpy(names, "mass\0", 5 * sizeof(char));
            strncpy(names + MAX_NAME_SIZE, "abundance\0", 10 * sizeof(char));

            for (ptrdiff_t  e = 0; e < element_amount; e++) {
                for (ptrdiff_t  ee = 0; ee < (elements + e)->iso_amount; ee++) {
                    strncpy(names + n_c * MAX_NAME_SIZE, ((elements + e)->isotopes +ee)->isotope, MAX_NAME_SIZE * sizeof(char));
                    n_c++;
                }
            }

            PROTECT(list_names = allocVector(STRSXP, iso_amount + 3));
            for(ptrdiff_t o = 0; o < n_c; o++){
                char tmp[ MAX_NAME_SIZE ];

                memcpy(tmp, names + o * MAX_NAME_SIZE, MAX_NAME_SIZE * sizeof(char));
                tmp[ MAX_NAME_SIZE - 1] = '\0';

                // #pragma GCC diagnostic push
                // #pragma GCC diagnostic ignored "-Wstringop-truncation"
                //     strncpy(tmp, names + o * MAX_NAME_SIZE, MAX_NAME_SIZE * sizeof(char));
                // #pragma GCC diagnostic pop

                SET_STRING_ELT(list_names, o,  mkChar(tmp));
            }
            SET_STRING_ELT(list_names, n_c,  mkChar("NAMES"));
            setAttrib(iso_pattern, R_NamesSymbol, list_names);
            SET_VECTOR_ELT(iso_pattern, n_c, list_names);

            free(elements);
            free(m_);
            free(a_);
            free(cc_);
            free(names);
            UNPROTECT(4);

#if SHOW_DETAILS == 1
            end = clock();
            cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
            Rprintf("\n    peak_amount: %zu, elapsed_time: %lf, ",c, cpu_time_used);
#endif
            return iso_pattern;
        }
    }


/*************************

 calculates isotopic pattern with algorithm 4

 sum:           character array of sum formula
 peak_limit:    maximum allowed peak amount
 threshold:     only peaks equal or above this threshold are returned
 iso_list:      character array of elements with their isotopes

 *************************/
SEXP iso_pattern_Call_4(SEXP sum, SEXP peak_limit, SEXP threshold, SEXP iso_list, SEXP rel_to_mono, SEXP return_iso_calc_amount) {

    char *s;
    char *i_l;
    int p_l;
    double t;
    int rtm;
    int rica;

    double* m_;
    double* a_;
    int* cc_;

    PROTECT(sum = AS_CHARACTER(sum));
    PROTECT(iso_list = AS_CHARACTER(iso_list));
    PROTECT(peak_limit = AS_INTEGER(peak_limit));
    PROTECT(threshold = AS_NUMERIC(threshold));
    PROTECT(rel_to_mono = AS_INTEGER(rel_to_mono));
    PROTECT(return_iso_calc_amount = AS_INTEGER(return_iso_calc_amount));

    s = R_alloc(strlen(CHARACTER_VALUE(sum)), sizeof(char));
    i_l = R_alloc(strlen(CHARACTER_VALUE(iso_list)), sizeof(char));
    p_l = INTEGER_VALUE(peak_limit);
    t = NUMERIC_VALUE(threshold);
    rtm = INTEGER_VALUE(rel_to_mono);
    rica = INTEGER_VALUE(return_iso_calc_amount);


    if (p_l >= INT_LIMIT || p_l < 1) {

        Rprintf("\ninvalid peak limit");
		UNPROTECT(6);
        return R_NilValue;
    }

    strcpy(s, CHARACTER_VALUE(sum));
    strcpy(i_l, CHARACTER_VALUE(iso_list));

    UNPROTECT(6);

#if SHOW_DETAILS == 1
    Rprintf("\n%s, algo4, ",s);
    clock_t start, end;
    double cpu_time_used = 0.0;
    start = clock();
#endif

    unsigned short element_amount = 0;
    unsigned short iso_amount = 0;
    size_t  peak_amount = 0;

    if ( rtm > 4 || rtm < 0) {
        Rprintf("\nERROR: wrong value for rtm");
        return R_NilValue;
    }

    Element* elements = (Element*)calloc(MAX_ELEMENTS, sizeof(Element));
    if(elements == NULL){
        Rprintf("\nERROR: cannot allocate memory for elements pointer");
        return R_NilValue;
    }

    if ( parse_sum_formula(elements, s, &element_amount, &iso_amount, i_l) ) {
        Rprintf("\nERROR: cannot parse sum formula with the given isolist");
        free(elements);
        return R_NilValue;
    }

    if (iso_amount >= MAX_ISO_SIZE
        || iso_amount <= 0
        || element_amount <= 0
        || element_amount >= MAX_ELEMENTS
        ) {
        free(elements);
        return R_NilValue;
    }

    // allocate pointer for peaks ///////////////////////////////////////////////////////////////
    //int msg_peaks = alloc_peaks(p_l, iso_amount);
    int msg_peaks = alloc_peaks(p_l, iso_amount, &m_, &a_, &cc_);
    if (msg_peaks) {
        Rprintf("\nERROR: pointer allocation, error code: %d", msg_peaks);
        free(elements);
        return R_NilValue;
    }

    SEXP list_names;
    char* l_n = (char*)calloc((iso_amount + 2) * MAX_NAME_SIZE, sizeof(char));
    if(l_n == NULL){
        Rprintf("\nERROR: cannot allocate memory column names");
        free(elements);
        free(m_);
        free(a_);
        free(cc_);
        return R_NilValue;
    }
    strncpy(l_n, "mass\0", 5 * sizeof(char));
    strncpy(l_n + MAX_NAME_SIZE, "abundance\0", 10 * sizeof(char));
    size_t index = 0;
    for (ptrdiff_t b = 0; b < element_amount; b++) {
        for (ptrdiff_t bb = 0; bb < (elements+b)->iso_amount; bb++) {
            strncpy(l_n + (index + 2) * MAX_NAME_SIZE, ((elements + b)->isotopes + bb)->isotope, MAX_NAME_SIZE*sizeof(char));
            index++;
        }
    }

    long double mono_abundance = -1.0;
    for (ptrdiff_t  i = 0; i < element_amount; i++) {
        mono_abundance *= pow((elements + i)->isotopes[0].abundance, (elements + i)->amount);
    }

    //int msg = calc_pattern_algo_3( elements, &peak_amount, t, iso_amount, element_amount, mono_abundance, p_l, l_n, rtm);
    long double max_a = 0;
    int msg = calc_pattern_algo_4(elements, &peak_amount, t, iso_amount, element_amount, mono_abundance, p_l, l_n, rtm, &m_, &a_, &cc_, &max_a);
     //Rprintf("\n peak amount %zu", peak_amount);
    if(msg || peak_amount == 0){
        Rprintf("\nERROR: cannot combine combinations, exit code: %d", msg);
        free(elements);
        free(m_);
        free(a_);
        free(cc_);
        free(l_n);
        return R_NilValue;

    }

    // return values to R ////////////////////////////////////////////////////////////////////
    if (rica) {

        SEXP list_names_elements;

        int container_amount = 0;
        for (int i = 0; i < element_amount; i++) {
            if ((elements+i)->all_iso_calc_amount > 0) {
                container_amount++;
            }
        }

        PROTECT(list_names_elements = allocVector(STRSXP, container_amount));
        SEXP out;
        out = PROTECT(allocVector(INTSXP, container_amount));

        int j = 0;
        char tmp[ MAX_NAME_SIZE ];
        for (ptrdiff_t i = 0; i < element_amount; i++) {

            if ((elements+i)->all_iso_calc_amount > 0) {
                if (j > 0) {
                    SET_STRING_ELT(list_names_elements, j - 1,  mkChar(tmp));
                }

                INTEGER(out)[j] = (elements+i)->all_iso_calc_amount;

                memcpy(tmp, (elements+i)->name, MAX_NAME_SIZE * sizeof(char));
                tmp[ MAX_NAME_SIZE - 1] = '\0';

                // #pragma GCC diagnostic push
                // #pragma GCC diagnostic ignored "-Wstringop-truncation"
                //     strncpy(tmp, (elements+i)->name, MAX_NAME_SIZE * sizeof(char));
                // #pragma GCC diagnostic pop

                j++;
            }else if(j > 0){
                strcat(tmp, (elements + i)->name);
            }

            if (i == element_amount - 1 && j > 0 && j <=     container_amount) {
                SET_STRING_ELT(list_names_elements, j - 1,  mkChar(tmp));
            }
        }

        setAttrib(out, R_NamesSymbol, list_names_elements);

        free(m_);
        free(a_);
        free(l_n);
        free(cc_);
        free(elements);
        UNPROTECT(2);
        return out;
    }else{

        SEXP iso_pattern = PROTECT(allocVector(VECSXP, 3 + iso_amount));
        SEXP mass_R;
        SEXP a_R;
        double *p_m_R;
        double *p_a_R;
        PROTECT(mass_R = NEW_NUMERIC((long)peak_amount));
        PROTECT(a_R = NEW_NUMERIC((long)peak_amount));
        p_m_R = NUMERIC_POINTER(mass_R);
        p_a_R = NUMERIC_POINTER(a_R);

        for (ptrdiff_t k = 0; k < peak_amount; k ++) {
            p_m_R[k] = *(m_ + k);
            p_a_R[k] = *(a_ + k);
        }

        SET_VECTOR_ELT(iso_pattern, 0, mass_R);
        SET_VECTOR_ELT(iso_pattern, 1, a_R);

        for (ptrdiff_t q = 0; q < iso_amount; q++) {
            SEXP compound_R;
            int *compound;
            PROTECT(compound_R = NEW_INTEGER((long)peak_amount));
            compound = INTEGER_POINTER(compound_R);

            for (ptrdiff_t r = 0; r < peak_amount; r++) {
                compound[r] = *(cc_ + r * iso_amount + q);
            }

            SET_VECTOR_ELT(iso_pattern, q + 2, compound_R);
            UNPROTECT(1);
        }

        PROTECT(list_names = allocVector(STRSXP, iso_amount + 3));
        for(ptrdiff_t o = 0; o < iso_amount + 2; o++){
            char tmp[ MAX_NAME_SIZE ];

            memcpy(tmp, l_n + o * MAX_NAME_SIZE, MAX_NAME_SIZE * sizeof(char));
            tmp[ MAX_NAME_SIZE - 1] = '\0';

            // #pragma GCC diagnostic push
            // #pragma GCC diagnostic ignored "-Wstringop-truncation"
            //     strncpy(tmp, l_n + o * MAX_NAME_SIZE, MAX_NAME_SIZE * sizeof(char));
            // #pragma GCC diagnostic pop

            SET_STRING_ELT(list_names, o,  mkChar(tmp));
        }
        SET_STRING_ELT(list_names, iso_amount + 2,  mkChar("NAMES"));
        setAttrib(iso_pattern, R_NamesSymbol, list_names);
        SET_VECTOR_ELT(iso_pattern, iso_amount + 2, list_names);

        free(m_);
        free(a_);
        free(l_n);
        free(cc_);
        free(elements);
        UNPROTECT(4);

#if SHOW_DETAILS == 1
        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
        Rprintf("\n    peak_amount: %zu, elapsed_time: %lf, ",peak_amount, cpu_time_used);
#endif
        return iso_pattern;
    }
}

    /*************************

     calculates isotopic pattern with algorithm 3

     sum:           character array of sum formula
     peak_limit:    maximum allowed peak amount
     threshold:     only peaks equal or above this threshold are returned
     iso_list:      character array of elements with their isotopes

     *************************/
    SEXP iso_pattern_Call_3(SEXP sum, SEXP peak_limit, SEXP threshold, SEXP iso_list, SEXP rel_to_mono, SEXP return_iso_calc_amount) {

        char *s;
        char *i_l;
        int p_l;
        double t;
        int rtm;
        int rica;

        double* m_;
        double* a_;
        int* cc_;

        PROTECT(sum = AS_CHARACTER(sum));
        PROTECT(iso_list = AS_CHARACTER(iso_list));
        PROTECT(peak_limit = AS_INTEGER(peak_limit));
        PROTECT(threshold = AS_NUMERIC(threshold));
        PROTECT(rel_to_mono = AS_INTEGER(rel_to_mono));
        PROTECT(return_iso_calc_amount = AS_INTEGER(return_iso_calc_amount));

        s = R_alloc(strlen(CHARACTER_VALUE(sum)), sizeof(char));
        i_l = R_alloc(strlen(CHARACTER_VALUE(iso_list)), sizeof(char));
        p_l = INTEGER_VALUE(peak_limit);
        t = NUMERIC_VALUE(threshold);
        rtm = INTEGER_VALUE(rel_to_mono);
        rica = INTEGER_VALUE(return_iso_calc_amount);

        if (p_l >= INT_LIMIT || p_l < 1) {

            Rprintf("\ninvalid peak limit");
			UNPROTECT(6);
            return R_NilValue;
        }

        strcpy(s, CHARACTER_VALUE(sum));
        strcpy(i_l, CHARACTER_VALUE(iso_list));

        UNPROTECT(6);

#if SHOW_DETAILS == 1
        Rprintf("\n%s, algo3, ",s);
        clock_t start, end;
        double cpu_time_used = 0.0;
        start = clock();
#endif

        unsigned short element_amount = 0;
        unsigned short iso_amount = 0;
        size_t  peak_amount = 0;

        if ( rtm > 4 || rtm < 0) {
            Rprintf("\nERROR: wrong value for rtm");
            return R_NilValue;
        }

        Element* elements = (Element*)calloc(MAX_ELEMENTS, sizeof(Element));
        if(elements == NULL){
            Rprintf("\nERROR: cannot allocate memory for elements pointer");
            return R_NilValue;
        }

        if ( parse_sum_formula(elements, s, &element_amount, &iso_amount, i_l) ) {
            Rprintf("\nERROR: cannot parse sum formula with the given isolist");
            free(elements);
            return R_NilValue;
        }

        if (iso_amount >= MAX_ISO_SIZE
            || iso_amount <= 0
            || element_amount <= 0
            || element_amount >= MAX_ELEMENTS
            ) {
            free(elements);
            return R_NilValue;
        }

        // allocate pointer for peaks ///////////////////////////////////////////////////////////////
        int msg_peaks = alloc_peaks(p_l, iso_amount, &m_, &a_, &cc_);
        if (msg_peaks) {
            Rprintf("\nERROR: pointer allocation, error code: %d", msg_peaks);
            free(elements);
            return R_NilValue;
        }

        SEXP list_names;
        char* l_n = (char*)calloc((iso_amount + 2) * MAX_NAME_SIZE, sizeof(char));
        if(l_n == NULL){
            Rprintf("\nERROR: cannot allocate memory column names");
            free(elements);
            free(m_);
            free(a_);
            free(cc_);
            return R_NilValue;
        }
        strncpy(l_n, "mass\0", 5 * sizeof(char));
        strncpy(l_n + MAX_NAME_SIZE, "abundance\0", 10 * sizeof(char));
        size_t index = 0;
        for (ptrdiff_t b = 0; b < element_amount; b++) {
            for (ptrdiff_t bb = 0; bb < (elements+b)->iso_amount; bb++) {
                strncpy(l_n + (index + 2) * MAX_NAME_SIZE, ((elements + b)->isotopes + bb)->isotope, MAX_NAME_SIZE*sizeof(char));
                index++;
            }
        }

//        long double mono_abundance = -1.0;
//        for (ptrdiff_t  i = 0; i < element_amount; i++) {
//            mono_abundance *= pow((elements + i)->isotopes[0].abundance, (elements + i)->amount);
//        }

        CompoundMulti* monoisotopic = (CompoundMulti*)calloc(1,sizeof(CompoundMulti));
        if(monoisotopic == NULL){
            Rprintf("\nERROR: cannot allocate memory for list name pointer");
            free(elements);
            free(m_);
            free(a_);
            free(cc_);
            free(l_n);
            return R_NilValue;
        }

        calc_monoisotopic(elements, element_amount, monoisotopic);
        long double mono_abundance = monoisotopic->abundance;
        free(monoisotopic);

        int msg = calc_pattern_algo_3( elements, &peak_amount, t, iso_amount, element_amount, mono_abundance, p_l, l_n, rtm, &m_, &a_, &cc_);
        if(msg || peak_amount == 0){
            Rprintf("\nERROR: cannot combine combinations, exit code: %d", msg);
            free(elements);
            free(m_);
            free(a_);
            free(cc_);
            free(l_n);
            return R_NilValue;

        }

        // return values to R ////////////////////////////////////////////////////////////////////
        if (rica) {

            SEXP list_names_elements;
            PROTECT(list_names_elements = allocVector(STRSXP, element_amount));
            SEXP out;
            out = PROTECT(allocVector(INTSXP, element_amount));

            for (ptrdiff_t i = 0; i< element_amount; i++) {
                INTEGER(out)[i] = (elements+i)->all_iso_calc_amount;

                char tmp[ MAX_NAME_SIZE ];

                memcpy(tmp, (elements+i)->name, MAX_NAME_SIZE * sizeof(char));
                tmp[ MAX_NAME_SIZE - 1] = '\0';

                // #pragma GCC diagnostic push
                // #pragma GCC diagnostic ignored "-Wstringop-truncation"
                //     strncpy(tmp, (elements+i)->name, MAX_NAME_SIZE * sizeof(char));
                // #pragma GCC diagnostic pop

                SET_STRING_ELT(list_names_elements, i,  mkChar(tmp));
            }

            setAttrib(out, R_NamesSymbol, list_names_elements);

            free(m_);
            free(a_);
            free(l_n);
            free(cc_);
            free(elements);
            UNPROTECT(2);
            return out;
        }else{

            SEXP iso_pattern = PROTECT(allocVector(VECSXP, 3 + iso_amount));
            SEXP mass_R;
            SEXP a_R;
            double *p_m_R;
            double *p_a_R;
            PROTECT(mass_R = NEW_NUMERIC((long)peak_amount));
            PROTECT(a_R = NEW_NUMERIC((long)peak_amount));
            p_m_R = NUMERIC_POINTER(mass_R);
            p_a_R = NUMERIC_POINTER(a_R);

            for (ptrdiff_t k = 0; k < peak_amount; k ++) {
                p_m_R[k] = *(m_ + k);
                p_a_R[k] = *(a_ + k);
            }

            SET_VECTOR_ELT(iso_pattern, 0, mass_R);
            SET_VECTOR_ELT(iso_pattern, 1, a_R);

            for (ptrdiff_t q = 0; q < iso_amount; q++) {
                SEXP compound_R;
                int *compound;
                PROTECT(compound_R = NEW_INTEGER((long)peak_amount));
                compound = INTEGER_POINTER(compound_R);

                for (ptrdiff_t r = 0; r < peak_amount; r++) {
                    compound[r] = *(cc_ + r * iso_amount + q);
                }

                SET_VECTOR_ELT(iso_pattern, q + 2, compound_R);
                UNPROTECT(1);
            }

            PROTECT(list_names = allocVector(STRSXP, iso_amount + 3));
            for(ptrdiff_t o = 0; o < iso_amount + 2; o++){
                char tmp[ MAX_NAME_SIZE ];

                memcpy(tmp, l_n + o * MAX_NAME_SIZE, MAX_NAME_SIZE * sizeof(char));
                tmp[ MAX_NAME_SIZE - 1] = '\0';

                // #pragma GCC diagnostic push
                // #pragma GCC diagnostic ignored "-Wstringop-truncation"
                //     strncpy(tmp, l_n + o * MAX_NAME_SIZE, MAX_NAME_SIZE * sizeof(char));
                // #pragma GCC diagnostic pop

                SET_STRING_ELT(list_names, o,  mkChar(tmp));
            }
            SET_STRING_ELT(list_names, iso_amount + 2,  mkChar("NAMES"));
            setAttrib(iso_pattern, R_NamesSymbol, list_names);
            SET_VECTOR_ELT(iso_pattern, iso_amount + 2, list_names);

            free(m_);
            free(a_);
            free(l_n);
            free(cc_);
            free(elements);
            UNPROTECT(4);

#if SHOW_DETAILS == 1
            end = clock();
            cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
            Rprintf("\n    peak_amount: %zu, elapsed_time: %lf, ",peak_amount, cpu_time_used);
#endif
            return iso_pattern;
        }
    }

    /*************************

     calculates isotopic pattern with algorithm 1

     sum:           character array of sum formula
     peak_limit:    maximum allowed peak amount
     threshold:     only peaks equal or above this threshold are returned
     iso_list:      character array of elements with their isotopes

     *************************/
    SEXP iso_pattern_Call(SEXP sum, SEXP peak_limit, SEXP threshold, SEXP iso_list, SEXP rel_to_mono, SEXP return_iso_calc_amount) {

        char *s;
        char *i_l;
        int p_l;
        double t;
        int rtm;
        int rica;

        double* m_;
        double* a_;
        int* cc_;

        PROTECT(sum = AS_CHARACTER(sum));
        PROTECT(iso_list = AS_CHARACTER(iso_list));
        PROTECT(peak_limit = AS_INTEGER(peak_limit));
        PROTECT(threshold = AS_NUMERIC(threshold));
        PROTECT(rel_to_mono = AS_INTEGER(rel_to_mono));
        PROTECT(return_iso_calc_amount = AS_INTEGER(return_iso_calc_amount));

        s = R_alloc(strlen(CHARACTER_VALUE(sum)), sizeof(char));
        i_l = R_alloc(strlen(CHARACTER_VALUE(iso_list)), sizeof(char));
        p_l = INTEGER_VALUE(peak_limit);
        t = NUMERIC_VALUE(threshold);
        rtm = INTEGER_VALUE(rel_to_mono);
        rica = INTEGER_VALUE(return_iso_calc_amount);


        if (p_l >= INT_LIMIT || p_l < 1) {

            Rprintf("\ninvalid peak limit");
			UNPROTECT(6);
            return R_NilValue;
        }

        strcpy(s, CHARACTER_VALUE(sum));
        strcpy(i_l, CHARACTER_VALUE(iso_list));

        UNPROTECT(6);

#if SHOW_DETAILS == 1
        Rprintf("\n%s, algo1, ",s);
        clock_t start, end;
        double cpu_time_used = 0.0;
        start = clock();
#endif

        unsigned short element_amount = 0;
        unsigned short iso_amount = 0;
        size_t  peak_amount = 0;

        if ( rtm > 4 || rtm < 0) {
            Rprintf("\nERROR: wrong value for rtm");
            return R_NilValue;
        }

        Element* elements = (Element*)calloc(MAX_ELEMENTS, sizeof(Element));
        if (elements == NULL) {
            Rprintf("\nERROR: cannot allocate memory for elements pointer");
            return R_NilValue;
        }

        if ( parse_sum_formula(elements, s, &element_amount, &iso_amount, i_l) ) {
            Rprintf("\nERROR: cannot parse sum formula with the given isolist");
            free(elements);
            return R_NilValue;
        }

        if (iso_amount >= MAX_ISO_SIZE
            || iso_amount <= 0
            || element_amount <= 0
            || element_amount >= MAX_ELEMENTS
            ) {
            free(elements);
            return R_NilValue;
        }

        // allocate pointer for peaks ///////////////////////////////////////////////////////////////
        int msg_peaks = alloc_peaks(p_l, iso_amount, &m_, &a_, &cc_);
        if (msg_peaks) {
            Rprintf("\nERROR: pointer allocation, error code: %d", msg_peaks);
            free(elements);
            return R_NilValue;
        }

        qsort(elements, element_amount, sizeof(Element), elements_sort_by_isoamount_inc);

        SEXP list_names;
        char* l_n = (char*)calloc((iso_amount + 2) * MAX_NAME_SIZE, sizeof(char));
        if(l_n == NULL){
            Rprintf("\nERROR: cannot allocate memory for list name pointer");
            free(elements);
            free(m_);
            free(a_);
            free(cc_);
            return R_NilValue;
        }
        strncpy(l_n, "mass\0", 5 * sizeof(char));
        strncpy(l_n + MAX_NAME_SIZE, "abundance\0", 10 * sizeof(char));
        size_t index = 0;
        for (ptrdiff_t b = 0; b < element_amount; b++) {
            for (ptrdiff_t bb = 0; bb < (elements+b)->iso_amount; bb++) {
                strncpy(l_n + (index + 2) * MAX_NAME_SIZE, ((elements + b)->isotopes + bb)->isotope, MAX_NAME_SIZE*sizeof(char));
                index++;
            }
        }

        CompoundMulti* monoisotopic = (CompoundMulti*)calloc(1,sizeof(CompoundMulti));
        if(monoisotopic == NULL){
            Rprintf("\nERROR: cannot allocate memory for list name pointer");
            free(elements);
            free(m_);
            free(a_);
            free(cc_);
            free(l_n);
            return R_NilValue;
        }

        calc_monoisotopic(elements, element_amount, monoisotopic);
        long double mono_abundance = monoisotopic->abundance;
        free(monoisotopic);

        int msg = calc_pattern_algo_1(elements, &peak_amount, t, iso_amount,element_amount, mono_abundance,p_l,rtm, &m_, &a_, &cc_);
        if(msg || peak_amount < 1){
            Rprintf("\nERROR: Unable to calculate isotope pattern, exit code %d",msg);
            free(cc_);
            free(l_n);
            free(m_);
            free(a_);
            free(elements);
            return R_NilValue;
        }

        const size_t name_size = MAX_NAME_SIZE * sizeof(char);
        if (rica) {

            SEXP list_names_elements;

            PROTECT(list_names_elements = allocVector(STRSXP, element_amount));

            SEXP out;
            out = PROTECT(allocVector(INTSXP, element_amount));

            for (ptrdiff_t i = 0; i< element_amount; i++) {
                INTEGER(out)[i] = (int)(elements+i)->all_iso_calc_amount;

                char tmp[ MAX_NAME_SIZE ];

                memcpy(tmp, (elements+i)->name, name_size);
                tmp[ MAX_NAME_SIZE - 1] = '\0';

                // #pragma GCC diagnostic push
                // #pragma GCC diagnostic ignored "-Wstringop-truncation"
                //     strncpy(tmp, (elements+i)->name, name_size);
                // #pragma GCC diagnostic pop

                SET_STRING_ELT(list_names_elements, i,  mkChar(tmp));
            }

            setAttrib(out, R_NamesSymbol, list_names_elements);

            free(elements);
            free(m_);
            free(a_);
            free(l_n);
            free(cc_);

            UNPROTECT(2);
            return out;
        }else{

            SEXP iso_pattern = PROTECT(allocVector(VECSXP, 3 + iso_amount));
            SEXP mass_R;
            SEXP a_R;
            double *p_m_R;
            double *p_a_R;
            PROTECT(mass_R = NEW_NUMERIC((long)peak_amount));
            PROTECT(a_R = NEW_NUMERIC((long)peak_amount));
            p_m_R = NUMERIC_POINTER(mass_R);
            p_a_R = NUMERIC_POINTER(a_R);

            for (ptrdiff_t k = 0; k < peak_amount; k ++) {
                p_m_R[k] = *(m_ + k);
                p_a_R[k] = *(a_ + k);
            }

            SET_VECTOR_ELT(iso_pattern, 0, mass_R);
            SET_VECTOR_ELT(iso_pattern, 1, a_R);

            for (ptrdiff_t q = 0; q < iso_amount; q++) {
                SEXP compound_R;
                int *compound;
                PROTECT(compound_R = NEW_INTEGER((long)peak_amount));
                compound = INTEGER_POINTER(compound_R);

                for (ptrdiff_t r = 0; r < peak_amount; r++) {
                    compound[r] = *(cc_ + r * iso_amount + q);
                }
                SET_VECTOR_ELT(iso_pattern, q + 2, compound_R);
                UNPROTECT(1);
            }

            PROTECT(list_names = allocVector(STRSXP, iso_amount + 3));

            for(ptrdiff_t o = 0; o < iso_amount + 2; o++){
                char tmp[ MAX_NAME_SIZE ];

                memcpy(tmp, l_n + o * MAX_NAME_SIZE, name_size);
                tmp[ MAX_NAME_SIZE - 1] = '\0';

                // #pragma GCC diagnostic push
                // #pragma GCC diagnostic ignored "-Wstringop-truncation"
                //     strncpy(tmp, l_n + o * MAX_NAME_SIZE, name_size);
                // #pragma GCC diagnostic pop

                SET_STRING_ELT(list_names, o,  mkChar(tmp));
            }
            SET_STRING_ELT(list_names, iso_amount + 2,  mkChar("NAMES"));
            setAttrib(iso_pattern, R_NamesSymbol, list_names);
            SET_VECTOR_ELT(iso_pattern, iso_amount + 2, list_names);

            free(elements);
            free(m_);
            free(a_);
            free(l_n);
            free(cc_);

            UNPROTECT(4);

#if SHOW_DETAILS == 1
            end = clock();
            cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
            Rprintf("\n    peak_amount: %zu, elapsed_time: %lf, ",peak_amount, cpu_time_used);
#endif
            return iso_pattern;
        }
    }

/*************************

 calculates isotopic pattern with algorithm 1

 sum:           character array of sum formula
 peak_limit:    maximum allowed peak amount
 threshold:     only peaks equal or above this threshold are returned
 iso_list:      character array of elements with their isotopes

 *************************/
SEXP iso_pattern(SEXP sum
                 , SEXP peak_limit
                 , SEXP threshold
                 , SEXP iso_list_elem
                 , SEXP iso_list_iso
                 , SEXP iso_list_mass
                 , SEXP iso_list_abu
                 , SEXP rel_to_mono
                 , SEXP return_iso_calc_amount
                 )
{
    char *s;

    int p_l;
    double t;
    int rtm;
    int rica;

    double* m_;
    double* a_;
    int* cc_;

    PROTECT(sum = AS_CHARACTER(sum));
    PROTECT(peak_limit = AS_INTEGER(peak_limit));
    PROTECT(threshold = AS_NUMERIC(threshold));
    PROTECT(rel_to_mono = AS_INTEGER(rel_to_mono));
    PROTECT(return_iso_calc_amount = AS_INTEGER(return_iso_calc_amount));

    s = R_alloc(strlen(CHARACTER_VALUE(sum)), sizeof(char));
    p_l = INTEGER_VALUE(peak_limit);
    t = NUMERIC_VALUE(threshold);
    rtm = INTEGER_VALUE(rel_to_mono);
    rica = INTEGER_VALUE(return_iso_calc_amount);

    strcpy(s, CHARACTER_VALUE(sum));

    UNPROTECT(5);

    // parse elements ///////////////////////////////////////////////////////////////

    PROTECT(iso_list_elem = AS_CHARACTER(iso_list_elem));
    PROTECT(iso_list_iso = AS_CHARACTER(iso_list_iso));
    PROTECT(iso_list_mass = AS_NUMERIC(iso_list_mass));
    PROTECT(iso_list_abu = AS_NUMERIC(iso_list_abu));

    if (LENGTH( iso_list_elem) != LENGTH(iso_list_iso)
        || LENGTH(iso_list_iso) != LENGTH(iso_list_mass)
        || LENGTH(iso_list_mass) != LENGTH(iso_list_abu)
        ) {
        Rprintf("\ninput vectors not of same dimensions");
        UNPROTECT(4);
        return R_NilValue;
    }

    char list[LENGTH(iso_list_elem) * MAX_NAME_SIZE];
    size_t size_list = sizeof(list);

    char* i_el_ = (char*)malloc(size_list);
    if (i_el_ == NULL) {
        Rprintf("\nERROR: cannot allocate memory for isotope information");
        UNPROTECT(4);
        return R_NilValue;
    }
    char* i_il_ = (char*)malloc(size_list);
    if (i_il_ == NULL) {
        Rprintf("\nERROR: cannot allocate memory for isotope information");
        free(i_el_);
        UNPROTECT(4);
        return R_NilValue;
    }

    for (int i = 0; i < LENGTH(iso_list_elem); i++) {
        strcpy(&i_el_[i*MAX_NAME_SIZE], CHARACTER_VALUE(STRING_ELT(iso_list_elem, i)));
        strcpy(&i_il_[i*MAX_NAME_SIZE], CHARACTER_VALUE(STRING_ELT(iso_list_iso, i)));
    }

    unsigned short element_amount = 0;
    unsigned short iso_amount = 0;
    unsigned short iso_amount_global = (unsigned short)LENGTH(iso_list_elem);

    Element* elements = (Element*)calloc(MAX_ELEMENTS, sizeof(Element));
    if (elements == NULL) {
        Rprintf("\nERROR: cannot allocate memory for elements pointer");
        UNPROTECT(4);
        free(i_el_);
        free(i_il_);
        return R_NilValue;
    }

    if ( parse_sum_formula_vector_R(
                                    elements,
                                    s,
                                    &element_amount,
                                    &iso_amount,
                                    iso_amount_global,
                                    i_el_,
                                    i_il_,
                                    REAL(iso_list_mass),
                                    REAL(iso_list_abu)
                                    )
        ) {
        Rprintf("\nERROR: cannot parse sum formula with the given isolist, ");
        free(elements);
        UNPROTECT(4);
        free(i_el_);
        free(i_il_);
        return R_NilValue;
    }
    UNPROTECT(4);
    free(i_el_);
    free(i_il_);
    qsort(elements, element_amount, sizeof(Element), elements_sort_by_isoamount_inc);

    // end parse elements ///////////////////////////////////////////////////////////////

    if (p_l >= INT_LIMIT || p_l < 1) {
        Rprintf("\ninvalid peak limit");
        return R_NilValue;
    }

#if SHOW_DETAILS == 1
    Rprintf("\n%s, algo1, ",s);
    clock_t start, end;
    double cpu_time_used = 0.0;
    start = clock();
#endif

    size_t  peak_amount = 0;

    if ( rtm > 4 || rtm < 0) {
        Rprintf("\nERROR: wrong value for rtm");
        return R_NilValue;
    }

    if (iso_amount >= MAX_ISO_SIZE
        || iso_amount <= 0
        || element_amount <= 0
        || element_amount >= MAX_ELEMENTS
        ) {
        free(elements);
        return R_NilValue;
    }

    // allocate pointer for peaks ///////////////////////////////////////////////////////////////
    int msg_peaks = alloc_peaks(p_l, iso_amount, &m_, &a_, &cc_);
    if (msg_peaks) {
        Rprintf("\nERROR: pointer allocation, error code: %d", msg_peaks);
        free(elements);
        return R_NilValue;
    }

    //SEXP list_names;
    char* l_n = (char*)calloc((iso_amount + 2) * MAX_NAME_SIZE, sizeof(char));
    if(l_n == NULL){
        Rprintf("\nERROR: cannot allocate memory for list name pointer");
        free(elements);
        free(m_);
        free(a_);
        free(cc_);
        return R_NilValue;
    }
    strncpy(l_n, "mass\0", 5 * sizeof(char));
    strncpy(l_n + MAX_NAME_SIZE, "abundance\0", 10 * sizeof(char));
    size_t index = 0;
    for (ptrdiff_t b = 0; b < element_amount; b++) {
        for (ptrdiff_t bb = 0; bb < (elements+b)->iso_amount; bb++) {
            strncpy(l_n + (index + 2) * MAX_NAME_SIZE, ((elements + b)->isotopes + bb)->isotope, MAX_NAME_SIZE*sizeof(char));
            index++;
        }
    }

    CompoundMulti* monoisotopic = (CompoundMulti*)calloc(1,sizeof(CompoundMulti));
    if(monoisotopic == NULL){
        Rprintf("\nERROR: cannot allocate memory for list name pointer");
        free(elements);
        free(m_);
        free(a_);
        free(cc_);
        free(l_n);
        return R_NilValue;
    }

    calc_monoisotopic(elements, element_amount, monoisotopic);
    long double mono_abundance = monoisotopic->abundance;
    free(monoisotopic);

    int msg = calc_pattern_algo_1(elements, &peak_amount, t, iso_amount,element_amount, mono_abundance,p_l,rtm, &m_, &a_, &cc_);
    if(msg || peak_amount < 1){
        Rprintf("\nERROR: Unable to calculate isotope pattern, exit code %d",msg);
        free(cc_);
        free(l_n);
        free(m_);
        free(a_);
        free(elements);
        return R_NilValue;
    }

    const size_t name_size = MAX_NAME_SIZE * sizeof(char);
    if (rica) {

        SEXP list_names_elements;

        PROTECT(list_names_elements = allocVector(STRSXP, element_amount));

        SEXP out;
        out = PROTECT(allocVector(INTSXP, element_amount));

        for (ptrdiff_t i = 0; i< element_amount; i++) {
            INTEGER(out)[i] = (int)(elements+i)->all_iso_calc_amount;

            char tmp[ MAX_NAME_SIZE ];

            memcpy(tmp, (elements+i)->name, name_size);
            tmp[ MAX_NAME_SIZE - 1] = '\0';

            // #pragma GCC diagnostic push
            // #pragma GCC diagnostic ignored "-Wstringop-truncation"
            //     strncpy(tmp, (elements+i)->name, name_size);
            // #pragma GCC diagnostic pop

            SET_STRING_ELT(list_names_elements, i,  mkChar(tmp));
        }

        setAttrib(out, R_NamesSymbol, list_names_elements);

        free(elements);
        free(m_);
        free(a_);
        free(l_n);
        free(cc_);

        UNPROTECT(2);
        return out;
    }else{

        //SEXP iso_pattern = PROTECT(allocVector(VECSXP, 3 + iso_amount));
        //SEXP mass_R;
        //SEXP a_R;
        //double *p_m_R;
        //double *p_a_R;
        //PROTECT(mass_R = NEW_NUMERIC((long)peak_amount));
        //PROTECT(a_R = NEW_NUMERIC((long)peak_amount));
        //p_m_R = NUMERIC_POINTER(mass_R);
        //p_a_R = NUMERIC_POINTER(a_R);

        //for (ptrdiff_t k = 0; k < peak_amount; k ++) {
            //p_m_R[k] = *(m_ + k);
            //p_a_R[k] = *(a_ + k);
        //}

        //SET_VECTOR_ELT(iso_pattern, 0, mass_R);
        //SET_VECTOR_ELT(iso_pattern, 1, a_R);

        //for (ptrdiff_t q = 0; q < iso_amount; q++) {
            //SEXP compound_R;
            //int *compound;
            //PROTECT(compound_R = NEW_INTEGER((long)peak_amount));
            //compound = INTEGER_POINTER(compound_R);

            //for (ptrdiff_t r = 0; r < peak_amount; r++) {
                //compound[r] = *(cc_ + r * iso_amount + q);
            //}

            //SET_VECTOR_ELT(iso_pattern, q + 2, compound_R);
            //UNPROTECT(1);
        //}

        //PROTECT(list_names = allocVector(STRSXP, iso_amount + 3));
        //for(ptrdiff_t o = 0; o < iso_amount + 2; o++){
            //char tmp[ MAX_NAME_SIZE ];
            //strncpy(tmp, l_n + o * MAX_NAME_SIZE, MAX_NAME_SIZE * sizeof(char));
            //SET_STRING_ELT(list_names, o,  mkChar(tmp));
        //}
        //SET_STRING_ELT(list_names, iso_amount + 2,  mkChar("NAMES"));
        //setAttrib(iso_pattern, R_NamesSymbol, list_names);
        //SET_VECTOR_ELT(iso_pattern, iso_amount + 2, list_names);

        //free(m_);
        //free(a_);
        //free(l_n);
        //free(cc_);
        //free(elements);
        //UNPROTECT(4);
        //return iso_pattern;


        int nx = (int)peak_amount, ny = iso_amount + 2;
        double *rans;
        //double tmp;
        SEXP ans;
        SEXP dim;
        SEXP dimnames;
        //SEXP dimx;
        SEXP dimy;

        PROTECT(ans = allocVector(REALSXP, nx*ny));
        //PROTECT(dimx = allocVector(STRSXP, nx));
        PROTECT(dimy = allocVector(STRSXP, ny));
        rans = REAL(ans);

        //char* str = (char*)malloc(MAX_NAME_SIZE * sizeof(char));

        for(int i = 0; i < nx; i++) {
            //tmp = i;

			//snprintf(str, MAX_NAME_SIZE, "%d", i + 1);
            //SET_STRING_ELT(dimx, i, mkChar(str));
            //SET_STRING_ELT(dimx, i, mkChar("x"));

            rans[i] = m_[i];
			rans[i + nx] = a_[i];

            for(int j = 2; j < ny; j++){
				if (i == 0) {
					SET_STRING_ELT(dimy, 0, mkChar("mass"));
					SET_STRING_ELT(dimy, 1, mkChar("abundance"));
					char tmp_s[ MAX_NAME_SIZE ];

                    memcpy(tmp_s, l_n + j* MAX_NAME_SIZE, MAX_NAME_SIZE * sizeof(char));
                    tmp_s[ MAX_NAME_SIZE - 1] = '\0';

                    // #pragma GCC diagnostic push
                    // #pragma GCC diagnostic ignored "-Wstringop-truncation"
                    //     strncpy(tmp_s, l_n + j* MAX_NAME_SIZE, MAX_NAME_SIZE * sizeof(char));
                    // #pragma GCC diagnostic pop

					SET_STRING_ELT(dimy, j,  mkChar(tmp_s));
				}
				rans[i + nx*j] = cc_[i * iso_amount + j - 2];
			}
        }

        PROTECT(dim = allocVector(INTSXP, 2));
        INTEGER(dim)[0] = nx;
        INTEGER(dim)[1] = ny;
        setAttrib(ans, R_DimSymbol, dim);

        PROTECT(dimnames = allocVector(VECSXP, 2));
        //SET_VECTOR_ELT(dimnames, 0, dimx);
        SET_VECTOR_ELT(dimnames, 1, dimy);
        setAttrib(ans, R_DimNamesSymbol, dimnames);

        UNPROTECT(4);

        free(m_);
        free(a_);
        free(cc_);
        free(l_n);
        free(elements);


#if SHOW_DETAILS == 1
        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
        Rprintf("\n    peak_amount: %zu, elapsed_time: %lf, ",peak_amount, cpu_time_used);
#endif

        return(ans);
    }
}

/*************************

 calculates isotopic pattern with algorithm 2

 sum:           character array of sum formula
 peak_limit:    maximum allowed peak amount
 thres:         only peaks equal or above this threshold are returned
 iso_list:      character array of elements with their isotopes

 *************************/
SEXP iso_pattern_2(SEXP sum
                 , SEXP peak_limit
                 , SEXP threshold
                 , SEXP iso_list_elem
                 , SEXP iso_list_iso
                 , SEXP iso_list_mass
                 , SEXP iso_list_abu
                 , SEXP rel_to_mono
                 , SEXP return_iso_calc_amount
                 )
{
    char *s;

    int p_l;
    double t;
    int rtm;
    int rica;

    double* m_;
    double* a_;
    int* cc_;

    PROTECT(sum = AS_CHARACTER(sum));
    PROTECT(peak_limit = AS_INTEGER(peak_limit));
    PROTECT(threshold = AS_NUMERIC(threshold));
    PROTECT(rel_to_mono = AS_INTEGER(rel_to_mono));
    PROTECT(return_iso_calc_amount = AS_INTEGER(return_iso_calc_amount));


    s = R_alloc(strlen(CHARACTER_VALUE(sum)), sizeof(char));
    p_l = INTEGER_VALUE(peak_limit);
    t = NUMERIC_VALUE(threshold);
    rtm = INTEGER_VALUE(rel_to_mono);
    rica = INTEGER_VALUE(return_iso_calc_amount);

    strcpy(s, CHARACTER_VALUE(sum));

    UNPROTECT(5);

    // parse elements ///////////////////////////////////////////////////////////////

    PROTECT(iso_list_elem = AS_CHARACTER(iso_list_elem));
    PROTECT(iso_list_iso = AS_CHARACTER(iso_list_iso));
    PROTECT(iso_list_mass = AS_NUMERIC(iso_list_mass));
    PROTECT(iso_list_abu = AS_NUMERIC(iso_list_abu));

    if (LENGTH( iso_list_elem) != LENGTH(iso_list_iso)
        || LENGTH(iso_list_iso) != LENGTH(iso_list_mass)
        || LENGTH(iso_list_mass) != LENGTH(iso_list_abu)
        ) {
        Rprintf("\ninput vectors not of same dimensions");
        UNPROTECT(4);
        return R_NilValue;
    }

    char list[LENGTH(iso_list_elem) * MAX_NAME_SIZE];
    size_t size_list = sizeof(list);

    char* i_el_ = (char*)malloc(size_list);
    if (i_el_ == NULL) {
        Rprintf("\nERROR: cannot allocate memory for isotope information");
        UNPROTECT(4);
        return R_NilValue;
    }
    char* i_il_ = (char*)malloc(size_list);
    if (i_il_ == NULL) {
        Rprintf("\nERROR: cannot allocate memory for isotope information");
        free(i_el_);
        UNPROTECT(4);
        return R_NilValue;
    }

    for (int i = 0; i < LENGTH(iso_list_elem); i++) {
        strcpy(&i_el_[i*MAX_NAME_SIZE], CHARACTER_VALUE(STRING_ELT(iso_list_elem, i)));
        strcpy(&i_il_[i*MAX_NAME_SIZE], CHARACTER_VALUE(STRING_ELT(iso_list_iso, i)));
    }

    unsigned short element_amount = 0;
    unsigned short iso_amount = 0;
    unsigned short iso_amount_global = (unsigned short)LENGTH(iso_list_elem);

    Element* elements = (Element*)calloc(MAX_ELEMENTS, sizeof(Element));
    if (elements == NULL) {
        Rprintf("\nERROR: cannot allocate memory for elements pointer");
        UNPROTECT(4);
        free(i_el_);
        free(i_il_);
        return R_NilValue;
    }

    if ( parse_sum_formula_vector_R(
                                    elements,
                                    s,
                                    &element_amount,
                                    &iso_amount,
                                    iso_amount_global,
                                    i_el_,
                                    i_il_,
                                    REAL(iso_list_mass),
                                    REAL(iso_list_abu)
                                    )
        ) {
        Rprintf("\nERROR: cannot parse sum formula with the given isolist, ");
        free(elements);
        UNPROTECT(4);
        free(i_el_);
        free(i_il_);
        return R_NilValue;
    }
    UNPROTECT(4);
    free(i_el_);
    free(i_il_);
    qsort(elements, element_amount, sizeof(Element), elements_sort_by_isoamount_inc);

    // end parse elements ///////////////////////////////////////////////////////////////


    if (p_l >= INT_LIMIT || p_l < 1) {

        Rprintf("\ninvalid peak limit");
        return R_NilValue;
    }

#if SHOW_DETAILS == 1
    Rprintf("\n%s, algo2, ",s);
    clock_t start, end;
    double cpu_time_used = 0.0;
    start = clock();
#endif

    long double max_a = 0.0;
    size_t  peak_amount = 0;
    int iso_count_stats = 0;

    if ( rtm > 4 || rtm < 0) {
        Rprintf("\nERROR:  wrong value for rtm");
        return R_NilValue;
    }

    if (iso_amount >= MAX_ISO_SIZE
        || iso_amount <= 0
        || element_amount <= 0
        || element_amount >= MAX_ELEMENTS
        ) {
        free(elements);
        return R_NilValue;
    }

    // allocate pointer for peaks ///////////////////////////////////////////////////////////////
    int msg_peaks = alloc_peaks(p_l, iso_amount, &m_, &a_, &cc_);
    if (msg_peaks) {
        Rprintf("\nERROR: pointer allocation, error code: %d", msg_peaks);
        free(elements);
        return R_NilValue;
    }

    // calc pattern algo 2 //////////////////////////////////////////////////////////////////////
    int msg = calc_pattern_algo_2(&max_a, elements, element_amount, t, &peak_amount, p_l, &iso_count_stats,rtm, &m_, &a_, &cc_);

    if(msg || peak_amount < 1){
        Rprintf("\nERROR: could not create combinations, error code: %d", msg);
        free(elements);
        free(m_);
        free(a_);
        free(cc_);
        return R_NilValue;
    }

    // return data to R ////////////////////////////////////////////////////////////////////////
    if (rica) {
        SEXP out = PROTECT(allocVector(INTSXP, 1));
        INTEGER(out)[0] = iso_count_stats;
        free(elements);
        free(m_);
        free(a_);
        free(cc_);
        UNPROTECT(1);
        return out;
    }else{

        SEXP iso_pattern = PROTECT(allocVector(VECSXP, 3 + iso_amount));
        SEXP mass_R;
        SEXP a_R;

        PROTECT(mass_R = NEW_NUMERIC((long)peak_amount));
        PROTECT(a_R = NEW_NUMERIC((long)peak_amount));
        double *p_m_R = NUMERIC_POINTER(mass_R);
        double *p_a_R = NUMERIC_POINTER(a_R);

        size_t c = 0;
        for (ptrdiff_t i = 0; i < peak_amount; i++) {
            double tmp = *(a_ + i);
            if (rtm != 2) {
                *(a_ + i) = (double)((100/ max_a) * tmp);
            }

            if (*(a_ + i) > t) {
                p_m_R[c] = *(m_ + i);

                if (rtm == 3 || rtm == 4) {
                    p_a_R[c] = tmp;
                }else{
                    p_a_R[c] = *(a_ + i);
                }

                c++;
            }
        }
        int r = 0;
        for (ptrdiff_t j = 0; j < iso_amount; j++) {
            SEXP compound_R;
            double *compound;
            PROTECT(compound_R = NEW_NUMERIC((long)c));
            compound = NUMERIC_POINTER(compound_R);
            r = 0;
            for (ptrdiff_t k = 0; k < peak_amount; k++) {
                if (*(a_ + k) > t) {
                    *(compound + r) = *(cc_ + k * iso_amount + j);
                    r++;
                }
            }
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
        if (names == NULL) {
            Rprintf("\nERROR: allocation names pointer failed");
            free(elements);
            free(m_);
            free(a_);
            free(cc_);
			UNPROTECT(3);
            return R_NilValue;
        }

        strncpy(names, "mass\0", 5 * sizeof(char));
        strncpy(names + MAX_NAME_SIZE, "abundance\0", 10 * sizeof(char));

        for (ptrdiff_t  e = 0; e < element_amount; e++) {
            for (ptrdiff_t  ee = 0; ee < (elements + e)->iso_amount; ee++) {
                strncpy(names + n_c * MAX_NAME_SIZE, ((elements + e)->isotopes +ee)->isotope, MAX_NAME_SIZE * sizeof(char));
                n_c++;
            }
        }

        PROTECT(list_names = allocVector(STRSXP, iso_amount + 3));
        for(ptrdiff_t o = 0; o < n_c; o++){
            char tmp[ MAX_NAME_SIZE ];

            memcpy(tmp, names + o * MAX_NAME_SIZE, MAX_NAME_SIZE * sizeof(char));
            tmp[ MAX_NAME_SIZE - 1] = '\0';

            // #pragma GCC diagnostic push
            // #pragma GCC diagnostic ignored "-Wstringop-truncation"
            //     strncpy(tmp, names + o * MAX_NAME_SIZE, MAX_NAME_SIZE * sizeof(char));
            // #pragma GCC diagnostic pop

            SET_STRING_ELT(list_names, o,  mkChar(tmp));
        }
        SET_STRING_ELT(list_names, n_c,  mkChar("NAMES"));
        setAttrib(iso_pattern, R_NamesSymbol, list_names);
        SET_VECTOR_ELT(iso_pattern, n_c, list_names);

        free(elements);
        free(m_);
        free(a_);
        free(cc_);
        free(names);
        
        UNPROTECT(4);

#if SHOW_DETAILS == 1
        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
        Rprintf("\n    peak_amount: %zu, elapsed_time: %lf, ",c, cpu_time_used);
#endif
        return iso_pattern;
    }
}


/*************************

 calculates isotopic pattern with algorithm 3

 sum:           character array of sum formula
 peak_limit:    maximum allowed peak amount
 threshold:     only peaks equal or above this threshold are returned
 iso_list:      character array of elements with their isotopes

 *************************/
SEXP iso_pattern_3(SEXP sum
                   , SEXP peak_limit
                   , SEXP threshold
                   , SEXP iso_list_elem
                   , SEXP iso_list_iso
                   , SEXP iso_list_mass
                   , SEXP iso_list_abu
                   , SEXP rel_to_mono
                   , SEXP return_iso_calc_amount
                   )
{
    char *s;

    int p_l;
    double t;
    int rtm;
    int rica;

    double* m_;
    double* a_;
    int* cc_;

    PROTECT(sum = AS_CHARACTER(sum));
    PROTECT(peak_limit = AS_INTEGER(peak_limit));
    PROTECT(threshold = AS_NUMERIC(threshold));
    PROTECT(rel_to_mono = AS_INTEGER(rel_to_mono));
    PROTECT(return_iso_calc_amount = AS_INTEGER(return_iso_calc_amount));


    s = R_alloc(strlen(CHARACTER_VALUE(sum)), sizeof(char));
    p_l = INTEGER_VALUE(peak_limit);
    t = NUMERIC_VALUE(threshold);
    rtm = INTEGER_VALUE(rel_to_mono);
    rica = INTEGER_VALUE(return_iso_calc_amount);

    strcpy(s, CHARACTER_VALUE(sum));

    UNPROTECT(5);

    // parse elements ///////////////////////////////////////////////////////////////
    PROTECT(iso_list_elem = AS_CHARACTER(iso_list_elem));
    PROTECT(iso_list_iso = AS_CHARACTER(iso_list_iso));
    PROTECT(iso_list_mass = AS_NUMERIC(iso_list_mass));
    PROTECT(iso_list_abu = AS_NUMERIC(iso_list_abu));

    if (LENGTH( iso_list_elem) != LENGTH(iso_list_iso)
        || LENGTH(iso_list_iso) != LENGTH(iso_list_mass)
        || LENGTH(iso_list_mass) != LENGTH(iso_list_abu)
        ) {
        Rprintf("\ninput vectors not of same dimensions");
        UNPROTECT(4);
        return R_NilValue;
    }

    char list[LENGTH(iso_list_elem) * MAX_NAME_SIZE];
    size_t size_list = sizeof(list);

    char* i_el_ = (char*)malloc(size_list);
    if (i_el_ == NULL) {
        Rprintf("\nERROR: cannot allocate memory for isotope information");
        UNPROTECT(4);
        return R_NilValue;
    }
    char* i_il_ = (char*)malloc(size_list);
    if (i_il_ == NULL) {
        Rprintf("\nERROR: cannot allocate memory for isotope information");
        free(i_el_);
        UNPROTECT(4);
        return R_NilValue;
    }

    for (int i = 0; i < LENGTH(iso_list_elem); i++) {
        strcpy(&i_el_[i*MAX_NAME_SIZE], CHARACTER_VALUE(STRING_ELT(iso_list_elem, i)));
        strcpy(&i_il_[i*MAX_NAME_SIZE], CHARACTER_VALUE(STRING_ELT(iso_list_iso, i)));
    }

    unsigned short element_amount = 0;
    unsigned short iso_amount = 0;
    unsigned short iso_amount_global = (unsigned short)LENGTH(iso_list_elem);

    Element* elements = (Element*)calloc(MAX_ELEMENTS, sizeof(Element));
    if (elements == NULL) {
        Rprintf("\nERROR: cannot allocate memory for elements pointer");
        UNPROTECT(4);
        free(i_el_);
        free(i_il_);
        return R_NilValue;
    }

    if ( parse_sum_formula_vector_R(
                                    elements,
                                    s,
                                    &element_amount,
                                    &iso_amount,
                                    iso_amount_global,
                                    i_el_,
                                    i_il_,
                                    REAL(iso_list_mass),
                                    REAL(iso_list_abu)
                                    )
        ) {
        Rprintf("\nERROR: cannot parse sum formula with the given isolist, ");
        free(elements);
        UNPROTECT(4);
        free(i_el_);
        free(i_il_);
        return R_NilValue;
    }
    UNPROTECT(4);
    free(i_el_);
    free(i_il_);
    qsort(elements, element_amount, sizeof(Element), elements_sort_by_isoamount_inc);

    // end parse elements ///////////////////////////////////////////////////////////////


    if (p_l >= INT_LIMIT || p_l < 1) {

        Rprintf("\ninvalid peak limit");
        return R_NilValue;
    }

#if SHOW_DETAILS == 1
    Rprintf("\n%s, algo3, ",s);
    clock_t start, end;
    double cpu_time_used = 0.0;
    start = clock();
#endif
    size_t  peak_amount = 0;

    if ( rtm > 4 || rtm < 0) {
        Rprintf("\nERROR: wrong value for rtm");
        return R_NilValue;
    }
    if (iso_amount >= MAX_ISO_SIZE
        || iso_amount <= 0
        || element_amount <= 0
        || element_amount >= MAX_ELEMENTS
        ) {
        free(elements);
        return R_NilValue;
    }

    // allocate pointer for peaks ///////////////////////////////////////////////////////////////
    int msg_peaks = alloc_peaks(p_l, iso_amount, &m_, &a_, &cc_);
    if (msg_peaks) {
        Rprintf("\nERROR: pointer allocation, error code: %d", msg_peaks);
        free(elements);
        return R_NilValue;
    }

    SEXP list_names;
    char* l_n = (char*)calloc((iso_amount + 2) * MAX_NAME_SIZE, sizeof(char));
    if(l_n == NULL){
        Rprintf("\nERROR: cannot allocate memory column names");
        free(elements);
        free(m_);
        free(a_);
        free(cc_);
        return R_NilValue;
    }
    strncpy(l_n, "mass\0", 5 * sizeof(char));
    strncpy(l_n + MAX_NAME_SIZE, "abundance\0", 10 * sizeof(char));
    size_t index = 0;
    for (ptrdiff_t b = 0; b < element_amount; b++) {
        for (ptrdiff_t bb = 0; bb < (elements+b)->iso_amount; bb++) {
            strncpy(l_n + (index + 2) * MAX_NAME_SIZE, ((elements + b)->isotopes + bb)->isotope, MAX_NAME_SIZE*sizeof(char));
            index++;
        }
    }

    long double mono_abundance = -1.0;
    for (ptrdiff_t  i = 0; i < element_amount; i++) {
        mono_abundance *= pow((elements + i)->isotopes[0].abundance, (elements + i)->amount);
    }

    int msg = calc_pattern_algo_3( elements, &peak_amount, t, iso_amount, element_amount, mono_abundance, p_l, l_n, rtm, &m_, &a_, &cc_);
    if(msg || peak_amount == 0){
        Rprintf("\nERROR: cannot combine combinations, exit code: %d", msg);
        free(elements);
        free(m_);
        free(a_);
        free(cc_);
        free(l_n);
        return R_NilValue;

    }

    // return values to R ////////////////////////////////////////////////////////////////////
    if (rica) {

        SEXP list_names_elements;
        PROTECT(list_names_elements = allocVector(STRSXP, element_amount));
        SEXP out;
        out = PROTECT(allocVector(INTSXP, element_amount));

        for (ptrdiff_t i = 0; i< element_amount; i++) {
            INTEGER(out)[i] = (elements+i)->all_iso_calc_amount;

            char tmp[ MAX_NAME_SIZE ];


            memcpy(tmp, (elements+i)->name, MAX_NAME_SIZE * sizeof(char));
            tmp[ MAX_NAME_SIZE - 1] = '\0';

            // #pragma GCC diagnostic push
            // #pragma GCC diagnostic ignored "-Wstringop-truncation"
            //     strncpy(tmp, (elements+i)->name, MAX_NAME_SIZE * sizeof(char));
            // #pragma GCC diagnostic pop

            SET_STRING_ELT(list_names_elements, i,  mkChar(tmp));
        }

        setAttrib(out, R_NamesSymbol, list_names_elements);

        free(m_);
        free(a_);
        free(l_n);
        free(cc_);
        free(elements);
        UNPROTECT(2);
        return out;
    }else{

        SEXP iso_pattern = PROTECT(allocVector(VECSXP, 3 + iso_amount));
        SEXP mass_R;
        SEXP a_R;
        double *p_m_R;
        double *p_a_R;
        PROTECT(mass_R = NEW_NUMERIC((long)peak_amount));
        PROTECT(a_R = NEW_NUMERIC((long)peak_amount));
        p_m_R = NUMERIC_POINTER(mass_R);
        p_a_R = NUMERIC_POINTER(a_R);

        for (ptrdiff_t k = 0; k < peak_amount; k ++) {
            p_m_R[k] = *(m_ + k);
            p_a_R[k] = *(a_ + k);
        }

        SET_VECTOR_ELT(iso_pattern, 0, mass_R);
        SET_VECTOR_ELT(iso_pattern, 1, a_R);

        for (ptrdiff_t q = 0; q < iso_amount; q++) {
            SEXP compound_R;
            int *compound;
            PROTECT(compound_R = NEW_INTEGER((long)peak_amount));
            compound = INTEGER_POINTER(compound_R);

            for (ptrdiff_t r = 0; r < peak_amount; r++) {
                compound[r] = *(cc_ + r * iso_amount + q);
            }

            SET_VECTOR_ELT(iso_pattern, q + 2, compound_R);
            UNPROTECT(1);
        }

        PROTECT(list_names = allocVector(STRSXP, iso_amount + 3));
        for(ptrdiff_t o = 0; o < iso_amount + 2; o++){
            char tmp[ MAX_NAME_SIZE ];


            memcpy(tmp, l_n + o * MAX_NAME_SIZE, MAX_NAME_SIZE * sizeof(char));
            tmp[ MAX_NAME_SIZE - 1] = '\0';

            // #pragma GCC diagnostic push
            // #pragma GCC diagnostic ignored "-Wstringop-truncation"
            //     strncpy(tmp, l_n + o * MAX_NAME_SIZE, MAX_NAME_SIZE * sizeof(char));
            // #pragma GCC diagnostic pop

            SET_STRING_ELT(list_names, o,  mkChar(tmp));
        }
        SET_STRING_ELT(list_names, iso_amount + 2,  mkChar("NAMES"));
        setAttrib(iso_pattern, R_NamesSymbol, list_names);
        SET_VECTOR_ELT(iso_pattern, iso_amount + 2, list_names);

        free(m_);
        free(a_);
        free(l_n);
        free(cc_);
        free(elements);
        UNPROTECT(4);

#if SHOW_DETAILS == 1
        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
        Rprintf("\n    peak_amount: %zu, elapsed_time: %lf, ",peak_amount, cpu_time_used);
#endif
        return iso_pattern;
    }
}


/*************************

 calculates isotopic pattern with algorithm 4

 sum:           character array of sum formula
 peak_limit:    maximum allowed peak amount
 threshold:     only peaks equal or above this threshold are returned
 iso_list:      character array of elements with their isotopes

 *************************/
SEXP iso_pattern_4(SEXP sum
                   , SEXP peak_limit
                   , SEXP threshold
                   , SEXP iso_list_elem
                   , SEXP iso_list_iso
                   , SEXP iso_list_mass
                   , SEXP iso_list_abu
                   , SEXP rel_to_mono
                   , SEXP return_iso_calc_amount
                   )
{
    char *s;

    int p_l;
    double t;
    int rtm;
    int rica;

    double* m_;
    double* a_;
    int* cc_;

    PROTECT(sum = AS_CHARACTER(sum));
    PROTECT(peak_limit = AS_INTEGER(peak_limit));
    PROTECT(threshold = AS_NUMERIC(threshold));
    PROTECT(rel_to_mono = AS_INTEGER(rel_to_mono));
    PROTECT(return_iso_calc_amount = AS_INTEGER(return_iso_calc_amount));


    s = R_alloc(strlen(CHARACTER_VALUE(sum)), sizeof(char));
    p_l = INTEGER_VALUE(peak_limit);
    t = NUMERIC_VALUE(threshold);
    rtm = INTEGER_VALUE(rel_to_mono);
    rica = INTEGER_VALUE(return_iso_calc_amount);

    strcpy(s, CHARACTER_VALUE(sum));

    UNPROTECT(5);

    // parse elements ///////////////////////////////////////////////////////////////
    PROTECT(iso_list_elem = AS_CHARACTER(iso_list_elem));
    PROTECT(iso_list_iso = AS_CHARACTER(iso_list_iso));
    PROTECT(iso_list_mass = AS_NUMERIC(iso_list_mass));
    PROTECT(iso_list_abu = AS_NUMERIC(iso_list_abu));

    if (LENGTH( iso_list_elem) != LENGTH(iso_list_iso)
        || LENGTH(iso_list_iso) != LENGTH(iso_list_mass)
        || LENGTH(iso_list_mass) != LENGTH(iso_list_abu)
        ) {
        Rprintf("\ninput vectors not of same dimensions");
        UNPROTECT(4);
        return R_NilValue;
    }

    char list[LENGTH(iso_list_elem) * MAX_NAME_SIZE];
    size_t size_list = sizeof(list);

    char* i_el_ = (char*)malloc(size_list);
    if (i_el_ == NULL) {
        Rprintf("\nERROR: cannot allocate memory for isotope information");
        UNPROTECT(4);
        return R_NilValue;
    }
    char* i_il_ = (char*)malloc(size_list);
    if (i_il_ == NULL) {
        Rprintf("\nERROR: cannot allocate memory for isotope information");
        free(i_el_);
        UNPROTECT(4);
        return R_NilValue;
    }

    for (int i = 0; i < LENGTH(iso_list_elem); i++) {
        strcpy(&i_el_[i*MAX_NAME_SIZE], CHARACTER_VALUE(STRING_ELT(iso_list_elem, i)));
        strcpy(&i_il_[i*MAX_NAME_SIZE], CHARACTER_VALUE(STRING_ELT(iso_list_iso, i)));
    }

    unsigned short element_amount = 0;
    unsigned short iso_amount = 0;
    unsigned short iso_amount_global = (unsigned short)LENGTH(iso_list_elem);

    Element* elements = (Element*)calloc(MAX_ELEMENTS, sizeof(Element));
    if (elements == NULL) {
        Rprintf("\nERROR: cannot allocate memory for elements pointer");
        UNPROTECT(4);
        free(i_el_);
        free(i_il_);
        return R_NilValue;
    }

    if ( parse_sum_formula_vector_R(
                                    elements,
                                    s,
                                    &element_amount,
                                    &iso_amount,
                                    iso_amount_global,
                                    i_el_,
                                    i_il_,
                                    REAL(iso_list_mass),
                                    REAL(iso_list_abu)
                                    )
        ) {
        Rprintf("\nERROR: cannot parse sum formula with the given isolist, ");
        free(elements);
        UNPROTECT(4);
        free(i_el_);
        free(i_il_);
        return R_NilValue;
    }
    UNPROTECT(4);
    free(i_el_);
    free(i_il_);
    qsort(elements, element_amount, sizeof(Element), elements_sort_by_isoamount_inc);

    // end parse elements ///////////////////////////////////////////////////////////////


    if (p_l >= INT_LIMIT || p_l < 1) {
        Rprintf("\ninvalid peak limit");
        return R_NilValue;
    }

#if SHOW_DETAILS == 1
    Rprintf("\n%s, algo4, ",s);
    clock_t start, end;
    double cpu_time_used = 0.0;
    start = clock();
#endif

    size_t  peak_amount = 0;

    if ( rtm > 4 || rtm < 0) {
        Rprintf("\nERROR: wrong value for rtm");
        return R_NilValue;
    }

    if (iso_amount >= MAX_ISO_SIZE
        || iso_amount <= 0
        || element_amount <= 0
        || element_amount >= MAX_ELEMENTS
        ) {
        free(elements);
        return R_NilValue;
    }

    // allocate pointer for peaks ///////////////////////////////////////////////////////////////
    int msg_peaks = alloc_peaks(p_l, iso_amount, &m_, &a_, &cc_);
    if (msg_peaks) {
        Rprintf("\nERROR: pointer allocation, error code: %d", msg_peaks);
        free(elements);
        return R_NilValue;
    }

    //SEXP list_names;
    char* l_n = (char*)calloc((iso_amount + 2) * MAX_NAME_SIZE, sizeof(char));
    if(l_n == NULL){
        Rprintf("\nERROR: cannot allocate memory column names");
        free(elements);
        free(m_);
        free(a_);
        free(cc_);
        return R_NilValue;
    }
    strncpy(l_n, "mass\0", 5 * sizeof(char));
    strncpy(l_n + MAX_NAME_SIZE, "abundance\0", 10 * sizeof(char));
    size_t index = 0;
    for (ptrdiff_t b = 0; b < element_amount; b++) {
        for (ptrdiff_t bb = 0; bb < (elements+b)->iso_amount; bb++) {
            strncpy(l_n + (index + 2) * MAX_NAME_SIZE, ((elements + b)->isotopes + bb)->isotope, MAX_NAME_SIZE*sizeof(char));
            index++;
        }
    }

    long double mono_abundance = -1.0;
    for (ptrdiff_t  i = 0; i < element_amount; i++) {
        mono_abundance *= pow((elements + i)->isotopes[0].abundance, (elements + i)->amount);
    }

    long double max_a = 0;
    int msg = calc_pattern_algo_4(elements, &peak_amount, t, iso_amount, element_amount, mono_abundance, p_l, l_n, rtm, &m_, &a_, &cc_, &max_a);
    if(msg || peak_amount == 0){
        Rprintf("\nERROR: cannot combine combinations, exit code: %d", msg);
        free(elements);
        free(m_);
        free(a_);
        free(cc_);
        free(l_n);
        return R_NilValue;

    }

    // return values to R ////////////////////////////////////////////////////////////////////
    if (rica) {

        SEXP list_names_elements;

        int container_amount = 0;
        for (int i = 0; i < element_amount; i++) {
            if ((elements+i)->all_iso_calc_amount > 0) {
                container_amount++;
            }
        }

        PROTECT(list_names_elements = allocVector(STRSXP, container_amount));
        SEXP out;
        out = PROTECT(allocVector(INTSXP, container_amount));

        int j = 0;
        char tmp[ MAX_NAME_SIZE ];
        for (ptrdiff_t i = 0; i < element_amount; i++) {

            if ((elements+i)->all_iso_calc_amount > 0) {
                if (j > 0) {
                    SET_STRING_ELT(list_names_elements, j - 1,  mkChar(tmp));
                }

                INTEGER(out)[j] = (elements+i)->all_iso_calc_amount;

                memcpy(tmp, (elements+i)->name, MAX_NAME_SIZE * sizeof(char));
                tmp[ MAX_NAME_SIZE - 1] = '\0';

                // #pragma GCC diagnostic push
                // #pragma GCC diagnostic ignored "-Wstringop-truncation"
                //     strncpy(tmp, (elements+i)->name, MAX_NAME_SIZE * sizeof(char));
                // #pragma GCC diagnostic pop

                j++;
            }else if(j > 0){
                strcat(tmp, (elements + i)->name);
            }

            if (i == element_amount - 1 && j > 0 && j <=     container_amount) {
                SET_STRING_ELT(list_names_elements, j - 1,  mkChar(tmp));
            }
        }

        setAttrib(out, R_NamesSymbol, list_names_elements);

        free(m_);
        free(a_);
        free(l_n);
        free(cc_);
        free(elements);
        UNPROTECT(2);
        return out;
    }else{

        //SEXP iso_pattern = PROTECT(allocVector(VECSXP, 3 + iso_amount));
        //SEXP mass_R;
        //SEXP a_R;
        //double *p_m_R;
        //double *p_a_R;
        //PROTECT(mass_R = NEW_NUMERIC((long)peak_amount));
        //PROTECT(a_R = NEW_NUMERIC((long)peak_amount));
        //p_m_R = NUMERIC_POINTER(mass_R);
        //p_a_R = NUMERIC_POINTER(a_R);

        //for (ptrdiff_t k = 0; k < peak_amount; k ++) {
            //p_m_R[k] = *(m_ + k);
            //p_a_R[k] = *(a_ + k);
        //}

        //SET_VECTOR_ELT(iso_pattern, 0, mass_R);
        //SET_VECTOR_ELT(iso_pattern, 1, a_R);

        //for (ptrdiff_t q = 0; q < iso_amount; q++) {
            //SEXP compound_R;
            //int *compound;
            //PROTECT(compound_R = NEW_INTEGER((long)peak_amount));
            //compound = INTEGER_POINTER(compound_R);

            //for (ptrdiff_t r = 0; r < peak_amount; r++) {
                //compound[r] = *(cc_ + r * iso_amount + q);
            //}

            //SET_VECTOR_ELT(iso_pattern, q + 2, compound_R);
            //UNPROTECT(1);
        //}

        //PROTECT(list_names = allocVector(STRSXP, iso_amount + 3));
        //for(ptrdiff_t o = 0; o < iso_amount + 2; o++){
            //char tmp[ MAX_NAME_SIZE ];
            //strncpy(tmp, l_n + o * MAX_NAME_SIZE, MAX_NAME_SIZE * sizeof(char));
            //SET_STRING_ELT(list_names, o,  mkChar(tmp));
        //}
        //SET_STRING_ELT(list_names, iso_amount + 2,  mkChar("NAMES"));
        //setAttrib(iso_pattern, R_NamesSymbol, list_names);
        //SET_VECTOR_ELT(iso_pattern, iso_amount + 2, list_names);

        //free(m_);
        //free(a_);
        //free(l_n);
        //free(cc_);
        //free(elements);
        //UNPROTECT(4);

        //return iso_pattern;


        int nx = (int)peak_amount, ny = iso_amount + 2;
        double *rans;
        //double tmp;
        SEXP ans;
        SEXP dim;
        SEXP dimnames;
        //SEXP dimx;
        SEXP dimy;

        PROTECT(ans = allocVector(REALSXP, nx*ny));
        //PROTECT(dimx = allocVector(STRSXP, nx));
        PROTECT(dimy = allocVector(STRSXP, ny));
        rans = REAL(ans);

        //char* str = (char*)malloc(MAX_NAME_SIZE * sizeof(char));

        for(int i = 0; i < nx; i++) {
            //tmp = i;

			//snprintf(str, MAX_NAME_SIZE, "%d", i + 1);
            //SET_STRING_ELT(dimx, i, mkChar("x"));

            rans[i] = m_[i];
			rans[i + nx] = a_[i];

            for(int j = 2; j < ny; j++){
					if (i == 0) {
						SET_STRING_ELT(dimy, 0, mkChar("mass"));
						SET_STRING_ELT(dimy, 1, mkChar("abundance"));
						char tmp_s[ MAX_NAME_SIZE ];

                        memcpy(tmp_s, l_n + j* MAX_NAME_SIZE, MAX_NAME_SIZE * sizeof(char));
                        tmp_s[ MAX_NAME_SIZE - 1] = '\0';

                        // #pragma GCC diagnostic push
                        // #pragma GCC diagnostic ignored "-Wstringop-truncation"
                        //     strncpy(tmp_s, l_n + j* MAX_NAME_SIZE, MAX_NAME_SIZE * sizeof(char));
                        // #pragma GCC diagnostic pop

						SET_STRING_ELT(dimy, j,  mkChar(tmp_s));
					}
					rans[i + nx*j] = cc_[i * iso_amount + j - 2];
				}
        }

        PROTECT(dim = allocVector(INTSXP, 2));
        INTEGER(dim)[0] = nx;
        INTEGER(dim)[1] = ny;
        setAttrib(ans, R_DimSymbol, dim);

        PROTECT(dimnames = allocVector(VECSXP, 2));
        //SET_VECTOR_ELT(dimnames, 0, dimx);
        SET_VECTOR_ELT(dimnames, 1, dimy);
        setAttrib(ans, R_DimNamesSymbol, dimnames);

        //setAttrib(ans, R_ColumnNamesSymbol, dimy);


        UNPROTECT(4);

        free(m_);
        free(a_);
        free(cc_);
        free(l_n);
        free(elements);

#if SHOW_DETAILS == 1
        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
        Rprintf("\n    peak_amount: %zu, elapsed_time: %lf, ",peak_amount, cpu_time_used);
#endif

        return(ans);
    }
}

    /*************************

     calculates points of interest; centroid,intensoid,valley with a given profile
     ATTENTION sort profile according to the mass before executing this function (increasing order)

     profile_mass:                      mass values array of the profile
     profile_abundance:                 abundance values array of the profile
     type:                              return type of point of interest
     centroid:                   		type = 0
     local maxima(intensoid):           type = 1
     valley:                     		type = 2

     *************************/
    SEXP iso_centroid_Call(SEXP profile_mass, SEXP profile_abundance, SEXP type) {

        PROTECT(profile_mass = AS_NUMERIC(profile_mass));
        PROTECT(profile_abundance = AS_NUMERIC(profile_abundance));
        PROTECT(type = AS_INTEGER(type));

        double *p_m = NUMERIC_POINTER(profile_mass);
        double *p_a = NUMERIC_POINTER(profile_abundance);
        int t = INTEGER_VALUE(type);
        int n = LENGTH(profile_mass);

        double cent_size[n];
        size_t cent_size_alloc = sizeof(cent_size);

        double* c_m = (double*)malloc(cent_size_alloc);
        double* c_a = (double*)malloc(cent_size_alloc);

        double centroid = 0.0;
        double sum_sticks = 0.0;

        double upper_sum_a = 0.0;
        double lower_sum_a = 0.0;
        double max_centroid = 0.0;
        double max_intensoid = 0.0;
        int centroid_count = 0;

        int j = 0;
        int  step = 1;

        for (ptrdiff_t  i = step; i < n - step; i++) {
            if (t == 1)
            {
                if (    (
                         *(p_a + i) > *(p_a + i + step) && *(p_a + i) > *(p_a + i - step)
                         )
                    &&
                    (
                     (*(p_m + i - step) < *(p_m + i)) && (*(p_m + i) < *(p_m + i + step))
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
                        if (sum_sticks > 0.0 && centroid_temp > 0.0 ) {
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
            for (ptrdiff_t l = 0; l < centroid_count; l++) {
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

        for (ptrdiff_t k = 0; k < j; k ++) {
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
                                     SEXP threshold,
                                     SEXP filter
                                     ) {
        int p_t;
        double *m;
        double *a;
        double *tr;
        int r;
        double t;
        int f = 0;

        PROTECT(profile_type = AS_INTEGER(profile_type));
        PROTECT(mass = AS_NUMERIC(mass));
        PROTECT(abundance = AS_NUMERIC(abundance));
        PROTECT(trace = AS_NUMERIC(trace));
        PROTECT(resolution = AS_INTEGER(resolution));
        PROTECT(threshold = AS_NUMERIC(threshold));
        PROTECT(filter = AS_INTEGER(filter));

        p_t = INTEGER_VALUE(profile_type);
        m = NUMERIC_POINTER(mass);
        a = NUMERIC_POINTER(abundance);
        tr = NUMERIC_POINTER(trace);
        r = INTEGER_VALUE(resolution);
        t = NUMERIC_VALUE(threshold);
        f = INTEGER_VALUE(filter);

        size_t tr_num = (size_t)LENGTH(trace);
        size_t n = (size_t)LENGTH(mass);

        if (tr_num >= INT_LIMIT || tr_num < 1) {
            Rprintf("\ninvalid amount of trace values");
			UNPROTECT(7);
            return R_NilValue;
        }
        if (n >= INT_LIMIT || n < 1) {
            Rprintf("\ninvalid peak limit");
			UNPROTECT(7);
            return R_NilValue;
        }

        double* p_m;
        double* p_a;

        p_m = (double*)malloc(sizeof(double[tr_num]));
        p_a = (double*)malloc(sizeof(double[tr_num]));


        unsigned int p_n = 0;

        int msg = calc_profile_with_trace(n, m, a, tr_num, tr, p_m, p_a, &p_n, r, p_t, t, f);
        if (msg) {
            Rprintf("\nCould not calculate profile");
			UNPROTECT(7);
            return R_NilValue;
        }

        UNPROTECT(7);

        SEXP profile = PROTECT(allocVector(VECSXP, 2));
        SEXP profile_mass_R;
        SEXP profile_a_R;
        double *p_m_R;
        double *p_a_R;

        PROTECT(profile_mass_R = NEW_NUMERIC(p_n));
        PROTECT(profile_a_R = NEW_NUMERIC(p_n));
        p_m_R = NUMERIC_POINTER(profile_mass_R);
        p_a_R = NUMERIC_POINTER(profile_a_R);

        for (ptrdiff_t k = 0; k < p_n; k ++) {
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

#endif
