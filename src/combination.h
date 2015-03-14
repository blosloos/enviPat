//
//  combination.h
//  CalcIsoStruct
//
//  Created by Christian Gerber on 11/29/12.
//  Copyright (c) 2012 EAWAG. All rights reserved.
//

#ifndef CalcIsoStruct_combination_h
#define CalcIsoStruct_combination_h


#include <stddef.h>

#include "isotope.h"
#include "element.h"
#include "preferences.h"

typedef struct Compound Compound;
typedef struct CompoundMulti CompoundMulti;
typedef struct Combination_1 Combination_1;
typedef struct Combination_3 Combination_3;
typedef struct Combination_4 Combination_4;
typedef struct CombinationMulti CombinationMulti;
typedef struct CombinationMulti_1_C CombinationMulti_1_C;
typedef struct CombinationMulti_2_C CombinationMulti_2_C;
typedef struct CombinationMulti_1_A CombinationMulti_1_A;
typedef struct CombinationMulti_2_A CombinationMulti_2_A;


struct Compound{
    unsigned int sum[ALGO4_BLOCK_SIZE * MAX_ISO_ELEM];
    long double abundance;
    double mass;
    unsigned int counter;
    unsigned short indicator_iso;
};

struct CompoundMulti{
    unsigned int sum[MAX_ISO_SIZE];
    unsigned int counter[MAX_ELEMENTS];
    long double abundance;
    double mass;
    unsigned short indicator_iso;
};

struct Combination_1 {
#if USE_REALLOC == 1
    Compound* compounds;
    Compound* a2_list;
#else
    Compound compounds[MAX_COMPOUNDS_1];
    Compound a2_list[MAX_COMPOUNDS_1_A2];
#endif
    Element element;
    long double max_abundance;
    size_t amount;
    size_t a2_amount;
#if SHOW_DETAILS == 1
    size_t maxA;
    size_t maxA2;
    int realloc_steps;
    int A_realloc_steps;
    int A2_realloc_steps;
#endif
};

struct Combination_3 {
#if USE_REALLOC == 1
    Compound* compounds;
#else
    Compound compounds[MAX_COMPOUNDS_3];
#endif
    Element element;
    long double max_abundance;
    size_t amount;
};

struct Combination_4 {
#if USE_REALLOC == 1
    Compound* compounds;
#else
    Compound compounds[MAX_COMPOUNDS_4];
#endif
    Element element;
    long double max_abundance;
    size_t amount;
};

struct CombinationMulti {
#if USE_REALLOC == 1
    CompoundMulti* compounds;
#else
    CompoundMulti compounds[MAX_COMPOUNDS_2_A2];
#endif
    long double max_abundance;
    size_t amount;
};

struct CombinationMulti_1_C {
    Compound compounds[MAX_ISO_SIZE];
    long double max_abundance;
    size_t amount;
};

struct CombinationMulti_2_C {
    CompoundMulti compounds[MAX_ISO_SIZE];
    long double max_abundance;
    size_t amount;
};

struct CombinationMulti_1_A {
#if USE_REALLOC == 1
    Compound* compounds;
#else
    Compound compounds[MAX_COMPOUNDS_1_A];
#endif
    long double max_abundance;
    size_t amount;
};

struct CombinationMulti_2_A {
#if USE_REALLOC == 1
    CompoundMulti* compounds;
#else
    CompoundMulti compounds[MAX_COMPOUNDS_2_A];
#endif
    long double max_abundance;
    size_t amount;
};

// algo 3 ////////////////////////////////////////////////////////////////////////////////
int create_combination_algo_3( Combination_3* combination,
                               Element *element,
                               double threshold);

int clean_combination_algo_3(Combination_3* combinations,
                              long double threshold,
                              size_t  comb_amount);

int clean_combination_algo_4( Combination_4* combinations, long double threshold, size_t element_amount);

int calc_pattern_algo_3(Element *elements,
                size_t  *peak_amount,
                double t,
                unsigned short iso_amount,
                unsigned short element_amount,
                long double a_monoisotopic,
                int p_l,
                char* l_n,
                int rtm,double** m_, double** a_, int** cc_);


int combine_combinations_algo_3(Combination_3* combinations,
                                double threshold,
                                unsigned short  element_amount,
                                size_t * peak_amount,
                                int  peak_limit,
                                unsigned short iso_amount,
                                long double max_abundance,
                                int rtm, double** m_, double** a_, int** cc_);


// algo 1 ////////////////////////////////////////////////////////////////////////////////
int calc_combination_max_abundance(Combination_1* combination,
                                   Element *element,
                                   Isotope2* isotopes,
                                   double threshold,
                                   CombinationMulti_1_A* A,
                                   long double mono_abundance,
                                   int rtm);

int create_combination_algo_1(   Combination_1 *combination,
                              Element *element,
                              Isotope2* isotopes,
                              long double clean_abundance,
                              long double threshold,
                              int peak_limit,
                              CombinationMulti_1_A* A);


int combine_combinations_algo_1(Combination_1* combinations,
                                double threshold,
                                unsigned short  element_amount,
                                size_t * peak_amount,
                                int  peak_limit,
                                unsigned short iso_amount,
                                long double max_abundance,
                                int rtm,double** m_, double** a_, int** cc_);

int calc_pattern_algo_1(Element *elements,
                size_t  *p_a,
                double threshold,
                unsigned short i_a,
                unsigned short e_a,
                long double mono_abundance,
                int  peak_limit,
                int rtm, double** m_, double** a_, int** cc_);


// algo 2 ////////////////////////////////////////////////////////////////////////////////
int calc_pattern_algo_2(
                        long double* max_a,
                        Element *elements,
                        unsigned short element_amount,
                        double threshold,
                        size_t * peak_amount,
                        int peak_limit,
                        int *iso_amount_stats,
                        int rtm, double **m_, double **a_, int **cc_);

// algo 4 ////////////////////////////////////////////////////////////////////////////////
int calc_pattern_algo_4(Element *elements,
                        size_t  *peak_amount,
                        double t,
                        unsigned short iso_amount,
                        unsigned short element_amount,
                        long double a_monoisotopic,
                        int p_l,
                        char* l_n,
                        int rtm,double** m_, double** a_, int** cc_, long double* max_a);

int create_combination_algo_4(
                            Combination_4* combinations,
                            Element *elements,
                            unsigned short from_element,
                            Isotope2 *isotopes,
                            CombinationMulti_2_A* A,
                            CombinationMulti* A2,
                            CombinationMulti_2_C* C,
                            double threshold,
                            int rtm,
                            long double mono_a);

int combine_combinations_algo_4(Combination_4* combinations,
                                double threshold,
                                unsigned short  element_amount,
                                size_t* peak_amount,
                                int  peak_limit,
                                unsigned short iso_amount,
                                long double max_abundance,
                                int rtm, double **m_, double **a_, int **cc_);

void print_compoundmulti(CompoundMulti cm);

int compound_sort_by_abundance_dec(const void *a, const void *b);
int compound_sort_by_abundance_inc(const void *a, const void *b);
int compoundmulti_sort_by_abundance_dec(const void *a, const void *b);
int compoundmulti_sort_by_abundance_inc(const void *a, const void *b);
int combinations_sort_by_amount_inc(const void *a, const void *b);
int combinations_sort_by_amount_dec(const void *a, const void *b);

void calc_monoisotopic_single(Element* element, Compound *monoisotopic);
void create_isotope_list_single(Element element, Isotope2 *isotopes);

void create_isotope_list(Element *elements, size_t element_amount, Isotope2 *isotopes, unsigned short *iso_c);
void calc_monoisotopic(Element* elements, size_t element_amount, CompoundMulti *monoisotopic);

#endif
