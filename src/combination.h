//
//  combination.h
//  CalcIsoStruct
//
//  Created by Christian Gerber on 11/29/12.
//  Copyright (c) 2012 EAWAG. All rights reserved.
//

#ifndef CalcIsoStruct_combination_h
#define CalcIsoStruct_combination_h

#include "isotope.h"
#include "element.h"
#include "preferences.h"

typedef struct Compound Compound;
typedef struct Combination Combination;
typedef struct Combination2 Combination2;
typedef struct CompoundMulti CompoundMulti;
typedef struct CombinationMulti CombinationMulti;
typedef struct CombinationMulti_C CombinationMulti_C;
typedef struct CombinationMulti_A CombinationMulti_A;

struct Compound{
    unsigned int sum[MAX_ISO_ELEM];
	double mass;
    double abundance;
    int counter;
};

struct Combination {
    Compound compounds[MAX_COMPOUNDS];
    Element element;
    double max_abundance;
    double max_mass;
    int amount;
};

struct CompoundMulti{	
	int sum[MAX_ISO_SIZE];
    unsigned short counter[MAX_ELEMENTS];
    double mass;
    double abundance;
    unsigned short indicator_iso;
};

struct CombinationMulti {
    CompoundMulti compounds[MAX_COMPOUNDS];
    double max_abundance;
    double max_mass;
    int amount;
};

struct Combination2 {
    Compound compounds[MAX_COMPOUNDS_2];
    CompoundMulti a2_list[MAX_COMPOUNDS_A2];
    Element element;
    double max_abundance;
    double max_mass;
    int amount;
    int a2_amount;
};

struct CombinationMulti_C {
    CompoundMulti compounds[MAX_ISO_SIZE];
    double max_abundance;
    double max_mass;
    int amount;
};

struct CombinationMulti_A {
    CompoundMulti compounds[MAX_COMPOUNDS_A];
    double max_abundance;
    double max_mass;
    int amount;
};



// algo 3 ////////////////////////////////////////////////////////////////////////////////
int create_combinations_algo_3( Combination* combination,
                                Element *element,
                                int n,
                                double threshold);

int calc_pattern_algo_3(Combination* combinations,
                         double threshold,
                         unsigned short element_amount,
                         double* mass,
                         double* a,
                         unsigned int* peak_amount,
                         unsigned int peak_limit);

int clean_combinations_algo_3(Combination* combinations,
                       double threshold,
                       unsigned short comb_amount);


// algo 1 ////////////////////////////////////////////////////////////////////////////////
int calc_combination_max_abundance(Combination2* combination,
                                   Element *element,
                                   double threshold,
                                   CombinationMulti_A* A,
                                   CombinationMulti_C* C);
                                   
int calc_combination_max_abundance_mono(Combination2* combination,
									Element *element,
									double threshold,
									CombinationMulti_A* A,
									CombinationMulti_C* C,
									double mono_abundance);

int create_combination_algo_1_mono(   Combination2 *combination,
                                 Element *element,
                                 double threshold,
                                 int peak_limit,
                                 CombinationMulti_A* A,
                                 CombinationMulti_C* C,
                                 double mono_abundance);

int create_combination_algo_1(   Combination2 *combination,
                                 Element *element,
                                 double clean_abundance,
                                 double threshold,
                                 int peak_limit,
                                 CombinationMulti_A* A,
                                 CombinationMulti_C* C);

int combine_combinations_algo_1(Combination2* combinations,
                                double threshold,
                                unsigned short element_amount,
                                double* m,
                                double* a,
                                int* cc,
                                unsigned int* peak_amount,
                                unsigned int peak_limit,
                                unsigned int iso_amount,
                                double max_abundance);

// algo 2 ////////////////////////////////////////////////////////////////////////////////
int calc_pattern_algo_2(double* m,
                         double* a,
                         int *cc,
                         double* max_a,
                         Element *elements,
                         int element_amount,
                         double threshold,
                         unsigned int* peak_amount,
                         int peak_limit);
                         
int calc_pattern_algo_2_mono(double* m,
                         double* a,
                         int *cc,
                         double* max_a,
                         Element *elements,
                         int element_amount,
                         double threshold,
                         unsigned int* peak_amount,
                         int peak_limit,
                         double mono_abundance);

int clean_combinations_2(Combination2* combinations,
                         double threshold,
                         unsigned short comb_amount);


void print_compoundmulti(CompoundMulti cm);

int compound_sort_by_abundance_dec(const void *a, const void *b);
int compound_sort_by_abundance_inc(const void *a, const void *b);
int compoundmulti_sort_by_abundance_dec(const void *a, const void *b);
int compoundmulti_sort_by_abundance_inc(const void *a, const void *b);
int combinations_sort_by_amount_inc(const void *a, const void *b);
int combinations_sort_by_amount_dec(const void *a, const void *b);

void calc_monoisotopic_single(Element* element, CompoundMulti *monoisotopic);
void create_isotope_list_single(Element *element, Isotope2 *isotopes, int *iso_c);

void create_isotope_list(Element *elements, int element_amount, Isotope2 *isotopes, int *iso_c);
void calc_monoisotopic(Element* elements, int element_amount, CompoundMulti *monoisotopic);

#endif
