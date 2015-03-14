//
//  preferences.h
//  CalcIsoStruct
//
//  Created by Christian Gerber on 1/30/13.
//  Copyright (c) 2013 EAWAG. All rights reserved.
//


#include <stddef.h>

#ifndef CalcIsoStruct_preferences_h
#define CalcIsoStruct_preferences_h

#define USE_REALLOC 1
#define SHOW_DETAILS 0
#define USE_IN_R 1
#define WITH_C_LIST 0

// allowed size of Element name, Isotope name
#define MAX_NAME_SIZE 10

// allowed size of all the isotopes over all occurring elements
#define MAX_ISO_SIZE 50

// allowed size of isotopes for a single element
#define MAX_ISO_ELEM 10

// allowed size of different elements within the sum formula
#define MAX_ELEMENTS 20

#define INT_LIMIT 2147483647


// algo_1 //////////////////////////////////////////////////////////////////////////////////////////////
#define MAX_COMPOUNDS_1 400000 // allowed size of different compounds within element container
#define MAX_COMPOUNDS_1_A 20000 // allowed size of different compounds within A list
#define MAX_COMPOUNDS_1_A2 200000 // allowed size of different compounds
#define ALLOC_START_1 10
#define ALLOC_FACTOR_1 2
#define ALLOC_START_1_A 2
#define ALLOC_FACTOR_1_A 2
#define ALLOC_START_1_A2 10
#define ALLOC_FACTOR_1_A2 2


// algo_2 //////////////////////////////////////////////////////////////////////////////////////////////
#define MAX_COMPOUNDS_2_A2 800000 // allowed size of different compounds within A2 list
#define MAX_COMPOUNDS_2_A 80000 // allowed size of different compounds within A list
#define ALLOC_START_2_A 10
#define ALLOC_FACTOR_2_A 2
#define ALLOC_START_2_A2 10
#define ALLOC_FACTOR_2_A2 2
  // mas iteration for algo2
#define MAX_ITERATION_2 4E8

// algo_3 //////////////////////////////////////////////////////////////////////////////////////////////
#define MAX_COMPOUNDS_3 200000 // allowed size of different compounds within single container
#define ALLOC_START_3 10
#define ALLOC_FACTOR_3 2


// algo_4 //////////////////////////////////////////////////////////////////////////////////////////////
#define ALGO4_BLOCK_SIZE 2
#define MAX_ITERATION_4 ALGO4_BLOCK_SIZE * 2E8
  // container constants
#define MAX_COMPOUNDS_4 ALGO4_BLOCK_SIZE * 400000 // allowed size of different compounds within multi-element container
#define ALLOC_START_4 10
#define ALLOC_FACTOR_4 2
  // A list constants
#define MAX_COMPOUNDS_4_A ALGO4_BLOCK_SIZE * 20000 // allowed size of different compounds within A list
#define ALLOC_START_4_A 10
#define ALLOC_FACTOR_4_A 2
  // A2 list constants
#define MAX_COMPOUNDS_4_A2 ALGO4_BLOCK_SIZE * 200000 // allowed size of different compounds within A2 list
#define ALLOC_START_4_A2 10
#define ALLOC_FACTOR_4_A2 2


// all algos ////////////////////////////////////////////////////////////////////////////////////////////
#define ALLOC_START_PL 100
#define ALLOC_FACTOR_PL 2

#endif


/* ERROR CODES
 * 0   : wrong theshold, failed to generate container elements
 * 22  : wrong value for relative to mono parameter
 * 90  : cannot alloc pointer for peaks mass values
 * 91  : cannot alloc pointer for peaks abundance values
 * 92  : cannot alloc pointer for peaks combination information values
 *
 * 1000: algo 1
 * * 1001 : A container is full in function calc_combination_max_abundance
 * * 1002 : A2 container is full in function calc_combination_max_abundance
 * * 1003 : element container is full in function calc_combination_max_abundance
 * 
 * * 1100 : empty container
 * * 1101 : A container is full in function create_combination_algo_1
 * * 1102 : A2 container is full in function create_combination_algo_1
 * * 1103 : element container is full in function create_combination_algo_1
 * 
 * * 1200 : amount of elements is zero
 * * 1201 : exeeded amount of peaks in combine_combinations_algo_1
 * * 1202 : exeeded amount of isotopes in combine_combinations_algo_1
 * 
 * * 109??: pointer allocation failure in calc_combination_max_abundance
 * * 119??: pointer allocation failure in create_combination_algo_1
 * 
 * 
 * 2000: algo 2
 * * 2001 : A container is full in calc_pattern_algo_2
 * * 2002 : A2 container is full in calc_pattern_algo_2
 * * 2003 : reached peak limit
 * * 2004 : isotope amount of elements is 0
 * * 2006 : reached maximum iteration
 * 
 * * 29?? : pointer allocation failure in calc_pattern_algo_2
 * 
 * 
 * 3000: algo 3
 * * 3001 : element amount is 0
 * * 3002 : exeeded peak limit
 * * 3003 : reached container limit
 * * 3202 : combination failure
 * 
 * * 39?? : pointer allocation failure
 * 
 * 4000: algo 4
 * * 4001 : A container is full
 * * 4002 : A2 container is full
 * * 4003 : reached container limit
 * * 4004 : amount of compounds within container is zero
 * * 4005 : empty isotope list for container
 * * 4006 : reached maximum interation while creating element container
 * * 4201 : reached peak limit while combining combinations
 * * 4202 : combination failure
 *
 * * 49?? : pointer allocation failure
 **/
 
