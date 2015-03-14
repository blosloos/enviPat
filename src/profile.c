//
//  profile.c
//  CalcIsoStruct
//
//  Created by Christian Gerber on 1/4/13.
//  Copyright (c) 2013 EAWAG. All rights reserved.
//


#include "profile.h"

int calc_profile_with_trace(    size_t n,
                                double* m,
                                double* a,
                                size_t tr_num,
                                double* trace,
                                double* profile_mass,
                                double* profile_a,
                                unsigned int *profile_n,
                                int res,
                                int profile_type,
                                double thres_profile,
                                int filter
                            ){
    
    Peak *peaks = (Peak*)malloc(sizeof(Peak[n]));
    double a_max = 0.0;

    for (ptrdiff_t k = 0; k < n; k++) {
        (peaks + k)->abundance = *(a+k);
        (peaks + k)->mass = *(m + k);
        
        if (*(a + k) > a_max) {
            a_max = *(a + k);
        }
    }
    
    
    qsort(peaks, n, sizeof(Peak), peak_sort_by_mass);
    qsort(trace, tr_num, sizeof(double), trace_sort_by_mass);
    
    unsigned int index_prev_m = 0;
    unsigned int c = 0;
    
    for (unsigned int i = 0; i < tr_num; i++) {
        
        double m = 0.0;
        m = *(trace + i);
        
        double value = 0.0;
        
        for (unsigned int j = index_prev_m; j < n; j++) {
            
            double mk = (peaks + j)->mass;
            double v = 0.0;
            
            if (profile_type == 0) {
                v = (peaks + j)->abundance * exp(-1.0 * (pow(m - mk, 2.0)*pow(res, 2.0)*log(256))/(2.0*pow(mk, 2.0)));
                if (thres_profile == 0.0) {
                    
                    if ( m < mk ) {
                        double threshold = a_max *exp(-1.0 * (pow(m - mk, 2.0)*pow(res, 2.0)*log(256))/(2.0*pow(mk, 2.0)));
                        if (threshold == 0.0) {
                            break;
                        }
                    }
                    
                    if (m > mk) {
                        double m_prev = (peaks + index_prev_m)->mass;
                        double threshold = a_max *exp(-1.0 * (pow(m - m_prev, 2.0)*pow(res, 2.0)*log(256))/(2.0*pow(m_prev, 2.0)));
                        if (threshold == 0.0) {
                            index_prev_m = j;
                        }
                    }
                }else if(filter == 0){
                    if (fabs(m - mk) > thres_profile) {
                        if (m < mk) {
                            break;
                        }else{
                            index_prev_m = j;
                        }
                    }
                    
                }
            }
            
            if (profile_type == 1) {
                v = ((peaks + j)->abundance * pow(mk, 2.0))/(pow(mk, 2.0) + 4.0 * pow(res, 2.0) * pow(m-mk, 2.0));
                if (thres_profile == 0.0) {
                    if ( m < mk ) {
                        double threshold = (a_max * pow(mk, 2.0))/(pow(mk, 2.0) + 4.0 * pow(res, 2.0) * pow(m-mk, 2.0));
                        if (threshold <= thres_profile) {
                            break;
                        }
                    }
                    
                    if (m > mk) {
                        double m_prev = (peaks + index_prev_m)->mass;
                        double threshold = (a_max * pow(m_prev, 2.0))/(pow(m_prev, 2.0) + 4.0 * pow(res, 2.0) * pow(m-m_prev, 2.0));
                        if (threshold <= thres_profile) {
                            index_prev_m = j;
                        }
                    }
                }else if(filter == 0){
                    if (fabs(m - mk) > thres_profile) {
                        if (m < mk) {
                            break;
                        }else{
                            index_prev_m = j;
                        }
                    }
                }
            }
            value += v;
        }
        if (filter == 0) {
            if (c > 0) {
                if (value > 0.0 || (value == 0.0 && *(profile_a + c - 1) > 0.0)) {
                    *(profile_a + c) = value;
                    *(profile_mass + c) = m;
                    c++;
                }
            }else{
                if (value > 0.0) {
                    *(profile_a + c) = value;
                    *(profile_mass + c) = m;
                    c++;
                }
            }
        }else{
            *(profile_a + c) = value;
            *(profile_mass + c) = m;
            c++;
        }
    }
    
    *profile_n = c;
    free(peaks);
    return 0;
}

int calc_profile(double* m, double* a, double* profile_mass, double* profile_a, unsigned long *num, int res, int profile_type, double thres_profile){

    Peak *peaks = (Peak*)malloc(*num * sizeof(Peak));
    double a_max = 0.0;
    
    for (ptrdiff_t k = 0; k < *num; k++) {
        (peaks + k)->abundance = *(a+k);
        (peaks + k)->mass = *(m + k);
        
        if (*(a + k) > a_max) {
            a_max = *(a + k);
        }
    }

    qsort(peaks, *num, sizeof(Peak), peak_sort_by_mass);
    
    double start = peaks->mass -  2 * peaks->mass / res;
    double end = (peaks + *num - 1)->mass +   2 * (peaks + *num - 1)->mass / res;
    
    double step = (end - start)/ res;
    
    unsigned int index_prev_m = 0;
    size_t c = 0;
    
    for (unsigned int i = 0; i < res; i++) {
        
        double m = 0.0;
        m = start + i* step;
        
        double value = 0.0;
        
        for (unsigned int j = index_prev_m; j < *num; j++) {
            
            double mk = (peaks + j)->mass;
            double v = 0.0;
            
            if (profile_type == 0) {
                v = (peaks + j)->abundance * exp(-1.0 * (pow(m - mk, 2.0)*pow(res, 2.0)*log(256))/(2.0*pow(mk, 2.0)));
                if (thres_profile == 0.0) {

                    if ( m < mk ) {
                        double threshold = a_max *exp(-1.0 * (pow(m - mk, 2.0)*pow(res, 2.0)*log(256))/(2.0*pow(mk, 2.0)));
                        if (threshold == 0.0) {
                            break;
                        }
                    }

                    if (m > mk) {
                        double m_prev = (peaks + index_prev_m)->mass;
                        double threshold = a_max *exp(-1.0 * (pow(m - m_prev, 2.0)*pow(res, 2.0)*log(256))/(2.0*pow(m_prev, 2.0)));
                        if (threshold == 0.0) {
                            index_prev_m = j;
                        } 
                    }
                }else{
                    if (fabs(m - mk) > thres_profile) {
                        if (m < mk) {
                            break;
                        }else{
                            index_prev_m = j;
                        }
                    }
                    
                }
            }
            
            if (profile_type == 1) {
                v = ((peaks + j)->abundance * pow(mk, 2.0))/(pow(mk, 2.0) + 4 * pow(res, 2.0) * pow(m-mk, 2.0));
                if (thres_profile == 0.0) {
                    if ( m < mk ) {
                        double threshold = (a_max * pow(mk, 2.0))/(pow(mk, 2.0) + 4 * pow(res, 2.0) * pow(m-mk, 2.0));
                        if (threshold <= thres_profile) {
                            break;
                        }
                    }
                    
                    if (m > mk) {
                        double m_prev = (peaks + index_prev_m)->mass;
                        double threshold = (a_max * pow(m_prev, 2.0))/(pow(m_prev, 2.0) + 4 * pow(res, 2.0) * pow(m-m_prev, 2.0));
                        if (threshold <= thres_profile) {
                            index_prev_m = j;
                        }
                    }
                }else{
                    if (fabs(m - mk) > thres_profile) {
                        if (m < mk) {
                            break;
                        }else{
                            index_prev_m = j;
                        }
                    }
                    
                }
            }
            value += v;
        }
        if (c > 0) {
            if (value > 0.0 || (value == 0.0 && *(profile_a + c - 1) > 0.0)) {
                *(profile_a + c) = value;
                *(profile_mass + c) = m;
                c++;
            }
        }else{
            if (value > 0.0) {
                *(profile_a + c) = value;
                *(profile_mass + c) = m;
                c++;
            }
        }
    }
    *num = c;
    free(peaks);
    return 0;
}
