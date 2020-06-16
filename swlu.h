/*************************************************************************
  > File Name: swlu.h
  > Author: Xu Kai
  > Created Time: Tue 03 Dec 2019 10:39:43 AM CST
 ************************************************************************/
#include <stdlib.h>                                                                                                                                                                      
static inline int env_swlu_prof_enabled(){
    //char *swlu_enable = getenv("SWLU_PROF");
    //if (swlu_enable != NULL && !strcmp(swlu_enable, "TRUE"))
    //    return 1;
    //return 0;
    return 1;
}
static int swlu_prof_enabled = 0;
void fem_swlu_prof_init(){
    if (env_swlu_prof_enabled()){
        swlu_prof_enabled = 1;
        swlu_prof_init_();
    }
}

void fem_swlu_debug_init(){
    swlu_debug_init_();
}

void fem_swlu_prof_print(){
    if (swlu_prof_enabled)
        swlu_prof_print_();
}

void fem_swlu_prof_start(){
    if (swlu_prof_enabled)
        swlu_prof_start_();
}

void fem_swlu_prof_stop(){
    if (swlu_prof_enabled)
        swlu_prof_stop_();
}

