#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void bdmcmc_for_multi_dim(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void bdmcmc_map_for_multi_dim(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_bdmcmc_mpl_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_bdmcmc_mpl_ma_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_bdmcmc_mpl_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_bdmcmc_mpl_map_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_global_mpl_hc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_rjmcmc_mpl_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_rjmcmc_mpl_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_bdmcmc_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_bdmcmc_ma_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_bdmcmc_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_bdmcmc_map_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_DMH_bdmcmc_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_DMH_bdmcmc_ma_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_DMH_bdmcmc_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_DMH_bdmcmc_map_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_DMH_rjmcmc_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_DMH_rjmcmc_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_rjmcmc_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_rjmcmc_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_bdmcmc_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_bdmcmc_ma_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_bdmcmc_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_bdmcmc_map_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_bdmcmc_mpl_ma(void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_bdmcmc_mpl_ma_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_bdmcmc_mpl_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_bdmcmc_mpl_map_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_DMH_bdmcmc_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_DMH_bdmcmc_ma_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_DMH_bdmcmc_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_DMH_bdmcmc_map_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_DMH_rjmcmc_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_DMH_rjmcmc_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_rjmcmc_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_rjmcmc_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_rjmcmc_mpl_ma(void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_rjmcmc_mpl_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void log_exp_mc(void *, void *, void *, void *, void *, void *, void *, void *);
extern void log_mpl_dis(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rcwish_c(void *, void *, void *, void *);
extern void rgcwish_c(void *, void *, void *, void *, void *);
extern void rgwish_c(void *, void *, void *, void *, void *);
extern void rwish_c(void *, void *, void *, void *);
extern void scale_free(void *, void *);

static const R_CMethodDef CEntries[] = {
    {"bdmcmc_for_multi_dim",             (DL_FUNC) &bdmcmc_for_multi_dim,             17},
    {"bdmcmc_map_for_multi_dim",         (DL_FUNC) &bdmcmc_map_for_multi_dim,         22},
    {"dgm_bdmcmc_mpl_ma",                (DL_FUNC) &dgm_bdmcmc_mpl_ma,                12},
    {"dgm_bdmcmc_mpl_ma_multi_update",   (DL_FUNC) &dgm_bdmcmc_mpl_ma_multi_update,   13},
    {"dgm_bdmcmc_mpl_map",               (DL_FUNC) &dgm_bdmcmc_mpl_map,               16},
    {"dgm_bdmcmc_mpl_map_multi_update",  (DL_FUNC) &dgm_bdmcmc_mpl_map_multi_update,  18},
    {"dgm_global_mpl_hc",                (DL_FUNC) &dgm_global_mpl_hc,                10},
    {"dgm_rjmcmc_mpl_ma",                (DL_FUNC) &dgm_rjmcmc_mpl_ma,                12},
    {"dgm_rjmcmc_mpl_map",               (DL_FUNC) &dgm_rjmcmc_mpl_map,               16},
    {"gcgm_bdmcmc_ma",                   (DL_FUNC) &gcgm_bdmcmc_ma,                   17},
    {"gcgm_bdmcmc_ma_multi_update",      (DL_FUNC) &gcgm_bdmcmc_ma_multi_update,      18},
    {"gcgm_bdmcmc_map",                  (DL_FUNC) &gcgm_bdmcmc_map,                  21},
    {"gcgm_bdmcmc_map_multi_update",     (DL_FUNC) &gcgm_bdmcmc_map_multi_update,     23},
    {"gcgm_DMH_bdmcmc_ma",               (DL_FUNC) &gcgm_DMH_bdmcmc_ma,               18},
    {"gcgm_DMH_bdmcmc_ma_multi_update",  (DL_FUNC) &gcgm_DMH_bdmcmc_ma_multi_update,  19},
    {"gcgm_DMH_bdmcmc_map",              (DL_FUNC) &gcgm_DMH_bdmcmc_map,              22},
    {"gcgm_DMH_bdmcmc_map_multi_update", (DL_FUNC) &gcgm_DMH_bdmcmc_map_multi_update, 24},
    {"gcgm_DMH_rjmcmc_ma",               (DL_FUNC) &gcgm_DMH_rjmcmc_ma,               18},
    {"gcgm_DMH_rjmcmc_map",              (DL_FUNC) &gcgm_DMH_rjmcmc_map,              22},
    {"gcgm_rjmcmc_ma",                   (DL_FUNC) &gcgm_rjmcmc_ma,                   17},
    {"gcgm_rjmcmc_map",                  (DL_FUNC) &gcgm_rjmcmc_map,                  21},
    {"ggm_bdmcmc_ma",                    (DL_FUNC) &ggm_bdmcmc_ma,                    12},
    {"ggm_bdmcmc_ma_multi_update",       (DL_FUNC) &ggm_bdmcmc_ma_multi_update,       13},
    {"ggm_bdmcmc_map",                   (DL_FUNC) &ggm_bdmcmc_map,                   16},
    {"ggm_bdmcmc_map_multi_update",      (DL_FUNC) &ggm_bdmcmc_map_multi_update,      18},
    {"ggm_bdmcmc_mpl_ma",                (DL_FUNC) &ggm_bdmcmc_mpl_ma,                 8},
    {"ggm_bdmcmc_mpl_ma_multi_update",   (DL_FUNC) &ggm_bdmcmc_mpl_ma_multi_update,    9},
    {"ggm_bdmcmc_mpl_map",               (DL_FUNC) &ggm_bdmcmc_mpl_map,               12},
    {"ggm_bdmcmc_mpl_map_multi_update",  (DL_FUNC) &ggm_bdmcmc_mpl_map_multi_update,  14},
    {"ggm_DMH_bdmcmc_ma",                (DL_FUNC) &ggm_DMH_bdmcmc_ma,                14},
    {"ggm_DMH_bdmcmc_ma_multi_update",   (DL_FUNC) &ggm_DMH_bdmcmc_ma_multi_update,   15},
    {"ggm_DMH_bdmcmc_map",               (DL_FUNC) &ggm_DMH_bdmcmc_map,               18},
    {"ggm_DMH_bdmcmc_map_multi_update",  (DL_FUNC) &ggm_DMH_bdmcmc_map_multi_update,  20},
    {"ggm_DMH_rjmcmc_ma",                (DL_FUNC) &ggm_DMH_rjmcmc_ma,                14},
    {"ggm_DMH_rjmcmc_map",               (DL_FUNC) &ggm_DMH_rjmcmc_map,               18},
    {"ggm_rjmcmc_ma",                    (DL_FUNC) &ggm_rjmcmc_ma,                    12},
    {"ggm_rjmcmc_map",                   (DL_FUNC) &ggm_rjmcmc_map,                   16},
    {"ggm_rjmcmc_mpl_ma",                (DL_FUNC) &ggm_rjmcmc_mpl_ma,                 8},
    {"ggm_rjmcmc_mpl_map",               (DL_FUNC) &ggm_rjmcmc_mpl_map,               12},
    {"log_exp_mc",                       (DL_FUNC) &log_exp_mc,                        8},
    {"log_mpl_dis",                      (DL_FUNC) &log_mpl_dis,                      11},
    {"rcwish_c",                         (DL_FUNC) &rcwish_c,                          4},
    {"rgcwish_c",                        (DL_FUNC) &rgcwish_c,                         5},
    {"rgwish_c",                         (DL_FUNC) &rgwish_c,                          5},
    {"rwish_c",                          (DL_FUNC) &rwish_c,                           4},
    {"scale_free",                       (DL_FUNC) &scale_free,                        2},
    {NULL, NULL, 0}
};

void R_init_BDgraph(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
