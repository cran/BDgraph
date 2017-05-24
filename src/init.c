#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void bdmcmc_for_multi_dim(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void bdmcmc_map_for_multi_dim(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_bdmcmc_mpl_binary_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_bdmcmc_mpl_binary_ma_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_bdmcmc_mpl_binary_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_bdmcmc_mpl_binary_map_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_bdmcmc_mpl_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_bdmcmc_mpl_ma_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_bdmcmc_mpl_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_bdmcmc_mpl_map_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_global_mpl_hc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_rjmcmc_mpl_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_rjmcmc_mpl_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_bdmcmc_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_bdmcmc_ma_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_bdmcmc_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_bdmcmc_map_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_DMH_bdmcmc_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_DMH_bdmcmc_ma_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_DMH_bdmcmc_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_DMH_bdmcmc_map_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_DMH_rjmcmc_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_DMH_rjmcmc_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_rjmcmc_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_rjmcmc_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_bdmcmc_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_bdmcmc_ma_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_bdmcmc_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_bdmcmc_map_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_bdmcmc_mpl_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_bdmcmc_mpl_ma_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_bdmcmc_mpl_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_bdmcmc_mpl_map_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_DMH_bdmcmc_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_DMH_bdmcmc_ma_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_DMH_bdmcmc_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_DMH_bdmcmc_map_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_DMH_rjmcmc_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_DMH_rjmcmc_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_rjmcmc_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_rjmcmc_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_rjmcmc_mpl_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_rjmcmc_mpl_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void log_exp_mc(void *, void *, void *, void *, void *, void *, void *, void *);
extern void log_mpl_dis(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rcwish_c(void *, void *, void *, void *);
extern void rgcwish_c(void *, void *, void *, void *, void *);
extern void rgwish_c(void *, void *, void *, void *, void *);
extern void rwish_c(void *, void *, void *, void *);
extern void scale_free(void *, void *);

static const R_CMethodDef CEntries[] = {
    {"bdmcmc_for_multi_dim",                   (DL_FUNC) &bdmcmc_for_multi_dim,                   17},
    {"bdmcmc_map_for_multi_dim",               (DL_FUNC) &bdmcmc_map_for_multi_dim,               22},
    {"dgm_bdmcmc_mpl_binary_ma",               (DL_FUNC) &dgm_bdmcmc_mpl_binary_ma,               13},
    {"dgm_bdmcmc_mpl_binary_ma_multi_update",  (DL_FUNC) &dgm_bdmcmc_mpl_binary_ma_multi_update,  14},
    {"dgm_bdmcmc_mpl_binary_map",              (DL_FUNC) &dgm_bdmcmc_mpl_binary_map,              17},
    {"dgm_bdmcmc_mpl_binary_map_multi_update", (DL_FUNC) &dgm_bdmcmc_mpl_binary_map_multi_update, 19},
    {"dgm_bdmcmc_mpl_ma",                      (DL_FUNC) &dgm_bdmcmc_mpl_ma,                      13},
    {"dgm_bdmcmc_mpl_ma_multi_update",         (DL_FUNC) &dgm_bdmcmc_mpl_ma_multi_update,         14},
    {"dgm_bdmcmc_mpl_map",                     (DL_FUNC) &dgm_bdmcmc_mpl_map,                     17},
    {"dgm_bdmcmc_mpl_map_multi_update",        (DL_FUNC) &dgm_bdmcmc_mpl_map_multi_update,        19},
    {"dgm_global_mpl_hc",                      (DL_FUNC) &dgm_global_mpl_hc,                      10},
    {"dgm_rjmcmc_mpl_ma",                      (DL_FUNC) &dgm_rjmcmc_mpl_ma,                      13},
    {"dgm_rjmcmc_mpl_map",                     (DL_FUNC) &dgm_rjmcmc_mpl_map,                     17},
    {"gcgm_bdmcmc_ma",                         (DL_FUNC) &gcgm_bdmcmc_ma,                         18},
    {"gcgm_bdmcmc_ma_multi_update",            (DL_FUNC) &gcgm_bdmcmc_ma_multi_update,            19},
    {"gcgm_bdmcmc_map",                        (DL_FUNC) &gcgm_bdmcmc_map,                        22},
    {"gcgm_bdmcmc_map_multi_update",           (DL_FUNC) &gcgm_bdmcmc_map_multi_update,           24},
    {"gcgm_DMH_bdmcmc_ma",                     (DL_FUNC) &gcgm_DMH_bdmcmc_ma,                     19},
    {"gcgm_DMH_bdmcmc_ma_multi_update",        (DL_FUNC) &gcgm_DMH_bdmcmc_ma_multi_update,        20},
    {"gcgm_DMH_bdmcmc_map",                    (DL_FUNC) &gcgm_DMH_bdmcmc_map,                    23},
    {"gcgm_DMH_bdmcmc_map_multi_update",       (DL_FUNC) &gcgm_DMH_bdmcmc_map_multi_update,       25},
    {"gcgm_DMH_rjmcmc_ma",                     (DL_FUNC) &gcgm_DMH_rjmcmc_ma,                     19},
    {"gcgm_DMH_rjmcmc_map",                    (DL_FUNC) &gcgm_DMH_rjmcmc_map,                    23},
    {"gcgm_rjmcmc_ma",                         (DL_FUNC) &gcgm_rjmcmc_ma,                         18},
    {"gcgm_rjmcmc_map",                        (DL_FUNC) &gcgm_rjmcmc_map,                        22},
    {"ggm_bdmcmc_ma",                          (DL_FUNC) &ggm_bdmcmc_ma,                          13},
    {"ggm_bdmcmc_ma_multi_update",             (DL_FUNC) &ggm_bdmcmc_ma_multi_update,             14},
    {"ggm_bdmcmc_map",                         (DL_FUNC) &ggm_bdmcmc_map,                         17},
    {"ggm_bdmcmc_map_multi_update",            (DL_FUNC) &ggm_bdmcmc_map_multi_update,            19},
    {"ggm_bdmcmc_mpl_ma",                      (DL_FUNC) &ggm_bdmcmc_mpl_ma,                       9},
    {"ggm_bdmcmc_mpl_ma_multi_update",         (DL_FUNC) &ggm_bdmcmc_mpl_ma_multi_update,         10},
    {"ggm_bdmcmc_mpl_map",                     (DL_FUNC) &ggm_bdmcmc_mpl_map,                     13},
    {"ggm_bdmcmc_mpl_map_multi_update",        (DL_FUNC) &ggm_bdmcmc_mpl_map_multi_update,        15},
    {"ggm_DMH_bdmcmc_ma",                      (DL_FUNC) &ggm_DMH_bdmcmc_ma,                      15},
    {"ggm_DMH_bdmcmc_ma_multi_update",         (DL_FUNC) &ggm_DMH_bdmcmc_ma_multi_update,         16},
    {"ggm_DMH_bdmcmc_map",                     (DL_FUNC) &ggm_DMH_bdmcmc_map,                     19},
    {"ggm_DMH_bdmcmc_map_multi_update",        (DL_FUNC) &ggm_DMH_bdmcmc_map_multi_update,        21},
    {"ggm_DMH_rjmcmc_ma",                      (DL_FUNC) &ggm_DMH_rjmcmc_ma,                      15},
    {"ggm_DMH_rjmcmc_map",                     (DL_FUNC) &ggm_DMH_rjmcmc_map,                     19},
    {"ggm_rjmcmc_ma",                          (DL_FUNC) &ggm_rjmcmc_ma,                          13},
    {"ggm_rjmcmc_map",                         (DL_FUNC) &ggm_rjmcmc_map,                         17},
    {"ggm_rjmcmc_mpl_ma",                      (DL_FUNC) &ggm_rjmcmc_mpl_ma,                       9},
    {"ggm_rjmcmc_mpl_map",                     (DL_FUNC) &ggm_rjmcmc_mpl_map,                     13},
    {"log_exp_mc",                             (DL_FUNC) &log_exp_mc,                              8},
    {"log_mpl_dis",                            (DL_FUNC) &log_mpl_dis,                            11},
    {"rcwish_c",                               (DL_FUNC) &rcwish_c,                                4},
    {"rgcwish_c",                              (DL_FUNC) &rgcwish_c,                               5},
    {"rgwish_c",                               (DL_FUNC) &rgwish_c,                                5},
    {"rwish_c",                                (DL_FUNC) &rwish_c,                                 4},
    {"scale_free",                             (DL_FUNC) &scale_free,                              2},
    {NULL, NULL, 0}
};

void R_init_BDgraph(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
