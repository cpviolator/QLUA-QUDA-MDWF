#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "qlayout.h"                                                 /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "latdirferm.h"                                              /* DEPS */
#include "qquda.h"                                                   /* DEPS */
#include "quda.h"
#include <string.h>

#define XQUDA_REAL double
#define XQUDA_Nc 3
#define XQUDA_Ns 4
#define XQUDA_DIM 4

static const char qudalib[] = "_quda";
static const char mtnGaugeParam[] = "_quda.mtGaugeParam";
static const char mtnInvertParam[] = "_quda.mtInvertParam";
static const char mtnMultigridParam[] = "_quda.mtMultigridParam";
static const char mtnMultigridInverter[] = "_quda.mtMultigridInverter";
// XXX? static const char mtnEigParam[] = "_quda.mtEigParam";

/* conversions between QUDA enums and Qlua strings */
#include "qq-enum.c"   /* DEPS */

#define GET_INT_VALUE(name) if (strcmp(fld, #name) == 0) { \
    lua_pushnumber(L, p->name); return 1; }
#define GET_DOUBLE_VALUE(name) if (strcmp(fld, #name) == 0) { \
    lua_pushnumber(L, p->name); return 1; }
#define GET_STRING_VALUE(name) if (strcmp(fld, #name) == 0) { \
    lua_pushstring(L, p->name); return 1; }
#define GET_NAMED_VALUE(fn,t,name)  if (strcmp(fld, #name) == 0) {	\
    lua_pushstring(L, xx2qlua_##t(L, #name, fn, p->name)); return 1; }
#define GET_PTR_VALUE(name) if (strcmp(fld,#name) == 0) { \
  if (p->name == NULL) {				     \
    lua_pushnil(L);					     \
  } else {						     \
    char msg[1024];					     \
    snprintf(msg, sizeof (msg) - 1, "<<%p>>", p->name);      \
    lua_pushstring(L, msg);				     \
  }							     \
  return 1; }
#define GET_INT_VEC_VALUE(name) if (strcmp(fld,#name) == 0) {	 \
  int len = (p->n_level < 0 || p->n_level > QUDA_MAX_MG_LEVEL) ? \
    QUDA_MAX_MG_LEVEL : p->n_level;				 \
  int i;							 \
  lua_createtable(L, len, 0);					 \
  for (i = 0; i < len; i++) {					 \
    lua_pushinteger(L, p->name[i]);				 \
    lua_rawseti(L, -2, i + 1);					 \
  }								 \
  return 1; }
#define GET_DOUBLE_VEC_VALUE(name) if (strcmp(fld,#name) == 0) {	\
    int len = (p->n_level < 0 || p->n_level > QUDA_MAX_MG_LEVEL) ?	\
      QUDA_MAX_MG_LEVEL : p->n_level;					\
    int i;								\
    lua_createtable(L, len, 0);						\
    for (i = 0; i < len; i++) {						\
      lua_pushnumber(L, p->name[i]);					\
      lua_rawseti(L, -2, i + 1);					\
    }									\
  return 1; }
//DMH
#define GET_DOUBLE_VEC_VALUE_COEFF(name) if (strcmp(fld,#name) == 0) {	\
    int len = p->Ls;							\
    int i;								\
    lua_createtable(L, len, 0);						\
    for (i = 0; i < len; i++) {						\
      lua_pushnumber(L, p->name[i]);					\
      lua_rawseti(L, -2, i + 1);					\
    }									\
  return 1; }
#define GET_NAMED_VEC_VALUE(fn,t,name) if (strcmp(fld,#name) == 0) {	\
    int len = (p->n_level < 0 || p->n_level > QUDA_MAX_MG_LEVEL) ?	\
      QUDA_MAX_MG_LEVEL : p->n_level;					\
    int i;								\
    lua_createtable(L, len, 0);						\
    for (i = 0; i < len; i++) {						\
      lua_pushstring(L, xx2qlua_##t(L, #name, fn, p->name[i]));		\
      lua_rawseti(L, -2, i + 1);					\
    }									\
  return 1; }
#define GET_INT_MAT_VALUE(name) if (strcmp(fld,#name) == 0) {		\
    int ilen = (p->n_level < 0 || p->n_level > QUDA_MAX_MG_LEVEL) ?	\
      QUDA_MAX_MG_LEVEL : p->n_level;					\
    int jlen = (XQUDA_DIM < QUDA_MAX_DIM)? XQUDA_DIM : QUDA_MAX_DIM;	\
    int i, j;								\
    lua_createtable(L, ilen, 0);					\
    for (i = 0; i < ilen; i++) {					\
      lua_createtable(L, jlen, 0);					\
      for (j = 0; j < jlen; j++) {					\
	lua_pushinteger(L, p->name[i][j]);				\
	lua_rawseti(L, -2, j + 1);					\
      }									\
      lua_rawseti(L, -2, i + 1);					\
    }									\
    return 1; }


#define PUT_INT_VALUE(name) if (strcmp(fld, #name) == 0) { \
    p->name = luaL_checkint(L, 3); return 0; }
#define PUT_DOUBLE_VALUE(name) if (strcmp(fld, #name) == 0) { \
    p->name = luaL_checknumber(L, 3); return 1; }
#define PUT_STRING_VALUE(name) if (strcmp(fld, #name) == 0) { \
    strcpy(p->name, lua_tostring(L, 3)); return 0; }
#define PUT_NAMED_VALUE(fn,t,name)  if (strcmp(fld, #name) == 0) {	\
    p->name = qlua2xx_##t(L, #name, fn, luaL_checkstring(L, 3)); return 0; }
#define PUT_INT_VEC_VALUE(name) if (strcmp(fld,#name) == 0) {	 \
  if (p->n_level < 0 || p->n_level > QUDA_MAX_MG_LEVEL) {	 \
    luaL_error(L, "bad value of n_level for setting " #name);	 \
  } else {							 \
    int *x = qlua_checkintarray(L, -1, p->n_level, NULL);	 \
    int i;							 \
    for (i = 0; i < p->n_level; i++) {				 \
      p->name[i] = x[i];					 \
    }								 \
    qlua_free(L, x);						 \
    lua_pop(L, 1);						 \
  }								 \
  return 0; }
#define PUT_DOUBLE_VEC_VALUE(name) if (strcmp(fld,#name) == 0) { \
  if (p->n_level < 0 || p->n_level > QUDA_MAX_MG_LEVEL) {	 \
    luaL_error(L, "bad value of n_level for setting " #name);	 \
  } else {							 \
    double *x = qlua_checknumberarray(L, -1, p->n_level, NULL);	 \
    int i;							 \
    for (i = 0; i < p->n_level; i++) {				 \
      p->name[i] = x[i];					 \
    }								 \
    qlua_free(L, x);						 \
    lua_pop(L, 1);						 \
  }								 \
  return 0; }

//DMH
#define PUT_DOUBLE_VEC_VALUE_COEFF(name) if (strcmp(fld,#name) == 0) { \
    double *x = qlua_checknumberarray(L, -1, 16, NULL);	       \
    int i;							       \
    for (i = 0; i < 16; i++) {				       \
      p->name[i] = x[i];					 \
    }								 \
    qlua_free(L, x);						 \
    lua_pop(L, 1);						 \
    								 \
    return 0; }

#define PUT_NAMED_VEC_VALUE(fn,t,name) if (strcmp(fld,#name) == 0) {	\
    if (p->n_level < 0 || p->n_level > QUDA_MAX_MG_LEVEL) {		\
      luaL_error(L, "bad value of n_level for setting " #name);		\
    } else {								\
      int i;								\
      const char *msg = "setting multigrid." #name;			\
      qlua_checktable(L, -1, msg);					\
      if (lua_objlen(L, -1) != p->n_level) {				\
	luaL_error(L, "bad array size for multigrid." #name);		\
      }									\
      for (i = 0; i < p->n_level; i++) {				\
	lua_pushinteger(L, i + 1);					\
	lua_gettable(L, -2);						\
	p->name[i] = qlua2xx_##t(L, #name, fn,				\
				 qlua_checkstring(L, -1, msg));		\
	lua_pop(L, 1);							\
      }									\
    }									\
    return 0; }
#define PUT_INT_MAT_VALUE(name) if (strcmp(fld,#name) == 0) {		\
    if (p->n_level < 0 || p->n_level > QUDA_MAX_MG_LEVEL) {		\
      luaL_error(L, "bad value of n_level for setting " #name);		\
    } else {								\
      int i, j;								\
      const char *msg = "setting multigrid." #name;			\
      qlua_checktable(L, -1, msg);					\
      if (lua_objlen(L, -1) != p->n_level) {				\
	luaL_error(L, "bad array size for multigrid." #name);		\
      }									\
      for (i = 0; i < p->n_level; i++) {				\
        int *x;								\
	lua_pushinteger(L, i + 1);					\
	lua_gettable(L, -2);						\
	x = qlua_checkintarray(L, -1, XQUDA_DIM, NULL);			\
	for (j = 0; j < XQUDA_DIM; j++) {				\
	  p->name[i][j] = x[j];						\
	}								\
	qlua_free(L, x);						\
	lua_pop(L, 1);							\
      }									\
    }									\
    return 0; }

/**** QudaGaugeParam */
static QudaGaugeParam *
qq_checkGaugeParam(lua_State *L, int idx)
{
  void *v = luaL_checkudata(L, idx, mtnGaugeParam);
  luaL_argcheck(L, v != 0, idx, "quda.GaugeParam expected");
  return v;
}

static int
qq_gp_fmt(lua_State *L)
{
  QudaGaugeParam *p = qq_checkGaugeParam(L, 1);
  lua_pushfstring(L, "quda.GaugeParam(%p)", p);
  return 1;
}

static int
qq_gp_get(lua_State *L)
{
  QudaGaugeParam *p = qq_checkGaugeParam(L, 1);
  const char *fld = luaL_checkstring(L, 2);

  GET_INT_VALUE(ga_pad);
  GET_INT_VALUE(llfat_ga_pad);
  GET_INT_VALUE(make_resident_gauge);
  GET_INT_VALUE(make_resident_mom);
  GET_INT_VALUE(mom_ga_pad);
  GET_INT_VALUE(overlap);
  GET_INT_VALUE(overwrite_mom);
  GET_INT_VALUE(return_result_gauge);
  GET_INT_VALUE(return_result_mom);
  GET_INT_VALUE(site_ga_pad);
  GET_INT_VALUE(staggered_phase_applied);
  GET_INT_VALUE(staple_pad);
  GET_INT_VALUE(use_resident_gauge);
  GET_INT_VALUE(use_resident_mom);
  GET_INT_VALUE(gauge_offset);
  GET_INT_VALUE(mom_offset);
  GET_INT_VALUE(site_size);

  GET_DOUBLE_VALUE(anisotropy);
  GET_DOUBLE_VALUE(i_mu);
  GET_DOUBLE_VALUE(scale);
  GET_DOUBLE_VALUE(tadpole_coeff);

  GET_NAMED_VALUE("GaugeParam", QudaFieldLocation, location);
  GET_NAMED_VALUE("GaugeParam", QudaGaugeFieldOrder, gauge_order);
  GET_NAMED_VALUE("GaugeParam", QudaGaugeFixed, gauge_fix);
  GET_NAMED_VALUE("GaugeParam", QudaLinkType, type);
  GET_NAMED_VALUE("GaugeParam", QudaPrecision, cpu_prec);
  GET_NAMED_VALUE("GaugeParam", QudaPrecision, cuda_prec);
  GET_NAMED_VALUE("GaugeParam", QudaPrecision, cuda_prec_precondition);
  GET_NAMED_VALUE("GaugeParam", QudaPrecision, cuda_prec_sloppy);
  GET_NAMED_VALUE("GaugeParam", QudaReconstructType, reconstruct);
  GET_NAMED_VALUE("GaugeParam", QudaReconstructType, reconstruct_precondition);
  GET_NAMED_VALUE("GaugeParam", QudaReconstructType, reconstruct_sloppy);
  GET_NAMED_VALUE("GaugeParam", QudaStaggeredPhase, staggered_phase_type);
  GET_NAMED_VALUE("GaugeParam", QudaTboundary, t_boundary);
  
  if (strcmp(fld, "X") == 0) {
    int i;
    lua_createtable(L, XQUDA_DIM, 0);
    for (i = 0; i < XQUDA_DIM; i++) {
      lua_pushinteger(L, p->X[i]);
      lua_rawseti(L, -2, i+1); /* lua indexing */
    }
    return 1;
  }
  return qlua_lookup(L, 2, mtnGaugeParam);
}

static int
qq_gp_put(lua_State *L)
{
  QudaGaugeParam *p = qq_checkGaugeParam(L, 1);
  const char *fld = luaL_checkstring(L, 2);

  PUT_INT_VALUE(ga_pad);
  PUT_INT_VALUE(llfat_ga_pad);
  PUT_INT_VALUE(make_resident_gauge);
  PUT_INT_VALUE(make_resident_mom);
  PUT_INT_VALUE(mom_ga_pad);
  PUT_INT_VALUE(overlap);
  PUT_INT_VALUE(overwrite_mom);
  PUT_INT_VALUE(return_result_gauge);
  PUT_INT_VALUE(return_result_mom);
  PUT_INT_VALUE(site_ga_pad);
  PUT_INT_VALUE(staggered_phase_applied);
  PUT_INT_VALUE(staple_pad);
  PUT_INT_VALUE(use_resident_gauge);
  PUT_INT_VALUE(use_resident_mom);
  PUT_INT_VALUE(gauge_offset);
  PUT_INT_VALUE(mom_offset);
  PUT_INT_VALUE(site_size);

  PUT_DOUBLE_VALUE(anisotropy);
  PUT_DOUBLE_VALUE(i_mu);
  PUT_DOUBLE_VALUE(scale);
  PUT_DOUBLE_VALUE(tadpole_coeff);

  PUT_NAMED_VALUE("GaugeParam", QudaFieldLocation, location);
  PUT_NAMED_VALUE("GaugeParam", QudaGaugeFieldOrder, gauge_order);
  PUT_NAMED_VALUE("GaugeParam", QudaGaugeFixed, gauge_fix);
  PUT_NAMED_VALUE("GaugeParam", QudaLinkType, type);
  PUT_NAMED_VALUE("GaugeParam", QudaPrecision, cpu_prec);
  PUT_NAMED_VALUE("GaugeParam", QudaPrecision, cuda_prec);
  PUT_NAMED_VALUE("GaugeParam", QudaPrecision, cuda_prec_precondition);
  PUT_NAMED_VALUE("GaugeParam", QudaPrecision, cuda_prec_sloppy);
  PUT_NAMED_VALUE("GaugeParam", QudaReconstructType, reconstruct);
  PUT_NAMED_VALUE("GaugeParam", QudaReconstructType, reconstruct_precondition);
  PUT_NAMED_VALUE("GaugeParam", QudaReconstructType, reconstruct_sloppy);
  PUT_NAMED_VALUE("GaugeParam", QudaStaggeredPhase, staggered_phase_type);
  PUT_NAMED_VALUE("GaugeParam", QudaTboundary, t_boundary);

  if (strcmp(fld, "X") == 0) {
    int *xx = qlua_checkintarray(L, -1, XQUDA_DIM, NULL);
    int i;
    for (i = 0; i < XQUDA_DIM; i++)
      p->X[i] = xx[i];
    qlua_free(L, xx);
    lua_pop(L, 1);
    return 0;
  }
  luaL_error(L, "invalid or unsettable GaugeParam element");
  return 0;
}

static int
qq_gp_print(lua_State *L)
{
  QudaGaugeParam *p = qq_checkGaugeParam(L, 1);
  printQudaGaugeParam(p);
  return 0;
}

static QudaGaugeParam *
qq_gp_create(lua_State *L)
{
  QudaGaugeParam *v = lua_newuserdata(L, sizeof (QudaGaugeParam));
  luaL_getmetatable(L, mtnGaugeParam);
  lua_setmetatable(L, -2);
  return v;
}

static int
qq_gp_copy(lua_State *L)
{
  QudaGaugeParam *p = qq_checkGaugeParam(L, 1);
  QudaGaugeParam *v = qq_gp_create(L);
  *v = *p;
  return 1;
}

static struct luaL_Reg mtGaugeParam[] = {
  { "__tostring", qq_gp_fmt   },
  { "__index",    qq_gp_get   },
  { "__newindex", qq_gp_put   },
  { "copy",       qq_gp_copy  },
  { "print",      qq_gp_print },
  { NULL,         NULL        }
};

static int
qq_gauge_param(lua_State *L)
{
  QudaGaugeParam *v = qq_gp_create(L);
  *v = newQudaGaugeParam();
  return 1;
}

/**** QudaInvertParam */
//DMH: remove static
QudaInvertParam *
qq_checkInvertParam(lua_State *L, int idx)
{
  void *v = luaL_checkudata(L, idx, mtnInvertParam);
  luaL_argcheck(L, v != 0, idx, "quda.InvertParam expected");
  return v;
}

static int
qq_ip_fmt(lua_State *L)
{
  QudaInvertParam *p = qq_checkInvertParam(L, 1);
  lua_pushfstring(L, "quda.InvertParam(%p)", p);
  return 1;
}

static int
qq_ip_get(lua_State *L)
{
  QudaInvertParam *p = qq_checkInvertParam(L, 1);
  const char *fld = luaL_checkstring(L, 2);

  GET_INT_VALUE(Ls);
  GET_INT_VALUE(Nsteps);
  GET_INT_VALUE(chrono_index);
  GET_INT_VALUE(chrono_make_resident);
  GET_INT_VALUE(chrono_max_dim);
  GET_INT_VALUE(chrono_replace_last);
  GET_INT_VALUE(chrono_use_resident);
  GET_INT_VALUE(cl_pad);
  GET_INT_VALUE(compute_action);
  GET_INT_VALUE(compute_clover);
  GET_INT_VALUE(compute_clover_inverse);
  GET_INT_VALUE(compute_clover_trlog);
  GET_INT_VALUE(compute_true_res);
  GET_INT_VALUE(deflation_grid);
  GET_INT_VALUE(eigcg_max_restarts);
  GET_INT_VALUE(gcrNkrylov);
  GET_INT_VALUE(heavy_quark_check);
  GET_INT_VALUE(iter);
  GET_INT_VALUE(make_resident_solution);
  GET_INT_VALUE(max_res_increase);
  GET_INT_VALUE(max_res_increase_total);
  GET_INT_VALUE(max_restart_num);
  GET_INT_VALUE(max_search_dim);
  GET_INT_VALUE(maxiter);
  GET_INT_VALUE(maxiter_precondition);
  GET_INT_VALUE(nev);
  GET_INT_VALUE(overlap);
  GET_INT_VALUE(pipeline);
  GET_INT_VALUE(precondition_cycle);
  GET_INT_VALUE(return_clover);
  GET_INT_VALUE(return_clover_inverse);
  GET_INT_VALUE(rhs_idx);
  GET_INT_VALUE(solution_accumulator_pipeline);
  GET_INT_VALUE(sp_pad);
  GET_INT_VALUE(use_resident_solution);
  GET_INT_VALUE(use_sloppy_partial_accumulator);

  GET_DOUBLE_VALUE(clover_coeff);
  GET_DOUBLE_VALUE(clover_rho);
  GET_DOUBLE_VALUE(eigenval_tol);
  GET_DOUBLE_VALUE(epsilon);
  GET_DOUBLE_VALUE(gflops);
  GET_DOUBLE_VALUE(inc_tol);
  GET_DOUBLE_VALUE(kappa);
  GET_DOUBLE_VALUE(m5);
  GET_DOUBLE_VALUE(mass);
  GET_DOUBLE_VALUE(mu);
  GET_DOUBLE_VALUE(omega);
  GET_DOUBLE_VALUE(reliable_delta);
  GET_DOUBLE_VALUE(secs);
  GET_DOUBLE_VALUE(tol);
  GET_DOUBLE_VALUE(tol_hq);
  GET_DOUBLE_VALUE(tol_precondition);
  GET_DOUBLE_VALUE(tol_restart);
  GET_DOUBLE_VALUE(true_res);
  GET_DOUBLE_VALUE(true_res_hq);

  GET_NAMED_VALUE("InvertParam", QudaCloverFieldOrder, clover_order);
  GET_NAMED_VALUE("InvertParam", QudaDagType, dagger);
  GET_NAMED_VALUE("InvertParam", QudaDiracFieldOrder, dirac_order);
  GET_NAMED_VALUE("InvertParam", QudaDslashType, dslash_type);
  GET_NAMED_VALUE("InvertParam", QudaDslashType, dslash_type_precondition);
  GET_NAMED_VALUE("InvertParam", QudaExtLibType, extlib_type);
  GET_NAMED_VALUE("InvertParam", QudaFieldLocation, clover_location);
  GET_NAMED_VALUE("InvertParam", QudaFieldLocation, input_location);
  GET_NAMED_VALUE("InvertParam", QudaFieldLocation, output_location);
  GET_NAMED_VALUE("InvertParam", QudaGammaBasis, gamma_basis);
  GET_NAMED_VALUE("InvertParam", QudaInverterType, inv_type);
  GET_NAMED_VALUE("InvertParam", QudaInverterType, inv_type_precondition);
  GET_NAMED_VALUE("InvertParam", QudaMassNormalization, mass_normalization);
  GET_NAMED_VALUE("InvertParam", QudaMatPCType, matpc_type);
  GET_NAMED_VALUE("InvertParam", QudaPrecision, clover_cpu_prec);
  GET_NAMED_VALUE("InvertParam", QudaPrecision, clover_cuda_prec);
  GET_NAMED_VALUE("InvertParam", QudaPrecision, clover_cuda_prec_precondition);
  GET_NAMED_VALUE("InvertParam", QudaPrecision, clover_cuda_prec_sloppy);
  GET_NAMED_VALUE("InvertParam", QudaPrecision, cpu_prec);
  GET_NAMED_VALUE("InvertParam", QudaPrecision, cuda_prec);
  GET_NAMED_VALUE("InvertParam", QudaPrecision, cuda_prec_precondition);
  GET_NAMED_VALUE("InvertParam", QudaPrecision, cuda_prec_ritz);
  GET_NAMED_VALUE("InvertParam", QudaPrecision, cuda_prec_sloppy);
  GET_NAMED_VALUE("InvertParam", QudaPreserveSource, preserve_source);
  GET_NAMED_VALUE("InvertParam", QudaResidualType, residual_type);
  GET_NAMED_VALUE("InvertParam", QudaSchwarzType, schwarz_type);
  GET_NAMED_VALUE("InvertParam", QudaSolutionType, solution_type);
  GET_NAMED_VALUE("InvertParam", QudaSolveType, solve_type);
  GET_NAMED_VALUE("InvertParam", QudaSolverNormalization, solver_normalization);
  GET_NAMED_VALUE("InvertParam", QudaTune, tune);
  GET_NAMED_VALUE("InvertParam", QudaTwistFlavorType, twist_flavor);
  GET_NAMED_VALUE("InvertParam", QudaUseInitGuess, use_init_guess);
  GET_NAMED_VALUE("InvertParam", QudaVerbosity, verbosity);
  GET_NAMED_VALUE("InvertParam", QudaVerbosity, verbosity_precondition);
  
  GET_PTR_VALUE(preconditioner);
  GET_PTR_VALUE(deflation_op);

  GET_DOUBLE_VEC_VALUE_COEFF(b_5);
  GET_DOUBLE_VEC_VALUE_COEFF(c_5);

  /*
  //DMH: Strip the b and c DW values
  if (strcmp(fld, "b_5") == 0) {
    printf("Enter b_5 get\n");   
    int i;
    lua_createtable(L, 16, 0);
    for (i = 0; i < 16; i++) {
      printf("b5 get %d, %f\n", i, p->b_5[i]);
      lua_pushnumber(L, p->b_5[i]);
      lua_rawseti(L, -2, i+1); // lua indexing 
    }
    return 1;
  }
  if (strcmp(fld, "c_5") == 0) {
    int i;
    lua_createtable(L, 16, 0);
    for (i = 0; i < 16; i++) {
      printf("c5 get %d, %f\n", i, p->c_5[i]);
      lua_pushnumber(L, p->c_5[i]);
      lua_rawseti(L, -2, i+1); // lua indexing 
    }
    return 1;
  }
  */
  
  return qlua_lookup(L, 2, mtnInvertParam);
}

static int
qq_ip_put(lua_State *L)
{
  QudaInvertParam *p = qq_checkInvertParam(L, 1);
  const char *fld = luaL_checkstring(L, 2);
  
  PUT_INT_VALUE(Ls);
  PUT_INT_VALUE(Nsteps);
  PUT_INT_VALUE(chrono_index);
  PUT_INT_VALUE(chrono_make_resident);
  PUT_INT_VALUE(chrono_max_dim);
  PUT_INT_VALUE(chrono_replace_last);
  PUT_INT_VALUE(chrono_use_resident);
  PUT_INT_VALUE(cl_pad);
  PUT_INT_VALUE(compute_action);
  PUT_INT_VALUE(compute_clover);
  PUT_INT_VALUE(compute_clover_inverse);
  PUT_INT_VALUE(compute_clover_trlog);
  PUT_INT_VALUE(compute_true_res);
  PUT_INT_VALUE(deflation_grid);
  PUT_INT_VALUE(eigcg_max_restarts);
  PUT_INT_VALUE(gcrNkrylov);
  PUT_INT_VALUE(heavy_quark_check);
  PUT_INT_VALUE(iter);
  PUT_INT_VALUE(make_resident_solution);
  PUT_INT_VALUE(max_res_increase);
  PUT_INT_VALUE(max_res_increase_total);
  PUT_INT_VALUE(max_restart_num);
  PUT_INT_VALUE(max_search_dim);
  PUT_INT_VALUE(maxiter);
  PUT_INT_VALUE(maxiter_precondition);
  PUT_INT_VALUE(nev);
  PUT_INT_VALUE(overlap);
  PUT_INT_VALUE(pipeline);
  PUT_INT_VALUE(precondition_cycle);
  PUT_INT_VALUE(return_clover);
  PUT_INT_VALUE(return_clover_inverse);
  PUT_INT_VALUE(rhs_idx);
  PUT_INT_VALUE(solution_accumulator_pipeline);
  PUT_INT_VALUE(sp_pad);
  PUT_INT_VALUE(use_resident_solution);
  PUT_INT_VALUE(use_sloppy_partial_accumulator);

  PUT_DOUBLE_VALUE(clover_coeff);
  PUT_DOUBLE_VALUE(clover_rho);
  PUT_DOUBLE_VALUE(eigenval_tol);
  PUT_DOUBLE_VALUE(epsilon);
  PUT_DOUBLE_VALUE(gflops);
  PUT_DOUBLE_VALUE(inc_tol);
  PUT_DOUBLE_VALUE(kappa);
  PUT_DOUBLE_VALUE(m5);
  PUT_DOUBLE_VALUE(mass);
  PUT_DOUBLE_VALUE(mu);
  PUT_DOUBLE_VALUE(omega);
  PUT_DOUBLE_VALUE(reliable_delta);
  PUT_DOUBLE_VALUE(secs);
  PUT_DOUBLE_VALUE(tol);
  PUT_DOUBLE_VALUE(tol_hq);
  PUT_DOUBLE_VALUE(tol_precondition);
  PUT_DOUBLE_VALUE(tol_restart);
  PUT_DOUBLE_VALUE(true_res);
  PUT_DOUBLE_VALUE(true_res_hq);

  PUT_NAMED_VALUE("InvertParam", QudaCloverFieldOrder, clover_order);
  PUT_NAMED_VALUE("InvertParam", QudaDagType, dagger);
  PUT_NAMED_VALUE("InvertParam", QudaDiracFieldOrder, dirac_order);
  PUT_NAMED_VALUE("InvertParam", QudaDslashType, dslash_type);
  PUT_NAMED_VALUE("InvertParam", QudaDslashType, dslash_type_precondition);
  PUT_NAMED_VALUE("InvertParam", QudaExtLibType, extlib_type);
  PUT_NAMED_VALUE("InvertParam", QudaFieldLocation, clover_location);
  PUT_NAMED_VALUE("InvertParam", QudaFieldLocation, input_location);
  PUT_NAMED_VALUE("InvertParam", QudaFieldLocation, output_location);
  PUT_NAMED_VALUE("InvertParam", QudaGammaBasis, gamma_basis);
  PUT_NAMED_VALUE("InvertParam", QudaInverterType, inv_type);
  PUT_NAMED_VALUE("InvertParam", QudaInverterType, inv_type_precondition);
  PUT_NAMED_VALUE("InvertParam", QudaMassNormalization, mass_normalization);
  PUT_NAMED_VALUE("InvertParam", QudaMatPCType, matpc_type);
  PUT_NAMED_VALUE("InvertParam", QudaPrecision, clover_cpu_prec);
  PUT_NAMED_VALUE("InvertParam", QudaPrecision, clover_cuda_prec);
  PUT_NAMED_VALUE("InvertParam", QudaPrecision, clover_cuda_prec_precondition);
  PUT_NAMED_VALUE("InvertParam", QudaPrecision, clover_cuda_prec_sloppy);
  PUT_NAMED_VALUE("InvertParam", QudaPrecision, cpu_prec);
  PUT_NAMED_VALUE("InvertParam", QudaPrecision, cuda_prec);
  PUT_NAMED_VALUE("InvertParam", QudaPrecision, cuda_prec_precondition);
  PUT_NAMED_VALUE("InvertParam", QudaPrecision, cuda_prec_ritz);
  PUT_NAMED_VALUE("InvertParam", QudaPrecision, cuda_prec_sloppy);
  PUT_NAMED_VALUE("InvertParam", QudaPreserveSource, preserve_source);
  PUT_NAMED_VALUE("InvertParam", QudaResidualType, residual_type);
  PUT_NAMED_VALUE("InvertParam", QudaSchwarzType, schwarz_type);
  PUT_NAMED_VALUE("InvertParam", QudaSolutionType, solution_type);
  PUT_NAMED_VALUE("InvertParam", QudaSolveType, solve_type);
  PUT_NAMED_VALUE("InvertParam", QudaSolverNormalization, solver_normalization);
  PUT_NAMED_VALUE("InvertParam", QudaTune, tune);
  PUT_NAMED_VALUE("InvertParam", QudaTwistFlavorType, twist_flavor);
  PUT_NAMED_VALUE("InvertParam", QudaUseInitGuess, use_init_guess);
  PUT_NAMED_VALUE("InvertParam", QudaVerbosity, verbosity);
  PUT_NAMED_VALUE("InvertParam", QudaVerbosity, verbosity_precondition);

  PUT_DOUBLE_VEC_VALUE_COEFF(b_5);
  PUT_DOUBLE_VEC_VALUE_COEFF(c_5);

  /*
  //DMH: Strip the b and c DW values  
  if (strcmp(fld, "b_5") == 0) {
    printf("Enter b5 put. Ls = %d\n", 16);   
    double *xx = qlua_checknumberarray(L, 3, 16, NULL);
    int i;
    for (i = 0; i < 16; i++) {
      printf("b5 put %d, %f\n", i, xx[i]);
      p->b_5[i] = xx[i];
    }
    qlua_free(L, xx);
    lua_pop(L, 1);
    return 0;
  }
  if (strcmp(fld, "c_5") == 0) {
    double *xx = qlua_checknumberarray(L, 3, 16, NULL);
    int i;
    for (i = 0; i < 16; i++) {
      printf("c5 put %d, %f\n", i, xx[i]);
      p->c_5[i] = xx[i];
    }
    qlua_free(L, xx);
    lua_pop(L, 1);
    return 0;
  }  
  */
  
  luaL_error(L, "invalid or unsettable InvertParam element");
  return 0;
}

static int
qq_ip_print(lua_State *L)
{
  QudaInvertParam *p = qq_checkInvertParam(L, 1);
  printQudaInvertParam(p);
  return 0;
}

static QudaInvertParam *
qq_ip_create(lua_State *L)
{
  QudaInvertParam *v = lua_newuserdata(L, sizeof (QudaInvertParam));
  luaL_getmetatable(L, mtnInvertParam);
  lua_setmetatable(L, -2);
  return v;
}

static int
qq_ip_copy(lua_State *L)
{
  QudaInvertParam *p = qq_checkInvertParam(L, 1);
  QudaInvertParam *v = qq_ip_create(L);
  *v = *p;
  return 1;
}

static struct luaL_Reg mtInvertParam[] = {
  { "__tostring", qq_ip_fmt    },
  { "__index",    qq_ip_get    },
  { "__newindex", qq_ip_put    },
  { "copy",       qq_ip_copy   },
  { "print",      qq_ip_print  },
  { NULL,         NULL         }
};

static int
qq_invert_param(lua_State *L)
{
  QudaInvertParam *v = qq_ip_create(L);
  *v = newQudaInvertParam();
  return 1;
}

/**** QudaMultigridParam */
static QudaMultigridParam *
qq_checkMultigridParam(lua_State *L, int idx)
{
  void *v = luaL_checkudata(L, idx, mtnMultigridParam);
  luaL_argcheck(L, v != 0, idx, "quda.MultigridtParam expected");
  return v;
}

static int
qq_mp_fmt(lua_State *L)
{
  QudaMultigridParam *p = qq_checkMultigridParam(L, 1);
  lua_pushfstring(L, "quda.MultigridParam(%p)", p);
  return 1;
}

static int
qq_mp_get(lua_State *L)
{
  QudaMultigridParam *p = qq_checkMultigridParam(L, 1);
  const char *fld = luaL_checkstring(L, 2);

  GET_INT_VALUE(n_level);
  GET_DOUBLE_VALUE(gflops);
  GET_DOUBLE_VALUE(secs);
  GET_STRING_VALUE(vec_infile);
  GET_STRING_VALUE(vec_outfile);
  GET_NAMED_VALUE("MultigridParam", QudaBoolean, generate_all_levels); 
  GET_NAMED_VALUE("MultigridParam", QudaBoolean, post_orthonormalize);
  GET_NAMED_VALUE("MultigridParam", QudaBoolean, pre_orthonormalize);
  GET_NAMED_VALUE("MultigridParam", QudaBoolean, run_verify);
  GET_NAMED_VALUE("MultigridParam", QudaComputeNullVector, compute_null_vector);
  GET_NAMED_VALUE("MultigridParam", QudaSetupType, setup_type);

  GET_INT_VEC_VALUE(n_vec);
  GET_INT_VEC_VALUE(nu_post);
  GET_INT_VEC_VALUE(nu_pre);
  GET_INT_VEC_VALUE(num_setup_iter);
  GET_INT_VEC_VALUE(setup_maxiter);
  GET_INT_VEC_VALUE(setup_maxiter_refresh);
  GET_INT_VEC_VALUE(smoother_schwarz_cycle);
  GET_INT_VEC_VALUE(spin_block_size);

  GET_INT_MAT_VALUE(geo_block_size);

  GET_DOUBLE_VEC_VALUE(coarse_solver_maxiter);
  GET_DOUBLE_VEC_VALUE(coarse_solver_tol);
  GET_DOUBLE_VEC_VALUE(mu_factor);
  GET_DOUBLE_VEC_VALUE(omega);
  GET_DOUBLE_VEC_VALUE(setup_tol);
  GET_DOUBLE_VEC_VALUE(smoother_tol);
  
  GET_NAMED_VEC_VALUE("MultigridParam", QudaBoolean, global_reduction);
  GET_NAMED_VEC_VALUE("MultigridParam", QudaFieldLocation, location);
  GET_NAMED_VEC_VALUE("MultigridParam", QudaFieldLocation, setup_location);
  GET_NAMED_VEC_VALUE("MultigridParam", QudaInverterType, coarse_solver);
  GET_NAMED_VEC_VALUE("MultigridParam", QudaInverterType, setup_inv_type);
  GET_NAMED_VEC_VALUE("MultigridParam", QudaInverterType, smoother);
  GET_NAMED_VEC_VALUE("MultigridParam", QudaMultigridCycleType, cycle_type);
  GET_NAMED_VEC_VALUE("MultigridParam", QudaPrecision, precision_null);
  GET_NAMED_VEC_VALUE("MultigridParam", QudaSchwarzType, smoother_schwarz_type);
  GET_NAMED_VEC_VALUE("MultigridParam", QudaSolutionType, coarse_grid_solution_type);
  GET_NAMED_VEC_VALUE("MultigridParam", QudaSolveType, smoother_solve_type);
  GET_NAMED_VEC_VALUE("MultigridParam", QudaVerbosity, verbosity);

  GET_PTR_VALUE(invert_param);

  return qlua_lookup(L, 2, mtnMultigridParam);
}

static int
qq_mp_put(lua_State *L)
{
  QudaMultigridParam *p = qq_checkMultigridParam(L, 1);
  const char *fld = luaL_checkstring(L, 2);

  PUT_INT_VALUE(n_level);
  PUT_DOUBLE_VALUE(gflops);
  PUT_DOUBLE_VALUE(secs);
  PUT_STRING_VALUE(vec_infile);
  PUT_STRING_VALUE(vec_outfile);
  PUT_NAMED_VALUE("MultigridParam", QudaBoolean, generate_all_levels); 
  PUT_NAMED_VALUE("MultigridParam", QudaBoolean, post_orthonormalize);
  PUT_NAMED_VALUE("MultigridParam", QudaBoolean, pre_orthonormalize);
  PUT_NAMED_VALUE("MultigridParam", QudaBoolean, run_verify);
  PUT_NAMED_VALUE("MultigridParam", QudaComputeNullVector, compute_null_vector);
  PUT_NAMED_VALUE("MultigridParam", QudaSetupType, setup_type);

  PUT_INT_VEC_VALUE(n_vec);
  PUT_INT_VEC_VALUE(nu_post);
  PUT_INT_VEC_VALUE(nu_pre);
  PUT_INT_VEC_VALUE(num_setup_iter);
  PUT_INT_VEC_VALUE(setup_maxiter);
  PUT_INT_VEC_VALUE(setup_maxiter_refresh);
  PUT_INT_VEC_VALUE(smoother_schwarz_cycle);
  PUT_INT_VEC_VALUE(spin_block_size);

  PUT_INT_MAT_VALUE(geo_block_size);
  
  PUT_DOUBLE_VEC_VALUE(coarse_solver_maxiter);
  PUT_DOUBLE_VEC_VALUE(coarse_solver_tol);
  PUT_DOUBLE_VEC_VALUE(mu_factor);
  PUT_DOUBLE_VEC_VALUE(omega);
  PUT_DOUBLE_VEC_VALUE(setup_tol);
  PUT_DOUBLE_VEC_VALUE(smoother_tol);

  PUT_NAMED_VEC_VALUE("MultigridParam", QudaBoolean, global_reduction);
  PUT_NAMED_VEC_VALUE("MultigridParam", QudaFieldLocation, location);
  PUT_NAMED_VEC_VALUE("MultigridParam", QudaFieldLocation, setup_location);
  PUT_NAMED_VEC_VALUE("MultigridParam", QudaInverterType, coarse_solver);
  PUT_NAMED_VEC_VALUE("MultigridParam", QudaInverterType, setup_inv_type);
  PUT_NAMED_VEC_VALUE("MultigridParam", QudaInverterType, smoother);
  PUT_NAMED_VEC_VALUE("MultigridParam", QudaMultigridCycleType, cycle_type);
  PUT_NAMED_VEC_VALUE("MultigridParam", QudaPrecision, precision_null);
  PUT_NAMED_VEC_VALUE("MultigridParam", QudaSchwarzType, smoother_schwarz_type);
  PUT_NAMED_VEC_VALUE("MultigridParam", QudaSolutionType, coarse_grid_solution_type);
  PUT_NAMED_VEC_VALUE("MultigridParam", QudaSolveType, smoother_solve_type);
  PUT_NAMED_VEC_VALUE("MultigridParam", QudaVerbosity, verbosity);

  luaL_error(L, "invalid or unsettable MultigridParam element");
  return 0;
}

static int
qq_mp_print(lua_State *L)
{
  QudaMultigridParam *p = qq_checkMultigridParam(L, 1);
  printQudaMultigridParam(p);
  return 0;
}

static QudaMultigridParam *
qq_mp_create(lua_State *L)
{
  QudaMultigridParam *v = lua_newuserdata(L, sizeof (QudaMultigridParam));
  luaL_getmetatable(L, mtnMultigridParam);
  lua_setmetatable(L, -2);
  return v;
}

static int
qq_mp_copy(lua_State *L)
{
  QudaMultigridParam *p = qq_checkMultigridParam(L, 1);
  QudaMultigridParam *v = qq_mp_create(L);
  *v = *p;
  return 1;
}

static struct luaL_Reg mtMultigridParam[] = {
  { "__tostring", qq_mp_fmt   },
  { "__index",    qq_mp_get   },
  { "__newindex", qq_mp_put   },
  { "copy",       qq_mp_copy  },
  { "print",      qq_mp_print },
  { NULL,         NULL        }
};

static int
qq_multigrid_param(lua_State *L)
{
  QudaMultigridParam *v = qq_mp_create(L);
  *v = newQudaMultigridParam();
  return 1;
}

/**** MultigridInverter */
typedef struct {
  QudaInvertParam master_invert;
  QudaMultigridParam multigrid_param;
  QudaInvertParam mg_inv;
  void *mg_preconditioner;
} MultigridInverter;

static MultigridInverter *
qq_checkMultigridInverter(lua_State *L, int idx)
{
  void *v = luaL_checkudata(L, idx, mtnMultigridInverter);
  luaL_argcheck(L, v != 0, idx, "quda.MultigridInverter expected");
  return v;
}

static int
qq_mi_fmt(lua_State *L)
{
  MultigridInverter *p = qq_checkMultigridInverter(L, 1);
  lua_pushfstring(L, "quda.MultigridInverter(%p)", p);
  return 1;
}


static int
qq_mi_gc(lua_State *L)
{
  MultigridInverter *p = qq_checkMultigridInverter(L, 1);
  if (p->mg_preconditioner)
    destroyMultigridQuda(p->mg_preconditioner);
  p->mg_preconditioner = NULL;
  return 0;
}

static int
qq_mi_close(lua_State *L)
{
  MultigridInverter *p = qq_checkMultigridInverter(L, 1);
  if (p->mg_preconditioner)
    destroyMultigridQuda(p->mg_preconditioner);
  p->mg_preconditioner = NULL;
  return 0;
}

static struct luaL_Reg mtMultigridInverter[] = {
  { "__tostring", qq_mi_fmt   },
  { "__gc",       qq_mi_gc    },
  { "close",      qq_mi_close },
  { NULL,         NULL        }
};

static MultigridInverter *
qq_mi_create(lua_State *L)
{
  MultigridInverter *v = lua_newuserdata(L, sizeof (MultigridInverter));
  luaL_getmetatable(L, mtnMultigridInverter);
  lua_setmetatable(L, -2);
  v->mg_preconditioner = NULL;
  return v;
}

/* A very conservative constructor of multigrid inverter */
static int
qq_make_multigrid(lua_State *L)
{
  QudaInvertParam *ip = qq_checkInvertParam(L, 1);
  QudaMultigridParam *mg = qq_checkMultigridParam(L, 2);
  QudaInvertParam *mg_ip = qq_checkInvertParam(L, 3);
  MultigridInverter *r = qq_mi_create(L);

  r->master_invert = *ip;
  r->multigrid_param = *mg;
  r->mg_inv = *mg_ip;

  r->multigrid_param.invert_param = &r->mg_inv;
  r->mg_preconditioner = newMultigridQuda(&r->multigrid_param);
  r->master_invert.preconditioner = r->mg_preconditioner;
  return 1;
}


/**** end of quda structures */
static int
qq_initQudaMemory(lua_State *L)
{
  initQudaMemory();
  return 0;
}

static int
qq_endQuda(lua_State *L)
{
  endQuda();
  return 0;
}

static int
qq_freeGaugeQuda(lua_State *L)
{
  freeGaugeQuda();
  return 0;
}

static int
qq_freeCloverQuda(lua_State *L)
{
  freeCloverQuda();
  return 0;
}

static int
qq_qChargeCuda(lua_State *L)
{
  double q = qChargeCuda();
  lua_pushnumber(L, q);
  return 1;
}

static int
qq_openMagma(lua_State *L)
{
  openMagma();
  return 0;
}

static int
qq_closeMagma(lua_State *L)
{
  closeMagma();
  return 0;
}

static int
qq_initQudaDevice(lua_State *L)
{
  int dev = lua_gettop(L) > 0? luaL_checkint(L, 1): -1;
  initQudaDevice(dev);

  return 0;
}

static int
qq_initQuda(lua_State *L)
{
  int dev = lua_gettop(L) > 0? luaL_checkint(L, 1): -1;
  initQuda(dev);

  return 0;
}

static int
qq_plaqQuda(lua_State *L)
{
  double plaq[3];
  plaqQuda(plaq);
  lua_pushnumber(L, plaq[0]);
  lua_pushnumber(L, plaq[1]);
  lua_pushnumber(L, plaq[2]);

  return 3;
}

static int
qq_setVerbosityQuda(lua_State *L)
{
  QudaVerbosity verb = qlua2xx_QudaVerbosity(L, "verbosity", "setQudaVerbosity", luaL_checkstring(L, 1));
  const char *prefix = lua_gettop(L) >= 2? luaL_checkstring(L, 2): "QUDA> ";
  setVerbosityQuda(verb, prefix, stdout);

  return 0;
}

static int
qq_initCommsGridQuda(lua_State *L)
{
  mLattice *S = qlua_checkLattice(L, 1);
  initCommsGridQuda(S->rank, S->net, qlua_comm_map, S);
  return 0;
}

//DMH: remove static
int
quda_index(const int x[], const int lo[], const int hi[])
{
  int d, k, p, v;
  for (d = XQUDA_DIM, v = 1, k = 0, p = 0; d--;) {
    k *= hi[d] - lo[d];
    k += x[d] - lo[d];
    p += x[d];
    v *= hi[d] - lo[d];
  }
  k /= 2;
  if (p & 1) k += v / 2;
  return k;
}

static void
get_gauge_field(XQUDA_REAL **q, QDP_D3_ColorMatrix **U, lua_State *L, int idx, int d)
{
  QLA_D3_ColorMatrix *Ux;
  mLattice *S = NULL;
  int lo[XQUDA_DIM];
  int hi[XQUDA_DIM];
  int subvol;
  int i;

  lua_pushnumber(L, d + 1); /* lua indexing */
  lua_gettable(L, idx);
  *U = qlua_checkLatColMat3(L, -1, NULL, XQUDA_Nc)->ptr;
  S = qlua_ObjLattice(L, -1);
  qlua_assert(S->rank == XQUDA_DIM, "expected rank 4 lattice");
  qlua_sublattice(lo, hi, S->node, S);
  for (i = 0, subvol = 1; i < XQUDA_DIM; i++) {
    subvol *= hi[i] - lo[i];
  }
  
  *q = qlua_malloc(L, subvol * 2 * XQUDA_Nc * XQUDA_Nc * sizeof (XQUDA_REAL));
  Ux = QDP_D3_expose_M(*U);
  for (i = 0; i < subvol; i++) {
    int a, b, ci, x[XQUDA_DIM];
    XQUDA_REAL *ptr;
    QDP_get_coords_L(S->lat, x, S->node, i);
    ci = quda_index(x, lo, hi);
    ptr = (*q) + ci * 2 * XQUDA_Nc * XQUDA_Nc;
    for (a = 0; a < XQUDA_Nc; a++) {
      for (b = 0; b < XQUDA_Nc; b++) {
	ptr[2 * XQUDA_Nc * a + 2 * b] = QLA_real(QLA_D3_elem_M(Ux[i], a, b));
	ptr[2 * XQUDA_Nc * a + 2 * b + 1] = QLA_imag(QLA_D3_elem_M(Ux[i], a, b));
      }
    }
  }
  QDP_D3_reset_M(*U);
  lua_pop(L, 2);
}

static void
free_gauge_field(XQUDA_REAL **q, QDP_D3_ColorMatrix **U, lua_State *L, int idx, int d)
{
  qlua_free(L, *q);
}

//DMH: remove static
void
get_fermion_field(XQUDA_REAL **q, QDP_D3_DiracFermion *f, lua_State *L, mLattice *S)
{
  QLA_D3_DiracFermion *fx;
  int lo[XQUDA_DIM];
  int hi[XQUDA_DIM];
  int subvol;
  int i;

  qlua_assert(S->rank == XQUDA_DIM, "expected rank 4 lattice");
  qlua_sublattice(lo, hi, S->node, S);
  for (i = 0, subvol = 1; i < XQUDA_DIM; i++) {
    subvol *= hi[i] - lo[i];
  }
  *q = qlua_malloc(L, subvol * 2 * XQUDA_Nc * XQUDA_Ns * sizeof (XQUDA_REAL));
  fx = QDP_D3_expose_D(f);
  for (i = 0; i < subvol; i++) {
    int c, d, ci, x[XQUDA_DIM];
    XQUDA_REAL *ptr;
    QDP_get_coords_L(S->lat, x, S->node, i);
    ci = quda_index(x, lo, hi);
    ptr = (*q) + ci * 2 * XQUDA_Ns * XQUDA_Nc;
    for (c = 0; c < XQUDA_Nc; c++) {
      for (d = 0; d < XQUDA_Ns; d++) {
	ptr[2 * XQUDA_Ns * c + 2 * d] = QLA_real(QLA_D3_elem_D(fx[i], c, d));
	ptr[2 * XQUDA_Ns * c + 2 * d + 1] = QLA_imag(QLA_D3_elem_D(fx[i], c, d));
      }
    }
  }
  QDP_D3_reset_D(f);
}

static void
put_fermion_field(QDP_D3_DiracFermion *f, XQUDA_REAL *q, lua_State *L, mLattice *S)
{
  QLA_D3_DiracFermion *fx;
  int lo[XQUDA_DIM];
  int hi[XQUDA_DIM];
  int subvol;
  int i;

  qlua_assert(S->rank == XQUDA_DIM, "expected rank 4 lattice");
  qlua_sublattice(lo, hi, S->node, S);
  for (i = 0, subvol = 1; i < XQUDA_DIM; i++) {
    subvol *= hi[i] - lo[i];
  }
  fx = QDP_D3_expose_D(f);
  for (i = 0; i < subvol; i++) {
    int c, d, ci, x[XQUDA_DIM];
    XQUDA_REAL *ptr;
    QDP_get_coords_L(S->lat, x, S->node, i);
    ci = quda_index(x, lo, hi);
    ptr = q + ci * 2 * XQUDA_Ns * XQUDA_Nc;
    for (c = 0; c < XQUDA_Nc; c++) {
      for (d = 0; d < XQUDA_Ns; d++) {
	double v_re = ptr[2 * XQUDA_Ns * c + 2 * d];
	double v_im = ptr[2 * XQUDA_Ns * c + 2 * d + 1];
	QLA_c_eq_r_plus_ir(QLA_D3_elem_D(fx[i], c, d), v_re, v_im);
      }
    }
  }
  QDP_D3_reset_D(f);
}

static void
free_fermion_field(XQUDA_REAL **q, QDP_D3_DiracFermion *f, lua_State *L, mLattice *S)
{
  qlua_free(L, *q);
}

static int
qq_loadGaugeQuda(lua_State *L)
{
  int i;
  QudaGaugeParam *p = qq_checkGaugeParam(L, 2);
  QDP_D3_ColorMatrix *U[XQUDA_DIM];
  XQUDA_REAL *qu[XQUDA_DIM];

  luaL_checktype(L, 1, LUA_TTABLE);
  CALL_QDP(L);
  for (i = 0; i < XQUDA_DIM; i++)
    get_gauge_field(&qu[i], &U[i], L, 1, i);
  loadGaugeQuda(qu, p);
  for (i = 0; i < XQUDA_DIM; i++)
    free_gauge_field(&qu[i], &U[i], L, 1, i);
  
  return 0;
}

static QudaInvertParam *
extract_inverter(lua_State *L, int idx)
{
  switch (qlua_qtype(L,idx)) {
  case qQudaInvertParam:
    return qq_checkInvertParam(L, idx);
  case qQudaMultigridInverter:
    return &qq_checkMultigridInverter(L,idx)->master_invert;
  default:
    luaL_error(L, "Expecting a Quda Inverter");
  }
  return NULL;
}

static int
qq_loadCloverQuda(lua_State *L)
{
  QudaInvertParam *p = extract_inverter(L, 1);
  CALL_QDP(L);
  loadCloverQuda(NULL, NULL, p);

  return 0;
}

static int
qq_invertQuda(lua_State *L)
{
  QDP_D3_DiracFermion *rhs = qlua_checkLatDirFerm3(L, 1, NULL, 3)->ptr;
  QudaInvertParam *p = extract_inverter(L, 2);
  mLattice *S = qlua_ObjLattice(L, 1);
  QDP_D3_DiracFermion *sol = qlua_newZeroLatDirFerm3(L, lua_gettop(L), 3)->ptr;
  XQUDA_REAL *q_rhs;
  XQUDA_REAL *q_sol;
  
  CALL_QDP(L);
  get_fermion_field(&q_rhs, rhs, L, S);
  get_fermion_field(&q_sol, sol, L, S);
  
  invertQuda(q_sol, q_rhs, p);
  put_fermion_field(sol, q_sol, L, S);
  free_fermion_field(&q_rhs, rhs, L, S);
  free_fermion_field(&q_sol, sol, L, S);

  return 1;
}

static int
qq_performAPEnStep(lua_State *L)
{
  int nSteps = luaL_checkint(L, 1);
  double alpha = luaL_checknumber(L, 2);
  performAPEnStep(nSteps, alpha);

  return 0;
}

static struct luaL_Reg fquda[] = {
  /* QUDA structures */
  {"GaugeParam",              qq_gauge_param           },
  {"InvertParam",             qq_invert_param          },
  {"MultigridParam",          qq_multigrid_param       },
  /* QUDA functions */
  {"initCommsGridQuda",       qq_initCommsGridQuda     },
  {"initQuda",                qq_initQuda              },
  {"initQudaDevice",          qq_initQudaDevice        },
  {"initQudaMemory",          qq_initQudaMemory        },
  {"endQuda",                 qq_endQuda               },
  {"setVerbosityQuda",        qq_setVerbosityQuda      },
  {"multigridQuda",           qq_make_multigrid        },
  {"freeCloverQuda",          qq_freeCloverQuda        },
  {"freeGaugeQuda",           qq_freeGaugeQuda         },
  {"invertQuda",              qq_invertQuda            },
  {"loadCloverQuda",          qq_loadCloverQuda        },
  {"loadGaugeQuda",           qq_loadGaugeQuda         },
  {"performAPEnStep",         qq_performAPEnStep       },
  {"plaqQuda",                qq_plaqQuda              },
  {"qChargeCuda",             qq_qChargeCuda           },
  {"openMagma",               qq_openMagma             },
  {"closeMagma",              qq_closeMagma            },
  { NULL,                     NULL                     }
};

int
init_quda(lua_State *L)
{
  luaL_register(L, qudalib, fquda);
  qlua_metatable(L, mtnGaugeParam, mtGaugeParam, qQudaGaugeParam);
  qlua_metatable(L, mtnInvertParam, mtInvertParam, qQudaInvertParam);
  qlua_metatable(L, mtnMultigridParam, mtMultigridParam, qQudaMultigridParam);
  qlua_metatable(L, mtnMultigridInverter, mtMultigridInverter, qQudaMultigridInverter);
  return 0;
}

void
fini_quda(void)
{
  /* free resources of GPU ?? */
}
