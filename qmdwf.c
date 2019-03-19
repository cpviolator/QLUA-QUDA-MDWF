#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "qmdwf.h"                                                   /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "qvector.h"                                                 /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "qlayout.h"                                                 /* DEPS */
#include "latreal.h"                                                 /* DEPS */
#include "latcomplex.h"                                              /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "latdirferm.h"                                              /* DEPS */
#include "latdirprop.h"                                              /* DEPS */
#include "crc32.h"                                                   /* DEPS */
#include "qend.h"                                                    /* DEPS */
#include "qlanczos.h"                                                /* DEPS */
#define QOP_MDWF_DEFAULT_PRECISION QDP_Precision
#include "qop-mdwf3.h"
#include "qmp.h"

//DMH
#include "qquda.h"                                                   /* DEPS */

/* NB: Code in this file relies on \gamma_5 = diag(1,1,-1,-1) */
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>

#if USE_Nc3
static const char mdwf_name[] = "MDWF";
static const char MDWFName[]         = "lattice.MDWF";
static const char MDWFDeflatorName[]       = "lattice.MDWF.Deflator";
static const char MDWFDeflatorStateName[]  = "lattice.MDWF.DeflatorState";

//DMH: debug tool
static void stackDump (lua_State *L) {
  int i;
  int top = lua_gettop(L);
  for (i = 1; i <= top; i++) {  /* repeat for each level */
    int t = lua_type(L, i);
    switch (t) {
      
    case LUA_TSTRING:  /* strings */
      printf("`%s'", lua_tostring(L, i));
      break;
      
    case LUA_TBOOLEAN:  /* booleans */
      printf(lua_toboolean(L, i) ? "true" : "false");
      break;
      
    case LUA_TNUMBER:  /* numbers */
      printf("%g", lua_tonumber(L, i));
      break;
      
    default:  /* other values */
      printf("%s", lua_typename(L, t));
      break;
      
    }
    printf("  ");  /* put a separator */
  }
  printf("\n");  /* end the listing */
}

typedef enum {
    DW_Shamir,
    DW_Borichi,
    DW_Chiu,
    DW_Moebius,
    DW_generic
} DW_type;

typedef int MDWFInverter(lua_State                         *L,
                         struct QOP_D3_MDWF_Fermion        *solution,
                         int                               *out_iters,
                         double                            *out_epsilon,
                         const struct QOP_D3_MDWF_Fermion  *rhs,
                         int                                log_level);

typedef struct {
    MDWFInverter    *proc;
    const char      *name;
} MDWFSolver;

typedef struct {
    struct QOP_MDWF_State        *state;
    struct QOP_MDWF_Parameters   *params;
    struct QOP_D3_MDWF_Gauge     *gauge;
    int                           Ls;
    const char                   *name;
    DW_type                       type;
} mMDWF;

typedef struct {
    int nev;
    int umax;
    int vmax;
    struct QOP_F3_MDWF_Deflator *deflator;
} mDeflatorState;


static mMDWF *qlua_newMDWF(lua_State *L, int Sidx);
static mMDWF *qlua_checkMDWF(lua_State *L,
                             int idx,
                             mLattice *S,
                             int live);

/* The generic solver */
typedef struct {
    QDP_Lattice *lat;
    int c, d;
    QLA_D3_DiracPropagator *in;
    QLA_D3_DiracPropagator *out;
    double s;
    int Ls;
} DW_P_env;

typedef struct {
    QDP_Lattice *lat;
    QLA_D3_DiracFermion *f;
    double s;
    int Ls;
} DW_D_env;

typedef struct {
    QDP_Lattice *lat;
    QLA_D3_DiracFermion **f;
    double s;           
} DW_5_env;

static void
q_DW_5_reader_scaled(double *v_re, double *v_im,
                     const int p[], int c, int d, void *e)
{
    DW_5_env *env = e;
    int i = QDP_index_L(env->lat, p);
    QLA_D3_DiracFermion *f = env->f[p[QOP_MDWF_DIM]];
    QLA_D_Complex z;
    double s = env->s;

    QLA_c_eq_c(z, QLA_elem_D(f[i], c, d));
    *v_re = s * QLA_real(z);
    *v_im = s * QLA_imag(z);
}

static void
q_DW_5_writer_scaled(const int p[], int c, int d,
                     double v_re, double v_im, void *e)
{
    DW_5_env *env = e;
    int i = QDP_index_L(env->lat, p);
    QLA_D3_DiracFermion *f = env->f[p[QOP_MDWF_DIM]];
    double s = env->s;

    QLA_real(QLA_elem_D(f[i], c, d)) = s * v_re;
    QLA_imag(QLA_elem_D(f[i], c, d)) = s * v_im;
}

static void
q_DW_D_reader_scaled(double *v_re, double *v_im,
                     const int p[], int c, int d, void *e)
{
    DW_D_env *env = e;
    int i = QDP_index_L(env->lat, p);
    QLA_D3_DiracFermion *f = env->f;
    double s = env->s;
    QLA_D_Complex z;

    QLA_c_eq_c(z, QLA_elem_D(f[i], c, d));
    *v_re = s * QLA_real(z);
    *v_im = s * QLA_imag(z);
}

static void
q_DW_D_writer_scaled(const int p[], int c, int d,
                     double v_re, double v_im, void *e)
{
    DW_D_env *env = e;
    int i = QDP_index_L(env->lat, p);
    QLA_D3_DiracFermion *f = env->f;
    double s = env->s;


    QLA_real(QLA_elem_D(f[i], c, d)) = v_re * s;
    QLA_imag(QLA_elem_D(f[i], c, d)) = v_im * s;
}

static void
q_DW_P_reader_scaled(double *v_re, double *v_im,
                     const int p[], int c, int d, void *e)
{
    DW_P_env *env = e;
    int i = QDP_index_L(env->lat, p);
    QLA_D_Complex zz;

    *v_re = 0;
    *v_im = 0;
    if (p[QOP_MDWF_DIM] == 0) {
        if (d < QOP_MDWF_FERMION_DIM / 2) {
            QLA_c_eq_c(zz, QLA_elem_P(env->in[i], c, d, env->c, env->d));
        } else {
            return;
        }
    } else if (p[QOP_MDWF_DIM] == env->Ls - 1) {
        if (d >= QOP_MDWF_FERMION_DIM / 2) {
            QLA_c_eq_c(zz, QLA_elem_P(env->in[i], c, d, env->c, env->d));
        } else {
            return;
        }
    } else {
        return;
    }

    *v_re = QLA_real(zz) * env->s;
    *v_im = QLA_imag(zz) * env->s;
}

static void
q_DW_P_writer_scaled(const int p[], int c, int d,
                     double val_re, double val_im, void *e)
{
    DW_P_env *env = e;
    int i = QDP_index_L(env->lat, p);

    val_re = val_re * env->s;
    val_im = val_im * env->s;
    if (p[QOP_MDWF_DIM] == 0) {
        if (d >= QOP_MDWF_FERMION_DIM / 2) {
            QLA_real(QLA_elem_P(env->out[i], c, d, env->c, env->d)) = val_re;
            QLA_imag(QLA_elem_P(env->out[i], c, d, env->c, env->d)) = val_im;
        }
    } else if (p[QOP_MDWF_DIM] == env->Ls - 1) {
        if (d < QOP_MDWF_FERMION_DIM / 2) {
            QLA_real(QLA_elem_P(env->out[i], c, d, env->c, env->d)) = val_re;
            QLA_imag(QLA_elem_P(env->out[i], c, d, env->c, env->d)) = val_im;
        }
    }
}

typedef struct {
    double scale;
    QDP_Lattice *lat;
    QLA_D_Real *dst[QOP_MDWF_DIM];
} DW_acurrent_env;

static void
q_DW_writer_axial_current(const int pos[QOP_MDWF_DIM],
                          int dir,
                          double value,
                          void *e)
{
    DW_acurrent_env *env = e;
    int i = QDP_index_L(env->lat, pos);

    env->dst[dir][i] += value * env->scale;
}

typedef struct {
    double scale;
    QDP_Lattice *lat;
    QLA_D_Real *dst;
} DW_midpoint_env;

static void
q_DW_writer_midpoint(const int pos[QOP_MDWF_DIM],
                     double value,
                     void *e)
{
    DW_midpoint_env *env = e;
    int i = QDP_index_L(env->lat, pos);

    env->dst[i] += value * env->scale;
}


static void
DW_5_reader(double *v_re, double *v_im,
            const int p[], int c, int d, void *e)
{
    DW_5_env *env = e;
    int i = QDP_index_L(env->lat, p);
    QLA_D3_DiracFermion **f = env->f;

    QLA_r_eq_Re_c(*v_re, QLA_elem_D(f[p[QOP_MDWF_DIM]][i], c, d));
    QLA_r_eq_Im_c(*v_im, QLA_elem_D(f[p[QOP_MDWF_DIM]][i], c, d));
}

static void
DW_5_writer(const int p[], int c, int d,
            double v_re, double v_im, void *e)
{
    DW_5_env *env = e;
    int i = QDP_index_L(env->lat, p);
    QLA_D3_DiracFermion **f = env->f;

    QLA_real(QLA_elem_D(f[p[QOP_MDWF_DIM]][i], c, d)) = v_re;
    QLA_imag(QLA_elem_D(f[p[QOP_MDWF_DIM]][i], c, d)) = v_im;
}


//DMH: Hack this function to use QUDA solver
static int q_dirac_solver(lua_State *L) {

  //DMH: QUDA addition
  //If 4th arg is true, use QUDA.
  QudaInvertParam *p;
  int USE_QUDA = 0;
  if ((lua_type(L, 4) == LUA_TBOOLEAN) && (lua_toboolean(L, 3) !=0)) {    
    p = qq_checkInvertParam(L, 5);
    USE_QUDA = 1;
  }
  
  MDWFSolver *solver = lua_touserdata(L, lua_upvalueindex(1));
  mMDWF *c = qlua_checkMDWF(L, lua_upvalueindex(2), NULL, 1);
  mLattice *S = qlua_ObjLattice(L, lua_upvalueindex(2));
  int Sidx = lua_gettop(L);
  int relaxed_p;
  long long fl1;
  double t1;
  double out_eps;
  int out_iters;
  int log_level;
  int i, index;

    switch (lua_type(L, 2)) {
    case LUA_TNONE:
    case LUA_TNIL:
        relaxed_p = 0;
        break;
    case LUA_TBOOLEAN:
        relaxed_p = lua_toboolean(L, 2);
        break;
    default:
        relaxed_p = 1;
        break;
    }

    if ((lua_type(L, 3) == LUA_TBOOLEAN) && (lua_toboolean(L, 3) != 0))
        log_level = (QOP_MDWF_LOG_CG_RESIDUAL |
                     QOP_MDWF_LOG_EIG_POSTAMBLE |
                     QOP_MDWF_LOG_EIG_UPDATE1);
    else
        log_level = 0;

    switch (qlua_qtype(L, 1)) {
    case qTable: {
      mLatDirFerm3 **psi = qlua_malloc(L, c->Ls * sizeof (mLatDirFerm3 *));
      mLatDirFerm3 **eta = qlua_malloc(L, c->Ls * sizeof (mLatDirFerm3 *));
      struct QOP_D3_MDWF_Fermion *c_psi;
      struct QOP_D3_MDWF_Fermion *c_eta;
      QLA_D3_DiracFermion **e_psi=qlua_malloc(L,c->Ls*sizeof (QLA_D3_DiracFermion *));
      QLA_D3_DiracFermion **e_eta=qlua_malloc(L,c->Ls*sizeof (QLA_D3_DiracFermion *));
      DW_5_env env;
      double norm5;
      int status;
      int i;
      const char *err_str = NULL;

        CALL_QDP(L);
        norm5 = 0;
        for (i = 0; i < c->Ls; i++) {
            QLA_D_Real normi;

            lua_pushnumber(L, i + 1); /* [sic] lua indexing */
            lua_gettable(L, 1);
            psi[i] = qlua_checkLatDirFerm3(L, -1, S, 3);
            QDP_D3_r_eq_norm2_D(&normi, psi[i]->ptr, S->all);
            norm5 += normi;
            lua_pop(L, 1);
        }
        if (norm5 == 0) {
            lua_createtable(L, c->Ls, 0);
            for (i = 0; i < c->Ls; i++) {
                qlua_newLatDirFerm3(L, Sidx, 3);
                lua_rawseti(L, -2, i + 1); /* [sic] Lua indexing */
            }
            lua_pushnumber(L, 0.0);
            lua_pushnumber(L, 0);
            return 3;
        }

        for (i = 0; i < c->Ls; i++) {
	  e_psi[i] = QDP_D3_expose_D(psi[i]->ptr);
        }
        norm5 = sqrt(norm5);
        env.lat = S->lat;
        env.f = e_psi;
        env.s = 1 / norm5;
        if (QOP_D3_MDWF_import_fermion(&c_psi, c->state, q_DW_5_reader_scaled, &env)) {
            err_str = "MDWF_import_fermion() failed";
                        goto err_done;
                }
        for (i = 0; i < c->Ls; i++) {
            QDP_D3_reset_D(psi[i]->ptr);
        }
        if (QOP_D3_MDWF_allocate_fermion(&c_eta, c->state)) {
                        err_str = "MDWF_allocate_fermion() failed";
                        goto err_done;
	}

	/*
	for(int i=0; i<c->Ls; i++) {
	  for(int j=0; j<256; j++) {	
	    for (int c = 0; c < 3; c++) {
	      for (int d = 0; d < 4; d++) {
		
		printf("CPU source: %d %d %d %d (%.12e,%.12e) %s\n", i, j, d, c,
		       QLA_real(QLA_D3_elem_D(e_psi[i][j], c, d)),
		       QLA_imag(QLA_D3_elem_D(e_psi[i][j], c, d)),
		       fabs(QLA_real(QLA_D3_elem_D(e_psi[i][j], c, d))) +
		       fabs(QLA_imag(QLA_D3_elem_D(e_psi[i][j], c, d))) > 0.0 ?
		       "<--------" : "");
	      }
	    }
	  }
	}
	*/
	
	//Here we split the function... crudely, but effectively.
	if(USE_QUDA == 1) {

	  int sp_sze_real = 2*12;
	  int subvol = 1;
	  int lo[4];
	  int hi[4];
	  
	  qlua_assert(S->rank == 4, "expected rank 4 lattice");
	  qlua_sublattice(lo, hi, S->node, S);
	  for (int k = 0; k < 4; k++) {
	    subvol *= hi[k] - lo[k];
	    printf("Subvol = %d\n", subvol);
	  }
	  
	  double *q_rhs = (double*)malloc(subvol*c->Ls*sp_sze_real*sizeof(double));
	  double *q_sol = (double*)malloc(subvol*c->Ls*sp_sze_real*sizeof(double));
	  
	  int vol = subvol; //DMH Careful....
	  
	  for(int i=0; i<c->Ls; i++) {
	    for(int j=0; j<subvol; j++) {
	      
	      int ci, x[4];
	      double *ptr;
	      QDP_get_coords_L(S->lat, x, S->node, j);
	      ci = quda_index(x, lo, hi);   // Should differ for different j...
	      
              index = i * vol * sp_sze_real + ci * sp_sze_real;
	      for (int c = 0; c < 3; c++) {
		for (int d = 0; d < 4; d++) {

		  q_rhs[index] = QLA_real(QLA_D3_elem_D(e_psi[i][j], c, d));
		  q_sol[index] = 0.0;
		  index++;
		  q_rhs[index] = QLA_imag(QLA_D3_elem_D(e_psi[i][j], c, d));
		  q_sol[index] = 0.0;
		  index++;

		  /*
		  printf("QUDA source: %d %d %d %d (%.12e,%.12e) %s\n", i, j, d, c,
			 QLA_real(QLA_D3_elem_D(e_psi[i][j], c, d)),
			 QLA_imag(QLA_D3_elem_D(e_psi[i][j], c, d)),
			 fabs(QLA_real(QLA_D3_elem_D(e_psi[i][j], c, d))) +
			 fabs(QLA_imag(QLA_D3_elem_D(e_psi[i][j], c, d))) > 0 ?
			 "<--------" : "");
		  */
		}
	      }
	    }
	  }
	  
	  invertQuda(q_sol, q_rhs, p);
	  
	  //Work out how to place QUDA sol data into QDP solution
	  //structure and safely pass back via the Lua stack.
	  
	  QOP_D3_MDWF_free_fermion(&c_psi);
	  lua_createtable(L, c->Ls, 0);
	  for (i = 0; i < c->Ls; i++) {
	    eta[i] = qlua_newLatDirFerm3(L, Sidx, 3);
	    e_eta[i] = QDP_D3_expose_D(eta[i]->ptr);
	    lua_rawseti(L, -2, i + 1); // [sic] Lua indexing
	  }
	  
	  env.lat = S->lat;
	  env.f = e_eta;
	  env.s = norm5;
	  QOP_D3_MDWF_export_fermion(q_DW_5_writer_scaled, &env, c_eta);
	  for (i = 0; i < c->Ls; i++) {
	    QDP_D3_reset_D(eta[i]->ptr);
	  }
	  
	  for(int i=0; i<c->Ls; i++) {
	    for(int j=0; j<subvol; j++) {
	      
	      int ci, x[4];
	      double *ptr;
	      QDP_get_coords_L(S->lat, x, S->node, j);
	      ci = quda_index(x, lo, hi);   // Should differ for different j...
	      
              index = i * vol * sp_sze_real + ci * sp_sze_real;
	      for (int c = 0; c < 3; c++) {
		for (int d = 0; d < 4; d++) {
		  
		  double re = q_sol[index];
                  index++;
		  double im = q_sol[index];
                  index++;
		  
		  QLA_c_eq_r_plus_ir(QLA_D3_elem_D(e_eta[i][j], c, d), re, im);
		  
		  //printf("QUDA sinkA(%d): %d %d(%d) %d %d (%.12e,%.12e)\n", vol, i, j, ci, d, c, re, im);
		  
		  /*
		  //if(i<2 && j<64){
		  printf("QUDA sink: %d %d %d %d (%.12e,%.12e) %s\n", i, j, d, c,
		  QLA_real(QLA_D3_elem_D(e_eta[i][j], c, d)),
			 QLA_imag(QLA_D3_elem_D(e_eta[i][j], c, d)),
			 fabs(QLA_real(QLA_D3_elem_D(e_eta[i][j], c, d))) +
			 fabs(QLA_imag(QLA_D3_elem_D(e_eta[i][j], c, d))) > 0.00 ?
			 "<--------" : "");
			 //}
			 */
		}
	      }
	    }
	  }
	  
	  for(int dex = 0; dex<c->Ls*subvol*sp_sze_real; dex+=2) {
	    //printf("QUDA sinkB: %d (%.12e,%.12e)\n", dex, q_sol[dex], q_sol[dex+1]);
	  }
	  /*
	    lua_createtable(L, c->Ls, 0);
	    for (i = 0; i < c->Ls; i++) {
            eta[i] = qlua_newLatDirFerm3(L, Sidx, 3);
            e_eta[i] = QDP_D3_expose_D(eta[i]->ptr);
            lua_rawseti(L, -2, i + 1); // [sic] Lua indexing 
	    }
	    
	    env.lat = S->lat;
	    env.f = e_eta;
	    env.s = norm5;
	    QOP_D3_MDWF_export_fermion(q_DW_5_writer_scaled, &env, c_eta);
	    for (i = 0; i < c->Ls; i++) {
            QDP_D3_reset_D(eta[i]->ptr);
	    }
	    
	  */
	  
	  QOP_D3_MDWF_free_fermion(&c_eta);
	  lua_pushnumber(L, out_eps * norm5);
	  lua_pushnumber(L, out_iters);
	  qlua_free(L, psi);
	  qlua_free(L, eta);
	  qlua_free(L, e_psi);
	  qlua_free(L, e_eta);
	  return 3;
	}
	else {
	  status = solver->proc(L, c_eta, &out_iters, &out_eps, c_psi, log_level);
	  
	  if (status)
            err_str = QOP_MDWF_error(c->state);
	  
	  QOP_MDWF_performance(&t1, &fl1, NULL, NULL, c->state);
	  if (t1 == 0)
            t1 = -1;
	  if (QDP_this_node == qlua_master_node)
            printf("MDWF %s solver: status = %d,"
                   " eps = %.4e, iters = %d, time = %.3f sec,"
                   " perf = %.2f MFlops/sec\n",
                   solver->name, status,
                   out_eps, out_iters, t1, fl1 * 1e-6 / t1);
	  QOP_D3_MDWF_free_fermion(&c_psi);
	  
	  lua_createtable(L, c->Ls, 0);
	  for (i = 0; i < c->Ls; i++) {
            eta[i] = qlua_newLatDirFerm3(L, Sidx, 3);
            e_eta[i] = QDP_D3_expose_D(eta[i]->ptr);
            lua_rawseti(L, -2, i + 1); /* [sic] Lua indexing */
	  }
	  
	  env.lat = S->lat;
	  env.f = e_eta;
	  env.s = norm5;
	  QOP_D3_MDWF_export_fermion(q_DW_5_writer_scaled, &env, c_eta);
	  for (i = 0; i < c->Ls; i++) {
            QDP_D3_reset_D(eta[i]->ptr);
	  }
	  QOP_D3_MDWF_free_fermion(&c_eta);

	  /*
	  for(int i=0; i<c->Ls; i++) {
	    for(int j=0; j<256; j++) {	
	      for (int c = 0; c < 3; c++) {
		for (int d = 0; d < 4; d++) {
		  
		  printf("CPU sink: %d %d %d %d (%.12e,%.12e) %s\n", i, j, d, c,
			 QLA_real(QLA_D3_elem_D(e_eta[i][j], c, d)),
			 QLA_imag(QLA_D3_elem_D(e_eta[i][j], c, d)),
			 fabs(QLA_real(QLA_D3_elem_D(e_eta[i][j], c, d))) +
			 fabs(QLA_imag(QLA_D3_elem_D(e_eta[i][j], c, d))) > 0.0 ?
			 "<--------" : "");
		}
	      }
	    }
	  }
	  */	  
	  
	  lua_pushnumber(L, out_eps * norm5);
	  lua_pushnumber(L, out_iters);
	  qlua_free(L, psi);
	  qlua_free(L, eta);
	  qlua_free(L, e_psi);
	  qlua_free(L, e_eta);
	  return 3;
	}
	
      err_done:
        qlua_free(L, psi);
        qlua_free(L, eta);
        qlua_free(L, e_psi);
        qlua_free(L, e_eta);
	return luaL_error(L, err_str);
    }
    case qLatDirFerm3: {
        mLatDirFerm3 *psi = qlua_checkLatDirFerm3(L, 1, S, 3);
        mLatDirFerm3 *eta = qlua_newZeroLatDirFerm3(L, Sidx, 3);
        struct QOP_D3_MDWF_Fermion *c_psi;
        struct QOP_D3_MDWF_Fermion *c_eta;
        DW_D_env env;
        QLA_D_Real rhs_norm2 = 0;
        double rhs_n;
        int status;
        const char *err_str;

        CALL_QDP(L);
        QDP_D3_r_eq_norm2_D(&rhs_norm2, psi->ptr, S->all);
        if (rhs_norm2 == 0) {
            lua_pushnumber(L, 0.0);
            lua_pushnumber(L, 0);
            if (c->type != DW_Shamir)
                return 3;
            lua_createtable(L, QOP_MDWF_DIM, 0);
            for (i = 0; i < QOP_MDWF_DIM; i++) {
                qlua_newLatReal(L, Sidx);
                lua_rawseti(L, -2, i + 1);
            }
            qlua_newLatReal(L, Sidx);
            return 5;
        }
        rhs_n = sqrt(rhs_norm2);
        env.lat = S->lat;
        env.Ls = c->Ls;
        env.f = QDP_D3_expose_D(psi->ptr);
        env.s = 1 / rhs_n;
        if (QOP_D3_MDWF_import_4d_fermion(&c_psi, c->state, q_DW_D_reader_scaled,
                                       &env))
            return luaL_error(L, "MDWF_import_fermion() failed");
        QDP_D3_reset_D(psi->ptr);

        if (QOP_D3_MDWF_allocate_fermion(&c_eta, c->state))
            return luaL_error(L, "MDWF_allocate_fermion() failed");

        status = solver->proc(L, c_eta, &out_iters, &out_eps, c_psi, log_level);

        if (status)
            err_str = QOP_MDWF_error(c->state);

        QOP_MDWF_performance(&t1, &fl1, NULL, NULL, c->state);
        if (t1 == 0)
            t1 = -1;
        if (QDP_this_node == qlua_master_node)
            printf("MDWF %s solver: status = %d,"
                   " eps = %.4e, iters = %d, time = %.3f sec,"
                   " perf = %.2f MFlops/sec\n",
                   solver->name, status,
                   out_eps, out_iters, t1, fl1 * 1e-6 / t1);

        QOP_D3_MDWF_free_fermion(&c_psi);

        if (status && !relaxed_p)
            return luaL_error(L, QOP_MDWF_error(c->state));

        env.lat = S->lat;
        env.f = QDP_D3_expose_D(eta->ptr);
        env.s = rhs_n;
        QOP_D3_MDWF_export_4d_fermion(q_DW_D_writer_scaled, &env, c_eta);
        QDP_D3_reset_D(eta->ptr);

        /* eta is on the stack already */
        lua_pushnumber(L, out_eps);
        lua_pushnumber(L, out_iters);

        int results;
        if (c->type == DW_Shamir) {
            QDP_D_Real *ca_ps[QOP_MDWF_DIM];
            DW_acurrent_env ac_env;

            ac_env.lat = S->lat;
            ac_env.scale = rhs_n * rhs_n;
            lua_createtable(L, QOP_MDWF_DIM, 0);
            for (i = 0; i < QOP_MDWF_DIM; i++) {
                ca_ps[i] = qlua_newZeroLatReal(L, Sidx)->ptr;
                ac_env.dst[i] = QDP_D_expose_R(ca_ps[i]);
                lua_rawseti(L, -2, i + 1);
            }

            if (QOP_D3_MDWF_axial_current(q_DW_writer_axial_current, &ac_env,
                                          c_eta, c->gauge))
                return luaL_error(L, "not enough memory");

            for (i = 0; i < QOP_MDWF_DIM; i++)
                QDP_D_reset_R(ca_ps[i]);

            QDP_D_Real *mp_ps = qlua_newZeroLatReal(L, Sidx)->ptr;
            DW_midpoint_env mp_env;

            mp_env.lat = S->lat;
            mp_env.dst = QDP_D_expose_R(mp_ps);
            mp_env.scale = rhs_n * rhs_n;

            if (QOP_D3_MDWF_midpoint_pseudo(q_DW_writer_midpoint, &mp_env,
                                            c_eta))
                return luaL_error(L, "not enough memory");

            QDP_D_reset_R(mp_ps);
            results = 5;
        } else {
            results = 3;
        }
        QOP_D3_MDWF_free_fermion(&c_eta);

        return results;
    }
    case qLatDirProp3: {
        mLatDirProp3 *psi = qlua_checkLatDirProp3(L, 1, S, 3);
        mLatDirProp3 *eta = qlua_newZeroLatDirProp3(L, Sidx, 3);
        struct QOP_D3_MDWF_Fermion *c_psi;
        struct QOP_D3_MDWF_Fermion *c_eta;
        DW_P_env env;
        QLA_D_Real rhs_norm2 = 0;
        int status;
        const char *err_str;

        lua_createtable(L, QOP_MDWF_COLORS, 0);  /* eps */
        lua_createtable(L, QOP_MDWF_COLORS, 0);  /* iters */
        CALL_QDP(L);
        if (QOP_D3_MDWF_allocate_fermion(&c_eta, c->state))
            return luaL_error(L, "MDWF_allocate_fermion() failed");

        QDP_D3_r_eq_norm2_P(&rhs_norm2, psi->ptr, S->all);
        if (rhs_norm2 == 0) {
            if (c->type != DW_Shamir)
                return 3;
            lua_createtable(L, QOP_MDWF_DIM, 0);
            for (i = 0; i < QOP_MDWF_DIM; i++) {
                qlua_newLatReal(L, Sidx);
                lua_rawseti(L, -2, i + 1);
            }
            qlua_newLatReal(L, Sidx);
            return 5;
        }
        double rhs_n = sqrt(rhs_norm2);
        env.Ls = c->Ls;
        env.lat = S->lat;
        env.in = QDP_D3_expose_P(psi->ptr);
        env.out = QDP_D3_expose_P(eta->ptr);

        QDP_D_Real *ca_ps[QOP_MDWF_DIM];
        DW_acurrent_env ac_env;
        ac_env.lat = S->lat;
        ac_env.scale = rhs_n * rhs_n;
        lua_createtable(L, QOP_MDWF_DIM, 0);
        for (i = 0; i < QOP_MDWF_DIM; i++) {
            ca_ps[i] = qlua_newZeroLatReal(L, Sidx)->ptr;
            ac_env.dst[i] = QDP_D_expose_R(ca_ps[i]);
            lua_rawseti(L, -2, i + 1);
        }
        QDP_D_Real *mp_ps = qlua_newZeroLatReal(L, Sidx)->ptr;
        DW_midpoint_env mp_env;

        mp_env.lat = S->lat;
        mp_env.dst = QDP_D_expose_R(mp_ps);
        mp_env.scale = rhs_n * rhs_n;

        for (env.c = 0; env.c < QOP_MDWF_COLORS; env.c++) {
            lua_createtable(L, QOP_MDWF_FERMION_DIM, 0); /* eps.c */
            lua_createtable(L, QOP_MDWF_FERMION_DIM, 0); /* iters.c */

            for (env.d = 0; env.d < QOP_MDWF_FERMION_DIM; env.d++) {
                env.s = 1 / rhs_n;
                if (QOP_D3_MDWF_import_fermion(&c_psi, c->state,
                                               q_DW_P_reader_scaled, &env))
                    return luaL_error(L, "MDWF_import_fermion() failed");
                status = solver->proc(L, c_eta, &out_iters, &out_eps, c_psi,
                                      log_level);

                if (status)
                    err_str = QOP_MDWF_error(c->state);

                QOP_MDWF_performance(&t1, &fl1, NULL, NULL, c->state);
                if (t1 == 0)
                    t1 = -1;
                if (QDP_this_node == qlua_master_node)
                    printf("MDWF %s solver: status = %d, c = %d, d = %d,"
                           " eps = %.4e, iters = %d, time = %.3f sec,"
                           " perf = %.2f MFlops/sec\n",
                           solver->name, status,
                           env.c, env.d, out_eps, out_iters, t1,
                           fl1 * 1e-6 / t1);
                QOP_D3_MDWF_free_fermion(&c_psi);
                if (status && !relaxed_p)
                    return luaL_error(L, QOP_MDWF_error(c->state));

                if (c->type == DW_Shamir) {
                    if (QOP_D3_MDWF_axial_current(q_DW_writer_axial_current,
                                                  &ac_env, c_eta, c->gauge))
                        return luaL_error(L, "not enough memory");
                    if (QOP_D3_MDWF_midpoint_pseudo(q_DW_writer_midpoint,
                                                    &mp_env, c_eta))
                        return luaL_error(L, "not enough memory");
                }

                env.s = rhs_n;
                QOP_D3_MDWF_export_fermion(q_DW_P_writer_scaled, &env, c_eta);
                lua_pushnumber(L, out_eps);
                lua_rawseti(L, -3, env.d + 1);
                lua_pushnumber(L, out_iters);
                lua_rawseti(L, -2, env.d + 1);
            }
            lua_rawseti(L, -5, env.c + 1);
            lua_rawseti(L, -5, env.c + 1);
        }
        QDP_D3_reset_P(psi->ptr);
        QDP_D3_reset_P(eta->ptr);
        QOP_D3_MDWF_free_fermion(&c_eta);

        for (i = 0; i < QOP_MDWF_DIM; i++)
            QDP_D_reset_R(ca_ps[i]);

        QDP_D_reset_R(mp_ps);

        if (c->type != DW_Shamir) {
            lua_pop(L, 2);
            return 3;
        }

        return 5;
    }
    default:
        break;
    }
    return luaL_error(L, "bad argument to MDWF solver");
}

/***** deflator state interface */
static mDeflatorState *q_checkDeflatorState(lua_State *L,
                                            int idx,
                                            mLattice *S,
                                            int live);

static int
q_DFS_gc(lua_State *L)
{
    mDeflatorState *d = q_checkDeflatorState(L, 1, NULL, 0);

    if (d->deflator)
        QOP_F3_MDWF_free_deflator(&d->deflator);
    d->deflator = 0;

    return 0;
}

static struct luaL_Reg mtMDWFDeflatorState[] = {
    { "__gc",         q_DFS_gc },
    { NULL, NULL},
};

static mDeflatorState *
q_newDeflatorState(lua_State *L, int Sidx)
{
    mDeflatorState *d= lua_newuserdata(L, sizeof (mDeflatorState));
    d->nev = 0;
    d->vmax = 0;
    d->umax = 0;
    d->deflator = 0;
    qlua_createLatticeTable(L, Sidx, mtMDWFDeflatorState, qMDWFDeflatorState,
                            MDWFDeflatorStateName);
    lua_setmetatable(L, -2);

    return d;
}

static mDeflatorState*
q_checkDeflatorState(lua_State *L, int idx, mLattice *S, int live)
{
    mDeflatorState *d = qlua_checkLatticeType(L, idx, qMDWFDeflatorState,
                                              MDWFDeflatorStateName);

    if (S) {
        mLattice *S1 = qlua_ObjLattice(L, idx);
        if (S1->id != S->id)
            luaL_error(L, "%s on a wrong lattice", MDWFDeflatorStateName);
        lua_pop(L, 1);
    }

    if (live && (d->deflator == 0))
        luaL_error(L, "Using closed MDWF.DeflatorState");

    return d;
}

/***** delfator interface */
static mMDWF *
q_Deflator_get_MDWF(lua_State *L, int idx, mLattice *S, int live)
{
    mMDWF *c;
    qlua_checktable(L, idx, "");

    lua_rawgeti(L, idx, 1);
    c = qlua_checkMDWF(L, -1, S, live);

    return c;
}

static mDeflatorState *
q_Deflator_get_State(lua_State *L, int idx, mLattice *S, int live)
{
    mDeflatorState *d;
    qlua_checktable(L, idx, "");

    lua_rawgeti(L, idx, 2);
    d = q_checkDeflatorState(L, -1, S, live);

    return d;
}

static int
q_DF_fmt(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 0);
    char fmt[72];

    if (d->deflator)
        sprintf(fmt, "Clover.Deflator(0x%p)", d->deflator);
    else
        sprintf(fmt, "Clover.Deflator(closed)");

    lua_pushstring(L, fmt);

    return 1;
}

static int
q_DF_close(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);

    QOP_F3_MDWF_free_deflator(&d->deflator);

    return 0;
}

static int
q_DF_reset(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);

    QOP_F3_MDWF_deflator_eigcg_reset(d->deflator);

    return 0;
}

static int
q_DF_stop(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);

    QOP_F3_MDWF_deflator_eigcg_stop(d->deflator);

    return 0;
}

static int
q_DF_resume(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);

    QOP_F3_MDWF_deflator_eigcg_resume(d->deflator);

    return 0;
}

static int
q_DF_eigenvalues(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);
    int df_cur_dim = QOP_F3_MDWF_deflator_current_dim(d->deflator);
    mVecReal *v = qlua_newVecReal(L, df_cur_dim);
    double *t = qlua_malloc(L, df_cur_dim * sizeof (double));

    CALL_QDP(L);
    int status = QOP_F3_MDWF_deflator_eigen(df_cur_dim, t, d->deflator);

    if (status == 0) {
        int i;
        for (i = 0 ; i < df_cur_dim ; i++)
            v->val[i] = t[i];
    }

    qlua_free(L, t);

    if (status == 0)
        return 1;
    else
        return 0;
}

/* return: current number of vectors in deflator space */
static int
q_DF_current_dim(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);
    lua_pushnumber(L, QOP_F3_MDWF_deflator_current_dim(d->deflator));
    return 1;
}
/* return: current number of vectors in deflator space */
static int
q_DF_start_load(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);
    if (QOP_F3_MDWF_deflator_start_load(d->deflator))
        return luaL_error(L, "MDWF_deflator_start_load() failed");
    lua_pushnumber(L, QOP_F3_MDWF_deflator_current_dim(d->deflator));
    return 1;
}
/* return: current number of vectors in deflator space */
static int
q_DF_stop_load(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);
    if (QOP_F3_MDWF_deflator_stop_load(d->deflator))
        return luaL_error(L, "MDWF_deflator_stop_load() failed");
    lua_pushnumber(L, QOP_F3_MDWF_deflator_current_dim(d->deflator));
    return 1;
}




/*****************************************************************************
 * evec rawdump hack with staggered save/load 
 *****************************************************************************/
#define EVEC_STAGGER_SAVE   1   /* do delay saves */
#define EVEC_STAGGER_LOAD   1   /* don't delay loads */
#define EVEC_BUFSIZE        1048576
typedef struct {
    int         ndim;
    const int   *subvol_dim;
    const int   *subvol_lo;
    int         L5;
    int         ncol, nspin;
    QDP_Lattice *lat;
    void        *buf;
} DW_general_env;

static size_t 
q_DW_evec_eopc_rawindex(const int loc_coord[], int s5, int ic, int is, const DW_general_env *e)
{
    int nspin_h = e->nspin / 2;
    int is_lo   = is % nspin_h,
        is_hi   = is / nspin_h;
    int eopc_x0 = is_hi + ((loc_coord[0] - e->subvol_lo[0]) / 2) * 2 ;
    size_t ind = ic + e->ncol * (is_lo + nspin_h * eopc_x0);
    size_t stride = e->ncol * nspin_h * e->subvol_dim[0];
    for (int i = 1 ; i < e->ndim ; i++) {
        ind     += stride * (loc_coord[i] - e->subvol_lo[i]);
        stride  *= e->subvol_dim[i];
    }
    return ind + stride * s5;
}
static void
q_DW_evec_eopc_rawdump_F(const int p[], int c, int d,
                  double v_re, double v_im, void *e_)
{
    DW_general_env *e = (DW_general_env *)e_;
    int ind = q_DW_evec_eopc_rawindex(p, p[QOP_MDWF_DIM], c, d, e);
    ((float complex *)(e->buf))[ind] = v_re + I * v_im;
}
static void
q_DW_evec_eopc_rawload_F(double *v_re, double *v_im, const int p[], int c, int d, void *e_)
{
    DW_general_env *e = (DW_general_env *)e_;
    size_t ind = q_DW_evec_eopc_rawindex(p, p[QOP_MDWF_DIM], c, d, e);
    float complex *cbuf = (float complex *)(e->buf);
    *v_re = creal(cbuf[ind]);
    *v_im = cimag(cbuf[ind]);
}

static int
q_DW_evec_eopc_rawdump_files(
        char **evecs_dir, char **evecs_file, char **evals_file, char **cksum_file,
        const char *dirname, const char *basename, 
        int rank, int dir_id, int group_id)
{
    *evecs_dir = *evecs_file = *evals_file = *cksum_file = NULL;
    int fname_len = strlen(dirname) + strlen(basename) + 64;
    if (       NULL == (*cksum_file = malloc(fname_len))
            || NULL == (*evals_file = malloc(fname_len))
            || NULL == (*evecs_file  = malloc(fname_len))
            || NULL == (*evecs_dir   = malloc(fname_len)))
        goto clearerr_1;

    snprintf(*evecs_dir,   fname_len, "%s/%04d", dirname, dir_id);
    snprintf(*evecs_file,  fname_len, "%s/%04d/%s.%07d", dirname, dir_id, basename, rank);
    snprintf(*evals_file, fname_len, "%s/%s.evals", dirname, basename);
    snprintf(*cksum_file, fname_len, "%s/%s.CRC32", dirname, basename);
    return 0;

clearerr_1:
    if (NULL != *evecs_dir)      { free(*evecs_dir);   *evecs_dir   = NULL; }
    if (NULL != *evecs_file)     { free(*evecs_file);  *evecs_file  = NULL; }
    if (NULL != *evals_file)    { free(*evals_file); *evals_file = NULL; }
    if (NULL != *cksum_file)    { free(*cksum_file); *cksum_file = NULL; }
    return 1;
}


#define FWRITE_SHORT_DELAY  3  /*sec*/
#define FWRITE_SHORT_COUNT  5
static int 
fwrite_shortproof(const void *ptr_, size_t size, size_t nmemb, FILE *stream, const char *msg)
{
    size_t count = 0, 
           count_tot = 0,
           i_retry  = 0,
           n_retry = FWRITE_SHORT_COUNT;
    const char *ptr = ptr_;
    n_retry = 0 < n_retry ? n_retry : 1;
    while (0 < nmemb && i_retry < n_retry) {
        count = fwrite(ptr, size, nmemb, stream);
        count_tot   += count;

        if ( (count < nmemb && 0 < count) || (0 == count && EINTR == errno)) {
            i_retry++;
            fprintf(stderr, "[%05d] %s: short fwrite %d/%d (size=%d); retry %d/%d in %ds\n", 
                    QDP_this_node, msg, (int)count, (int)nmemb, (int)size, (int)i_retry, 
                    (int)n_retry, (int)FWRITE_SHORT_DELAY);
            ptr     += size * count;
            nmemb   -= count;
            sleep(FWRITE_SHORT_DELAY);
        } else 
            break;
    }
    return count_tot;
}
#define FREAD_SHORT_DELAY  3  /*sec*/
#define FREAD_SHORT_COUNT  5
static int 
fread_shortproof(void *ptr_, size_t size, size_t nmemb, FILE *stream, const char *msg)
{
    size_t count = 0, 
           count_tot = 0,
           i_retry  = 0,
           n_retry = FREAD_SHORT_COUNT;
    char *ptr = ptr_;
    n_retry = 0 < n_retry ? n_retry : 1;
    while (0 < nmemb && i_retry < n_retry) {
        count = fread(ptr, size, nmemb, stream);
        count_tot   += count;

        if ( (count < nmemb && 0 < count) || (0 == count && EINTR == errno)) {
            i_retry++;
            fprintf(stderr, "[%05d] %s: short fread %d/%d (size=%d); retry %d/%d in %ds\n", 
                    QDP_this_node, msg, (int)count, (int)nmemb, (int)size, (int)i_retry, 
                    (int)n_retry, (int)FREAD_SHORT_DELAY);
            ptr     += size * count;
            nmemb   -= count;
            sleep(FREAD_SHORT_DELAY);
        } else 
            break;
    }
    return count_tot;
}

/* (deflator, filename, n_evec, opttable) 
    if 0 < n_evec and n_evec < df->dim, save only the first n_evec
 */
static int 
q_DF_evecs_rawdump(lua_State *L)
{
    int status = 0;

    const int evec_ncol     = 3;
    const int evec_nspin    = 4;

    struct QOP_F3_MDWF_HalfFermion *c_evec = NULL;
    void *evecs_buf = NULL;
    double *evals_buf = NULL;
    char *filename = NULL;
    char 
         *evecs_dir   = NULL,
         *evecs_file  = NULL,
         *evals_file = NULL,
         *cksum_file = NULL;
    unsigned long *crc32_buf = NULL;

    int test_bigendian = 1;
    int mach_bigendian = (0 == *(char *)&test_bigendian);
    int file_bigendian = 1;     /* use as global default */
    int wordsize       = 4;     /* single precision : default for now*/
    
    
    /* parse options */
    mDeflatorState *d   = q_Deflator_get_State(L, 1, NULL, 1);  /* [-0,+1,v] */
    mMDWF *c            = q_Deflator_get_MDWF(L, 1, NULL, 1);   /* [-0,+1,v] */
    mLattice *S         = qlua_ObjLattice(L, -1);               /* [-0,+1,v] */
    filename            = qlua_strdup(L, luaL_checkstring(L, 2));/* [-0,+0,v] */
    if (NULL == filename)
        luaL_error(L, "invalid filename or strdup error");
    int n_evec          = qlua_checkint(L, 3, "number of eigenvectors to save");
    int rank_stride     = 1;
    int verbose         = 0;
    if (qlua_checkopt_paramtable(L, 4)) {
        rank_stride = qlua_tabkey_intopt(L, 4, "rank_stride",  1);
        verbose     = qlua_tabkey_intopt(L, 4, "verbose",      0);
    }                                                           /* [-0,+0,-] */
    int df_dim  = QOP_F3_MDWF_deflator_current_dim(d->deflator);
    if (n_evec <= 0 || df_dim < n_evec)
        n_evec = df_dim;
    
    /* get mesh */
    if (4 != S->rank)
        luaL_error(L, "dim!=4 is not supported");
    int n_nodes = QMP_get_number_of_nodes();
    int subvol_lo[QLUA_MAX_LATTICE_RANK],
        subvol_hi[QLUA_MAX_LATTICE_RANK],
        subvol_dim[QLUA_MAX_LATTICE_RANK];
    qlua_sublattice(subvol_lo, subvol_hi, S->node, (void*)S);
    for (int i = 0 ; i < S->rank ; i++)
        subvol_dim[i] = subvol_hi[i] - subvol_lo[i];
    if (0 != subvol_dim[0] % 2)
        luaL_error(L, "subvol_dim[0] must be even");

    
    CALL_QDP(L);
    
    /* create file names */
    int n_groups = rank_stride;
    int my_i_group  = S->node % n_groups,     /* sequence in writing */
        my_i_dir    = S->node / n_groups,     /* directory to save to */
        my_rank     = S->node;
    char *dirname = ".", 
         *basename = filename,
         *k;
    if (NULL != (k = strrchr(filename, '/'))) {
        *k = '\0';          /* separate dirname, basename */
        dirname = filename;
        basename= k + 1;
    }
    if (q_DW_evec_eopc_rawdump_files(&evecs_dir, &evecs_file, &evals_file, &cksum_file,
            dirname, basename, my_rank, my_i_dir, my_i_group)) 
        luaL_error(L, "not enough memory");
    
    if (QOP_F3_MDWF_allocate_half_fermion(&c_evec, c->state))
        luaL_error(L, "MDWF_allocate_half_fermion() failed");

    /* alloc memory */    
    size_t evec_size_real = QDP_sites_on_node_L(S->lat) / 2 * c->Ls * evec_ncol * evec_nspin * 2;
    size_t evec_size_byte = evec_size_real * wordsize;
    evecs_buf   = qlua_malloc(L, evec_size_byte);
    evals_buf   = qlua_malloc(L, df_dim * sizeof(double));
    crc32_buf   = qlua_malloc(L, n_nodes * sizeof(crc32_buf[0]));
    if (NULL == evecs_buf || NULL == evals_buf || NULL == crc32_buf)
        luaL_error(L, "not enough memory"); 
    

    DW_general_env w_env;
    w_env.lat   = S->lat;
    w_env.buf   = evecs_buf;
    w_env.nspin = evec_nspin;   /* avoid QDP defaults */
    w_env.ncol  = evec_ncol;    /* avoid QDP defaults */
    w_env.ndim  = S->rank;
    w_env.subvol_dim = subvol_dim;
    w_env.subvol_lo  = subvol_lo;
    w_env.L5    = c->Ls;
    uint32_t crc32 = 0x0;

    if(mkdir(evecs_dir, 0777))
        if (EEXIST != errno)
            luaL_error(L, "mkdir %s: %s", evecs_dir, strerror(errno));

#if EVEC_STAGGER_SAVE
    for (int i_group = 0 ; i_group < n_groups ; i_group++) { if (i_group == my_i_group) {
#endif/*EVEC_STAGGER_SAVE*/
    FILE *f_evec_out = fopen(evecs_file, "w");
    if (NULL == f_evec_out)
        luaL_error(L, "[%05d]fopen %s: %s", QDP_this_node, evecs_file, strerror(errno));
    if (setvbuf(f_evec_out, NULL, _IOFBF, EVEC_BUFSIZE))
        printf("[%4d] cannot set bufsize=%d\n", S->node, EVEC_BUFSIZE);

    for (int i_evec = 0 ; i_evec < n_evec ; i_evec++) {
        if (QOP_F3_MDWF_deflator_extract_vector(c_evec, d->deflator, i_evec))
            luaL_error(L, "QOP_F3_MDWF_deflator_extract_vector: %s", QOP_MDWF_error(c->state));
        if (QOP_F3_MDWF_export_half_fermion(q_DW_evec_eopc_rawdump_F, &w_env, c_evec)) 
            luaL_error(L, "MDWF_export_half_fermion() failed");
        /* convert machine->file and update crc */
        if (   (mach_bigendian && !file_bigendian) 
            || (!mach_bigendian && file_bigendian) )
            swap_endian(evecs_buf, wordsize, evec_size_real);
        crc32 = crc32_fast(evecs_buf, evec_size_byte, crc32);

        if (1 != fwrite_shortproof(evecs_buf, evec_size_byte, 1, f_evec_out, evecs_file))
            luaL_error(L, "[%05d]fwrite %s: %s", QDP_this_node, evecs_file, strerror(errno));
    }
    fflush(f_evec_out);
    fsync(fileno(f_evec_out));
    if (fclose(f_evec_out))
        luaL_error(L, "[%05d]fclose %s: %s", QDP_this_node, evecs_file, strerror(errno));
#if EVEC_STAGGER_SAVE
    }  QMP_barrier(); }
#endif/*EVEC_STAGGER_SAVE*/

    /* collect crc32 from all nodes and write to a file */
    for (int i = 0 ; i < n_nodes ; i++) 
        crc32_buf[i] = 0;
    crc32_buf[S->node] = crc32;
    MPI_Allreduce(MPI_IN_PLACE, crc32_buf, n_nodes, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    
    if (S->node == qlua_master_node) {
        FILE *f_cksum_out = NULL;
        if (NULL == (f_cksum_out = fopen(cksum_file, "w")))
            luaL_error(L, "%s: %s", cksum_file, strerror(errno));
        for (int i = 0 ; i < n_nodes ; i++)
            fprintf(f_cksum_out, "%lX\n", (unsigned long)(crc32_buf[i]));
        if (fclose(f_cksum_out))
            luaL_error(L, "%s: %s", cksum_file, strerror(errno));
    }
    /* get e.values and write to a file */
    if (S->node == qlua_master_node) {
        if (0 != (status = QOP_F3_MDWF_deflator_eigen(df_dim, evals_buf, d->deflator)))
            luaL_error(L, "MDWF_deflator_eigen() failed");

        FILE *f_eval_out = NULL;
        if (NULL == (f_eval_out = fopen(evals_file, "w")))
            luaL_error(L, "%s: %s", evals_file, strerror(errno));
        for (int i = 0 ; i < df_dim ; i++)
            fprintf(f_eval_out, "%.15e\n", evals_buf[i]);
        if (fclose(f_eval_out))
            luaL_error(L, "%s: %s", evals_file, strerror(errno));
    }

    lua_pushnumber(L, n_evec);  /* return the number of saved e.vectors */
    
    /* cleanup */
#define FREE_NONNULL(p) do { if (NULL != (p)) { free(p) ; (p) = NULL; } } while(0)

    FREE_NONNULL(evecs_buf);
    FREE_NONNULL(evals_buf);
    FREE_NONNULL(crc32_buf);
    FREE_NONNULL(filename);
    FREE_NONNULL(evecs_dir);
    FREE_NONNULL(evecs_file);
    FREE_NONNULL(evals_file);
    FREE_NONNULL(cksum_file);
#undef FREE_NONNULL
    if (NULL != c_evec)         QOP_F3_MDWF_free_half_fermion(&c_evec);

    return 1;
}

/* (deflator, filename, n_evec, opttable) 
 */
static int q_DF_evecs_rawload(lua_State *L)
{
    char strbuf[1024];

    const int evec_ncol     = 3;
    const int evec_nspin    = 4;

    struct QOP_F3_MDWF_HalfFermion *c_evec = NULL;
    struct QOP_F3_MDWF_Gauge *gaugeF = NULL;
    void *evecs_buf = NULL;
    char *filename = NULL;
    char 
         *evecs_dir   = NULL,
         *evecs_file  = NULL,
         *evals_file = NULL,
         *cksum_file = NULL;
    unsigned long *crc32_buf = NULL;

    int test_bigendian = 1;
    int mach_bigendian = (0 == *(char *)&test_bigendian);
    int file_bigendian = 1;     /* use as global default */
    int wordsize       = 4;     /* single precision : default for now*/
    

    /* parse options */
    mDeflatorState *d   = q_Deflator_get_State(L, 1, NULL, 1);  /* [-0,+1,v] */
    mMDWF *c            = q_Deflator_get_MDWF(L, 1, NULL, 1);   /* [-0,+1,v] */
    mLattice *S         = qlua_ObjLattice(L, -1);               /* [-0,+1,v] */
    filename            = qlua_strdup(L, luaL_checkstring(L, 2));       /* [-0,+0,v] */
    if (NULL == filename)
        luaL_error(L, "invalid filename or strdup error");
    
    int n_evec          = qlua_checkint(L, 3, "number of eigenvectors to load");
    int rank_stride     = 1;
    int verbose         = 0;
    if (qlua_checkopt_paramtable(L, 4)) {
        rank_stride = qlua_tabkey_intopt(L, 4, "rank_stride",  1);
        verbose     = qlua_tabkey_intopt(L, 4, "verbose",      0);
    }                                                           /* [-0,+0,-] */

    /* get mesh */
    if (4 != S->rank)
        luaL_error(L, "dim!=4 is not supported");
    int n_nodes = QMP_get_number_of_nodes();
    int subvol_lo[QLUA_MAX_LATTICE_RANK],
        subvol_hi[QLUA_MAX_LATTICE_RANK],
        subvol_dim[QLUA_MAX_LATTICE_RANK];
    qlua_sublattice(subvol_lo, subvol_hi, S->node, (void*)S);
    for (int i = 0 ; i < S->rank ; i++)
        subvol_dim[i] = subvol_hi[i] - subvol_lo[i];
    if (0 != subvol_dim[0] % 2) 
        luaL_error(L, "subvol_dim[0] must be even");


    CALL_QDP(L);
    
    /* create file names */
    int n_groups = rank_stride;
    int my_i_group  = S->node % n_groups,     /* sequence in writing */
        my_i_dir    = S->node / n_groups,     /* directory to save to */
        my_rank     = S->node;
    char *dirname = ".", 
         *basename = filename,
         *k;
    if (NULL != (k = strrchr(filename, '/'))) {
        *k = '\0';          /* separate dirname, basename */
        dirname = filename;
        basename= k + 1;
    }
    if (q_DW_evec_eopc_rawdump_files(&evecs_dir, &evecs_file, &evals_file, &cksum_file,
            dirname, basename, my_rank, my_i_dir, my_i_group))
        luaL_error(L, "not enough memory");
    
    /* alloc memory */    
    size_t evec_size_real = QDP_sites_on_node_L(S->lat) / 2 * c->Ls * evec_ncol * evec_nspin * 2;
    size_t evec_size_byte = evec_size_real * wordsize;
    evecs_buf   = qlua_malloc(L, evec_size_byte);
    crc32_buf   = qlua_malloc(L, n_nodes * sizeof(crc32_buf[0]));
    if (NULL == evecs_buf || NULL == crc32_buf)
        luaL_error(L, "not enough memory"); 
   
    DW_general_env w_env;
    w_env.lat   = S->lat;
    w_env.buf   = evecs_buf;
    w_env.nspin = evec_nspin;   /* avoid QDP defaults */
    w_env.ncol  = evec_ncol;    /* avoid QDP defaults */
    w_env.ndim  = S->rank;
    w_env.subvol_dim = subvol_dim;
    w_env.subvol_lo  = subvol_lo;
    w_env.L5    = c->Ls;
    uint32_t crc32 = 0x0;

    if (QOP_MDWF_gauge_float_from_double(&gaugeF, c->gauge))
        return luaL_error(L, "MDWF_gauge_float_from_double() failed");

    if (QOP_F3_MDWF_deflator_start_load(d->deflator))
        return luaL_error(L, "MDWF_deflator_stop_load() failed");
    FILE *f_evec_in = NULL;
    if (NULL == (f_evec_in = fopen(evecs_file, "r")))
        luaL_error(L, "%s: %s", evecs_file, strerror(errno));
    if (setvbuf(f_evec_in, NULL, _IOFBF, EVEC_BUFSIZE))
        printf("[%4d] cannot set bufsize=%d\n", S->node, EVEC_BUFSIZE);

    for (int i_evec = 0 ; i_evec < n_evec ; i_evec++) {
#if EVEC_STAGGER_LOAD 
        for (int i_group = 0 ; i_group < n_groups ; i_group++) { if (i_group == my_i_group) {
#endif
        if (1 != fread_shortproof(evecs_buf, evec_size_byte, 1, f_evec_in, evecs_file))
            luaL_error(L, "%s: %s", evecs_file, strerror(errno));
#if EVEC_STAGGER_LOAD 
        }  QMP_barrier();  }
#endif

        /* update crc and convert file->machine */
        crc32 = crc32_fast(evecs_buf, evec_size_byte, crc32);
        if (   (mach_bigendian && !file_bigendian) 
            || (!mach_bigendian && file_bigendian) )
            swap_endian(evecs_buf, wordsize, evec_size_real);
        
        if (QOP_F3_MDWF_import_half_fermion(&c_evec, c->state, q_DW_evec_eopc_rawload_F, &w_env)) 
            luaL_error(L, "MDWF_import_half_fermion() failed");

        if (QOP_F3_MDWF_deflator_add_vector(c->params, gaugeF, d->deflator, c_evec))
            luaL_error(L, QOP_MDWF_error(c->state));

        if (NULL != c_evec)         
            QOP_F3_MDWF_free_half_fermion(&c_evec);
        if (0 == QDP_this_node) {
            time_t tv = time(NULL);
            printf("DF_evecs_rawload: loaded evec[%05d] %s", 
                    i_evec, ctime_r(&tv, strbuf));
        }
    }
    fclose(f_evec_in);
    if (QOP_F3_MDWF_deflator_stop_load(d->deflator))
        return luaL_error(L, "MDWF_deflator_start_load() failed");

    QOP_F3_MDWF_deflator_eigcg_stop(d->deflator); /* do not run eigcg by default */

    /* read crc32 from file and check */
    FILE *f_cksum_in = NULL;
    if (NULL == (f_cksum_in = fopen(cksum_file, "r")))
        luaL_error(L, "%s: %s", cksum_file, strerror(errno));
    for (int i = 0 ; i < n_nodes ; i++)
        fscanf(f_cksum_in, "%lX", crc32_buf + i);
    fclose(f_cksum_in);

    if (crc32_buf[S->node] != crc32) {
        snprintf(strbuf, sizeof(strbuf), 
                "[%4d] %s: crc32 mismatch: %lX(data) != %lX(expected)", 
                S->node, evecs_file, (unsigned long)crc32, 
                (unsigned long)(crc32_buf[S->node]));
        fprintf(stderr, "%s\n", strbuf);
        luaL_error(L, strbuf);
    }
    lua_pushnumber(L, n_evec);  /* return the number of loaded e.vectors */
    
    /* cleanup */
#define FREE_NONNULL(p) do { if (NULL != (p)) { free(p) ; (p) = NULL; } } while(0)
    FREE_NONNULL(evecs_buf);
    FREE_NONNULL(crc32_buf);
    FREE_NONNULL(filename);
    FREE_NONNULL(evecs_dir);
    FREE_NONNULL(evecs_file);
    FREE_NONNULL(evals_file);
    FREE_NONNULL(cksum_file);
#undef FREE_NONNULL
    if (NULL != gaugeF)         QOP_F3_MDWF_free_gauge(&gaugeF);

    return 1;


    lua_pushnumber(L, n_evec);  /* return the number of loaded e.vectors */
    return 1;
}

/* call:
   * DeflatorObj
   * {DiracFermion, ...} (Ls total)
   return current number of vectors */
static int
q_DF_add_vector(lua_State *L)
{
    struct QOP_F3_MDWF_Gauge *gaugeF;

    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1); /* push1 */
    mMDWF *c = q_Deflator_get_MDWF(L, 1, NULL, 1); /* push1 */
    mLattice *S = qlua_ObjLattice(L, -1);
    /*int Sidx = lua_gettop(L);*/

    if (!qlua_checkopt_table(L, 2))
        return luaL_error(L, "expect table, got %s",
                          lua_typename(L, lua_type(L, 2)));
    if (lua_objlen(L, 2) != c->Ls)
        return luaL_error(L, "expect Ls=%d DiracFermions",
                          c->Ls);

    int n_vec = QOP_F3_MDWF_deflator_current_dim(d->deflator);
    if (d->umax <= n_vec)
        return luaL_error(L, "vector space is full");

    /* import table {DiracFermion, ...} into MDWF HalfFermion */
    struct QOP_F3_MDWF_HalfFermion *c_psi;
    mLatDirFerm3 **psi = qlua_malloc(L, c->Ls * sizeof (mLatDirFerm3 *));
    QLA_D3_DiracFermion **e_psi = qlua_malloc(L, c->Ls * sizeof (QLA_D3_DiracFermion *));
    double normi, norm5 = 0.;
    int i;
    CALL_QDP(L);
    for (i = 0 ; i < c->Ls ; i++) {
        lua_pushnumber(L, 1 + i); /* Lua index = 1..Ls */
        lua_gettable(L, 2);
        psi[i] = qlua_checkLatDirFerm3(L, -1, S, 3); /* push ; if fails, leak here? */
        QDP_D3_r_eq_norm2_D(&normi, psi[i]->ptr, S->all);
        norm5 += normi;
        lua_pop(L, 1);
    }
    if (norm5 <= 0.)
      return luaL_error(L, "zero-norm vector");
    norm5 = sqrt(norm5);

    for (i = 0 ; i < c->Ls ; i++)
        e_psi[i] = QDP_D3_expose_D(psi[i]->ptr);

    DW_5_env    r_env;
    r_env.lat   = S->lat;
    r_env.f     = e_psi;
    r_env.s     = 1. / norm5;
    if (QOP_F3_MDWF_import_half_fermion(&c_psi, c->state, q_DW_5_reader_scaled, &r_env))
      return luaL_error(L, "MDWF_import_fermion() failed");
    if (QOP_MDWF_gauge_float_from_double(&gaugeF, c->gauge))
      return luaL_error(L, "MDWF_gauge_float_from_double() failed");
    if (QOP_F3_MDWF_deflator_add_vector(c->params, gaugeF, d->deflator, c_psi))
      return luaL_error(L, QOP_MDWF_error(c->state));

    /* cleanup */
    QOP_F3_MDWF_free_gauge(&gaugeF);
    QOP_F3_MDWF_free_half_fermion(&c_psi);

    for (i = 0 ; i < c->Ls ; i++)
        QDP_D3_reset_D(psi[i]->ptr);

    qlua_free(L, psi);
    qlua_free(L, e_psi);

    /* normal return */
    lua_pushnumber(L, QOP_F3_MDWF_deflator_current_dim(d->deflator));
    return 1;
}

static int
q_DF_get_vector(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);
    mMDWF *c = q_Deflator_get_MDWF(L, 1, NULL, 1);
    mLattice *S = qlua_ObjLattice(L, -1);
    int Sidx = lua_gettop(L);

    int num_vec = QOP_F3_MDWF_deflator_current_dim(d->deflator);
    int idx_vec = qlua_checkint(L, 2, "expect vector index");
    if (idx_vec < 0 || num_vec <= idx_vec)
        return luaL_error(L, "expect vector index 0 <= i < dim");

    /* resulting table */
    lua_createtable(L, c->Ls, 0);
    int Tidx = lua_gettop(L);

    /* export (c->Ls) DiracFermion's */
    struct QOP_F3_MDWF_HalfFermion *c_psi;
    mLatDirFerm3 **psi = qlua_malloc(L, c->Ls * sizeof (mLatDirFerm3 *));
    QLA_D3_DiracFermion **e_psi = qlua_malloc(L, c->Ls * sizeof (QLA_D3_DiracFermion *));
    int i;

    CALL_QDP(L);
    for (i = 0 ; i < c->Ls ; i++) {
        psi[i]  = qlua_newLatDirFerm3(L, Sidx, 3);
        QDP_D3_D_eq_zero(psi[i]->ptr, S->all);
        e_psi[i]= QDP_D3_expose_D(psi[i]->ptr);
    }

    DW_5_env w_env;
    w_env.lat = S->lat;
    w_env.f   = e_psi;
    w_env.s   = 1.;
    if (QOP_F3_MDWF_allocate_half_fermion(&c_psi, c->state))
      return luaL_error(L, "MDWF_allocate_half_fermion() failed");
    if (QOP_F3_MDWF_deflator_extract_vector(c_psi, d->deflator, idx_vec))
      return luaL_error(L, QOP_MDWF_error(c->state));
    if (QOP_F3_MDWF_export_half_fermion(q_DW_5_writer_scaled, &w_env, c_psi))
      return luaL_error(L, "MDWF_export_half_fermion() failed");

    for (i = 0 ; i < c->Ls ; i++) {
        QDP_D3_reset_D(psi[i]->ptr);
        lua_rawseti(L, Tidx, c->Ls - i); /* reverse order; qlua index = 1..Ls */
    }

    /* clean up */
    QOP_F3_MDWF_free_half_fermion(&c_psi);
    qlua_free(L, psi);
    qlua_free(L, e_psi);
    /* normal return */
    return 1;
}


static int
q_DF_deflated_mixed_solver(lua_State *L,
                           struct QOP_D3_MDWF_Fermion *solution,
                           int *out_iters,
                           double *out_epsilon,
                           const struct QOP_D3_MDWF_Fermion *rhs,
                           int log_level)
{
    mMDWF          *c = qlua_checkMDWF(L, lua_upvalueindex(2), NULL, 1);
    mLattice       *S = qlua_ObjLattice(L, lua_upvalueindex(2));
    mDeflatorState *d = q_checkDeflatorState(L, lua_upvalueindex(3), S, 1);
    double      f_eps = luaL_checknumber(L, lua_upvalueindex(4));
    int   inner_iters = luaL_checkint(L, lua_upvalueindex(5));
    double        eps = luaL_checknumber(L, lua_upvalueindex(6));
    int     max_iters = luaL_checkint(L, lua_upvalueindex(7));

    lua_pop(L, 1);
    return QOP_MDWF_deflated_mixed_D_CG(solution, out_iters, out_epsilon,
                                        c->params, rhs, c->gauge, rhs,
                                        d->deflator,
                                        inner_iters, f_eps,
                                        max_iters, eps,
                                        log_level);
}

static MDWFSolver deflated_mixed_solver = {
    q_DF_deflated_mixed_solver, "eigCG"
};

static int
q_DF_make_mixed_solver(lua_State *L)
{
    double inner_eps = luaL_checknumber(L, 2);
    int inner_iter = luaL_checkint(L, 3);
    double eps = luaL_checknumber(L, 4);
    int max_iter = luaL_checkint(L, 5);

    lua_pushlightuserdata(L, &deflated_mixed_solver); /* clo[1]: solver */
    q_Deflator_get_MDWF(L, 1, NULL, 1);               /* clo[2]: MDWF */
    q_Deflator_get_State(L, 1, NULL, 1);              /* clo[3]: Deflator */
    lua_pushnumber(L, inner_eps);                     /* clo[4]: inner_eps */
    lua_pushnumber(L, inner_iter);                    /* clo[5]: inner_iter */
    lua_pushnumber(L, eps);                           /* clo[6]: epsilon */
    lua_pushnumber(L, max_iter);                      /* clo[7]: max_iter */
    lua_pushcclosure(L, q_dirac_solver, 7);

    return 1;
}

static struct luaL_Reg mtMDWFDeflator[] = {
    { "__tostring",         q_DF_fmt               },
    { "__newindex",         qlua_nowrite           },
    { "mixed_solver",       q_DF_make_mixed_solver },
    { "close",              q_DF_close             },
    { "reset",              q_DF_reset             },
    { "stop",               q_DF_stop              },
    { "resume",             q_DF_resume            },
    { "eigenvalues",        q_DF_eigenvalues       },
    { "current_dim",        q_DF_current_dim       },
    { "start_load",         q_DF_start_load        },
    { "stop_load",          q_DF_stop_load         },
    { "add_vector",         q_DF_add_vector        },
    { "get_vector",         q_DF_get_vector        },
    { "evecs_rawdump",      q_DF_evecs_rawdump     },
    { "evecs_rawload",      q_DF_evecs_rawload     },
#if 0 /* XXX deflator extra methods */
    { "truncate",           q_DF_truncate          },
    { "get_counter",        q_DF_get_counter       },
    { "put_counter",        q_DF_put_counter       },
    { "write_eigspace",     q_DF_write             },
    { "read_eigspace",      q_DF_read              },
#endif /* XXX deflator extra methods */
    { NULL,                 NULL                   }
};

static int
q_DW_make_deflator(lua_State *L)
{
    mMDWF *c = qlua_checkMDWF(L, 1, NULL, 1);
    qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    int vmax = luaL_checkint(L, 2);
    int nev = luaL_checkint(L, 3);
    double eps = luaL_checknumber(L, 4);
    int umax = luaL_checkint(L, 5);

    if (nev <= 0 || vmax <= 0 || umax <= 0)
        return luaL_error(L, "bad eigenspace size");
    if (2 * nev >= vmax)
        return luaL_error(L, "eigcg VMAX: must satisfy VMAX > 2*NEV");

    if ((c->state == 0) || (c->gauge == 0))
        return luaL_error(L, "closed MDWF used");

    lua_createtable(L, 2, 0);
    lua_pushvalue(L, 1);
    lua_rawseti(L, -2, 1);
    mDeflatorState *d = q_newDeflatorState(L, Sidx);
    d->nev = nev;
    d->vmax = vmax;
    d->umax = umax;

    CALL_QDP(L);
    if (QOP_F3_MDWF_create_deflator(&d->deflator, c->state,
                                    vmax, nev, eps, umax))
        return luaL_error(L, "MDWF_create_deflator() failed");

    lua_rawseti(L, -2, 2);
    qlua_createLatticeTable(L, Sidx, mtMDWFDeflator, qMDWFDeflator,
                            MDWFDeflatorName);
    lua_setmetatable(L, -2);

    return 1;
}

#ifdef HAS_ARPACK
/* operator for lanczos : function and oblique arg */
typedef struct {
    struct QOP_MDWF_State           *mdwf_state;
    struct QOP_MDWF_Parameters      *mdwf_params;
    struct QOP_F3_MDWF_Gauge        *mdwf_gauge;

    /* workspace: must be allocated before calling mdwf_eoprec_op */
    struct QOP_F3_MDWF_HalfFermion  *x,
                                    *y;

    /* polynomial acc parameters :
       compute orthogonal polynomial of 'poly_n'-degree of the Op
       using 3-term recursion
       p_{i+1}(Op).x := [(poly_a[i] + poly_b[i]*Op) . p_{i}(Op)
                      + poly_c[i]*p_{i-1}(Op)] . x
       where n = 0 .. {poly_n-1}, p_0(Op).x := x, p_{-1}.x := 0
     */
    int                             poly_n; /* degree */
    double                          *poly_a,
                                    *poly_b,
                                    *poly_c;
} op_MDWF_F3_eoprec_MdagM_arg_s;

void
op_MDWF_F3_eoprec_MdagM_op(
        int loc_dim,
        float complex *x,
        float complex *y,
        void *op_arg) /* x<-op(y) */
{
    long long fl1, fl2;
    double t1, t2;

    op_MDWF_F3_eoprec_MdagM_arg_s *a = (op_MDWF_F3_eoprec_MdagM_arg_s *)op_arg;
    QLUA_ASSERT(2 * loc_dim == QOP_MDWF_half_fermion_size(a->mdwf_state));

    QOP_F3_MDWF_half_fermion_from_blas(a->y, (float *)y, 2 * loc_dim);
    if (0 < a->poly_n) {
        QOP_F3_MDWF_MxM_poly(
                a->x, NULL, a->mdwf_params, a->mdwf_gauge, a->y,
                a->poly_n, a->poly_a, a->poly_b, a->poly_c);
        QOP_MDWF_performance(&t1, &fl1, NULL, NULL, a->mdwf_state);

        QOP_F3_MDWF_blas_from_half_fermion((float *)x, 2 * loc_dim, a->x);
    } else {
        QOP_F3_MDWF_M_operator(
                a->x, a->mdwf_params, a->mdwf_gauge, a->y);
        QOP_MDWF_performance(&t1, &fl1, NULL, NULL, a->mdwf_state);

        QOP_F3_MDWF_M_operator_conjugated(
                a->y, a->mdwf_params, a->mdwf_gauge, a->x);
        QOP_MDWF_performance(&t2, &fl2, NULL, NULL, a->mdwf_state);
        t1  += t2;
        fl1 += fl2;

        QOP_F3_MDWF_blas_from_half_fermion((float *)x, 2 * loc_dim, a->y);
    }

    if (t1 == 0)
        t1 = -1;
    if (QDP_this_node == qlua_master_node)
        printf("MDWF MdagM(poly_n=%d): time = %.3f sec,"
               " perf = %.2f MFlops/sec\n",
               (0 < a->poly_n ? a->poly_n : 1),
               t1, fl1 * 1e-6 / t1);
}

/* double precision operator */
typedef struct {
    struct QOP_MDWF_State           *mdwf_state;
    struct QOP_MDWF_Parameters      *mdwf_params;
    struct QOP_D3_MDWF_Gauge        *mdwf_gauge;

    /* workspace: must be allocated before calling mdwf_eoprec_op */
    struct QOP_D3_MDWF_HalfFermion  *x,
                                    *y;

    /* polynomial acc parameters :
       compute orthogonal polynomial of 'poly_n'-degree of the Op
       using 3-term recursion
       p_{i+1}(Op).x := [(poly_a[i] + poly_b[i]*Op) . p_{i}(Op)
                      + poly_c[i]*p_{i-1}(Op)] . x
       where n = 0 .. {poly_n-1}, p_0(Op).x := x, p_{-1}.x := 0
     */
    int                             poly_n; /* degree */
    double                          *poly_a,
                                    *poly_b,
                                    *poly_c;
} op_MDWF_D3_eoprec_MdagM_arg_s;

void
op_MDWF_F3_eoprec_MdagM_double_op(
        int loc_dim,
        float complex *x,
        float complex *y,
        void *op_arg) /* x<-op(y) */
{
    long long fl1, fl2;
    double t1, t2;
    int i;
    double complex *d_buf = NULL;

    op_MDWF_D3_eoprec_MdagM_arg_s *a = (op_MDWF_D3_eoprec_MdagM_arg_s *)op_arg;
    QLUA_ASSERT(2 * loc_dim == QOP_MDWF_half_fermion_size(a->mdwf_state));

    d_buf = malloc(sizeof(d_buf[0]) * loc_dim);
    QLUA_ASSERT(NULL != d_buf);

    for (i = 0 ; i < loc_dim ; i++)
        d_buf[i] = y[i];
    QOP_D3_MDWF_half_fermion_from_blas(a->y, (double *)d_buf, 2 * loc_dim);

    if (0 < a->poly_n) {
        QOP_D3_MDWF_MxM_poly(
                a->x, NULL, a->mdwf_params, a->mdwf_gauge, a->y,
                a->poly_n, a->poly_a, a->poly_b, a->poly_c);
        QOP_MDWF_performance(&t1, &fl1, NULL, NULL, a->mdwf_state);

        QOP_D3_MDWF_blas_from_half_fermion((double *)d_buf, 2 * loc_dim, a->x);
        for (i = 0 ; i < loc_dim ; i++)
            x[i] = d_buf[i];
    } else {
        QOP_D3_MDWF_M_operator(
                a->x, a->mdwf_params, a->mdwf_gauge, a->y);
        QOP_MDWF_performance(&t1, &fl1, NULL, NULL, a->mdwf_state);

        QOP_D3_MDWF_M_operator_conjugated(
                a->y, a->mdwf_params, a->mdwf_gauge, a->x);
        QOP_MDWF_performance(&t2, &fl2, NULL, NULL, a->mdwf_state);
        t1  += t2;
        fl1 += fl2;

        QOP_D3_MDWF_blas_from_half_fermion((double *)d_buf, 2 * loc_dim, a->y);
        for (i = 0 ; i < loc_dim ; i++)
            x[i] = d_buf[i];
    }

    free(d_buf);

    if (t1 == 0)
        t1 = -1;
    if (QDP_this_node == qlua_master_node)
        printf("MDWF MdagM(poly_n=%d): time = %.3f sec,"
               " perf = %.2f MFlops/sec\n",
               (0 < a->poly_n ? a->poly_n : 1),
               t1, fl1 * 1e-6 / t1);
}

/* MDWF:eig_deflator_lanczos(nev, ncv, max_iter, tol, [param])
   param = {
     [cheb_accel]={n, a, b},
     [eigcg] = {vmax, nev, eps, umax}
     ... }
   return: deflator object, n_converged, n_iter

   notes:
     * first, Lanczos/Arnoldi is used to compute eigenvectors;
     * converged vectors are saved into a newly created EigCG deflator,
       which is set to frozen state;
     * only min(eigcg_umax, #converged) vectors will be saved;
     * the user may specify eigcg_umax, eigcg_vmax, eigcg_nev values
       to continue EigCG after Lanczos (which is perhaps pointless though);
     * if any of (eigcg_umax, eigcg_vmax, eigcg_nev) are set to invalid
       values, the code will figure out the best values to continue operation
 */

static int
q_DW_make_deflator_lanczos(lua_State *L)
{
    const char *err_str = NULL;
#define LANCZOS_MXM_DOUBLE  1
    /* by default, search for ev with smallest real part */
    const char *lanczos_which= "SR";
    const char *arpack_logfile = NULL;
    struct QOP_F3_MDWF_HalfFermionMat *hfm = NULL;
    /* operator parameters, init to empty */
    struct QOP_F3_MDWF_Gauge *gaugeF = NULL;

    /* v0 import, all stay NULL unless v0 is specified */
    mLatDirFerm3 **qlua_v0 = NULL;
    QLA_D3_DiracFermion **qla_v0 = NULL;
    struct QOP_F3_MDWF_HalfFermion *lib_v0 = NULL;
    float complex *blas_v0   = NULL;

#ifdef LANCZOS_MXM_DOUBLE
    op_MDWF_D3_eoprec_MdagM_arg_s op_arg;
#else
    op_MDWF_F3_eoprec_MdagM_arg_s op_arg;
#endif/*LANCZOS_MXM_DOUBLE*/
    op_arg.poly_n = -1;
    op_arg.poly_a = op_arg.poly_b = op_arg.poly_c = NULL;
    op_arg.x = op_arg.y = NULL;

    int status;
    int do_lanczos_inplace = 0;
    int loc_dim;
    int n_iters, nconv;
    int n_evecs;
    int i;
    int Sidx;
    float complex *evec = NULL,
                  *eval = NULL;

    mLattice *S = qlua_ObjLattice(L, 1);
    Sidx = lua_gettop(L);
    mMDWF *c = qlua_checkMDWF(L, 1, NULL, 1);
    if (NULL == c->state || NULL == c->gauge)
      return luaL_error(L, "closed MDWF used");

    mDeflatorState *d = NULL;

    /* parse parameters */
    int nev, ncv, max_iter;
    double tol;
    nev     = qlua_checkint(L, 2, "expect NEV at #2");
    ncv     = qlua_checkint(L, 3, "expect NCV at #3");
    max_iter= qlua_checkint(L, 4, "expect MAX_ITER at #4");
    tol     = luaL_checknumber(L, 5);

    /* parse optional parameters */
    int eigcg_vmax  = 0,
        eigcg_umax  = 0,
        eigcg_nev   = 0;
    double eigcg_eps= 0.;

    if (qlua_checkopt_paramtable(L, 6)) {
        if (qlua_tabpushopt_key(L, 6, "cheb_accel")) {
            /* Chebyshev acceleration parameters */
#ifdef LANCZOS_MXM_DOUBLE
            op_MDWF_D3_eoprec_MdagM_arg_s *a = &op_arg;
#else
            op_MDWF_F3_eoprec_MdagM_arg_s *a = &op_arg;
#endif
            int cheb_n = -1;
            double cheb_a = 0.,
                   cheb_b = 0.,
                   cheb_x0= 1.;
            int do_norm = 0;

            if (0 <= op_arg.poly_n)
              return luaL_error(L, "more than one poly.accel. parameter");

            cheb_n  = qlua_tabidx_int(L, -1, 1);
            if (cheb_n < 0)
              return luaL_error(L, "poly.degree must be positive");

            cheb_a  = qlua_tabidx_double(L, -1, 2);
            cheb_b  = qlua_tabidx_double(L, -1, 3);
            if (cheb_a == cheb_b)
              return luaL_error(L, "invalid segment [a;b]");

            if (qlua_tabpushopt_idx(L, -1, 4)) {
                do_norm = 1;
                cheb_x0 = luaL_checknumber(L, -1);
                lua_pop(L, 1);
            }

            /* set up poly parameters for P_n(x) = T_n(x'), with
                x:[a;b] -> x':[-1;+1] */
            a->poly_a = a->poly_b = a->poly_c = NULL;
            a->poly_a = qlua_malloc(L, cheb_n * sizeof(a->poly_a[0]));
            a->poly_b = qlua_malloc(L, cheb_n * sizeof(a->poly_b[0]));
            a->poly_c = qlua_malloc(L, cheb_n * sizeof(a->poly_c[0]));
            if (NULL == a->poly_a
                    || NULL == a->poly_b
                    || NULL == a->poly_c)
              return luaL_error(L, "not enough memory");

            a->poly_a[0] = (-cheb_b - cheb_a) / (cheb_b - cheb_a);
            a->poly_b[0] = 2. / (cheb_b - cheb_a);
            a->poly_c[0] = 1;
            for (i = 1; i < cheb_n ; i++) {
                a->poly_a[i] = 2 * a->poly_a[0];
                a->poly_b[i] = 2 * a->poly_b[0];
                a->poly_c[i] = -1;
            }

            if (do_norm)
                QOP_MDWF_poly_normalize(cheb_n,
                        a->poly_a, a->poly_b, a->poly_c,
                        cheb_x0, 1e-8);

            a->poly_n = cheb_n;

            lua_pop(L, 1);
        }

        if (qlua_tabpushopt_key(L, 6, "which")) {
            lanczos_which = luaL_checkstring(L, -1);
            if (NULL == lanczos_which ||
                    (  strcmp("SR", lanczos_which)
                    && strcmp("LR", lanczos_which)
                    && strcmp("SI", lanczos_which)
                    && strcmp("LI", lanczos_which)
                    && strcmp("SM", lanczos_which)
                    && strcmp("LM", lanczos_which)))
              return luaL_error(L, "invalid value for which='%s'",
                                NULL == lanczos_which ? "null" : lanczos_which);
            lua_pop(L, 1);
        }

        if (qlua_tabpushopt_key(L, 6, "eigcg")) {
            eigcg_vmax  = qlua_tabidx_int(L, -1, 1);
            eigcg_nev   = qlua_tabidx_int(L, -1, 2);
            eigcg_eps   = qlua_tabidx_double(L, -1, 3);
            eigcg_umax  = qlua_tabidx_int(L, -1, 4);
            lua_pop(L, 1);
        }

        if (qlua_tabpushopt_key(L, 6, "arpack_logfile")) {
            arpack_logfile = luaL_checkstring(L, -1);
            lua_pop(L, 1);
        }

        if (qlua_tabpushopt_key(L, 6, "inplace")) {
            if (lua_isboolean(L, -1))
              do_lanczos_inplace = lua_toboolean(L, -1);
            else
              return luaL_error(L, "'inplace' : expect boolean");
            lua_pop(L, 1);
        }
        if (qlua_tabpushopt_key(L, 6, "v0")) {
            /* parse D5 fermion v0 */
            qlua_v0 = qlua_malloc(L, c->Ls * sizeof(mLatDirFerm3 *));
            qla_v0  = qlua_malloc(L, c->Ls * sizeof(QLA_D3_DiracFermion *));
            if (NULL == qlua_v0 || NULL == qla_v0)
                return luaL_error(L, "not enough memory");
            for (i = 0; i < c->Ls ; i++) {  /* FIXME:QLOPT */
                lua_pushnumber(L, i + 1);   /* [sic] lua indexing */
                lua_gettable(L, -2);        /* stack top: tab, index */
                qlua_v0[i] = qlua_checkLatDirFerm3(L, -1, S, 3);
                lua_pop(L, 1);
            }
            lua_pop(L, 1);
        }
    }

#if 0    // no auto-setting of eigcg
    if (eigcg_vmax <= 0)
        eigcg_vmax = ncv;
    if (eigcg_nev <= 0) {
        eigcg_nev = 2; /* quite arbitrary, honestly */
        if (eigcg_vmax < 2 * eigcg_nev)
            eigcg_nev = eigcg_vmax / 2;
    }
    if (eigcg_vmax < 2 * eigcg_nev)
      return luaL_error(L, "eigcg VMAX: must satisfy VMAX > 2*NEV");
#endif


    /* create Qlua deflator object {MDWF, DeflatorState} and set its META */
    /* FIXME should be refactored into a separate function?;
       this code duplicates parts of q_DW_make_deflator(qmdwf.c) */
    qlua_ObjLattice(L, 1);
    lua_createtable(L, 2, 0);
    lua_pushvalue(L, 1);
    lua_rawseti(L, -2, 1);
    if (NULL == (d = q_newDeflatorState(L, Sidx)))
      return luaL_error(L, "cannot create deflator state");
    lua_rawseti(L, -2, 2);
    qlua_createLatticeTable(L, Sidx, mtMDWFDeflator, qMDWFDeflator,
            MDWFDeflatorName);
    lua_setmetatable(L, -2);



    CALL_QDP(L);

    gaugeF = NULL;
    if (QOP_MDWF_gauge_float_from_double(&gaugeF, c->gauge))
      return luaL_error(L, "QOP_MDWF_gauge_float_from_double() failed");

#ifdef LANCZOS_MXM_DOUBLE
    op_arg.mdwf_state = c->state;
    op_arg.mdwf_params= c->params;
    op_arg.mdwf_gauge = c->gauge;

    op_arg.x = op_arg.y = NULL;
    if (QOP_D3_MDWF_allocate_half_fermion(&op_arg.x, c->state)
            || QOP_D3_MDWF_allocate_half_fermion(&op_arg.y, c->state))
      return luaL_error(L, "cannot allocate HalfFermion");
#else
    op_arg.mdwf_state   = c->state;
    op_arg.mdwf_params  = c->params;
    if (QOP_MDWF_gauge_float_from_double(&(op_arg.mdwf_gauge), c->gauge))
      return luaL_error(L, "MDWF_gauge_float_from_double() failed");

    op_arg.x = op_arg.y = NULL;
    if (QOP_F3_MDWF_allocate_half_fermion(&op_arg.x, c->state)
            || QOP_F3_MDWF_allocate_half_fermion(&op_arg.y, c->state))
      return luaL_error(L, "cannot allocate HalfFermion");
#endif/*LANCZOS_MXM_DOUBLE*/

    MPI_Comm mpi_comm = MPI_COMM_WORLD; /* FIXME any better choice? */

    QLUA_ASSERT(0 == QOP_MDWF_half_fermion_size(c->state) % 2);
    loc_dim = QOP_MDWF_half_fermion_size(c->state) / 2; /* sic! 
                loc_dim = number of cplx = number of real / 2 */

    /* import v0 : qdp -> qla -> lib -> blas */
    if (NULL != qlua_v0) {
        qla_v0  = malloc(c->Ls * sizeof(qla_v0[0]));
        blas_v0 = malloc(loc_dim * sizeof(blas_v0[0]));
        if (NULL == qla_v0  || NULL == blas_v0) {
            err_str = "not enough memory";
            goto err_end;
        };
        double norm2i, v0_norm5_sq = 0.;
        for (i = 0; i < c->Ls ; i++) {  
            QDP_D3_r_eq_norm2_D(&norm2i, qlua_v0[i]->ptr, S->all);
            v0_norm5_sq += norm2i;
            qla_v0[i] = QDP_D3_expose_D(qlua_v0[i]->ptr);
        }
        if (v0_norm5_sq <= 0.) {
            err_str = "zero initial vector";
            goto err_end;
        }
        DW_5_env env;
        env.lat = S->lat;
        env.f   = qla_v0;
//        env.s   = 1. / sqrt(v0_norm5_sq);
        env.s   = 1.;   /* no rescaling : eigensolver will do it */
        if (QOP_F3_MDWF_import_half_fermion(&lib_v0, c->state, q_DW_5_reader_scaled, &env)) {
            err_str = "MDWF_import_half_fermion() failed";
            goto err_end;
        }
        for (i = 0; i < c->Ls ; i++)
            QDP_D3_reset_D(qlua_v0[i]->ptr);
        QOP_F3_MDWF_blas_from_half_fermion((float *)blas_v0, 2 * loc_dim, lib_v0);

        free(qla_v0);  qla_v0 = NULL;
        QOP_F3_MDWF_free_half_fermion(&lib_v0);  lib_v0 = NULL;
    }


    /* run Arnoldi/Lanczos iterations */
    n_iters = nconv = 0;
    evec = eval = NULL;
    if (do_lanczos_inplace) {
        int inplace_umax = ncv < eigcg_umax ? eigcg_umax : ncv;
        float *hfm_blas_ptr = NULL;
        int hfm_nrow_loc = 0,
            hfm_ncol     = 0,
            hfm_ld       = 0;

        if (NULL == (eval = qlua_malloc(L, sizeof(eval[0]) * nev)))
          return luaL_error(L, "not enough memory");
        if (QOP_F3_MDWF_alloc_half_fermion_matrix(&hfm, c->state, inplace_umax))
          return luaL_error(L, "MDWF_alloc_half_fermion_matrix failed");
        if (QOP_F3_MDWF_blas_view_half_fermion_matrix(hfm, &hfm_nrow_loc, &hfm_ncol, &hfm_blas_ptr, &hfm_ld))
          return luaL_error(L, QOP_MDWF_error(c->state));

        if (0 != (status = lanczos_inplace_float(
                L, mpi_comm,
#ifdef LANCZOS_MXM_DOUBLE
                op_MDWF_F3_eoprec_MdagM_double_op,
#else
                op_MDWF_F3_eoprec_MdagM_op,
#endif/*LANCZOS_MXM_DOUBLE*/
                &op_arg,
                lanczos_which, loc_dim, nev, ncv, max_iter, tol, blas_v0,
                eval, (float complex *)hfm_blas_ptr, hfm_ld, hfm_ncol,
                &n_iters, &nconv, arpack_logfile))) {
            QOP_F3_MDWF_free_half_fermion_matrix(&hfm);
            return luaL_error(L, "lanczos_float_inplace returned %d", status);
        }
        if (QOP_F3_MDWF_create_deflator_inplace(
                    &(d->deflator), c->params, gaugeF, &hfm,
                    nconv, eigcg_vmax, eigcg_nev, eigcg_eps, eigcg_umax))
          return luaL_error(L, QOP_MDWF_error(c->state));

    } else {
        if (0 != (status = lanczos_float(
                L, mpi_comm,
#ifdef LANCZOS_MXM_DOUBLE
                op_MDWF_F3_eoprec_MdagM_double_op,
#else
                op_MDWF_F3_eoprec_MdagM_op,
#endif/*LANCZOS_MXM_DOUBLE*/
                &op_arg,
                lanczos_which, loc_dim, nev, ncv, max_iter, tol, blas_v0,
                &eval, &evec, &n_iters, &nconv, arpack_logfile)))
          return luaL_error(L, "lanczos_float returned %d", status);
        /* FIXME rewrite with clear explanation for the choice of
           default values */
        if (eigcg_umax <= 0 || eigcg_umax <= nconv)
            eigcg_umax = nconv;


        if (QOP_F3_MDWF_create_deflator(&d->deflator, c->state, eigcg_vmax, eigcg_nev, eigcg_eps, eigcg_umax))
          return luaL_error(L, "MDWF_create_deflator() failed");

        /* fill deflator with e.vecs */
        if (QOP_F3_MDWF_deflator_start_load(d->deflator))
          return luaL_error(L, "MDWF_deflator_start_load() failed");
        /* (have space for eigcg_umax) */
        n_evecs = (nconv <= eigcg_umax ? nconv : eigcg_umax);

        struct QOP_F3_MDWF_HalfFermion *hf_buf = NULL;
        if (QOP_F3_MDWF_allocate_half_fermion(&hf_buf, c->state))
          return luaL_error(L, "cannot allocate HalfFermion");

        for (i = 0 ; i < n_evecs ; i++) {
            QOP_F3_MDWF_half_fermion_from_blas(hf_buf,
                    (float *)(evec + i * loc_dim), 2 * loc_dim);
            if (QOP_F3_MDWF_deflator_add_vector(c->params, gaugeF,
                                                d->deflator, hf_buf)) {
                return luaL_error(L, "MDWF_deflator_add_vector() failed");
            }
        }

        if (QOP_F3_MDWF_deflator_stop_load(d->deflator))
            return luaL_error(L, "MDWF_deflator_end_load() failed");
    }

    /* initialize MDWF side of deflator */
    d->nev  = eigcg_nev;
    d->vmax = eigcg_vmax;
    d->umax = eigcg_umax;

    CALL_QDP(L);

    QOP_F3_MDWF_deflator_eigcg_stop(d->deflator);

err_end:
    /* cleanup */
    if (NULL != qlua_v0) qlua_free(L, qlua_v0);
    if (NULL != qla_v0)  free(qla_v0);
    if (NULL != lib_v0)  QOP_F3_MDWF_free_half_fermion(&lib_v0);
    if (NULL != blas_v0) free(blas_v0);

    if (NULL != evec) free(evec);
    if (NULL != eval) free(eval);
    if (NULL != gaugeF) QOP_F3_MDWF_free_gauge(&(gaugeF));
#ifdef LANCZOS_MXM_DOUBLE
    if (NULL != op_arg.y) QOP_D3_MDWF_free_half_fermion(&op_arg.y);
    if (NULL != op_arg.x) QOP_D3_MDWF_free_half_fermion(&op_arg.x);
#else
    if (NULL != op_arg.mdwf_gauge) QOP_F3_MDWF_free_gauge(&(op_arg.mdwf_gauge));
    if (NULL != op_arg.y) QOP_F3_MDWF_free_half_fermion(&op_arg.y);
    if (NULL != op_arg.x) QOP_F3_MDWF_free_half_fermion(&op_arg.x);
#endif/*LANCZOS_MXM_DOUBLE*/

    if (NULL != op_arg.poly_a) qlua_free(L, op_arg.poly_a);
    if (NULL != op_arg.poly_b) qlua_free(L, op_arg.poly_b);
    if (NULL != op_arg.poly_c) qlua_free(L, op_arg.poly_c);

    if (NULL != hfm)  QOP_F3_MDWF_free_half_fermion_matrix(&hfm);

    /* return error, if any */
    if (NULL != err_str)
        luaL_error(L, err_str);

    /* successful exit */
    lua_pushnumber(L, nconv);
    lua_pushnumber(L, n_iters);
    return 3; /* deflator object, n_converged, n_iter */
}

#endif /* HAS_ARPACK */


/***** clover interface */
static int
q_DW_fmt(lua_State *L)
{
    char fmt[72];
    mMDWF *c = qlua_checkMDWF(L, 1, NULL, 0);

    if (c->state)
        sprintf(fmt, "MDWF[%s]", c->name);
    else
        sprintf(fmt, "MDWF(closed)");

    lua_pushstring(L, fmt);

    return 1;
}

static int
q_DW_gc(lua_State *L)
{
    mMDWF *c = qlua_checkMDWF(L, 1, NULL, 0);

    if (c->gauge)
        QOP_D3_MDWF_free_gauge(&c->gauge);
    c->gauge = 0;

    if (c->params)
        QOP_MDWF_free_parameters(&c->params);
    c->params = 0;

    if (c->state)
        QOP_MDWF_fini(&c->state);
    c->state = 0;

    return 0;
}

static int
q_DW_close(lua_State *L)
{
    mMDWF *c = qlua_checkMDWF(L, 1, NULL, 1);

    if (c->gauge)
        QOP_D3_MDWF_free_gauge(&c->gauge);
    c->gauge = 0;

    if (c->params)
        QOP_MDWF_free_parameters(&c->params);
    c->params = 0;

    if (c->state)
        QOP_MDWF_fini(&c->state);
    c->state = 0;

    return 0;
}

static int
q_DW_operator(lua_State *L,
              const char *name,
              int (*op)(struct QOP_D3_MDWF_Fermion *result,
                        const struct QOP_MDWF_Parameters *params,
                        const struct QOP_D3_MDWF_Gauge *gauge,
                        const struct QOP_D3_MDWF_Fermion *source))
{
    mMDWF *c = qlua_checkMDWF(L, 1, NULL, 1);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);

    switch (qlua_qtype(L, 2)) {
    case qTable: {
        mLatDirFerm3 **psi = qlua_malloc(L, c->Ls * sizeof (mLatDirFerm3 *));
        mLatDirFerm3 **eta = qlua_malloc(L, c->Ls * sizeof (mLatDirFerm3 *));
        struct QOP_D3_MDWF_Fermion *c_psi;
        struct QOP_D3_MDWF_Fermion *c_eta;
        QLA_D3_DiracFermion **e_psi = qlua_malloc(L, c->Ls * sizeof (QLA_D3_DiracFermion *));
        QLA_D3_DiracFermion **e_eta = qlua_malloc(L, c->Ls * sizeof (QLA_D3_DiracFermion *));
        DW_5_env env;
                char *err_str = NULL;
        int i;

        CALL_QDP(L);
        for (i = 0; i < c->Ls; i++) {
            lua_pushnumber(L, i + 1); /* [sic] lua indexing */
            lua_gettable(L, 2);
            psi[i] = qlua_checkLatDirFerm3(L, -1, S, 3);
            e_psi[i] = QDP_D3_expose_D(psi[i]->ptr);
            lua_pop(L, 1);
        }

        env.lat = S->lat;
        env.f = e_psi;
        if (QOP_D3_MDWF_import_fermion(&c_psi, c->state, DW_5_reader, &env)) {
            err_str = "MDWF_import_fermion() failed";
                        goto err_end;
                }
        for (i = 0; i < c->Ls; i++)
            QDP_D3_reset_D(psi[i]->ptr);

        if (QOP_D3_MDWF_allocate_fermion(&c_eta, c->state)) {
            err_str = "MDWF_create_fermion() failed";
                        goto err_end;
                }

        (*op)(c_eta, c->params, c->gauge, c_psi);

        QOP_D3_MDWF_free_fermion(&c_psi);

        lua_createtable(L, c->Ls, 0);
        for (i = 0; i < c->Ls; i++) {
            eta[i] = qlua_newLatDirFerm3(L, Sidx, 3);
            e_eta[i] = QDP_D3_expose_D(eta[i]->ptr);
            lua_rawseti(L, -2, i + 1);
        }

        env.lat = S->lat;
        env.f = e_eta;
        if (QOP_D3_MDWF_export_fermion(DW_5_writer, &env, c_eta)) {
            err_str = "MDWF_export_fermion() failed";
                        goto err_end;
                }
        for (i = 0; i < c->Ls; i++) {
            QDP_D3_reset_D(eta[i]->ptr);
        }

        QOP_D3_MDWF_free_fermion(&c_eta);

                qlua_free(L, psi);
                qlua_free(L, eta);
                qlua_free(L, e_psi);
                qlua_free(L, e_eta);
        return 1;
                err_end:
                qlua_free(L, psi);
                qlua_free(L, eta);
                qlua_free(L, e_psi);
                qlua_free(L, e_eta);
                return luaL_error(L, err_str);
    }
    default:
        break;
    }

        return luaL_error(L, "bad arguments in MDWF:%s", name);
}

static int
q_DW_D(lua_State *L)
{
    return q_DW_operator(L, "D", QOP_D3_MDWF_DDW_operator);
}

static int
q_DW_Dx(lua_State *L)
{
    return q_DW_operator(L, "Dx", QOP_D3_MDWF_DDW_operator_conjugated);
}

/* the standard clover solver */
static int
q_DW_std_solver(lua_State *L,
                struct QOP_D3_MDWF_Fermion *solution,
                int *out_iters,
                double *out_epsilon,
                const struct QOP_D3_MDWF_Fermion *rhs,
                int log_level)
{
    mMDWF *c = qlua_checkMDWF(L, lua_upvalueindex(2), NULL, 1);
    double eps = luaL_checknumber(L, lua_upvalueindex(3));
    int max_iters = luaL_checkint(L, lua_upvalueindex(4));

    return QOP_D3_MDWF_DDW_CG(solution, out_iters, out_epsilon,
                           c->params, rhs, c->gauge, rhs, max_iters, eps,
                           log_level);
}

static MDWFSolver std_solver = { q_DW_std_solver, "CG" };

static int
q_DW_make_solver(lua_State *L)
{
    qlua_checkMDWF(L, 1, NULL, 1);     /* mClover */
    luaL_checknumber(L, 2);            /* double epsilon */
    (void)luaL_checkint(L, 3);         /* int max_iter */

    lua_pushlightuserdata(L, &std_solver);  /* cl[1]: solver */
    lua_pushvalue(L, 1);                    /* cl[2]: mClover */
    lua_pushvalue(L, 2);                    /* cl[3]: epsilon */
    lua_pushvalue(L, 3);                    /* cl[4]: max_iters */
    lua_pushcclosure(L, q_dirac_solver, 4);

    return 1;
}

/* the mixed clover solver */
static int
q_DW_mixed_solver(lua_State *L,
                  struct QOP_D3_MDWF_Fermion *solution,
                  int *out_iters,
                  double *out_epsilon,
                  const struct QOP_D3_MDWF_Fermion *rhs,
                  int log_level)
{
    mMDWF *c = qlua_checkMDWF(L, lua_upvalueindex(2), NULL, 1);
    double f_eps    = luaL_checknumber(L, lua_upvalueindex(3));
    int inner_iters = luaL_checkint(L, lua_upvalueindex(4));
    double eps      = luaL_checknumber(L, lua_upvalueindex(5));
    int max_iters   = luaL_checkint(L, lua_upvalueindex(6));

    return QOP_MDWF_mixed_DDW_CG(solution, out_iters, out_epsilon,
                                 c->params, rhs, c->gauge, rhs,
                                 inner_iters, f_eps,
                                 max_iters, eps,
                                 log_level);
}

static MDWFSolver mixed_solver = { q_DW_mixed_solver, "mixedCG" };

static int
q_DW_debugmesilly(lua_State *L)
{
    /* stack: 1:mMDWF, : debugmesilly ( 2:"op_name", 3:{Fermion} ) */
    char *err_str = NULL;
    mMDWF *c    = qlua_checkMDWF(L, 1, NULL, 1); 
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    const char *op_name    = luaL_checkstring(L, 2);

//    fermion x   = ...
//    QLA_D3_DiracFermion **x = qlua_malloc(L, c->Ls * sizeof (QLA_D3_DiracFermion *));
    mLatDirFerm3 **qlua_x = qlua_malloc(L, c->Ls * sizeof (mLatDirFerm3 *));
    mLatDirFerm3 **qlua_y = qlua_malloc(L, c->Ls * sizeof (mLatDirFerm3 *));
    struct QOP_D3_MDWF_Fermion *c_x;
    struct QOP_D3_MDWF_Fermion *c_y;
    QLA_D3_DiracFermion **e_x = qlua_malloc(L, c->Ls * sizeof (QLA_D3_DiracFermion *));
    QLA_D3_DiracFermion **e_y = qlua_malloc(L, c->Ls * sizeof (QLA_D3_DiracFermion *));
    DW_5_env env;
    //    int status;
    int i;
    CALL_QDP(L);

    for (i = 0; i < c->Ls; i++) {
        lua_pushnumber(L, i + 1); /* [sic] lua indexing */
        lua_gettable(L, 3);
        qlua_x[i] = qlua_checkLatDirFerm3(L, -1, S, 3);
        e_x[i]  = QDP_D3_expose_D(qlua_x[i]->ptr);
        lua_pop(L, 1);
    }

    env.lat = S->lat;
    env.f = e_x;

    if (QOP_D3_MDWF_import_fermion(&c_x, c->state, DW_5_reader, &env)) {
        err_str = "MDWF_import_fermion() failed";
        goto err_end;
    }
    for (i = 0; i < c->Ls; i++)
        QDP_D3_reset_D(qlua_x[i]->ptr);

    if (QOP_D3_MDWF_allocate_fermion(&c_y, c->state)) {
        err_str = "MDWF_create_fermion() failed";
        goto err_end;
    }

    QOP_D3_MDWF_debugmesilly(c_y, c->params, c->gauge, op_name, c_x);

    QOP_D3_MDWF_free_fermion(&c_x);

    lua_createtable(L, c->Ls, 0);
    for (i = 0; i < c->Ls; i++) {
        qlua_y[i] = qlua_newLatDirFerm3(L, Sidx, 3);
        e_y[i] = QDP_D3_expose_D(qlua_y[i]->ptr);
        lua_rawseti(L, -2, i + 1);
    }

    env.lat = S->lat;
    env.f = e_y;
    if (QOP_D3_MDWF_export_fermion(DW_5_writer, &env, c_y)) {
        err_str = "MDWF_export_fermion() failed";
        goto err_end;
    }
    for (i = 0; i < c->Ls; i++) 
        QDP_D3_reset_D(qlua_y[i]->ptr);
    

    QOP_D3_MDWF_free_fermion(&c_y);

    qlua_free(L, qlua_x);
    qlua_free(L, qlua_y);
    qlua_free(L, e_x);
    qlua_free(L, e_y);
    return 1;

err_end:
    qlua_free(L, qlua_x);
    qlua_free(L, qlua_y);
    qlua_free(L, e_x);
    qlua_free(L, e_y);
    return luaL_error(L, err_str);
}

static int
q_DW_make_mixed_solver(lua_State *L)
{
    qlua_checkMDWF(L, 1, NULL, 1);     /* mClover */
    luaL_checknumber(L, 2);            /* double inner_epsilon */
    (void)luaL_checkint(L, 3);         /* int inner_iter */
    luaL_checknumber(L, 4);            /* double epsilon */
    (void)luaL_checkint(L, 5);         /* int max_iter */

    lua_pushlightuserdata(L, &mixed_solver);  /* cl[1]: solver */
    lua_pushvalue(L, 1);                      /* cl[2]: mMDWF */
    lua_pushvalue(L, 2);                      /* cl[3]: inner_epsilon */
    lua_pushvalue(L, 3);                      /* cl[4]: inner_iters */
    lua_pushvalue(L, 4);                      /* cl[5]: epsilon */
    lua_pushvalue(L, 5);                      /* cl[6]: max_iters */
    lua_pushcclosure(L, q_dirac_solver, 6);

    return 1;
}

typedef struct {
    QDP_Lattice *lat;
    int lattice[QOP_MDWF_DIM];
    QLA_D_Complex bf[QOP_MDWF_DIM];
    QLA_D3_ColorMatrix *uf[QOP_MDWF_DIM];
} QCArgs;

static void
q_DW_u_reader(double *v_re, double *v_im,
              int d, const int p[], int a, int b, void *env)
{
    QLA_D_Complex z;
    QCArgs *args = env;
    int i = QDP_index_L(args->lat, p);

    if (p[d] == (args->lattice[d] - 1)) {
        QLA_c_eq_c_times_c(z, args->bf[d], QLA_elem_M(args->uf[d][i], a, b));
    } else {
        QLA_c_eq_c(z, QLA_elem_M(args->uf[d][i], a, b));
    }

    *v_re = QLA_real(z);
    *v_im = QLA_imag(z);
}

/* MDWF constructors helper:
 *    Lua(U[4], bc[4], Ls, ..., [opt_idx_:{}]) => mMDWF
 *    opt_idx_ > 0
 */
static mMDWF *
q_mdwf(lua_State *L, int opt_idx_)
{
    int i;
    int Ls = luaL_checkint(L, 3);
    
    enum QOP_MDWF_eopc_type mdwf_pctype = QOP_MDWF_PC_DEFAULT;
    int mdwf_parity = 0;
    if (0 < opt_idx_ && qlua_checkopt_paramtable(L, opt_idx_)) {
        mdwf_parity = qlua_tabkey_intopt(L, opt_idx_, "parity", 0);
        const char *pctype_str = qlua_tabkey_stringopt(L, opt_idx_, "pctype", "eopc2");
        if (!strcmp("eopc2", pctype_str)) mdwf_pctype = QOP_MDWF_EOPC2;
        else if (!strcmp("eopc2p", pctype_str)) mdwf_pctype = QOP_MDWF_EOPC2PRIME;
        else luaL_error(L, "unknown MDWF preconditioner '%s'", pctype_str);
    }

    luaL_checktype(L, 1, LUA_TTABLE);
    lua_pushnumber(L, 1);
    lua_gettable(L, 1);
    qlua_checkLatColMat3(L, -1, NULL, 3);
    mLattice *S = qlua_ObjLattice(L, -1);
    int Sidx = lua_gettop(L);
    mMDWF *c = qlua_newMDWF(L, Sidx);

    if (Ls < 1)
        luaL_error(L, "Bad value of Ls");
    if (S->rank != QOP_MDWF_DIM)
        luaL_error(L, "MDWF is not implemented for #L=%d", S->rank);
    if (QDP_Ns != QOP_MDWF_FERMION_DIM)
        luaL_error(L, "MDWF does not support Ns=%d", QDP_Ns);

    c->Ls = Ls;
    QCArgs args;
    luaL_checktype(L, 2, LUA_TTABLE);
    for (i = 0; i < QOP_MDWF_DIM; i++) {
        lua_pushnumber(L, i + 1);
        lua_gettable(L, 2);
        switch (qlua_qtype(L, -1)) {
        case qReal:
            QLA_c_eq_r_plus_ir(args.bf[i], lua_tonumber(L, -1), 0);
            break;
        case qComplex:
            QLA_c_eq_c(args.bf[i], *qlua_checkComplex(L, -1));
            break;
        default:
            luaL_error(L, "bad MDWF boundary condition type");
        }
        lua_pop(L, 1);
    }

    QDP_D3_ColorMatrix *UF[QOP_MDWF_DIM];

    luaL_checktype(L, 1, LUA_TTABLE);
    CALL_QDP(L);

    /* extract U from the arguments */
    for (i = 0; i < QOP_MDWF_DIM; i++) {
        lua_pushnumber(L, i + 1); /* [sic] lua indexing */
        lua_gettable(L, 1);
                UF[i] = qlua_checkLatColMat3(L, -1, S, 3)->ptr;
        args.uf[i] = QDP_D3_expose_M(UF[i]);
        lua_pop(L, 1);
    }
    args.lat = S->lat;
    QDP_latsize_L(S->lat, args.lattice);

    struct QOP_MDWF_Config cc;
    cc.self = S->node;
    cc.master_p = QMP_is_primary_node();
    cc.rank = S->rank;
    cc.lat = S->dim;
    cc.ls = Ls;
    cc.net = S->net;
    cc.neighbor_up = S->neighbor_up;
    cc.neighbor_down = S->neighbor_down;
    cc.sublattice = qlua_sublattice;
    cc.env = S;
    cc.parity = mdwf_parity;
    cc.eopc_type = mdwf_pctype;
    if (QOP_MDWF_init(&c->state, &cc))
        luaL_error(L, "MDWF_init() failed");

    if (QOP_D3_MDWF_import_gauge(&c->gauge, c->state, q_DW_u_reader, &args))
        luaL_error(L, "MDWF_import_gauge() failed");

    for (i = 0; i < QOP_MDWF_DIM; i++) {
        QDP_D3_reset_M(UF[i]);
    }
    return c;
}



static void
q_mdwf_generic_check_complex_or_number(double *r, double *r_im, 
        lua_State *L, int pos)
{
    if (lua_type(L, pos) == LUA_TNUMBER) {
        *r = luaL_checknumber(L, pos);
        *r_im = 0.;
    } else {
        QLA_D_Complex *z = qlua_checkComplex(L, pos);
        *r = QLA_real(*z);
        *r_im = QLA_imag(*z);
    }
}
/* parse coefficients (b5 or c5) that can be a (table of) real or 
 * complex numbers
 * `r' an `r_im' must be pre-allocated to hold at least `n' numbers each 
 * the function either returns successfully or crashes in a Lua error
 */
static void
q_mdwf_generic_check_coeff_table(int n, double *r, double *r_im, 
        lua_State *L, int pos)
{
    int i;
    double a, a_im;
    if (lua_type(L, pos) == LUA_TTABLE) {
        for (i = 0 ; i < n ; i++) {
            lua_pushnumber(L, i + 1); /* [sic] lua indexing */
            lua_gettable(L, pos);
            q_mdwf_generic_check_complex_or_number(r + i, r_im + i, L, -1);
            lua_pop(L, 1);
        }
    } else {
        q_mdwf_generic_check_complex_or_number(&a, &a_im, L, pos);
        for (i = 0 ; i < n ; i++) {
            r[i] = a;
            r_im[i] = a_im;
        }
    }
}

/*
 * qcd.MDWF.generic(U[4],          -- [1] Gauge, gauge field
 *                  ubc[4],        -- [2] double[4], gauge boundary conditions
 *                  Ls,            -- [3] int, flavor dimension size
 QOP_MDWF_*                  M5,            -- [4] double
 *                  mf,            -- [5] double
 *                  b5[Ls],        -- [6] double[Ls]
 *                  c5[Ls])        -- [7] double[Ls]
 */
static int
q_mdwf_generic(lua_State *L)
{
    mMDWF *M = q_mdwf(L, 8);
    int Ls = luaL_checkint(L, 3);
    double M5 = luaL_checknumber(L, 4);
    double mf = luaL_checknumber(L, 5);
    double *b5 = qlua_malloc(L, Ls * sizeof (double)); /* [6] */
    double *c5 = qlua_malloc(L, Ls * sizeof (double)); /* [7] */
    double *b5_im = qlua_malloc(L, Ls * sizeof (double)); /* [6] */
    double *c5_im = qlua_malloc(L, Ls * sizeof (double)); /* [7] */
    int i;
    int have_cplx, st;
    /* TODO support case with b5,c5=double/complex/table */
    q_mdwf_generic_check_coeff_table(Ls, b5, b5_im, L, 6);
    q_mdwf_generic_check_coeff_table(Ls, c5, c5_im, L, 7);
    
    for (i = 0, have_cplx = 0; i < Ls; i++) {
        if (0. != b5_im[i] || 0. != c5_im[i])
            have_cplx = 1;
        /*printf("[%3d]\t(%+13.8e+j*%+13.8e)\t(%+13.8e+j*%+13.8e)\n", 
                i, b5[i], b5_im[i], c5[i], c5_im[i]);*/
    }

    M->name = "generic";
    M->type = DW_generic;
    if (have_cplx)
        st = QOP_MDWF_set_complex(&M->params, M->state, b5, b5_im, c5, c5_im, -M5, mf);

    else
        st = QOP_MDWF_set_generic(&M->params, M->state, b5, c5, -M5, mf);
    if (st) {
        qlua_free(L, b5);
        qlua_free(L, c5);
        qlua_free(L, b5_im);
        qlua_free(L, c5_im);
        return luaL_error(L, "Not enough space");
    }

    qlua_free(L, b5);
    qlua_free(L, c5);
    qlua_free(L, b5_im);
    qlua_free(L, c5_im);
    return 1;
}

/*
 * qcd.MDWF.Moebius(U[4],          -- [1] Gauge, gauge field
 *                  ubc[4],        -- [2] double[4], gauge boundary conditions
 *                  Ls,            -- [3] int, flavor dimension size
 *                  M5,            -- [4] double
 *                  mf,            -- [5] double
 *                  b5[Ls],        -- [6] double[Ls]
 *                  kappa)         -- [7] double
 */
static int
q_mdwf_Moebius(lua_State *L)
{
    mMDWF *M = q_mdwf(L, 8);
    int Ls = luaL_checkint(L, 3);
    double M5 = luaL_checknumber(L, 4);
    double mf = luaL_checknumber(L, 5);
    double *b5 = qlua_malloc(L, Ls * sizeof (double)); /* [6] */
    double kappa = luaL_checknumber(L, 7);
    int i;

    for (i = 0; i < Ls; i++) {
        lua_pushnumber(L, i + 1); /* [sic] lua indexing */
        lua_gettable(L, 6);
        b5[i] = luaL_checknumber(L, -1);
        lua_pop(L, 1);
    }

    M->name = "Moebius";
    M->type = DW_Moebius;
    if (QOP_MDWF_set_Moebius(&M->params, M->state, b5, kappa, -M5, mf)) {
                qlua_free(L, b5);
        return luaL_error(L, "Not enough space");
        }
        qlua_free(L, b5);
    return 1;
}

/*
 * qcd.MDWF.Shamir(U[4],           -- [1] Gauge, gauge field
 *                  ubc[4],        -- [2] double[4], gauge boundary conditions
 *                  Ls,            -- [3] int, flavor dimension size
 *                  M5,            -- [4] double
 *                  mf,            -- [5] double
 *                  a5)            -- [6] double
 */
static int
q_mdwf_Shamir(lua_State *L)
{
    mMDWF *M = q_mdwf(L, 7);
    double M5 = luaL_checknumber(L, 4);
    double mf = luaL_checknumber(L, 5);
    double a5 = luaL_checknumber(L, 6);

    M->name = "Shamir";
    M->type = DW_Shamir;
    if (QOP_MDWF_set_Shamir(&M->params, M->state, a5, -M5, mf))
        return luaL_error(L, "Not enough space");

    return 1;
}

/*
 * qcd.MDWF.Borichi(U[4],          -- [1] Gauge, gauge field
 *                  ubc[4],        -- [2] double[4], gauge boundary conditions
 *                  Ls,            -- [3] int, flavor dimension size
 *                  M5,            -- [4] double
 *                  mf,            -- [5] double
 *                  a5)            -- [6] double
 */
static int
q_mdwf_Borichi(lua_State *L)
{
    mMDWF *M = q_mdwf(L, 7);
    double M5 = luaL_checknumber(L, 4);
    double mf = luaL_checknumber(L, 5);
    double a5 = luaL_checknumber(L, 6);

    M->name = "Borichi";
    M->type = DW_Borichi;
    if (QOP_MDWF_set_Borichi(&M->params, M->state, a5, -M5, mf))
        return luaL_error(L, "Not enough space");

    return 1;
}

/*
 * qcd.MDWF.Chiu(U[4],          -- [1] Gauge, gauge field
 *               ubc[4],        -- [2] double[4], gauge boundary conditions
 *               Ls,            -- [3] int, flavor dimension size
 *               M5,            -- [4] double
 *               mf,            -- [5] double
 *               a5[Ls])        -- [6] double[Ls]
 */
static int
q_mdwf_Chiu(lua_State *L)
{
    mMDWF *M = q_mdwf(L, 7);
    int Ls = luaL_checkint(L, 3);
    double M5 = luaL_checknumber(L, 4);
    double mf = luaL_checknumber(L, 5);
    double *a5 = qlua_malloc(L, Ls * sizeof (double)); /* [6] */
    int i;

    for (i = 0; i < Ls; i++) {
        lua_pushnumber(L, i + 1); /* [sic] lua indexing */
        lua_gettable(L, 6);
        a5[i] = luaL_checknumber(L, -1);
        lua_pop(L, 1);
    }

    M->name = "Chiu";
    M->type = DW_Chiu;
    if (QOP_MDWF_set_Chiu(&M->params, M->state, a5, -M5, mf)) {
                qlua_free(L, a5);
        return luaL_error(L, "Not enough space");
        }

        qlua_free(L, a5);
    return 1;
}


static struct luaL_Reg mtMDWF[] = {
    { "__tostring",                 q_DW_fmt                        },
    { "__gc",                       q_DW_gc                         },
    { "close",                      q_DW_close                      },
    { "D",                          q_DW_D                          },
    { "Dx",                         q_DW_Dx                         },
    { "solver",                     q_DW_make_solver                },
    { "mixed_solver",               q_DW_make_mixed_solver          },
    { "debugmesilly",               q_DW_debugmesilly               },
    { "eig_deflator",               q_DW_make_deflator              },
#ifdef HAS_ARPACK
    { "eig_deflator_lanczos",       q_DW_make_deflator_lanczos      },
#endif /* HAS_ARPACK */
    { NULL,                         NULL                            }
};

static mMDWF *
qlua_newMDWF(lua_State *L, int Sidx)
{
    mMDWF *c = lua_newuserdata(L, sizeof (mMDWF));

    c->state = 0;
    c->gauge = 0;
    c->params = 0;
    qlua_createLatticeTable(L, Sidx, mtMDWF, qMDWF, MDWFName);
    lua_setmetatable(L, -2);

    return c;
}

static mMDWF *
qlua_checkMDWF(lua_State *L, int idx, mLattice *S, int live)
{
    mMDWF *c = qlua_checkLatticeType(L, idx, qMDWF, MDWFName);

    if (S) {
        mLattice *S1 = qlua_ObjLattice(L, idx);
        if (S1->id != S->id)
            luaL_error(L, "%s on a wrong lattice", MDWFName);
        lua_pop(L, 1);
    }

    if (live && (c->state == 0 || c->gauge == 0))
        luaL_error(L, "using closed qcd.MDWF");

    return c;
}

static struct luaL_Reg fMDWF[] = {
    { "generic",       q_mdwf_generic },
    { "Moebius",       q_mdwf_Moebius },
    { "Shamir",        q_mdwf_Shamir  },
    { "Borichi",       q_mdwf_Borichi },
    { "Chiu",          q_mdwf_Chiu    },
    { NULL,            NULL           }
};

int
init_mdwf(lua_State *L)
{
    lua_getglobal(L, qcdlib);
    lua_newtable(L);
    luaL_register(L, NULL, fMDWF);
    lua_setfield(L, -2, mdwf_name);
    lua_pop(L, 1);
    return 0;
}
#else /* USE_Nc3 */
int
init_mdwf(lua_State *L)
{
    return 0;
}
#endif /* USE_Nc3 */

void
fini_mdwf(void)
{
}
