/* Created by Language version: 6.2.0 */
/* NOT VECTORIZED */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib.h"
#undef PI
 
#include "md1redef.h"
#include "section.h"
#include "md2redef.h"

#if METHOD3
extern int _method3;
#endif

#undef exp
#define exp hoc_Exp
extern double hoc_Exp();
 
#define _threadargscomma_ /**/
#define _threadargs_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gbar _p[0]
#define ar2 _p[1]
#define m _p[2]
#define h _p[3]
#define s _p[4]
#define ena _p[5]
#define ina _p[6]
#define thegna _p[7]
#define Dm _p[8]
#define Dh _p[9]
#define Ds _p[10]
#define _g _p[11]
#define _ion_ena	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static int _hoc_alps();
 static int _hoc_alpv();
 static int _hoc_bets();
 static int _hoc_trap0();
 static int _hoc_trates();
 static int _mechtype;
extern int nrn_get_mechtype();
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range();
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 ret(1.);
}
 /* connect user functions to hoc names */
 static IntFunc hoc_intfunc[] = {
 "setdata_Naf_b1", _hoc_setdata,
 "alps_Naf_b1", _hoc_alps,
 "alpv_Naf_b1", _hoc_alpv,
 "bets_Naf_b1", _hoc_bets,
 "trap0_Naf_b1", _hoc_trap0,
 "trates_Naf_b1", _hoc_trates,
 0, 0
};
#define alps alps_Naf_b1
#define alpv alpv_Naf_b1
#define bets bets_Naf_b1
#define trap0 trap0_Naf_b1
 extern double alps();
 extern double alpv();
 extern double bets();
 extern double trap0();
 /* declare global and static user variables */
#define Rd Rd_Naf_b1
 double Rd = 0.03;
#define Rg Rg_Naf_b1
 double Rg = 0.01;
#define Rb Rb_Naf_b1
 double Rb = 0.124;
#define Ra Ra_Naf_b1
 double Ra = 0.4;
#define a0s a0s_Naf_b1
 double a0s = 0.0003;
#define gms gms_Naf_b1
 double gms = 0.2;
#define hmin hmin_Naf_b1
 double hmin = 0.5;
#define htau htau_Naf_b1
 double htau = 0;
#define hinf hinf_Naf_b1
 double hinf = 0;
#define mmin mmin_Naf_b1
 double mmin = 0.02;
#define mtau mtau_Naf_b1
 double mtau = 0;
#define minf minf_Naf_b1
 double minf = 0;
#define qq qq_Naf_b1
 double qq = 10;
#define q10 q10_Naf_b1
 double q10 = 2;
#define qg qg_Naf_b1
 double qg = 1.5;
#define qd qd_Naf_b1
 double qd = 1.5;
#define qa qa_Naf_b1
 double qa = 7.2;
#define qinf qinf_Naf_b1
 double qinf = 4;
#define smax smax_Naf_b1
 double smax = 10;
#define sinf sinf_Naf_b1
 double sinf = 0;
#define tq tq_Naf_b1
 double tq = -55;
#define thi2 thi2_Naf_b1
 double thi2 = -45;
#define thi1 thi1_Naf_b1
 double thi1 = -45;
#define tha tha_Naf_b1
 double tha = -30;
#define thinf thinf_Naf_b1
 double thinf = -50;
#define taus taus_Naf_b1
 double taus = 0;
#define vvs vvs_Naf_b1
 double vvs = 2;
#define vvh vvh_Naf_b1
 double vvh = -58;
#define vhalfs vhalfs_Naf_b1
 double vhalfs = -60;
#define zetas zetas_Naf_b1
 double zetas = 12;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "tha_Naf_b1", "mV",
 "qa_Naf_b1", "mV",
 "Ra_Naf_b1", "/ms",
 "Rb_Naf_b1", "/ms",
 "thi1_Naf_b1", "mV",
 "thi2_Naf_b1", "mV",
 "qd_Naf_b1", "mV",
 "qg_Naf_b1", "mV",
 "Rg_Naf_b1", "/ms",
 "Rd_Naf_b1", "/ms",
 "qq_Naf_b1", "mV",
 "tq_Naf_b1", "mV",
 "thinf_Naf_b1", "mV",
 "qinf_Naf_b1", "mV",
 "vhalfs_Naf_b1", "mV",
 "a0s_Naf_b1", "ms",
 "zetas_Naf_b1", "1",
 "gms_Naf_b1", "1",
 "smax_Naf_b1", "ms",
 "vvh_Naf_b1", "mV",
 "vvs_Naf_b1", "mV",
 "mtau_Naf_b1", "ms",
 "htau_Naf_b1", "ms",
 "sinf_Naf_b1", "ms",
 "taus_Naf_b1", "ms",
 "gbar_Naf_b1", "mho/cm2",
 "ar2_Naf_b1", "1",
 0,0
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 static double s0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "tha_Naf_b1", &tha_Naf_b1,
 "qa_Naf_b1", &qa_Naf_b1,
 "Ra_Naf_b1", &Ra_Naf_b1,
 "Rb_Naf_b1", &Rb_Naf_b1,
 "thi1_Naf_b1", &thi1_Naf_b1,
 "thi2_Naf_b1", &thi2_Naf_b1,
 "qd_Naf_b1", &qd_Naf_b1,
 "qg_Naf_b1", &qg_Naf_b1,
 "mmin_Naf_b1", &mmin_Naf_b1,
 "hmin_Naf_b1", &hmin_Naf_b1,
 "q10_Naf_b1", &q10_Naf_b1,
 "Rg_Naf_b1", &Rg_Naf_b1,
 "Rd_Naf_b1", &Rd_Naf_b1,
 "qq_Naf_b1", &qq_Naf_b1,
 "tq_Naf_b1", &tq_Naf_b1,
 "thinf_Naf_b1", &thinf_Naf_b1,
 "qinf_Naf_b1", &qinf_Naf_b1,
 "vhalfs_Naf_b1", &vhalfs_Naf_b1,
 "a0s_Naf_b1", &a0s_Naf_b1,
 "zetas_Naf_b1", &zetas_Naf_b1,
 "gms_Naf_b1", &gms_Naf_b1,
 "smax_Naf_b1", &smax_Naf_b1,
 "vvh_Naf_b1", &vvh_Naf_b1,
 "vvs_Naf_b1", &vvs_Naf_b1,
 "minf_Naf_b1", &minf_Naf_b1,
 "hinf_Naf_b1", &hinf_Naf_b1,
 "mtau_Naf_b1", &mtau_Naf_b1,
 "htau_Naf_b1", &htau_Naf_b1,
 "sinf_Naf_b1", &sinf_Naf_b1,
 "taus_Naf_b1", &taus_Naf_b1,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(), nrn_init(), nrn_state();
 static void nrn_cur(), nrn_jacob();
 
static int _ode_count(), _ode_map(), _ode_spec(), _ode_matsol();
 
#define _cvode_ieq _ppvar[3]._i
 /* connect range variables in _p that hoc is supposed to know about */
 static char *_mechanism[] = {
 "6.2.0",
"Naf_b1",
 "gbar_Naf_b1",
 "ar2_Naf_b1",
 0,
 0,
 "m_Naf_b1",
 "h_Naf_b1",
 "s_Naf_b1",
 0,
 0};
 static Symbol* _na_sym;
 
static void nrn_alloc(_prop)
	Prop *_prop;
{
	Prop *prop_ion, *need_memb();
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 12, _prop);
 	/*initialize range parameters*/
 	gbar = 0.01;
 	ar2 = 1;
 	_prop->param = _p;
 	_prop->param_size = 12;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 
}
 static _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 _Naf_b1_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("na", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
  hoc_register_dparam_size(_mechtype, 4);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 Naf_b1 /home/ximing/Documents/from_axon/ca1n1-mod/i686/Naf_b1.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double _zmexp , _zhexp , _zsexp ;
static int _reset;
static char *modelname = "Sodium fast channel";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static _modl_cleanup(){ _match_recurse=1;}
static trates();
 
static int _ode_spec1(), _ode_matsol1();
 static int _slist1[3], _dlist1[3];
 static int states();
 
double alpv (  _lv )  
	double _lv ;
 {
   double _lalpv;
 _lalpv = 1.0 / ( 1.0 + exp ( ( _lv - vvh ) / vvs ) ) ;
   
return _lalpv;
 }
 
static int _hoc_alpv() {
  double _r;
   _r =  alpv (  *getarg(1) ) ;
 ret(_r);
}
 
double alps (  _lv )  
	double _lv ;
 {
   double _lalps;
 _lalps = exp ( 1.e-3 * zetas * ( _lv - vhalfs ) * 9.648e4 / ( 8.315 * ( 273.16 + celsius ) ) ) ;
   
return _lalps;
 }
 
static int _hoc_alps() {
  double _r;
   _r =  alps (  *getarg(1) ) ;
 ret(_r);
}
 
double bets (  _lv )  
	double _lv ;
 {
   double _lbets;
 _lbets = exp ( 1.e-3 * zetas * gms * ( _lv - vhalfs ) * 9.648e4 / ( 8.315 * ( 273.16 + celsius ) ) ) ;
   
return _lbets;
 }
 
static int _hoc_bets() {
  double _r;
   _r =  bets (  *getarg(1) ) ;
 ret(_r);
}
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   trates ( _threadargscomma_ v , ar2 ) ;
   Dm = ( minf - m ) / mtau ;
   Dh = ( hinf - h ) / htau ;
   Ds = ( sinf - s ) / taus ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 trates ( _threadargscomma_ v , ar2 ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / mtau )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / htau )) ;
 Ds = Ds  / (1. - dt*( ( ( ( - 1.0 ) ) ) / taus )) ;
}
 /*END CVODE*/
 static int states () {_reset=0;
 {
   trates ( _threadargscomma_ v , ar2 ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / mtau)))*(- ( ( ( minf ) ) / mtau ) / ( ( ( ( - 1.0) ) ) / mtau ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / htau)))*(- ( ( ( hinf ) ) / htau ) / ( ( ( ( - 1.0) ) ) / htau ) - h) ;
    s = s + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / taus)))*(- ( ( ( sinf ) ) / taus ) / ( ( ( ( - 1.0) ) ) / taus ) - s) ;
   }
  return 0;
}
 
static int  trates (  _lvm , _la2 )  
	double _lvm , _la2 ;
 {
   double _la , _lb , _lc , _lqt ;
 _lqt = pow( q10 , ( ( celsius - 24.0 ) / 10.0 ) ) ;
   _la = trap0 ( _threadargscomma_ _lvm , tha , Ra , qa ) ;
   _lb = trap0 ( _threadargscomma_ - _lvm , - tha , Rb , qa ) ;
   mtau = 1.0 / ( _la + _lb ) / _lqt ;
   if ( mtau < mmin ) {
     mtau = mmin ;
     }
   minf = _la / ( _la + _lb ) ;
   _la = trap0 ( _threadargscomma_ _lvm , thi1 , Rd , qd ) ;
   _lb = trap0 ( _threadargscomma_ - _lvm , - thi2 , Rg , qg ) ;
   htau = 1.0 / ( _la + _lb ) / _lqt ;
   if ( htau < hmin ) {
     htau = hmin ;
     }
   hinf = 1.0 / ( 1.0 + exp ( ( _lvm - thinf ) / qinf ) ) ;
   _lc = alpv ( _threadargscomma_ _lvm ) ;
   sinf = _lc + _la2 * ( 1.0 - _lc ) ;
   taus = bets ( _threadargscomma_ _lvm ) / ( a0s * ( 1.0 + alps ( _threadargscomma_ _lvm ) ) ) ;
   if ( taus < smax ) {
     taus = smax ;
     }
    return 0; }
 
static int _hoc_trates() {
  double _r;
   _r = 1.;
 trates (  *getarg(1) , *getarg(2) ) ;
 ret(_r);
}
 
double trap0 (  _lv , _lth , _la , _lq )  
	double _lv , _lth , _la , _lq ;
 {
   double _ltrap0;
 if ( fabs ( _lv - _lth ) > 1e-6 ) {
     _ltrap0 = _la * ( _lv - _lth ) / ( 1.0 - exp ( - ( _lv - _lth ) / _lq ) ) ;
     }
   else {
     _ltrap0 = _la * _lq ;
     }
   
return _ltrap0;
 }
 
static int _hoc_trap0() {
  double _r;
   _r =  trap0 (  *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) ) ;
 ret(_r);
}
 
static int _ode_count(_type) int _type;{ return 3;}
 
static int _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
     _ode_spec1 ();
  }}
 
static int _ode_map(_ieq, _pv, _pvdot, _pp, _ppd, _atol, _type) int _ieq, _type; double** _pv, **_pvdot, *_pp, *_atol; Datum* _ppd; { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 3; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static int _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
 _ode_matsol1 ();
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h = h0;
  m = m0;
  s = s0;
 {
   trates ( _threadargscomma_ v , ar2 ) ;
   m = minf ;
   h = hinf ;
   s = sinf ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ena = _ion_ena;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   thegna = gbar * m * m * m * h * s ;
   ina = thegna * ( v - ena ) ;
   }
 _current += ina;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ena = _ion_ena;
 _g = _nrn_current(_v + .001);
 	{ double _dina;
  _dina = ina;
 _rhs = _nrn_current(_v);
  _ion_dinadv += (_dina - ina)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
 double _break, _save;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _break = t + .5*dt; _save = t;
 v=_v;
{
  ena = _ion_ena;
 { {
 for (; t < _break; t += dt) {
 error =  states();
 if(error){fprintf(stderr,"at line 66 in file Naf_b1.mod:\n        SOLVE states METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 
}}
 t = _save;
 } }}

}

static terminal(){}

static _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
 _slist1[2] = &(s) - _p;  _dlist1[2] = &(Ds) - _p;
_first = 0;
}