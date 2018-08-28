/// NDBD.hh ///
/// This class has been designed to handle Neutrinoless Double-Beta Decay (NDBD).
/// Its main purpose is to be used in imsrg_util.hh/cc for the M0nu Operator TBMEs.
/// In particlar: M0nuGTF_adpt_Op, M0nuT_adpt_Op, M0nu_PrintIntegrand (so far).
/// NOTE: to change the NDBD parameters, change the "#define" macros below, *and* those in Constants.hh
/// by: Charlie Payne (CP), 2016/2017/2018
/// research supervisor: Jason Holt, postdoc: Ragnar Stroberg
/// The work has been done primarily for my MSc RA at UBC/TRIUMF.
/// CP = cgpayne@triumf.ca
/// JH = jholt@triumf.ca
/// RS = sstroberg@triumf.ca


#ifndef NDBD_hh
#define NDBD_hh 1

#include "Constants.hh"
#include "ModelSpace.hh"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/factorials.hpp>


/*
#ifndef PI
  #define PI 3.14159265359
#endif
*/

// PSH style NDBD parameters
#define NDBD_MPRO   938.2720814 // [MeV]
#define NDBD_MPION  139.57018   // [MeV]
#define NDBD_MAGMOM 4.706       // [\mu_N]
#define NDBD_G0V    1.0         // [unitless]
#define NDBD_G0A    1.27        // [unitless]
#define NDBD_CUTV   850.0       // [MeV]
#define NDBD_CUTA   1086.0      // [MeV]
#define NDBD_R0     1.2         // [fm]
//
/* MH style NDBD parameters
#define NDBD_MPRO   938.2720814 // [MeV]
#define NDBD_MPION  139.57018   // [MeV]
#define NDBD_MAGMOM 3.7         // [\mu_N]
#define NDBD_G0V    1.0         // [unitless]
#define NDBD_G0A    1.25        // [unitless]
#define NDBD_CUTV   850.0       // [MeV]
#define NDBD_CUTA   1086.0      // [MeV]
#define NDBD_R0     1.2         // [fm]
*/
/* JE style NDBD parameters (should also change M_NUCLEON in Constants.hh)
#define NDBD_MPRO   938.27231 // [MeV]
#define NDBD_MPION  139.57    // [MeV]
#define NDBD_MAGMOM 3.706     // [\mu_N]
#define NDBD_G0V    1.0       // [unitless]
#define NDBD_G0A    1.27      // [unitless]
#define NDBD_CUTV   850.0     // [MeV]
#define NDBD_CUTA   1086.0    // [MeV]
#define NDBD_R0     1.2       // [fm]
*/

#define NDBD_R_LIMIT 10000  // max number of subintervals for integration wrt dr
#define NDBD_R_ABS   1e-7   // absolute error tolerance for integration wrt dr:  "converged" if abs[x(i) - x(i-1)] < epsabs...
#define NDBD_R_REL   1e-5   // relative " " " " " dr:                            ...or "converged" if abs[1 - x(i)/x(i-1)] < epsrel
#define NDBD_R_RB    1e+4   // when using QAG, we integrate over [RA,RB] ~= [0,Inf), so make RB [fm] reasonably large (reasonable value = ????)
#define NDBD_Q_LIMIT 10000  // max number of subintervals for integration wrt dq
#define NDBD_Q_ABS   1e-7   // absolute error tolerance for integration wrt dq:  "converged" if abs[x(i) - x(i-1)] < epsabs...
#define NDBD_Q_REL   1e-4   // relative " " " " " dq:                            ...or "converged" if abs[1 - x(i)/x(i-1)] < epsrel
#define NDBD_Q_QB    2.5e+3 // when using QAG, we integrate over [QA,QB] ~= [0,Inf), so make QB [MeV] reasonably large (reasonable value = ????)

#define NDBD_GT  "GT" // Gamow-Teller
#define NDBD_F   "F"  // Fermi
#define NDBD_T   "T"  // Tensor
#define NDBD_PSH "PSH" // Payne-Storberg-Holt method
//#define NDBD_JE  "JE"  // Jonh Engel method (intended to code this at some point for benchmarking... never did...)

#define NDBD_SRC_AV18  "AV18"           // AKA: Argonne V18 potential
#define NDBD_SRC_CDB   "CD-Bonn"        // AKA: charge-dependent Bonn potential
#define NDBD_SRC_MS    "Miller-Spencer" // see the paper: Ann.Phys.100.562(1976)
#define NDBD_SRC_DEBUG "debug"
#define NDBD_SRC_ON    1
#define NDBD_SRC_OFF   0

#define NDBD_INT_ANAL  1 // use an analytic integration formula
#define NDBD_INT_QAGIU 2 // use GSL qagiu integration
#define NDBD_INT_QAG   3 // use GSL qag integration


// class non-specific function prototypes
 inline int cpPhase(int i) {return (i%2)==0 ? 1 : -1;}; // calculates the phase of an integer
 string cpFormat(double val); // my personal output format for a number
 int decimalgen(int a, int b, int c, int d, int maxb, int maxc, int maxd); // converts (a,b,c,d) in base (maxa+1,maxb+1,maxc+1,maxd+1) to an ordered decimal integer
 void cpGSLerror(string errstr, int status, size_t limit, double epsabs, double epsrel, double abserr, double result, int n, int l, int np, int lp); // my personal GSL integration error handling
 double rbmeNONE(double hw, int rho, double x, int n, int l, int np, int lp); // analytical version of the relative Bessel's matrix elements (RBMEs), with *no* SRCs
 double frRBME(double r, void *params); // gsl_function handle for integration wrt dr (RBME)
 double rbmeNONEqagiu(double hw, int rho, double x, int n, int l, int np, int lp); // gsl_integration_qagiu version of the RBMEs, with *no* SRCs
 double rbmeSRC(double hw, int rho, double x, int n, int l, int np, int lp, double a, double b, double c); // analytical version of the RBMEs, with SRCs
 double frRBMEjt(double r, void *params); // gsl_function handle for integration wrt dr (RBME), including Jastrow-type correlation function for the SRCs
 double rbmeSRCqagiu(double hw, int rho, double x, int n, int l, int np, int lp, double a, double b, double c); // gsl_integration_qagiu version of the RBMEs, with SRCs
 //double rbmeNONErjrs(double hw, int rho, double x, int n, int l, int np, int lp); // right-justified Riemann's sum integration version (for JE) of the RBMEs, with *no* SRCs
 double hF(double q, double naughtV, double lambdaV); // Fermi (F) type form factor
 double hGT(double q, double Mpro, double Mpion, double MagMom, double naughtV, double naughtA, double lambdaV, double lambdaA); // Gamow-Teller (GT) type form factor
 double hT(double q, double Mpro, double Mpion, double MagMom, double naughtV, double naughtA, double lambdaV, double lambdaA); // Tensor (T) type form factor

 double rbmeGrab(int srchit, int drswitch, double hw, int rho, double q, int n, int l, int np, int lp, double a, double b, double c); // will be called by fqF/GT/T below
 double fqF(double q, void *params); // gsl_function handle for integration wrt dq (F)
 double fqGT(double q, void *params); // " " " " " " (GT)
 double fqT(double q, void *params); // " " " " " " (T)


// class prototype
 class NDBD
 {
   public:
   // members
   double hw; // oscillator basis frequency [MeV]
   int e2max; // 2*emax from modelspace [unitless]
   double mpro; // the proton mass [MeV] for the g-factors
   double mpion; // the pion mass [MeV] for the g-factors
   double magmom; // the difference between the (anomolous?) magnetic moment of a proton and neutron [\mu_N]
   double g0V; // the vector g-factor at zero momentum [unitless]
   double g0A; // the axial-vector g-factor at zero momentum [unitless]
   double cutoffV; // the vector finite-size cutoff regularization parameter [MeV]
   double cutoffA; // the axial-vector finite-size cutoff regularization parameter [MeV]
   double r0; // from the equation R = r0*A^(1/3) [fm]
   string type; // F, GT, or T
   double Ebar; // the closure energy [MeV]
   double srcparams[3]; // the SRC parameters
   int srchit; // 0 = none, 1 = o/w
   string method; // "PSH" (formerly known as "CP") or "JE", etc
   // constructors/destructors
   NDBD(double HW, int etwomax, string decaytype, double Ec, string src, string themethod);
   ~NDBD();
   // member functions
   void SetPSH(int n, int l, int np, int lp, int &dqs, int &drs); // set two-tiered integral switches (on/off)
   void SetIparams(int n, int l, int np, int lp, int drs, double Iparams[18]); // set gsl_function handle parameters, for integration wrt dq
   double dqdrPSH(int n, int l, int np, int lp); // computes the relevant integral wrt dq and dr of the neutrino potential
   //inline long int IntHash(int n, int l, int np) {return ((n + np)*(n + np + 1) + 2*(n + np + l))*((n + np)*(n + np + 1) + 2*(n + np + l + 1))/8 + n*np + l;}; // via modified pairing function, but...
   //void IntUnHash(long int key, int &n, int &l, int &np); // ...can't figure out how to do this
   uint64_t IntHash(uint64_t n, uint64_t l, uint64_t np, uint64_t lp); // Hashtag for IntList cacheing
   void IntUnHash(uint64_t key, uint64_t &n, uint64_t &l, uint64_t &np, uint64_t &lp); // the inverse of the above IntHash
   void PreCalcIntegrals(); // caches the relevant integrals wrt dq and dr
   double GetIntegral(int n, int l, int np, int lp); // get an integral from the IntList cache or calculate it (parallelization dependent)
   uint64_t T6jHash(uint64_t l1, uint64_t L1, uint64_t R, uint64_t L2, uint64_t l2); // Hashtag for T6jList cacheing
   void T6jUnHash(uint64_t key, uint64_t &l1, uint64_t &L1, uint64_t &R, uint64_t &L2, uint64_t &l2); // the inverse of the above T6jHash
   void PreCalcT6j(); // caches the relevant 6j-symbols for the Tensor M0nu component
   double GetT6j(int l1, int L1, int R, int L2, int l2); // get a 6j from the T6jList cache or calculate it (parallelization dependent), will use AngMom::SixJ(l1,L1,R,L2,l2,2) - notice the 2
   // static members
   IMSRGProfiler profiler;
   //static unordered_map<long int,double> IntList; // "Int" = "Integral", cache set by PreCalcIntegrals above
   static unordered_map<uint64_t,double> IntList; // "Int" = "Integral", cache set by PreCalcIntegrals above
   static unordered_map<uint64_t,double> T6jList; // "T6j" = "Tensor 6j-symbol", cache set by PreCalcT6j above
 };


// M0nu Operator specific function prototypes
 inline double cpNorm(int i, int j) {return i==j ? ISQRT2 : 1.0;}; // for anti-symmetrization normalization of TBMEs
 string cpBarcode(string Type, double hw, double emax); // returns a string of five random letters and a time stamp, to uniquely identify an M0nu run
 void M0nuHeader(NDBD zvbb, string dirname, string reduced, int Anuc, string src, double Rnuc, double prefact, string barcode); // creates a parameter header file for M0nu
 inline int PairFN(int a, int b) {return ((a + b)*(a + b + 1)/2) + b;}; // a Cantor pairing function used for cacheing CG in M0nuT


#endif
