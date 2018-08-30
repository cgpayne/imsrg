
#include "NDBD.hh"
#include "imsrg_util.hh"
#include "AngMom.hh"


using namespace imsrg_util;


/// My personal output format for a number
  string cpFormat(double val)
  {
    int E = floor(log10(val));
    double N = val/pow(10,E);
    stringstream Vss;
    Vss<<setprecision(3)<<N<<"e"<<E;
    return Vss.str();
  }


/// Converts (a,b,c,d) in base (maxa+1,maxb+1,maxc+1,maxd+1) to an ordered decimal integer
/// eg: for (maxa,maxb,maxc,maxd) = (1,1,1,1) decimalgen would convert (0,1,1,0) to 6, ie: a binary to decimal converter
/// NOTE: to make comparisons between such decimals, keep (maxa,maxb,maxc,maxd) consistent
/// PS: I tots made this thing up, did some testing and it seemed to work...
///     ... hopefully it's good, but could be a source of weird bugs if anything comes up with missing integrals wrt dq, see GetIntegral(...)
  int decimalgen(int a, int b, int c, int d, int maxb, int maxc, int maxd)
  {
    int coeff1 = maxd + 1;
    int coeff2 = (maxc + 1)*coeff1;
    int coeff3 = (maxb + 1)*coeff2;
    return coeff3*a + coeff2*b + coeff1*c + d; // eg: (0,1,1,0) -> (2*2*2)*0 + (2*2)*1 + (2)*1 + 0 = 6
  }


/// My personal GSL integration error handling
/// NOTE: must put "#pragma omp critical" calls around this function when calling it in the code
/// NOTE: in order to use this, one must set "gsl_set_error_handler_off();" beforehand
  void cpGSLerror(string errstr, int status, size_t limit, double epsabs, double epsrel, double abserr, double result, int n, int l, int np, int lp)
  {
    if (status)
    {
      cout<<"WARNING "<<errstr;
      if (status == GSL_EMAXITER)
      {
        cout<<" the maximum number of subdivisions was exceeded."<<endl;
      }
      else if (status == GSL_EROUND)
      {
        cout<<" cannot reach tolerance because of roundoff error, or roundoff error was detected in the extrapolation table."<<endl;
      }
      else if (status == GSL_ESING)
      {
        cout<<" a non-integrable singularity or other bad integrand behavior was found in the integration interval."<<endl;
      }
      else if (status == GSL_EDIVERGE)
      {
        cout<<" the integral is divergent, or too slowly convergent to be integrated numerically."<<endl;
      }
      else
      {
        cout<<"apparently I'm missing a GSL error status handle..."<<endl;
      }
      cout<<"limit  = "<<limit<<endl;
      cout<<"epsabs = "<<cpFormat(epsabs)<<endl;
      cout<<"epsrel = "<<cpFormat(epsrel)<<endl;
      cout<<"abserr = "<<cpFormat(abserr)<<endl;
      cout<<"result = "<<setprecision(12)<<result<<endl;
      cout<<"n = "<<n<<", l = "<<l<<", np = "<<np<<", lp = "<<lp<<endl;
      double errat = abs(abserr/result);
      if (errat <= epsrel) // this if-else is a last resort, but only seems to rarely trigger for Tensor with high emax for a few n,l,np,lp
      {
        cout<<"moving forward since abs(abserr/result) = "<<cpFormat(errat)<<" <= epsrel..."<<endl;
      }
      else if (abserr <= epsabs)
      {
        cout<<"moving forward since abserr <= epsabs..."<<endl;
      }
      else
      {
        cout<<"ERROR 37728: my GSL integration failed! :("<<endl;
        cout<<"we have a GSL error status (see the warning above)"<<endl;
        cout<<"and abs(abserr/result) = "<<cpFormat(errat)<<" > epsrel..."<<endl;
        cout<<"and abserr > epsabs..."<<endl;
        cout<<"therefore exiting..."<<endl;
        exit(1);
      }
    }
  }


/// This calculates the analytical integral wrt dr of the RBMEs
/// NONE => SRCs are *not* employed here
/// ie) RBME = < n l | j_\rho(\sqrt{2}qr) | n' l' >, via Equation (4.64) of my thesis
/// NOTE: the \sqrt{2} with the q is to stay consistent with r_{rel} = (r_1 - r_2)/SQRT2 in Moshinksy brackets
/// NOTE: rho must be divisible by 2! >:|
  double rbmeNONE(double hw, int rho, double x, int n, int l, int np, int lp)
  {
    double rbme = 0; // initialize the double sum over k, k'
    double kaptestnn = ((l + lp - rho)/2.0); // the tgamma_ratio below will only take non-negative...
    double kaptestint = abs(kaptestnn - floor(kaptestnn)); // ...integers as arguments
    if ((kaptestnn >= 0) and (kaptestint == 0))
    {
      double bosc = 1.0/sqrt(M_NUCLEON*hw); // the oscillator length [MeV] via Equation (2.4) of my thesis
      double yolo = x*x*bosc*bosc/2.0; // commonly reoccuring argument in the RBMEs
      double normy = (1.0/pow(2,rho+1))*sqrt(PI*boost::math::tgamma_ratio(n + 1, n + l + 1.5)*boost::math::tgamma_ratio(np + 1, np + lp + 1.5)); // the normalization factor
      double rholag = rho + 0.5; // for the Laguerre polynomial in the sum
      for (int k=0; k<=n; k++)
      {
        for (int kp=0; kp<=np; kp++)
        {
          int kappa = kaptestnn + k + kp;
          double coeff = cpPhase(k + kp)*gsl_sf_gammainv(k + 1)*gsl_sf_gammainv(kp + 1); // the sum coefficient
          double bci = boost::math::tgamma_ratio(n + l + 1.5, n - k + 1); // this is only part of the generalized binomial coeff...
          double bcj = boost::math::tgamma_ratio(np + lp + 1.5, np - kp + 1); // ...we put the other part within kapfact below
          double bcij = bci*bcj;
          double kapfact = boost::math::tgamma_ratio(kappa + 1, max(l + k, lp + kp) + 1.5)*gsl_sf_gammainv(min(l + k, lp + kp) + 1.5);
          double Lag = gsl_sf_laguerre_n(kappa,rholag,yolo); // from GSL -- L^a_n(x) = gsl_sf_laguerre_n(const int n, const double a, const double x)
          rbme += coeff*bcij*kapfact*Lag; // perform the summation
        }
      }
      rbme *= normy*pow(SQRT2*x*bosc,rho)*exp(-yolo); // see Equation (4.64) in my thesis
    }
    return rbme;
  }


/// gsl_function handle for integration wrt dr (RBME)
/// this will be called by rbmeNONEqagiu below
  double frRBME(double r, void *params)
  {
    double hw = ((double*) params)[0];
    int rho = ((double*) params)[1];
    double x = ((double*) params)[2];
    int n = ((double*) params)[3];
    int l = ((double*) params)[4];
    int np = ((double*) params)[5];
    int lp = ((double*) params)[6];
    x /= HBARC; // put it into [fm^-1] since r is in [fm]
    x *= SQRT2; // this is to stay consistent with r_{rel} = (r_1 - r_2)/SQRT2 in Moshinksy brackets
    return r*r*HO_Radial_psi(n,l,hw,r)*gsl_sf_bessel_jl(rho,x*r)*HO_Radial_psi(np,lp,hw,r);
  }


/// This calculates the gsl_integration_qagiu wrt dr of the RBMEs
/// NONE => SRCs are *not* employed here
  double rbmeNONEqagiu(double hw, int rho, double x, int n, int l, int np, int lp)
  {
    gsl_set_error_handler_off(); // set standard GSL error handling off, so I can do my own
    int statusr; // this will hold the GSL error status
    string errstrr="7828877-1: qagiu dr:"; // a unique GSL error string
    size_t limitr = NDBD_R_LIMIT; // the number of double precision intervals that can be held by the workspace below
    gsl_integration_workspace * workspacer = gsl_integration_workspace_alloc(limitr); // GSL: "One workspace may be used multiple times...", "...workspaces should be allocated on a per-thread basis."
    double frparams[] = {hw, (double)rho, x, (double)n, (double)l, (double)np, (double)lp};
    gsl_function Fr;
    Fr.function = &frRBME; // frRBME is defined above
    Fr.params = &(frparams[0]); // just need address of the first element of array, because the rest are sequential
    double r1 = 0; // integrate over the range [r1,Inf)
    double epsabsr = NDBD_R_ABS; // absolute error tolerance for the integration
    double epsrelr = NDBD_R_REL; // relative error tolerance for the adaptive integration
    double rbme = 0; // qagiu integration result
    double abserrr; // qagiu integration error estimate
    /*
    double r2 = 10000; // (for debugging) integrate over [r1,r2] ~= [0,Inf), so make r2 reasonably large
    double GKkeyr = 6; // (for debugging) 6 => 61 point Gauss-Kronrod quadrature
    statusr = gsl_integration_qag(&Fr,r1,r2,epsabsr,epsrelr,limitr,GKkeyr,workspacer,&rbme,&abserrr); // (for debugging) perform the qag integration (over r)
    */
    statusr = gsl_integration_qagiu(&Fr,r1,epsabsr,epsrelr,limitr,workspacer,&rbme,&abserrr); // perform the qagiu integration (over r)
    #pragma omp critical
    {
      cpGSLerror(errstrr,statusr,limitr,epsabsr,epsrelr,abserrr,rbme,n,l,np,lp); // check if the GSL integration worked
    }
    gsl_integration_workspace_free(workspacer); // free the allocated memory for the integration workspace, since GSL is written in C
    return rbme;
  }


/// This calculates the analytical integral wrt dr of the RBMEs, with SRCs
/// SRC => SRCs are employed here
/// ie) RBME = < n l | j_\rho(\sqrt{2}qr) | n' l' >, via Equation (4.70) and (4.62) of my thesis
/// NOTE: the \sqrt{2} with the q (and subsequently a 2 with the "a" and "b" SRC parameters) is to stay consistent with r_{rel} = (r_1 - r_2)/SQRT2 in Moshinksy brackets
/// NOTE: rho must be divisible by 2! >:|
  double rbmeSRC(double hw, int rho, double x, int n, int l, int np, int lp, double a, double b, double c)
  {
    double sum1 = 0; // initialize the double sums over k, k'
    double sum2 = 0;
    double sum3 = 0;
    double sum4 = 0;
    double sum5 = 0;
    double sum6 = 0;
    double rbme = 0; // the final RBME value
    double kaptestnn = ((l + lp - rho)/2.0); // the tgamma_ratio below will only take non-negative...
    double kaptestint = abs(kaptestnn - floor(kaptestnn)); // ...integers as arguments
    if ((kaptestnn >= 0) and (kaptestint == 0))
    {
      double bosc = 1.0/sqrt(M_NUCLEON*hw); // the oscillator length [MeV^-1], from Equation (2.4) of my thesis
      double boscsq = bosc*bosc;
      double sqb = SQRT2*x*bosc;
      double sqbsq = sqb*sqb;
      double owbsq23 = 1 + 2*a*boscsq; // 1 + wb^2, NOTE: owbsq1 = 1
      double owbsq456 = 1 + 4*a*boscsq;
      double yolo1 = sqbsq/4.0; // commonly reoccuring arguments in the RBMEs
      double yolo23 = sqbsq/(4.0*owbsq23);
      double yolo456 = sqbsq/(4.0*owbsq456);
      double normy = (1.0/pow(2,rho+1))*sqrt(PI*boost::math::tgamma_ratio(n + 1, n + l + 1.5)*boost::math::tgamma_ratio(np + 1, np + lp + 1.5)); // the normalization factor
      double rholag = rho + 0.5; // for the Laguerre polynomial in the sum
      for (int k=0; k<=n; k++)
      {
        for (int kp=0; kp<=np; kp++)
        {
          double coeff = cpPhase(k + kp)*gsl_sf_gammainv(k + 1)*gsl_sf_gammainv(kp + 1); // the sum coefficient
          double bci = boost::math::tgamma_ratio(n + l + 1.5, n - k + 1); // this is only part of the generalized binomial coeff...
          double bcj = boost::math::tgamma_ratio(np + lp + 1.5, np - kp + 1); // ...we put the other part within kapfacts below
          double bcij = bci*bcj;
          double CBC = coeff*bcij;
          int kappa124 = kaptestnn + k + kp; // kappas
          int kappa35 = kappa124 + 1;
          int kappa6 = kappa124 + 2;
          double tempmax = max(l + k, lp + kp) + 1.5; // just for the kapfacts below
          double tempgi = gsl_sf_gammainv(min(l + k, lp + kp) + 1.5); // " " " " " "
          double kapfact124 = boost::math::tgamma_ratio(kappa124 + 1, tempmax)*tempgi; // kappa factorials
          double kapfact35 = boost::math::tgamma_ratio(kappa35 + 1, tempmax)*tempgi;
          double kapfact6 = boost::math::tgamma_ratio(kappa6 + 1, tempmax)*tempgi;
          double tempks124 = kappa124 + rho + 1.5; // just for the pre-factors below
          double tempks35 = kappa35 + rho + 1.5; // " " " " "
          double pre2 = pow(owbsq23,-tempks124); // pre-factors, NOTE: pre1 = 1
          double pre3 = pow(owbsq23,-tempks35);
          double pre4 = pow(owbsq456,-tempks124);
          double pre5 = pow(owbsq456,-tempks35);
          double pre6 = pow(owbsq456,-(kappa6 + rho + 1.5));
          double Lag1 = gsl_sf_laguerre_n(kappa124,rholag,yolo1); // from GSL -- L^a_n(x) = gsl_sf_laguerre_n(const int n, const double a, const double x)
          double Lag2 = gsl_sf_laguerre_n(kappa124,rholag,yolo23);
          double Lag3 = gsl_sf_laguerre_n(kappa35,rholag,yolo23);
          double Lag4 = gsl_sf_laguerre_n(kappa124,rholag,yolo456);
          double Lag5 = gsl_sf_laguerre_n(kappa35,rholag,yolo456);
          double Lag6 = gsl_sf_laguerre_n(kappa6,rholag,yolo456);
          sum1 += CBC*kapfact124*Lag1; // perform the summations
          sum2 += CBC*pre2*kapfact124*Lag2;
          sum3 += CBC*pre3*kapfact35*Lag3;
          sum4 += CBC*pre4*kapfact124*Lag4;
          sum5 += CBC*pre5*kapfact35*Lag5;
          sum6 += CBC*pre6*kapfact6*Lag6;
        }
      }
      double finpre124 = normy*pow(sqb,rho);
      double finpre35 = finpre124*boscsq;
      double finpre6 = finpre35*boscsq;
      double exp1 = exp(-yolo1);
      double exp23 = exp(-yolo23);
      double exp456 = exp(-yolo456);
      double rbme1 = finpre124*exp1*sum1; // the individual RBMEs which are summed below
      double rbme2 = finpre124*exp23*sum2;
      double rbme3 = finpre35*exp23*sum3;
      double rbme4 = finpre124*exp456*sum4;
      double rbme5 = finpre35*exp456*sum5;
      double rbme6 = finpre6*exp456*sum6;
      rbme = rbme1 - 2*c*rbme2 + 4*b*c*rbme3 + c*c*rbme4 - 4*b*c*c*rbme5 + 4*b*b*c*c*rbme6; // see Equation (4.70) in my thesis
    }
    return rbme;
  }


/// gsl_function handle for integration wrt dr (RBME), including Jastrow-type correlation function for the SRCs
/// this will be called by rbmeSRCqagiu below
  double frRBMEjt(double r, void *params)
  {
    double hw = ((double*) params)[0];
    int rho = ((double*) params)[1];
    double x = ((double*) params)[2];
    int n = ((double*) params)[3];
    int l = ((double*) params)[4];
    int np = ((double*) params)[5];
    int lp = ((double*) params)[6];
    double a = ((double*) params)[7]/HBARC/HBARC; // rescale to [fm^-2] since r is in [fm]
    double b = ((double*) params)[8]/HBARC/HBARC; // rescale to [fm^-2] since r is in [fm]
    double c = ((double*) params)[9];
    x /= HBARC; // put it into [fm^-1] since r is in [fm]
    x *= SQRT2; // this is to stay consistent with r_{rel} = (r_1 - r_2)/SQRT2 in Moshinksy brackets
    a *= 2.0; // " " " " " " " " " " " "
    b *= 2.0; // " " " " " " " " " " " "
    double rsq = r*r;
    double JTF = -c*exp(-a*rsq)*(1 - b*rsq); // the Jastrow-type correlation function
    return rsq*HO_Radial_psi(n,l,hw,r)*gsl_sf_bessel_jl(rho,x*r)*HO_Radial_psi(np,lp,hw,r)*(1 + JTF)*(1 + JTF); // see Equation (6) of PRC.86.067304(2012), for example
  }


/// This calculates the gsl_integration_qagiu wrt dr of the RBMEs, with SRCs
/// SRC => SRCs are employed here
  double rbmeSRCqagiu(double hw, int rho, double x, int n, int l, int np, int lp, double a, double b, double c)
  {
    gsl_set_error_handler_off(); // set standard GSL error handling off, so I can do my own
    int statusr; // this will hold the GSL error status
    string errstrr="7828877-2: qagiu dr:"; // a unique GSL error string
    size_t limitr = NDBD_R_LIMIT; // the number of double precision intervals that can be held by the workspace below
    gsl_integration_workspace * workspacer = gsl_integration_workspace_alloc(limitr); // GSL: "One workspace may be used multiple times...", "...workspaces should be allocated on a per-thread basis."
    double frparams[] = {hw, (double)rho, x, (double)n, (double)l, (double)np, (double)lp, (double)a, (double)b, (double)c};
    gsl_function Fr;
    Fr.function = &frRBMEjt; // frRBMEjt is defined above
    Fr.params = &(frparams[0]); // just need address of the first element of array, because the rest are sequential
    double r1 = 0; // integrate over the range [r1,Inf)
    double epsabsr = NDBD_R_ABS; // absolute error tolerance for the integration
    double epsrelr = NDBD_R_REL; // relative error tolerance for the adaptive integration
    double rbme = 0; // qagiu integration result
    double abserrr; // qagiu integration error estimate
    /*
    double r2 = 10000; // (for debugging) integrate over [r1,r2] ~= [0,Inf), so make r2 reasonably large
    double GKkeyr = 6; // (for debugging) 6 => 61 point Gauss-Kronrod quadrature
    statusr = gsl_integration_qag(&Fr,r1,r2,epsabsr,epsrelr,limitr,GKkeyr,workspacer,&rbme,&abserrr); // (for debugging) perform the qag integration (over r)
    */
    statusr = gsl_integration_qagiu(&Fr,r1,epsabsr,epsrelr,limitr,workspacer,&rbme,&abserrr); // perform the qagiu integration (over r)
    #pragma omp critical
    {
      cpGSLerror(errstrr,statusr,limitr,epsabsr,epsrelr,abserrr,rbme,n,l,np,lp); // check if the GSL integration worked
    }
    gsl_integration_workspace_free(workspacer); // free the allocated memory for the integration workspace, since GSL is written in C
    return rbme;
  }


/// Fermi type form factor
/// q should be in units of [MeV]
  double hF(double q, double naughtV, double lambdaV)
  {
    double gV = naughtV/pow((1.0 + ((q*q)/(lambdaV*lambdaV))),2); // from Equation (4.14) of my thesis
    return (gV*gV)/(naughtV*naughtV); // from Equation (4.11) of my thesis
  }


/// Gamow-Teller type form factor
/// q should be in units of [MeV]
  double hGT(double q, double Mpro, double Mpion, double MagMom, double naughtV, double naughtA, double lambdaV, double lambdaA)
  {
    double qsq = q*q; // q squared [MeV^2]
    double Mprosq = Mpro*Mpro; // the proton mass squared [MeV^2]
    double Mpionsq = Mpion*Mpion; // the pion mass squared [MeV^2]
    double gV = naughtV/pow((1.0 + (qsq/(lambdaV*lambdaV))),2); // from Equation (4.14) of my thesis
    double gA = naughtA/pow((1.0 + (qsq/(lambdaA*lambdaA))),2); // " " " " " "
    double gP = (2*Mpro*gA)/(qsq + Mpionsq); // " " " " " "
    double gM = MagMom*gV; // " " " " " "
    return ((gA*gA) - ((gA*gP*qsq)/(3*Mpro)) + (pow(gP*qsq,2)/(12*Mprosq)) + ((gM*gM*qsq)/(6*Mprosq)))/(naughtA*naughtA); // from Equation (4.12) of my thesis
  }


/// Tensor type form factor
/// q should be in units of [MeV]
  double hT(double q, double Mpro, double Mpion, double MagMom, double naughtV, double naughtA, double lambdaV, double lambdaA)
  {
    double qsq = q*q; // q squared [MeV^2]
    double Mprosq = Mpro*Mpro; // the proton mass squared [MeV^2]
    double Mpionsq = Mpion*Mpion; // the pion mass squared [MeV^2]
    double gV = naughtV/pow((1.0 + (qsq/(lambdaV*lambdaV))),2); // from Equation (4.14) of my thesis
    double gA = naughtA/pow((1.0 + (qsq/(lambdaA*lambdaA))),2); // " " " " " "
    double gP = (2*Mpro*gA)/(qsq + Mpionsq); // " " " " " "
    double gM = MagMom*gV; // " " " " " "
    return (((gA*gP*qsq)/(3*Mpro)) - (pow(gP*qsq,2)/(12*Mprosq)) + ((gM*gM*qsq)/(12*Mprosq)))/(naughtA*naughtA); // from Equation (4.13) of my thesis
  }


/// This is an essential step of "dr" in the dqdrPSH
/// it grabs the relevant rbme... function from the many defined above
/// it will be called by fqF/GT/T below
  double rbmeGrab(int srchit, int drswitch, double hw, int rho, double q, int n, int l, int np, int lp, double a, double b, double c)
  {
    double RBMEval = 0; // relative Bessel's matrix element value
    if (srchit == NDBD_SRC_ON)
    {
      if (drswitch == NDBD_INT_ANAL)
      {
        RBMEval = rbmeSRC(hw,rho,q,n,l,np,lp,a,b,c); // integrate wrt dr analytically, with SRCs
      }
      else if (drswitch == NDBD_INT_QAGIU)
      {
        RBMEval = rbmeSRCqagiu(hw,rho,q,n,l,np,lp,a,b,c); // integrate wrt dr via qagiu, with SRCs
      }
      else
      {
        cout<<"ERROR 37794824-1: dqrswitch has gone out of bounds!?"<<endl;
        cout<<"drswitch = "<<drswitch<<endl;
        exit(1);
      }
    }
    else
    {
      if (drswitch == NDBD_INT_ANAL)
      {
        RBMEval = rbmeNONE(hw,rho,q,n,l,np,lp); // integrate wrt dr analytically, with *no* SRCs
      }
      else if (drswitch == NDBD_INT_QAGIU)
      {
        RBMEval = rbmeNONEqagiu(hw,rho,q,n,l,np,lp); // integrate wrt dr via qagiu, with *no* SRCs
      }
      else
      {
        cout<<"ERROR 37794824-2: dqrswitch has gone out of bounds!?"<<endl;
        cout<<"drswitch = "<<drswitch<<endl;
        exit(1);
      }
    }
    return RBMEval;
  }


/// gsl_function handle for integration wrt dq (F)
/// notice that "lp" is *not* needed here
/// this will be called by dqdrPSH for the Fermi type
  double fqF(double q, void *params)
  {
    int n = ((double*) params)[0];
    int l = ((double*) params)[1];
    int np = ((double*) params)[2];
    //int lp = ((double*) params)[3];
    int drswitch = ((double*) params)[4];
    double hw = ((double*) params)[5];
    double naughtV = ((double*) params)[6];
    double lambdaV = ((double*) params)[7];
    double Ebar = ((double*) params)[8];
    int srchit = ((double*) params)[9];
    double a = ((double*) params)[10];
    double b = ((double*) params)[11];
    double c = ((double*) params)[12];
    double ff = hF(q,naughtV,lambdaV); // call the F form factor
    int rho = 0; // order of the Bessel's function is zero for F
    double RBME = rbmeGrab(srchit,drswitch,hw,rho,q,n,l,np,l,a,b,c); // grab the appropriate RBME
    return (q/(q + Ebar))*ff*RBME; // return the neutrino potential
  }


/// gsl_function handle for integration wrt dq (GT)
/// notice that "lp" is *not* needed here
/// this will be called by dqdrPSH for the Gamow-Teller type
  double fqGT(double q, void *params)
  {
    int n = ((double*) params)[0];
    int l = ((double*) params)[1];
    int np = ((double*) params)[2];
    //int lp = ((double*) params)[3];
    int drswitch = ((double*) params)[4];
    double hw = ((double*) params)[5];
    double naughtV = ((double*) params)[6];
    double lambdaV = ((double*) params)[7];
    double Ebar = ((double*) params)[8];
    int srchit = ((double*) params)[9];
    double a = ((double*) params)[10];
    double b = ((double*) params)[11];
    double c = ((double*) params)[12];
    double Mpro = ((double*) params)[13];
    double Mpion = ((double*) params)[14];
    double MagMom = ((double*) params)[15];
    double naughtA = ((double*) params)[16];
    double lambdaA = ((double*) params)[17];
    double ff = hGT(q,Mpro,Mpion,MagMom,naughtV,naughtA,lambdaV,lambdaA); // call the GT form factor
    int rho = 0; // order of the Bessel's function is zero for GT
    double RBME = rbmeGrab(srchit,drswitch,hw,rho,q,n,l,np,l,a,b,c); // grab the appropriate RBME
    return (q/(q + Ebar))*ff*RBME; // return the neutrino potential
  }


/// gsl_function handle for integration wrt dq (T)
/// notice that "lp" is *indeed* needed here
/// this will be called by dqdrPSH for the Tensor type
  double fqT(double q, void *params)
  {
    int n = ((double*) params)[0];
    int l = ((double*) params)[1];
    int np = ((double*) params)[2];
    int lp = ((double*) params)[3];
    int drswitch = ((double*) params)[4];
    double hw = ((double*) params)[5];
    double naughtV = ((double*) params)[6];
    double lambdaV = ((double*) params)[7];
    double Ebar = ((double*) params)[8];
    int srchit = ((double*) params)[9];
    double a = ((double*) params)[10];
    double b = ((double*) params)[11];
    double c = ((double*) params)[12];
    double Mpro = ((double*) params)[13];
    double Mpion = ((double*) params)[14];
    double MagMom = ((double*) params)[15];
    double naughtA = ((double*) params)[16];
    double lambdaA = ((double*) params)[17];
    double ff = hT(q,Mpro,Mpion,MagMom,naughtV,naughtA,lambdaV,lambdaA); // call the T form factor
    int rho = 2; // order of the Bessel's function is two for T
    double RBME = rbmeGrab(srchit,drswitch,hw,rho,q,n,l,np,lp,a,b,c); // grab the appropriate RBME
    return (q/(q + Ebar))*ff*RBME; // return the neutrino potential
  }


/// Constructor for the NDBD (neutrinoless double-beta decay) class
/// I've made this class to keep all versions of M0nu consistent, via the constants set below
  //unordered_map<long int,double> NDBD::IntList; // global static member decleration
  unordered_map<uint64_t,double> NDBD::IntList; // global static member decleration
  unordered_map<uint64_t,double> NDBD::T6jList; // " " " "
  NDBD::NDBD(double HW, int etwomax, string decaytype, double Ec, string src, string themethod)
  {
    cout<<"Constructing NDBD..."<<endl;
    hw = HW; // construction set by parameter from M0nu
    e2max = etwomax; // " " " " " "
    ////// construct the NDBD parameters //////
    mpro = NDBD_MPRO; // [MeV]
    mpion = NDBD_MPION; // [MeV]
    magmom = NDBD_MAGMOM; // [\mu_N]
    g0V = NDBD_G0V; // [unitless]
    g0A = NDBD_G0A; // [unitless]
    cutoffV = NDBD_CUTV; // [MeV]
    cutoffA = NDBD_CUTA; // [MeV]
    r0 = NDBD_R0; // [fm]
    ///////////////////////////////////////////
    type = decaytype; // construction set by parameter from M0nu, premptive error protection below
    if (type == NDBD_F) {}
    else if (type == NDBD_GT) {}
    else if (type == NDBD_T) {}
    else
    {
      cout<<"ERROR 8973-1: with M0nu type = "<<decaytype<<endl;
      cout<<"this M0nu type has not been added to NDBD yet!"<<endl;
      exit(1);
    }
    cout<<"M0nu type               =  "<<type<<endl;
    Ebar = Ec; // construction set by parameter from M0nu
    cout<<"     closure energy     =  "<<Ebar<<" MeV"<<endl;
    // construct the SRC parameters set by parameter from M0nu
    cout<<"     SRC choice         =  ";
    if (src == NDBD_SRC_AV18)
    {
      srcparams[0] = 1.59*HBARC*HBARC; // [MeV^2]
      srcparams[1] = 1.45*HBARC*HBARC; // [MeV^2]
      srcparams[2] = 0.92; // [unitless]
      srchit = NDBD_SRC_ON; // turn *on* SRCs
      cout<<NDBD_SRC_AV18<<endl;
    }
    else if (src == NDBD_SRC_CDB)
    {
      srcparams[0] = 1.52*HBARC*HBARC; // [MeV^2]
      srcparams[1] = 1.88*HBARC*HBARC; // [MeV^2]
      srcparams[2] = 0.46; // [unitless]
      srchit = NDBD_SRC_ON; // turn *on* SRCs
      cout<<NDBD_SRC_CDB<<endl;
    }
    else if (src == NDBD_SRC_MS)
    {
      srcparams[0] = 1.10*HBARC*HBARC; // [MeV^2] (0.80, 0.90, 1.00, 1.10, 1.20, 1.30) <- [fm^-2] paired with values below, standard is 1.10
      srcparams[1] = 0.68*HBARC*HBARC; // [MeV^2] (0.49, 0.55, 0.62, 0.68, 0.74, 0.80) <- [fm^-2] " " " above, " " 0.68
      srcparams[2] = 1.0; // [unitless]
      srchit = NDBD_SRC_ON; // turn *on* SRCs
      cout<<NDBD_SRC_MS<<endl;
    }
    else if (src == NDBD_SRC_DEBUG)
    {
      srcparams[0] = 0; // [MeV^2]
      srcparams[1] = 0; // [MeV^2]
      srcparams[2] = 0; // [unitless]
      srchit = NDBD_SRC_ON; // turn *on* SRCs
      cout<<NDBD_SRC_DEBUG<<endl;
    }
    else // src == "none" or otherwise
    {
      for (int i=0; i<3; i++)
      {
        srcparams[i] = 0; // set them to zero to avoid variable warnings
      }
      srchit = NDBD_SRC_OFF; // turn *off* SRCs
      cout<<"none"<<endl;
    }
    method = themethod; // construction set by parameter from M0nu, premtive error protection below
    if (method == NDBD_PSH) {}
    //else if (method == NDBD_JE) {}
    else
    {
      cout<<"ERROR 638463-1: with M0nu method = "<<themethod<<endl;
      cout<<"this integration method has not been added to NDBD yet!"<<endl;
      exit(1);
    }
    cout<<"     integration method =  "<<method<<endl;
  }


/// Destructor for the NDBD class
  NDBD::~NDBD()
  {
    cout<<"...Destructing NDBD("<<type<<","<<Ebar<<","<<srchit<<","<<method<<")"<<endl<<endl;
  }


/// These "second tier" low and medium values were found via tests using M0nu_PrintIntegrand
/// for the switches, "dqs" and "drs" we use:        NDBD_INT_ANAL = analytic, NDBD_INT_QAGIU = qagiu, NDBD_INT_QAG = qag
/// it was found that the qagiu over dq breaks for:  n_low,l_low,np_low > 9+,9+,9+
/// " " " " " analytic formula for dr breaks for:    n_med,l_med,np_med > 8,19+,8 or 10,10+,10 or 11,4+,11
/// hence we've set n_low,l_low,np_low = 4,8,4 and n_med,l_med,np_med = 10,10,10 below
/// NOTE: I also found that if drs=NDBD_INT_QAGIU, then we must have dqs=NDBD_INT_QAG or the program chokes
/// NOTE: I also found that setting high mode on for all orbits gives consistent (but slower) results to about machine precision (good sign)
  void NDBD::SetPSH(int n, int l, int np, int lp, int &dqs, int &drs)
  {
    /*
    int n_low  = 100;
    int l_low  = 100;
    int np_low = 100;
    int lp_low = 100;
    int n_med  = 200;
    int l_med  = 200;
    int np_med = 200;
    int lp_med = 200;
    */
    /*
    int n_low  = -1;
    int l_low  = -1;
    int np_low = -1;
    int lp_low = -1;
    int n_med  = 100;
    int l_med  = 100;
    int np_med = 100;
    int lp_med = 100;
    */
    /*
    int n_low  = -1;
    int l_low  = -1;
    int np_low = -1;
    int lp_low = -1;
    int n_med  = -1;
    int l_med  = -1;
    int np_med = -1;
    int lp_med = -1;
    */
    //
    int n_low  = 4;
    int l_low  = 8;
    int np_low = 4;
    int lp_low = 8;
    int n_med  = 10;
    int l_med  = 10;
    int np_med = 10;
    int lp_med = 10;
    //
    dqs = 0; // this will set the integration style over dq, set to zero for precaution
    drs = 0; // " " " " " " " dr, " " " " "
    if (n<=n_low and l<=l_low and np<=np_low and lp<=lp_low) // qagiu for dq, and analytic for dr
    {
//cout<<"low"<<endl;
      dqs = NDBD_INT_QAGIU;
      drs = NDBD_INT_ANAL;
    }
    else if (n<=n_med and l<=l_med and np<=np_med and lp<=lp_med) // qag for dq, and analytic for dr
    {
//cout<<"medium"<<endl;
      dqs = NDBD_INT_QAG;
      drs = NDBD_INT_ANAL;
    }
    else // qag for dq, and qagiu for dr
    {
//cout<<"high"<<endl;
      dqs = NDBD_INT_QAG;
      drs = NDBD_INT_QAGIU;
    }
  }


/// Set gsl_function handle parameters, for integration wrt dq
  void NDBD::SetIparams(int n, int l, int np, int lp, int drs, double Iparams[18])
  {
    Iparams[0] = n;
    Iparams[1] = l;
    Iparams[2] = np;
    Iparams[3] = lp;
    Iparams[4] = (double)drs;
    Iparams[5] = hw;
    Iparams[6] = g0V;
    Iparams[7] = cutoffV;
    Iparams[8] = Ebar;
    Iparams[9] = (double)srchit;
    Iparams[10] = srcparams[0];
    Iparams[11] = srcparams[1];
    Iparams[12] = srcparams[2];
    Iparams[13] = mpro; // these params are only...
    Iparams[14] = mpion;
    Iparams[15] = magmom;
    Iparams[16] = g0A;
    Iparams[17] = cutoffA; // ...used in GT and T
  }


/// Perform the integration over dr and then dq
/// This method is what I mean by "two-tiered adaptive"
/// ie) Note the difference between when dqswitch == NDBD_INT_QAG or NDBD_INT_QAGIU
  double NDBD::dqdrPSH(int n, int l, int np, int lp)
  {
    int dqswitch,drswitch; // dq and dr switches, see NDBD::SetPSH
    gsl_function Fq; // this function will point to fqF/GT/T, NOTE: adding this to the NDBD constructor will break parallelization in PreCalcIntegrals()
    double fqparams[18]; // parameters for Fq
    // set up the switches and gsl_function handle
    SetPSH(n,l,np,lp,dqswitch,drswitch); // set the integration methods based on relative orbit height
    SetIparams(n,l,np,lp,drswitch,fqparams); // fqparams = Iparams
    Fq.params = &(fqparams[0]); // just need address of the first element of array, because the rest are sequential
    if (type == NDBD_GT)
    {
      Fq.function = &fqGT; // point to the Gamow-Teller type
    }
    else if (type == NDBD_F)
    {
      Fq.function = &fqF; // point to the Fermi type
    }
    else if (type == NDBD_T)
    {
      Fq.function = &fqT; // point to the Tensor type
    }
    else
    {
      cout<<"ERROR 8973-2: with M0nu type = "<<type<<endl;
      cout<<"this M0nu type has not been added to NDBD.dqdrPSH yet!"<<endl;
      exit(1);
    }
    // set up stuff needed for gsl integration
    gsl_set_error_handler_off(); // set standard GSL error handling off, so I can do my own
    int statusq; // this will hold the GSL error status
    string errstrq; // a unique GSL error string
    size_t limitq = omp_get_num_threads()*NDBD_Q_LIMIT; // the number of double precision intervals that can be held by the workspace below
    gsl_integration_workspace * workspaceq = gsl_integration_workspace_alloc(limitq); // GSL: "One workspace may be used multiple times...", "...workspaces should be allocated on a per-thread basis."
    double q1 = 0; // starting q for the integration is zero by definition, don't change this
    double epsabsq = NDBD_Q_ABS; // absolute error tolerance for the integration
    double epsrelq = NDBD_Q_REL; // relative error tolerance for the adaptive integration
    double resultq = 0; // qag(iu) integration result
    double abserrq; // qag(iu) integration error estimate
    if (dqswitch == NDBD_INT_QAGIU)
    {
      errstrq = "7828877-3: qagiu dq:";
      statusq = gsl_integration_qagiu(&Fq,q1,epsabsq,epsrelq,limitq,workspaceq,&resultq,&abserrq); // perform the qagiu integration (over q)
    }
    else if (dqswitch == NDBD_INT_QAG)
    {
      errstrq = "7828877-4: qag dq:";
      double q2 = 2500; // integrate over [q1,q2] ~= [0,Inf), so make q2 reasonably large
      double GKkeyq = 6; // 6 => 61 point Gauss-Kronrod quadrature
      statusq = gsl_integration_qag(&Fq,q1,q2,epsabsq,epsrelq,limitq,GKkeyq,workspaceq,&resultq,&abserrq); // perform the qag integration (over q)
    }
    else
    {
      cout<<"ERROR 37794824-3: dqswitch has gone out of bounds!?"<<endl;
      cout<<"dqswitch = "<<dqswitch<<endl;
      exit(1);
    }
    #pragma omp critical
    {
      cpGSLerror(errstrq,statusq,limitq,epsabsq,epsrelq,abserrq,resultq,n,l,np,lp); // check if the GSL integration worked
    }
    gsl_integration_workspace_free(workspaceq); // free the allocated memory for the integration workspace, since GSL is written in C
    /*
    cout<<"limitq = "<<limitq<<endl; // debugging...
    cout<<"epsabsq = "<<epsabsq<<endl;
    cout<<"epsrelq = "<<epsrelq<<endl;
    cout<<"abserrq = "<<abserrq<<endl; // ...debugging
    */
    cout<<"in dqdrPSH: n="<<n<<", l="<<l<<", np="<<np<<", lp="<<lp<<"; resultq = "<<resultq<<endl; // debugging (wtf)
    return resultq;
  }


/// Hashtag for the IntList cache
  uint64_t NDBD::IntHash(uint64_t n, uint64_t l, uint64_t np, uint64_t lp)
  {
    return   (n  << 21)
           + (l  << 15)
           + (np <<  8)
           +  lp;
  }


/// Inverse Hashtag for the IntList cache
  void NDBD::IntUnHash(uint64_t key, uint64_t &n, uint64_t &l, uint64_t &np, uint64_t &lp)
  {
    n  = (key >> 21) & 0x7FL;
    l  = (key >> 15) & 0x3FL;
    np = (key >>  8) & 0x7FL;
    lp = (key      ) & 0xFFL;
  }


/// NDBD member function which computes and cache's the relevant integrals wrt dq and dr
/// method = NDBD_PSH => for n,l,n' = 0,0,0 up to nlow,llow,nplow we use:            qagiu for dq, and analytic for dr
///                   for n,l,n' = nlow,llow,nplow up to nmed,lmed,npmed we use:  qag for dq,   and analytic for dr
///                   for n,l,n' = nmed++,lmed++,npmed++ and above we use:        qag for dq,   and qagiu for dr
/// method = "JE" => modified Simpson's rule for dq, and a right-justified Riemann's sum for dr (not coded up yet)
  void NDBD::PreCalcIntegrals()
  {
    double t_start_pci = omp_get_wtime(); // profiling (s)
    cout<<"calculating integrals wrt dq and dr..."<<endl;
    int maxn = e2max/2;
    int maxl = e2max;
    int maxnp = e2max/2;
    int maxlp = e2max;
    // fill up the KEYS (via NDBD::IntHash) first
    //vector<long int> KEYS; // if I ever get that version working...
    vector<uint64_t> KEYS;
    if (type == NDBD_F or type == NDBD_GT)
    {
      for (int n=0; n<=maxn; n++)
      {
        for (int l=0; l<=maxl; l++)
        {
          int tempminnp = n; // NOTE: need not start from 'int np=0' since IntHash(n,l,np,l) = IntHash(np,l,n,l), by construction
          for (int np=tempminnp; np<=maxnp; np++)
          {
            //long int key = IntHash(n,l,np);
            uint64_t key = IntHash(n,l,np,l);
            KEYS.push_back(key);
            IntList[key] = 0.0; // "Make sure eveything's in there to avoid a rehash in the parallel loop" (RS)
          }
        }
      }
    }
    else if (type == NDBD_T)
    {
      for (int n=0; n<=maxn; n++)
      {
        for (int l=0; l<=maxl; l++)
        {
          int tempminnp = n; // NOTE: need not start from 'int np=0' since IntHash(n,l,np,lp) = IntHash(np,lp,n,l), by construction
          for (int np=tempminnp; np<=maxnp; np++)
          {
            int tempminlp = (n==np ? l : 0); // NOTE: need not start from 'int lp=0' since IntHash(n,l,np,lp) = IntHash(np,lp,n,l), by construction
            for (int lp=tempminlp; lp<=maxlp; lp++)
            {
              //long int key = IntHash(n,l,np,lp);
              uint64_t key = IntHash(n,l,np,lp);
              KEYS.push_back(key);
              IntList[key] = 0.0; // "Make sure eveything's in there to avoid a rehash in the parallel loop" (RS)
            }
          }
        }
      }
    }
    else
    {
      cout<<"ERROR 8973-3: with M0nu type = "<<type<<endl;
      cout<<"this M0nu type has not been added to NDBD.PreCalcIntegrals() yet!"<<endl;
      exit(1);
    }
    // now we'll fill up the static member NDBD.IntList, declared above the constructor
    if (method == NDBD_PSH)
    {
      #pragma omp parallel for schedule(dynamic,1) // this works as long as the gsl_function handle is within this for-loop
      for (size_t i=0; i<KEYS.size(); i++)
      {
        //long int key = KEYS[i];
        //int n,l,np,lp; // get these from IntUnHash below
        uint64_t key = KEYS[i];
        uint64_t n,l,np,lp; // get these from IntUnHash below
        IntUnHash(key,n,l,np,lp);
        IntList[key] = dqdrPSH(n,l,np,lp); // these have been ordered by the above loops such that we take the "lowest" value of decimalgen(n,l,np,lp,maxl,maxnp,maxlp), see GetIntegral(...)
        cout<<"in PreCalcIntegrals: n="<<n<<", l="<<l<<", np="<<np<<", lp="<<lp<<"; I = "<<IntList[key]<<endl; // debugging
      }
    }
    /*
    if (method == "JE")
    {
      //
    }
    */
    else
    {
      cout<<"ERROR 638463-3: with M0nu method = "<<method<<endl;
      cout<<"this integration method has not been added to NDBD.PreCalcIntegrals() yet!"<<endl;
      exit(1);
    }
    cout<<"...done calculating the integrals"<<endl;
    cout<<"IntList has "<<IntList.bucket_count()<<" buckets and a load factor "<<IntList.load_factor()
      <<", estimated storage ~= "<<((IntList.bucket_count() + IntList.size())*(sizeof(size_t) + sizeof(void*)))/(1024.0*1024.0*1024.0)<<" GB"<<endl; // copied from (RS)
    profiler.timer["PreCalcIntegrals"] += omp_get_wtime() - t_start_pci; // profiling (r)
  }


/// Get an integral from the IntList cache or calculate it (parallelization dependent)
  double NDBD::GetIntegral(int n, int l, int np, int lp)
  {
    //int maxn = e2max/2;
    int maxl = e2max;
    int maxnp = e2max/2;
    int maxlp = e2max;
    int order1 = decimalgen(n,l,np,lp,maxl,maxnp,maxlp);
    int order2 = decimalgen(np,lp,n,l,maxl,maxnp,maxlp); // notice I was careful here with the order of maxl,maxnp,maxlp to make proper comparison
    if (order1 > order2)
    {
      swap(n,np); // using symmetry IntHash(n,l,np,lp) = IntHash(np,lp,n,l)
      swap(l,lp); // " " " " "
    }
    //long int key = IntHash(n,l,np,lp); // if I ever get that version working...
    uint64_t key = IntHash(n,l,np,lp);
    auto it = IntList.find(key);
    if (it != IntList.end()) // return what we've found
    {
      return it -> second;
    }
    else // if we didn't find it, calculate it and add it to the list!
    {
      double integral;
      if (method == NDBD_PSH)
      {
        integral = dqdrPSH(n,l,np,lp);
      }
      /*
      else if (method == "JE")
      {
        //
      }
      */
      else
      {
        cout<<"ERROR 638463-2: with M0nu method = "<<method<<endl;
        cout<<"this integration method has not been added to NDBD.GetIntegral(...) yet!"<<endl;
        exit(1);
      }
      if (omp_get_num_threads() >= 2)
      {
        printf("DANGER!!!!!!!  Updating IntList inside a parellel loop breaks thread safety!\n");
        printf("   I shouldn't be here in GetIntegral(%d, %d, %d, %d):   key =%lx   integral=%f\n",n,l,np,lp,key,integral);
        exit(EXIT_FAILURE);
      }
      IntList[key] = integral;
      return integral;
    }
  }


/// Hashtag for the T6jList cache
  uint64_t NDBD::T6jHash(uint64_t l1, uint64_t L1, uint64_t R, uint64_t L2, uint64_t l2)
  {
    return   (l1 << 40)
           + (L1 << 30)
           + (R  << 20)
           + (L2 << 10)
           +  l2;
  }


/// Inverse Hashtag for the T6jList cache
  void NDBD::T6jUnHash(uint64_t key, uint64_t &l1, uint64_t &L1, uint64_t &R, uint64_t &L2, uint64_t &l2)
  {
    l1 = (key >> 40) & 0x3FFL;
    L1 = (key >> 30) & 0x3FFL;
    R  = (key >> 20) & 0x3FFL;
    L2 = (key >> 10) & 0x3FFL;
    l2 = (key      ) & 0x3FFL;
  }


/// NDBD member function which caches the relevant 6j-symbols for the Tensor M0nu component
/// { l1 L1 R }
/// { L2 l2 2 } NOTE: the "2" in the bottom-right entry
  void NDBD::PreCalcT6j()
  {
    double t_start_pctsj = omp_get_wtime(); // profiling (s)
    cout<<"calculating 6j-symbols for the Tensor component..."<<endl;
    vector<uint64_t> KEYS;
    for (int l1=0; l1<=e2max; l1++)
    {
      for (int l2=0; l2<=e2max; l2++)
      {
        for (int L1=0; L1<=e2max; L1++)
        {
          for (int L2=0; L2<=e2max; L2++)
          {
            for (int R=0; R<=e2max; R++)
            {
              uint64_t key = T6jHash(l1,L1,R,L2,l2);
              KEYS.push_back(key);
              T6jList[key] = 0.0; // "Make sure eveything's in there to avoid a rehash in the parallel loop" (RS)
            }
          }
        }
      }
    }
    #pragma omp parallel for schedule(dynamic,1)
    for (size_t i=0; i<KEYS.size(); i++)
    {
      uint64_t l1,l2,L1,L2,R;
      uint64_t key = KEYS[i];
      T6jUnHash(key,l1,L1,R,L2,l2);
      T6jList[key] = AngMom::SixJ(l1,L1,R,L2,l2,2);
    }
    cout<<"...done calculating the 6j-symbols"<<endl;
    cout<<"T6jList has "<<T6jList.bucket_count()<<" buckets and a load factor "<<T6jList.load_factor()
      <<", estimated storage ~= "<<((T6jList.bucket_count() + T6jList.size())*(sizeof(size_t) + sizeof(void*)))/(1024.0*1024.0*1024.0)<<" GB"<<endl; // copied from (RS)
    profiler.timer["PreCalcT6j"] += omp_get_wtime() - t_start_pctsj; // profiling (r)
  }


/// Get a 6j from the T6jList cache or calculate it (parallelization dependent)
/// { l1 L1 R }
/// { L2 l2 2 } NOTE: the "2" in the bottom-right entry
  double NDBD::GetT6j(int l1, int L1, int R, int L2, int l2)
  {
    uint64_t key = T6jHash(l1,L1,R,L2,l2);
    auto it = T6jList.find(key);
    if (it != T6jList.end()) // return what we've found
    {
      return it -> second;
    }
    else // if we didn't find it, calculate it and add it to the list!
    {
      double sixj = AngMom::SixJ(l1,L1,R,L2,l2,2);
      if (omp_get_num_threads() >= 2)
      {
        printf("DANGER!!!!!!!  Updating T6jList inside a parellel loop breaks thread safety!\n");
        printf("   I shouldn't be here in GetT6j(%d, %d, %d, %d, %d):   key =%lx   sixj=%f\n",l1,L1,R,L2,l2,key,sixj);
        exit(EXIT_FAILURE);
      }
      T6jList[key] = sixj;
      return sixj;
    }
  }


/// returns a string of five random letters and a time stamp, to uniquely identify an M0nu run
  string cpBarcode(string Type, double hw, double emax)
  {
    double typenum = 0; // to help distinguish barcodes between different decay types in M0nu...
    if (Type == NDBD_GT)
    {
      typenum = 2389; // ...it's arbitrary, as long as it's different...
    }
    else if (Type == NDBD_F)
    {
      typenum = 1278; // ...from this...
    }
    else
    {
      typenum = 3018; // ...and this
    }
    double seed = omp_get_wtime() + typenum + 10*(hw + emax); // seed rand with time, typenum, hw, and emax
    srand(seed); // clear rand upon execution
    string alpha = "abcdefghijklmnopqrstuvwxyz";
    int v1 = rand()%26;
    int v2 = rand()%26;
    int v3 = rand()%26;
    int v4 = rand()%26;
    int v5 = rand()%26;
    time_t now = time(0);
    struct tm tstruct;
    char timestamp[80];
    tstruct = *localtime(&now);
    strftime(timestamp, sizeof(timestamp), "%y%m%d%H%M", &tstruct);
    stringstream rrrrr;
    rrrrr<<alpha[v1]<<alpha[v2]<<alpha[v3]<<alpha[v4]<<alpha[v5]<<timestamp;
    string barcode = rrrrr.str(); // a string of five random letters and a time stamp
    return barcode;
  }


/// This outputs an M0nu header file to the directory "dirname"
/// if you're intending to output the M0nu TBME, it *absolutely* should be called for any version of M0nu created
/// and then WriteM0nu should be used
  void M0nuHeader(NDBD zvbb, string dirname, string reduced, int Anuc, string src, double Rnuc, double prefact, string barcode)
  {
    ofstream header(dirname+"M0nu_header_"+barcode+".txt");
    if (!header)
    {
      cout<<"ERROR 432337: with M0nu header"<<endl;
      cout<<"a likely cause is with dirname = "<<dirname<<endl;
      exit(1);
    }
    //header<<" debugging = running with new optimization of KEYS ordering with decimalgen..."<<endl; // debugging message...
    //header<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl; // ...debugging message
    header<<" decay type             =  "<<zvbb.type<<" (lit)"<<endl;
    header<<" A                      =  "<<Anuc<<endl;
    header<<" E_closure              =  "<<setprecision(12)<<zvbb.Ebar<<" [MeV]"<<endl;
    header<<" emax                   =  "<<zvbb.e2max/2<<endl;
    header<<" hw                     =  "<<setprecision(12)<<zvbb.hw<<" [MeV]"<<endl;
    header<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
    header<<" integration method     =  "<<zvbb.method<<endl;
    header<<" dr |  GSL limit        =  "<<NDBD_R_LIMIT<<endl;
    header<<" dr |  abstol           =  "<<cpFormat(NDBD_R_ABS)<<endl;
    header<<" dr |  reltol           =  "<<cpFormat(NDBD_R_REL)<<endl;
    header<<" dq || GSL limit        =  "<<NDBD_Q_LIMIT<<endl;
    header<<" dq || abstol           =  "<<cpFormat(NDBD_Q_ABS)<<endl;
    header<<" dq || reltol           =  "<<cpFormat(NDBD_Q_REL)<<endl;
    header<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
    header<<" reduced?               =  ";
    if (reduced == "NR")
    {
      header<<"No"<<endl;
    }
    else
    {
      header<<"Yes"<<endl;
    }
    double precheck = (2*Rnuc)/PI; // this is the common prefactor of M0nu used in the literature
    header<<" global pre-factor      =  "<<setprecision(12)<<prefact/precheck<<endl;
    header<<" R = r_0*A^(1/3)        =  "<<setprecision(12)<<Rnuc<<" [fm]"<<endl;
    header<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
    header<<" srchit                 =  "<<zvbb.srchit<<endl;
    header<<" are SRCs on?           =  ";
    if (zvbb.srchit == NDBD_SRC_ON)
    {
      header<<"Yes"<<endl;
    }
    else
    {
      header<<"No"<<endl;
    }
    header<<" SRCs                   =  "<<src<<endl;
    header<<" a                      =  "<<setprecision(12)<<zvbb.srcparams[0]/HBARC/HBARC<<" [fm^-2]"<<endl;
    header<<" b                      =  "<<setprecision(12)<<zvbb.srcparams[1]/HBARC/HBARC<<" [fm^-2]"<<endl;
    header<<" c                      =  "<<setprecision(12)<<zvbb.srcparams[2]<<" [unitless]"<<endl;
    header<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
    header<<" m_proton               =  "<<setprecision(12)<<zvbb.mpro<<" [MeV]"<<endl;
    header<<" m_pion                 =  "<<setprecision(12)<<zvbb.mpion<<" [MeV]"<<endl;
    header<<" \\mu_p - \\mu_n          =  "<<setprecision(12)<<zvbb.magmom<<" [\\mu_N]"<<endl;
    header<<" g_{V,0}                =  "<<setprecision(12)<<zvbb.g0V<<" [unitless]"<<endl;
    header<<" g_{A,0}                =  "<<setprecision(12)<<zvbb.g0A<<" [unitless]"<<endl;
    header<<" cutoff_V               =  "<<setprecision(12)<<zvbb.cutoffV<<" [MeV]"<<endl;
    header<<" cutoff_A               =  "<<setprecision(12)<<zvbb.cutoffA<<" [MeV]"<<endl;
    header<<" r_0 in R               =  "<<setprecision(12)<<zvbb.r0<<" [fm]"<<endl;
    header<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
    header<<" HBARC                  =  "<<setprecision(12)<<HBARC<<" [MeV*fm]"<<endl;
    header<<" M_NUCLEON              =  "<<setprecision(12)<<M_NUCLEON<<" [MeV]"<<endl;
    header.close();
  }


/// FIN ///
