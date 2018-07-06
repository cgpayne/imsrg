
#include "imsrg_util.hh"
#include "AngMom.hh"
#include "NDBD.hh"
#include <gsl/gsl_integration.h>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <vector>

using namespace AngMom;

/// imsrg_util namespace. Used to define some helpful functions.
namespace imsrg_util
{

 std::vector<string> split_string(string s, string delimiter)
 {
  vector<string> out;
   size_t last = 0;
   size_t next = 0;
   while ((next = s.find(delimiter, last)) != string::npos)
   {
     out.push_back( s.substr(last, next-last) );
     last = next + 1;
   }
   out.push_back( s.substr(last) );
   return out;
 }

 Operator OperatorFromString(ModelSpace& modelspace, string &opname)
 {
      vector<string> opnamesplit = split_string( opname, "_" );  // split string on _ into a vector of string so that, e.g. "R2_p1"  =>  {"R2", "p1"}

           if (opname == "R2_p1")         return R2_1body_Op(modelspace,"proton") ;
      else if (opname == "R2_p2")         return R2_2body_Op(modelspace,"proton") ;
      else if (opname == "R2_n1")         return R2_1body_Op(modelspace,"neutron") ;
      else if (opname == "R2_n2")         return R2_2body_Op(modelspace,"neutron") ;
      else if (opname == "Rp2")           return Rp2_corrected_Op(modelspace,modelspace.GetTargetMass(),modelspace.GetTargetZ()) ;
      else if (opname == "Rn2")           return Rn2_corrected_Op(modelspace,modelspace.GetTargetMass(),modelspace.GetTargetZ()) ;
      else if (opname == "Rm2")           return Rm2_corrected_Op(modelspace,modelspace.GetTargetMass(),modelspace.GetTargetZ()) ;
      else if (opname == "E1")            return ElectricMultipoleOp(modelspace,1) ;
      else if (opname == "E2")            return ElectricMultipoleOp(modelspace,2) ;
      else if (opname == "E3")            return ElectricMultipoleOp(modelspace,3) ;
      else if (opname == "E4")            return ElectricMultipoleOp(modelspace,4) ;
      else if (opname == "E5")            return ElectricMultipoleOp(modelspace,5) ;
      else if (opname == "E6")            return ElectricMultipoleOp(modelspace,6) ;
      else if (opname == "E2int")         return IntrinsicElectricMultipoleOp(modelspace,2) ; // Untested
      else if (opname == "nE2")           return NeutronElectricMultipoleOp(modelspace,2) ;
      else if (opname == "M1")            return MagneticMultipoleOp(modelspace,1) ;
      else if (opname == "M2")            return MagneticMultipoleOp(modelspace,2) ;
      else if (opname == "M3")            return MagneticMultipoleOp(modelspace,3) ;
      else if (opname == "M4")            return MagneticMultipoleOp(modelspace,4) ;
      else if (opname == "M5")            return MagneticMultipoleOp(modelspace,5) ;
      else if (opname == "M1p")           return MagneticMultipoleOp_pn(modelspace,1,"proton") ;
      else if (opname == "M1n")           return MagneticMultipoleOp_pn(modelspace,1,"neutron") ;
      else if (opname == "Fermi")         return AllowedFermi_Op(modelspace) ;
      else if (opname == "GamowTeller")   return AllowedGamowTeller_Op(modelspace) ;
      else if (opname == "Iso2")          return Isospin2_Op(modelspace) ;
      else if (opname == "R2CM")          return R2CM_Op(modelspace) ;
      else if (opname == "TCM")           return TCM_Op(modelspace) ;
      else if (opname == "Rso")           return RpSpinOrbitCorrection(modelspace) ;
      else if (opname == "RadialOverlap") return RadialOverlap(modelspace); // Untested...
      else if (opname == "Sigma")         return Sigma_Op(modelspace);
      else if (opname == "Sigma_p")       return Sigma_Op_pn(modelspace,"proton");
      else if (opname == "Sigma_n")       return Sigma_Op_pn(modelspace,"neutron");
      else if (opname == "L2rel")         return L2rel_Op(modelspace); // Untested...
      else if (opname == "DGT")           return DGT_Op(modelspace); // Double-Gamow-Teller TBMEs (CP), Untested...
      else if (opnamesplit[0] =="HCM")
      {
         if ( opnamesplit.size() == 1 ) return HCM_Op(modelspace);
         double hw_HCM; // frequency of trapping potential
         istringstream( opnamesplit[1] ) >> hw_HCM;
         int A = modelspace.GetTargetMass();
         return TCM_Op(modelspace) + 0.5*A*M_NUCLEON*hw_HCM*hw_HCM/HBARC/HBARC*R2CM_Op(modelspace);
      }
      else if (opnamesplit[0] == "VCM") // GetHCM with a different frequency, ie HCM_24 for hw=24
      {
         double hw_VCM; // frequency of trapping potential
         istringstream(opnamesplit[1]) >> hw_VCM;
         int A = modelspace.GetTargetMass();
         return 0.5*A*M_NUCLEON*hw_VCM*hw_VCM/HBARC/HBARC*R2CM_Op(modelspace); 
      }
      else if (opnamesplit[0] == "Rp2Z") // Get point proton radius for specified Z, e.g. Rp2Z_10 for neon
      {
        int Z_rp;
        istringstream(opnamesplit[1]) >> Z_rp;
        return Rp2_corrected_Op(modelspace,modelspace.GetTargetMass(),Z_rp) ;
      }
      else if (opnamesplit[0] == "Rp2AZ") // Get point proton radius for specified A and Z, e.g. Rp2AZ_20_10 for neon
      {
        int A_rp;
        int Z_rp;
        istringstream(opnamesplit[1]) >> A_rp;
        istringstream(opnamesplit[2]) >> Z_rp;
        return Rp2_corrected_Op(modelspace,A_rp,Z_rp) ;
      }
      else if (opnamesplit[0] == "Rn2Z") // Get point neutron radius for specified Z
      {
        int Z_rp;
        istringstream(opnamesplit[1]) >> Z_rp;
        return Rn2_corrected_Op(modelspace,modelspace.GetTargetMass(),Z_rp) ;
      }
      else if (opnamesplit[0] == "rhop") // point proton  density at position r, e.g. rhop_1.25
      {
        double rr;
        istringstream(opnamesplit[1]) >> rr;
        return ProtonDensityAtR(modelspace,rr);
      }
      else if (opnamesplit[0] == "rhon") // point neutron density at position r
      {
        double rr;
        istringstream(opnamesplit[1]) >> rr;
        NeutronDensityAtR(modelspace,rr);
      }
      else if (opnamesplit[0] == "OneOcc") // Get occupation of specified orbit, e.g. OneOccp_1p3
      {
         index_t ind = modelspace.String2Index( {  opnamesplit[1] } )[0];
         Orbit& oi = modelspace.GetOrbit(ind);
         return NumberOp(modelspace,oi.n,oi.l,oi.j2,oi.tz2) ;
      }
      else if (opnamesplit[0]== "AllOcc") // Get occupation of orbit, summed over all values of radial quantum number n, e.g. AllOccpp3
      {
         index_t ind = modelspace.String2Index( { "0"+ opnamesplit[1] } )[0];
         Orbit& oi = modelspace.GetOrbit(ind);
         return NumberOpAlln(modelspace,oi.l,oi.j2,oi.tz2) ;
      }
      else if (opnamesplit[0] == "protonFBC") // Fourier bessel coefficient of order nu
      {
         int nu;
         istringstream(opnamesplit[1]) >> nu;
         return FourierBesselCoeff( modelspace, nu, 8.0, modelspace.proton_orbits);
      }
      else if (opnamesplit[0] == "neutronFBC") // Fourier bessel coefficient of order nu
      {
         int nu;
         istringstream(opnamesplit[1]) >> nu;
         return FourierBesselCoeff( modelspace, nu, 8.0, modelspace.neutron_orbits) ;
      }
      /*
      else if (opnamesplit[0] == "M0nu" and opnamesplit[1] == "adpt") // 0\nu\beta\beta decay TBME, M0nu_adpt_${dirout}_${Decay}_${Reduced}_${Ec}_${SRC} (CP)
      {
         string dirname; // the directory to hold M0nu_header_[barcode][date].txt
         string Type; // the decay type (F = "Fermi", GT = "Gamow-Teller", T = "Tensor")
         string reduced; // to reduced the matrix elements or not (NR = "not reduced", R = "reduced")
         double Ediff; // the closure energy [MeV]
         string src; // chosen SRC parameters (none, AV18, CD-Bonn, Miller-Spencer, debug)
         istringstream(opnamesplit[2]) >> dirname; // ${dirout}
cout<<"dirname = "<<dirname<<endl; // debug
         istringstream(opnamesplit[3]) >> Type; // ${Decay}
cout<<"Type = "<<Type<<endl; // debug
         istringstream(opnamesplit[4]) >> reduced; // ${Reduced}
cout<<"reduced = "<<reduced<<endl; // debug
         istringstream(opnamesplit[5]) >> Ediff; // ${Ec}
cout<<"Ediff = "<<Ediff<<endl; // debug
         istringstream(opnamesplit[6]) >> src; // ${SRC}
cout<<"src = "<<src<<endl; // debug
         double hw = modelspace.GetHbarOmega();
         int emax = modelspace.GetEmax();
         string barcode = cpBarcode(Type,hw,emax); // cpBarcode() found in NDBD.cc
         cout<<endl<<"barcode = "<<barcode<<endl<<endl;
         opname = opnamesplit[0] +
           '_' + opnamesplit[1] +
           '_' + opnamesplit[3] +
           '_' + opnamesplit[4] +
           '_' + opnamesplit[5] +
           '_' + opnamesplit[6] + barcode; // reset opname to not have dirname in it (skipped [2], makes the opname too long), but add the barcode for imsrg++.cc
cout<<"(new) opname = "<<opname<<endl; // debug
exit(1); // debug
         if (Type == "GT" or Type == "F")
         {
           return M0nuGTF_adpt_Op(modelspace,dirname,Type,reduced,Ediff,src,barcode);
         }
         else if (Type == "T")
         {
           return M0nuT_adpt_Op(modelspace,dirname,reduced,Ediff,src,barcode);
         }
         else
         {
           cout<<"ERROR 8973-4: with M0nu Type = "<<Type<<endl;
           cout<<"this M0nu type has not been coded as an M0nu Operator yet!"<<endl;
           exit(1);
         }
      }
      */
      // NOTE: opnamesplit does not work for the following operator, because a directory is an argument, and said directory may contain the delimiter "_" within (see attempt above)
      else if (opname.substr(0,10) == "M0nu_adpt_") // 0\nu\beta\beta decay TBME, M0nu_adpt_${dirout}_${Decay}_${Reduced}_${Ec}_${SRC} (CP)
      {
         string dirname; // the directory to hold M0nu_header_[barcode][date].txt
         string Type; // the decay type (F = "Fermi", GT = "Gamow-Teller", T = "Tensor")
         string reduced; // to reduced the matrix elements or not (NR = "not reduced", R = "reduced")
         double Ediff; // the closure energy [MeV]
         string src; // chosen SRC parameters (none, AV18, CD-Bonn, Miller-Spencer, debug)
         string segment;
         stringstream ssopnames1(opname);
         vector<string> seglist1;
         while(getline(ssopnames1,segment,'/'))
         {
           seglist1.push_back(segment);
         }
         for (unsigned int i=1; i<seglist1.size()-1; i++)
         {
           dirname = dirname + "/" + seglist1.at(i);
         }
         dirname = dirname + "/"; // ${dirout}
         string tempname = opname;
         string dirsub = "_" + dirname;
         string::size_type nd = dirsub.length();
         for (string::size_type i = tempname.find(dirsub); i != string::npos; i = tempname.find(dirsub))
         {
           tempname.erase(i,nd);
         }
         stringstream ssopnames2(tempname);
         vector<string> seglist2;
         while(getline(ssopnames2,segment,'_'))
         {
           seglist2.push_back(segment);
         }
         Type = seglist2.at(2); // ${Decay}
         reduced = seglist2.at(3); // ${Reduced}
         Ediff = stod(seglist2.at(4)); // ${Ec}
         src = seglist2.at(5); // ${SRC}
         double hw = modelspace.GetHbarOmega();
         int emax = modelspace.GetEmax();
         string barcode = cpBarcode(Type,hw,emax); // cpBarcode() found in NDBD.cc
         cout<<endl<<"barcode = "<<barcode<<endl<<endl;
         opname = tempname + '_' + barcode; // reset opname to not have dirname in it (makes the opname too long), but add the barcode for imsrg++.cc
         if (Type == "GT" or Type == "F")
         {
           return M0nuGTF_adpt_Op(modelspace,dirname,Type,reduced,Ediff,src,barcode);
         }
         else if (Type == "T")
         {
           return M0nuT_adpt_Op(modelspace,dirname,reduced,Ediff,src,barcode);
         }
         else
         {
           cout<<"ERROR 8973-4: with M0nu Type = "<<Type<<endl;
           cout<<"this M0nu type has not been coded as an M0nu Operator yet!"<<endl;
           exit(1);
         }
      }
      // NOTE: opnamesplit does not work for the following operator, because a directory is an argument, and said directory may contain the delimiter "_" within
      else if (opname.substr(0,20) == "M0nu_PrintIntegrand_") // 0\nu\beta\beta decay integrand printer, M0nu_PrintIntegrand_${dirout}_${Mode}_${Ec}_${SRC}_${nn}_${ll}_${np}_${lp}_${qa}_${qb} (CP)
      {
         string dirdat; // the directory to hold output data
         int mode; // chosen mode (0 = output to file first, 1 = cout integration result first)
         double Ediff; // the closure energy [MeV]
         string src; // chosen SRC parameters (none, AV18, CD-Bonn, Miller-Spencer)
         int nr,lr,npr,lpr; // the integrand parameters for the RBMEs
         double qmin,qmax; // print the integrand over the interval [qmin,qmax] in units of [MeV]
         string segment;
         stringstream ssopnames1(opname);
         vector<string> seglist1;
         while(getline(ssopnames1,segment,'/'))
         {
           seglist1.push_back(segment);
         }
         for (unsigned int i=1; i<seglist1.size()-1; i++)
         {
           dirdat = dirdat + "/" + seglist1.at(i);
         }
         dirdat = dirdat + "/"; // ${dirout}
         string tempname = opname;
         string dirsub = "_" + dirdat;
         string::size_type nd = dirsub.length();
         for (string::size_type i = tempname.find(dirsub); i != string::npos; i = tempname.find(dirsub))
         {
           tempname.erase(i,nd);
         }
         stringstream ssopnames2(tempname);
         vector<string> seglist2;
         while(getline(ssopnames2,segment,'_'))
         {
           seglist2.push_back(segment);
         }
         mode = stoi(seglist2.at(2)); // ${Mode}
         Ediff = stod(seglist2.at(3)); // ${Ec}
         src = seglist2.at(4); // ${SRC}
         nr = stoi(seglist2.at(5)); // ${nn}
         lr = stoi(seglist2.at(6)); // ${ll}
         npr = stoi(seglist2.at(7)); // ${np}
         lpr = stoi(seglist2.at(8)); // ${lp}
         qmin = stod(seglist2.at(9)); // ${qa}
         qmax = stod(seglist2.at(10)); // ${qb}
         M0nu_PrintIntegrand(modelspace,dirdat,mode,Ediff,src,nr,lr,npr,lpr,qmin,qmax);
      }
//      else if (opnamesplit[0] == "DMNREFT") // point radius density at position r, e.g. rhop1.25
//      {
//        double q;
//        int J;
//        string dmopname = opnamesplit[1];
//        istringstream(opnamesplit[2]) >> q;
//        istringstream(opnamesplit[3]) >> J;
//
//        std::map<string, Operator (*)(ModelSpace&, int, double) > dmop = { {"M",       &DM_NREFT::M},
//                                                                           {"Sigma",   &DM_NREFT::Sigma},
//                                                                           {"Sigmap",  &DM_NREFT::Sigmap},
//                                                                           {"Sigmapp", &DM_NREFT::Sigmapp},
//                                                                           {"Delta",   &DM_NREFT::Delta},
//                                                                           {"Deltap",  &DM_NREFT::Deltap},
//                                                                           {"Phipp",   &DM_NREFT::Phipp},
//                                                                           {"Phitp",   &DM_NREFT::Phitp},
//                                                                           {"Omega",   &DM_NREFT::Omega},
//                                                                         };
//        if ( dmop.find(dmopname) != dmop.end() )
//        {
//        return dmop[dmopname](modelspace, J, q );
//        }
//      }
      else //need to remove from the list
      {
         cout << "Unknown operator: " << opname << endl;
      }
      return Operator();
 
 }

 Operator NumberOp(ModelSpace& modelspace, int n, int l, int j2, int tz2)
 {
   Operator NumOp = Operator(modelspace);
   int indx = modelspace.Index1(n,l,j2,tz2);
   NumOp.ZeroBody = 0;
   NumOp.EraseOneBody();
   NumOp.EraseTwoBody();
   NumOp.OneBody(indx,indx) = 1;
   return NumOp;
 }

 Operator NumberOpAlln(ModelSpace& modelspace, int l, int j2, int tz2)
 {
   Operator NumOp = Operator(modelspace);
   for (int n=0;n<=(modelspace.GetEmax()-l)/2;++n)
   {
     int indx = modelspace.Index1(n,l,j2,tz2);
     NumOp.OneBody(indx,indx) = 1;
   }
   return NumOp;
 }

 double HO_density(int n, int l, double hw, double r)
 {
    double v = M_NUCLEON * hw / (HBARC*HBARC);
    double Norm = pow(v/2.,1.5+l) * M_SQRT2/M_SQRTPI * pow(2,n+2*l+3) * gsl_sf_fact(n) / gsl_sf_doublefact(2*n + 2*l + 1);
    double L = gsl_sf_laguerre_n(n, l+0.5, v*r*r);
    double rho = Norm * pow(r,2*l) * exp(-v * r*r) * L * L;
    return rho;
 }

double HO_Radial_psi(int n, int l, double hw, double r)
{
   double b = sqrt( (HBARC*HBARC) / (hw * M_NUCLEON) );
   double x = r/b;
   double Norm = 2*sqrt( gsl_sf_fact(n) * pow(2,n+l) / M_SQRTPI / gsl_sf_doublefact(2*n+2*l+1) / pow(b,3.0) );
   double L = gsl_sf_laguerre_n(n,l+0.5,x*x);
   double psi = Norm * pow(x,l) * exp(-x*x*0.5) * L;
   return psi;
}

 // Just do the HF transformation
 vector<double> GetOccupationsHF(HartreeFock& hf)
 {
    ModelSpace* modelspace = hf.Hbare.modelspace;
    int norb = modelspace->GetNumberOrbits();
    vector<double> occupation(norb);

    for (int i=0; i<norb; ++i)
    {
      Orbit & oi = modelspace->GetOrbit(i);
      // Get the number operator for orbit i
      Operator N_bare = NumberOp(*modelspace,oi.n,oi.l,oi.j2,oi.tz2);
      // Transform it to the normal-ordered HF basis
      Operator N_NO = hf.TransformToHFBasis(N_bare).DoNormalOrdering();
      occupation[i] = N_NO.ZeroBody;
      cout << oi.n << " " << oi.l << " " << oi.j2 << "/2 " << occupation[i] << endl;
    }
    return occupation;
 }

 // Do the full IMSRG transformation
 vector<double> GetOccupations(HartreeFock& hf, IMSRGSolver& imsrgsolver)
 {
    ModelSpace* modelspace = imsrgsolver.modelspace;
    int norb = modelspace->GetNumberOrbits();
    vector<double> occupation(norb,0);

    for (int i=0; i<norb; ++i)
    {
      Orbit & oi = modelspace->GetOrbit(i);
      // Get the number operator for orbit i
      Operator N_bare = NumberOp(*modelspace,oi.n,oi.l,oi.j2,oi.tz2);
      // Transform it to the normal-ordered HF basis
      Operator N_NO = hf.TransformToHFBasis(N_bare).DoNormalOrdering();
      // Transform to the imsrg evolved basis
      Operator N_final = imsrgsolver.Transform(N_NO);
      occupation[i] = N_final.ZeroBody;
    }
    return occupation;
 }

 vector<double> GetDensity( vector<double>& occupation, vector<double>& R, vector<int>& orbits, ModelSpace& modelspace )
 {
     int nr_steps = R.size();
     double hw = modelspace.GetHbarOmega();
     vector<double> dens(nr_steps,0);
     for (int& i : orbits)
     {
       Orbit & oi = modelspace.GetOrbit(i);
       for (int ir=0;ir<nr_steps;++ir)
       {
          dens[ir] += HO_density(oi.n,oi.l,hw,R[ir]) * occupation[i];
       }
     }
     return dens;
 }

 Operator Single_Ref_1B_Density_Matrix(ModelSpace& modelspace)
 {
    Operator DM(modelspace,0,0,0,2);
//    for (index_t a : modelspace.holes)
    for (auto& a : modelspace.holes)
    {
//       index_t a = it_a.first;
        Orbit& oa = modelspace.GetOrbit(a);
//       DM.OneBody(a,a) = 1.0;
       DM.OneBody(a,a) = oa.j2+1.0;
    }
    return DM;
 }


 double Get_Charge_Density( Operator& DM, double r)
 {
   ModelSpace* modelspace = DM.GetModelSpace();
   double hw = modelspace->GetHbarOmega();
   double rho=0;
   for (index_t i : modelspace->proton_orbits )
   {
      Orbit& oi = modelspace->GetOrbit(i);
      rho += (oi.j2+1) * DM.OneBody(i,i) * HO_density(oi.n, oi.l, hw, r);
   }
   return rho;
 }

/*
 double Get_Charge_Density( Operator& DM, double r)
 {
   ModelSpace* modelspace = DM.GetModelSpace();
   double hw = modelspace->GetHbarOmega();
   double rho=0;
   for (index_t i : modelspace->proton_orbits )
   {
      Orbit& oi = modelspace->GetOrbit(i);
      if (abs(DM.OneBody(i,i))<1e-7) continue;
//      cout << i << " " << (oi.j2+1)*DM.OneBody(i,i) << endl;
//      rho += (oi.j2+1) * DM.OneBody(i,i) * HO_density(oi.n, oi.l, hw, r);
      rho += (1) * DM.OneBody(i,i) * HO_density(oi.n, oi.l, hw, r);
   }
   return rho;
 }
*/


// double Get_Charge_Density( Operator& DM, double r)
// {
//   ModelSpace* modelspace = DM.GetModelSpace();
//   double hw = modelspace->GetHbarOmega();
//   double rho=0;
//   for (index_t i : modelspace->proton_orbits )
//   {
//      Orbit& oi = modelspace->GetOrbit(i);
//      for ( index_t j : DM.OneBodyChannels[{oi.l,oi.j2,oi.tz2}] )
//      {
//        if (abs(DM.OneBody(i,j))<1e-7) continue;
//        Orbit& oj = modelspace->GetOrbit(j);
//        rho += DM.OneBody(i,j) * HO_Radial_psi(oi.n, oi.l, hw, r) * HO_Radial_psi(oj.n, oj.l, hw, r);
//      }
//   }
//   return rho;
// }

Operator KineticEnergy_Op(ModelSpace& modelspace)
{
   Operator T(modelspace);
   int norbits = modelspace.GetNumberOrbits();
   double hw = modelspace.GetHbarOmega();
   for (int a=0;a<norbits;++a)
   {
      Orbit & oa = modelspace.GetOrbit(a);
      T.OneBody(a,a) = 0.5 * hw * (2*oa.n + oa.l +3./2); 
      for ( int b : T.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
      {
         if (b<=a) continue;
         Orbit & ob = modelspace.GetOrbit(b);
         if (oa.n == ob.n+1)
            T.OneBody(a,b) = 0.5 * hw * sqrt( (oa.n)*(oa.n + oa.l +1./2));
         else if (oa.n == ob.n-1)
            T.OneBody(a,b) = 0.5 * hw * sqrt( (ob.n)*(ob.n + ob.l +1./2));
         T.OneBody(b,a) = T.OneBody(a,b);
      }
   }
   return T;
}

/// Relative kinetic energy, includingthe hw/A factor
/// \f[
/// T_{rel} = T - T_{CM}
/// \f]
 Operator Trel_Op(ModelSpace& modelspace)
 {
   Operator Trel( KineticEnergy_Op(modelspace) );
   Trel -= TCM_Op(modelspace);
   return Trel;
 }


/// Center of mass kinetic energy, including the hw/A factor
/// \f[
/// T = \frac{\hbar\omega}{A}\sum_{ij} t_{ij} a^{\dagger}_{i} a_{j} + \frac{\hbar\omega}{A}\frac{1}{4}\sum_{ijkl} t_{ijkl} a^{\dagger}_{i}a^{\dagger}_{j}a_{l}a_{k}
/// \f]
/// with a one-body piece
/// \f[
/// t_{ij} = \frac{1}{\hbar\omega} \left\langle i | T_{12} | j \right\rangle = \frac{1}{2}(2n_i+\ell_i+3/2) \delta_{ij} + \frac{1}{2}\sqrt{n_j(n_j+\ell_j+\frac{1}{2})} \delta_{n_i,n_j-1}\delta_{k_i k_j}
/// \f]
/// where \f$k\f$ labels all quantum numbers other than \f$n\f$ and a two-body piece
/// \f[  
/// t_{ijkl} = \frac{1}{\hbar\omega} \left\langle ij | (T^{CM}_{12} - T^{rel}_{12}) | kl \right\rangle
/// \f]
 Operator TCM_Op(ModelSpace& modelspace)
 {
   double t_start = omp_get_wtime();
   int E2max = modelspace.GetE2max();
   double hw = modelspace.GetHbarOmega();
   int A = modelspace.GetTargetMass();
   Operator TcmOp = Operator(modelspace);
   TcmOp.SetHermitian();
   // One body piece = p**2/(2mA)
   int norb = modelspace.GetNumberOrbits();
   for (int i=0; i<norb; ++i)
   {
      Orbit & oi = modelspace.GetOrbit(i);
      for (int j : TcmOp.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
      {
         Orbit & oj = modelspace.GetOrbit(j);
         if (j<i) continue;
         double tij = 0;
         if (oi.n == oj.n) tij = 0.5*(2*oi.n+oi.l + 1.5) * hw/A;
         else if (oi.n == oj.n-1) tij = 0.5*sqrt(oj.n*(oj.n+oj.l + 0.5)) * hw/A;
         TcmOp.OneBody(i,j) = tij;
         TcmOp.OneBody(j,i) = tij;
      }
   }

   // Two body piece = 2*p1*p2/(2mA) = (Tcm-Trel)/A
   int nchan = modelspace.GetNumberTwoBodyChannels();
   modelspace.PreCalculateMoshinsky();
   #pragma omp parallel for schedule(dynamic,1)  // In order to make this parallel, need to precompute the Moshinsky brackets.
   for (int ch=0; ch<nchan; ++ch)
   {
      TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
      int nkets = tbc.GetNumberKets();
      for (int ibra=0;ibra<nkets;++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         Orbit & oi = modelspace.GetOrbit(bra.p);
         Orbit & oj = modelspace.GetOrbit(bra.q);
         if ( 2*(oi.n+oj.n)+oi.l+oj.l > E2max) continue;
         for (int iket=ibra;iket<nkets;++iket)
         {
            
            Ket & ket = tbc.GetKet(iket);
            Orbit & ok = modelspace.GetOrbit(ket.p);
            Orbit & ol = modelspace.GetOrbit(ket.q);
            if ( 2*(ok.n+ol.n)+ok.l+ol.l > E2max) continue;
//            double p1p2 = Calculate_p1p2(modelspace,bra,ket,tbc.J) * hw/A;
            double p1p2 = Calculate_p1p2(modelspace,bra,ket,tbc.J) / A;
            if (abs(p1p2)>1e-7)
            {
              TcmOp.TwoBody.SetTBME(ch,ibra,iket,p1p2);
            }
         }
      }
   }
   TcmOp.profiler.timer["TCM_Op"] += omp_get_wtime() - t_start;
   return TcmOp;
 }


 // evaluate <bra| p1*p2/2 | ket> 
/// This returns the antisymmetrized J-coupled two body matrix element of \f$ \vec{p}_1 \cdot \vec{p}_2 / (m) \f$.
/// The formula is
/// \f{eqnarray*}{
/// \frac{1}{m}
/// \left \langle a b \right| \vec{p}_1 \cdot \vec{p}_2 \left| c d \right \rangle_J
/// = \frac{1}{\sqrt{(1+\delta_{ab})(1+\delta_{cd})}} &\sum\limits_{LS}
/// \left[ \begin{array}{ccc}
///  \ell_a & s_a & j_a \\
///  \ell_b & s_b & j_b \\
///  L      & S   & J
/// \end{array} \right]
/// \left[ \begin{array}{ccc}
///  \ell_c & s_c & j_c \\
///  \ell_d & s_d & j_d \\
///  L      & S   & J
/// \end{array} \right] & \\
/// & \times \sum\limits_{\substack{N_{ab}N_{cd} \Lambda \\ n_{ab} n_{cd} \lambda}} \mathcal{A}_{abcd}^{\lambda S} \times
/// \left\langle N_{ab}\Lambda n_{ab} \lambda | n_{a} \ell_{a} n_{b} \ell_{b} \right\rangle_{L} &
/// \left\langle N_{cd}\Lambda n_{cd} \lambda | n_{c} \ell_{c} n_{d} \ell_{d} \right\rangle_{L} \\ & \times
/// \left( \left\langle N_{ab}\Lambda | t_{cm} | N_{cd} \Lambda \right \rangle
/// -\left\langle n_{ab}\lambda | t_{rel} | n_{cd} \lambda \right \rangle \right)
/// \f}
/// The antisymmetrization factor \f$ \mathcal{A}_{abcd}^{\lambda S} \f$ ensures that
/// the relative wave function is antisymmetrized.
/// It is given by \f$ \mathcal{A}_{abcd}^{\lambda S} = \left|t_{za}+t_{zc}\right| + \left| t_{za} + t_{zd} \right| (-1)^{\lambda + S + \left|T_z\right|}\f$ .
///
/// The center-of-mass and relative kinetic energies can be found by the same equation as used in the one-body piece of TCM_Op()
 double Calculate_p1p2(ModelSpace& modelspace, Ket & bra, Ket & ket, int J)
 {
   Orbit & oa = modelspace.GetOrbit(bra.p);
   Orbit & ob = modelspace.GetOrbit(bra.q);
   Orbit & oc = modelspace.GetOrbit(ket.p);
   Orbit & od = modelspace.GetOrbit(ket.q);

   int na = oa.n;
   int nb = ob.n;
   int nc = oc.n;
   int nd = od.n;

   int la = oa.l;
   int lb = ob.l;
   int lc = oc.l;
   int ld = od.l;

   double ja = oa.j2/2.0;
   double jb = ob.j2/2.0;
   double jc = oc.j2/2.0;
   double jd = od.j2/2.0;

   int fab = 2*na + 2*nb + la + lb;
   int fcd = 2*nc + 2*nd + lc + ld;
   // p1*p2 only connects kets with delta N = 0,1 ==> delta E = 0,2
   if (abs(fab-fcd)>2 or abs(fab-fcd)%2 >0 ) return 0; 

   double sa,sb,sc,sd;
   sa=sb=sc=sd=0.5;

   double p1p2=0;

   // First, transform to LS coupling using 9j coefficients
   for (int Lab=abs(la-lb); Lab<= la+lb; ++Lab)
   {
     for (int Sab=0; Sab<=1; ++Sab)
     {
       if ( abs(Lab-Sab)>J or Lab+Sab<J) continue;

       double njab = NormNineJ(la,sa,ja, lb,sb,jb, Lab,Sab,J);
       if (njab == 0) continue;
       int Scd = Sab;
       int Lcd = Lab;
       double njcd = NormNineJ(lc,sc,jc, ld,sd,jd, Lcd,Scd,J);
       if (njcd == 0) continue;

       // Next, transform to rel / com coordinates with Moshinsky tranformation
       for (int N_ab=0; N_ab<=fab/2; ++N_ab)  // N_ab = CoM n for a,b
       {
         for (int Lam_ab=0; Lam_ab<= fab-2*N_ab; ++Lam_ab) // Lam_ab = CoM l for a,b
         {
           int Lam_cd = Lam_ab; // tcm and trel conserve lam and Lam, ie relative and com orbital angular momentum
           for (int lam_ab=(fab-2*N_ab-Lam_ab)%2; lam_ab<= (fab-2*N_ab-Lam_ab); lam_ab+=2) // lam_ab = relative l for a,b
           {
              if (Lab<abs(Lam_ab-lam_ab) or Lab>(Lam_ab+lam_ab) ) continue;
              // factor to account for antisymmetrization

              int asymm_factor = (abs(bra.op->tz2+ket.op->tz2) + abs(bra.op->tz2+ket.oq->tz2)*modelspace.phase( lam_ab + Sab ))/ 2;
              if ( asymm_factor ==0 ) continue;

              int lam_cd = lam_ab; // tcm and trel conserve lam and Lam
              int n_ab = (fab - 2*N_ab-Lam_ab-lam_ab)/2; // n_ab is determined by energy conservation

              double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);

              if (abs(mosh_ab)<1e-8) continue;

              for (int N_cd=max(0,N_ab-1); N_cd<=N_ab+1; ++N_cd) // N_cd = CoM n for c,d
              {
                int n_cd = (fcd - 2*N_cd-Lam_cd-lam_cd)/2; // n_cd is determined by energy conservation
                if (n_cd < 0) continue;
                if  (n_ab != n_cd and N_ab != N_cd) continue;

                double mosh_cd = modelspace.GetMoshinsky(N_cd,Lam_cd,n_cd,lam_cd,nc,lc,nd,ld,Lcd);
                if (abs(mosh_cd)<1e-8) continue;

                double tcm = 0;
                double trel = 0;
                if (n_ab == n_cd)
                {
                  if      (N_ab == N_cd)   tcm = (2*N_ab+Lam_ab+1.5);
                  else if (N_ab == N_cd+1) tcm = sqrt(N_ab*( N_ab+Lam_ab+0.5));
                  else if (N_ab == N_cd-1) tcm = sqrt(N_cd*( N_cd+Lam_ab+0.5));
                }
                if (N_ab == N_cd)
                {
                  if      (n_ab == n_cd)   trel = (2*n_ab+lam_ab+1.5);
                  else if (n_ab == n_cd+1) trel = sqrt(n_ab*( n_ab+lam_ab+0.5));
                  else if (n_ab == n_cd-1) trel = sqrt(n_cd*( n_cd+lam_cd+0.5));
                }
                double prefactor = njab * njcd * mosh_ab * mosh_cd * asymm_factor;
                p1p2 += (tcm-trel) * prefactor ;

              } // N_cd
           } // lam_ab
         } // Lam_ab
       } // N_ab

     } // Sab
   } // Lab

   // normalize. The 0.5 comes from t ~ 0.5 * (N+3/2) hw
   p1p2 *= 0.5*modelspace.GetHbarOmega() / sqrt((1.0+bra.delta_pq())*(1.0+ket.delta_pq()));
   return p1p2 ;

 }



////////////////////////////////////////////////////////////////////////////
/////////////  IN PROGRESS...  doesn't work yet, and for now it's slower.  /
////////////////////////////////////////////////////////////////////////////

 void Calculate_p1p2_all(Operator& OpIn)
 {
   ModelSpace* modelspace = OpIn.GetModelSpace();
//   modelspace->PreCalculateMoshinsky();
   for ( int ch : modelspace->SortedTwoBodyChannels )
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      int parity = tbc.parity;
//      int Tz = tbc.Tz;
      arma::mat& MatJJ = OpIn.TwoBody.GetMatrix(ch);
      int nkets_JJ = tbc.GetNumberKets();

      // Find the maximum oscillator energy for this channel
      Ket& ketlast = tbc.GetKet( tbc.GetNumberKets()-1 );
      int emax_ket = 2*ketlast.op->n + 2*ketlast.oq->n + ketlast.op->l + ketlast.oq->l;

      vector<array<int,6>> JacobiBasis;  // L,S,N,Lambda,n,lambda
      for (int L=max(J-1,0); L<=J+1; ++L)
      {
       for ( int S=abs(J-L); S<=1; ++S)
       {
        for ( int N=0; N<=emax_ket/2; ++N )
        {
         for ( int Lambda=0; Lambda<=(emax_ket-2*N); ++Lambda)
         {
          for ( int lambda=abs(L-Lambda)+(L+parity)%2; lambda<=min(Lambda+L,emax_ket-2*N-Lambda); lambda+=2)
          {
           for ( int n =0; n<=(emax_ket-2*N-Lambda-lambda)/2; ++n)
           {
             JacobiBasis.push_back({L,S,N,Lambda,n,lambda});
           }
          }
         }
        }
       }
      }


      int nkets_Jacobi = JacobiBasis.size();
      arma::mat MatJacobi(nkets_Jacobi,nkets_Jacobi);
      arma::mat Trans(nkets_Jacobi,nkets_JJ);
      int n_nonzero = 0;
      for (int iJJ=0; iJJ<nkets_JJ; ++iJJ)
      {
        Ket & ket = tbc.GetKet(iJJ);
        int la = ket.op->l;
        int lb = ket.oq->l;
        int na = ket.op->n;
        int nb = ket.oq->n;
        double ja = ket.op->j2*0.5;
        double jb = ket.oq->j2*0.5;
//        int ta = ket.op->n;
//        int tb = ket.oq->n;
        for (int iJac=0; iJac<nkets_Jacobi; ++iJac)
        {
          int L      = JacobiBasis[iJac][0];
          int S      = JacobiBasis[iJac][1];
          int N      = JacobiBasis[iJac][2];
          int Lambda = JacobiBasis[iJac][3];
          int n      = JacobiBasis[iJac][4];
          int lambda = JacobiBasis[iJac][5];

//          int Asym = 1; // Fix this...
          double ninej = NormNineJ(la,0.5,ja,lb,0.5,jb,L,S,J);
          if (abs(ninej)<1e-6) continue;
          double mosh = modelspace->GetMoshinsky(N,Lambda,n,lambda,na,la,nb,lb,L);
          if (abs(mosh)<1e-6) continue;
          Trans(iJac,iJJ) = ninej * mosh;
          n_nonzero += 1;
           
        }
      }
      for (int i=0; i<nkets_Jacobi; ++i)
      {
        int Li      = JacobiBasis[i][0];
        int Si      = JacobiBasis[i][1];
        int Ni      = JacobiBasis[i][2];
        int Lambdai = JacobiBasis[i][3];
        int ni      = JacobiBasis[i][4];
        int lambdai = JacobiBasis[i][5];
        for (int j=0; j<nkets_Jacobi; ++j)
        {
          int Lj      = JacobiBasis[j][0];
          int Sj      = JacobiBasis[j][1];
          int Nj      = JacobiBasis[j][2];
          int Lambdaj = JacobiBasis[j][3];
          int nj      = JacobiBasis[j][4];
          int lambdaj = JacobiBasis[j][5];
          if(Li!=Lj or Si!=Sj or Lambdai!=Lambdaj or lambdai!=lambdaj) continue;
          double tcm = 0;
          double trel = 0;
          if (ni == nj)
          {
            if      (Ni == Nj)   tcm = (2*Ni+Lambdai+1.5);
            else if (Ni == Nj+1) tcm = sqrt(Ni*( Ni+Lambdai+0.5));
            else if (Ni == Nj-1) tcm = sqrt(Nj*( Nj+Lambdai+0.5));
          }
          if (Ni == Nj)
          {
            if      (ni == nj)   tcm = (2*ni+lambdai+1.5);
            else if (ni == nj+1) tcm = sqrt(ni*( ni+lambdai+0.5));
            else if (ni == nj-1) tcm = sqrt(nj*( nj+lambdai+0.5));
          }
          MatJacobi(i,j) = tcm - trel;
        }
      }


      MatJJ = Trans.t() * MatJacobi * Trans;
      cout << "ch = " << ch << "   size of JJ basis = " << nkets_JJ << "  size of Jacobi Basis = " << nkets_Jacobi << "   nonzero matrix elements = " << n_nonzero << endl;

   }
 }

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////





// Center of mass R^2, in units of fm^2
/// Returns
/// \f[ 
/// R^{2}_{CM} = \left( \frac{1}{A}\sum_{i}\vec{r}_{i}\right)^2 =
/// \frac{1}{A^2} \left( \sum_{i}r_{i}^{2} + 2\sum_{i<j}\vec{r}_i\cdot\vec{r}_j  \right)
/// \f]
/// evaluated in the oscillator basis.


 Operator R2CM_Op(ModelSpace& modelspace)
 {
//   Operator R2cmOp = Operator(modelspace);
   Operator R2cmOp = RSquaredOp(modelspace);

//   unsigned int norb = modelspace.GetNumberOrbits();
//   double oscillator_b2 = (HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega());
//   for (unsigned int i=0; i<norb; ++i)
//   {
//      Orbit & oi = modelspace.GetOrbit(i);
//      for (auto j : R2cmOp.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
//      {
//         if (j<i) continue;
//         Orbit & oj = modelspace.GetOrbit(j);
//         double rij = 0;
//         if (oi.n == oj.n)        rij = (2*oi.n+oi.l + 1.5);
//         else if (oi.n == oj.n-1) rij = -sqrt(oj.n*(oj.n+oj.l + 0.5));
////         R2cmOp.OneBody(i,j) = rij;
////         R2cmOp.OneBody(j,i) = rij;
//         R2cmOp.OneBody(i,j) = rij * oscillator_b2;
//         R2cmOp.OneBody(j,i) = rij * oscillator_b2;
//      }
//   }

   int nchan = modelspace.GetNumberTwoBodyChannels();
   modelspace.PreCalculateMoshinsky();
   #pragma omp parallel for schedule(dynamic,1) 
   for (int ch=0; ch<nchan; ++ch)
   {
      TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
      int nkets = tbc.GetNumberKets();
      for (int ibra=0;ibra<nkets;++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         for (int iket=ibra;iket<nkets;++iket)
         {
            Ket & ket = tbc.GetKet(iket);
            // factor of 2 comes from limiting sum to i<j. Otherwise it would be r1*r2 + r2*r1.
            double mat_el = 2*Calculate_r1r2(modelspace,bra,ket,tbc.J); 
             
            R2cmOp.TwoBody.SetTBME(ch,ibra,iket,mat_el);
            R2cmOp.TwoBody.SetTBME(ch,iket,ibra,mat_el);
         }
      }
   }
//   double hw = modelspace.GetHbarOmega();
   int A = modelspace.GetTargetMass();
//   return R2cmOp * (HBARC*HBARC/M_NUCLEON/hw)/(A*A);
   return R2cmOp /(A*A);
 }



// Intrinsic point proton radius squared
/// Returns
/// \f[ 
/// R_p^{2} = \frac{1}{Z} \sum_{p}\left(\vec{r}_{p}-\vec{R}_{CM}\right)^2 =
/// R^2_{CM} + \frac{A-2}{AZ} \sum_{p}r_{p}^{2} - \frac{4}{AZ}\sum_{i<j}\vec{r}_i\cdot\vec{r}_j  \right)
/// \f]
/// evaluated in the oscillator basis.
 Operator Rp2_corrected_Op(ModelSpace& modelspace, int A, int Z)
 {
   if (Z==0) return 0.0*KineticEnergy_Op(modelspace);
   return R2CM_Op(modelspace) + (A-2.0)/(A*Z)*R2_1body_Op(modelspace,"proton")
                                   - 4./(A*Z)*R2_2body_Op(modelspace,"proton")
                                   + 1./Z * RpSpinOrbitCorrection(modelspace);
 }

 Operator Rn2_corrected_Op(ModelSpace& modelspace, int A, int Z)
 {
   if (Z==A) return 0.0*KineticEnergy_Op(modelspace);
   return R2CM_Op(modelspace) + (A-2.0)/(A*(A-Z))*R2_1body_Op(modelspace,"neutron")
                                   - 4./(A*(A-Z))*R2_2body_Op(modelspace,"neutron");
 }

 Operator Rm2_corrected_Op(ModelSpace& modelspace, int A, int Z)
 {
   return (1./A)*RSquaredOp(modelspace) - R2CM_Op(modelspace)  ;
 }



 // Evaluate <bra | r1*r2 | ket>, omitting the factor (hbar * omega) /(m * omega^2)
/// Returns the normalized, anti-symmetrized, J-coupled, two-body matrix element of \f$  \vec{r}_1\cdot\vec{r}_2 \f$.
/// Calculational details are similar to Calculate_p1p2().
 double Calculate_r1r2(ModelSpace& modelspace, Ket & bra, Ket & ket, int J)
 {
   Orbit & oa = modelspace.GetOrbit(bra.p);
   Orbit & ob = modelspace.GetOrbit(bra.q);
   Orbit & oc = modelspace.GetOrbit(ket.p);
   Orbit & od = modelspace.GetOrbit(ket.q);

   int na = oa.n;
   int nb = ob.n;
   int nc = oc.n;
   int nd = od.n;

   int la = oa.l;
   int lb = ob.l;
   int lc = oc.l;
   int ld = od.l;

   double ja = oa.j2*0.5;
   double jb = ob.j2*0.5;
   double jc = oc.j2*0.5;
   double jd = od.j2*0.5;

   int fab = 2*na + 2*nb + la + lb;
   int fcd = 2*nc + 2*nd + lc + ld;
   if (abs(fab-fcd)%2 >0) return 0; // p1*p2 only connects kets with delta N = 0,1
   if (abs(fab-fcd)>2) return 0; // p1*p2 only connects kets with delta N = 0,1

   double sa,sb,sc,sd;
   sa=sb=sc=sd=0.5;

   double r1r2=0;

   // First, transform to LS coupling using 9j coefficients
   for (int Lab=abs(la-lb); Lab<= la+lb; ++Lab)
   {
     for (int Sab=0; Sab<=1; ++Sab)
     {
       if ( abs(Lab-Sab)>J or Lab+Sab<J) continue;

       double njab = NormNineJ(la,sa,ja, lb,sb,jb, Lab,Sab,J);
       if (njab == 0) continue;
       int Scd = Sab;
       int Lcd = Lab;
       double njcd = NormNineJ(lc,sc,jc, ld,sd,jd, Lcd,Scd,J);
       if (njcd == 0) continue;

       // Next, transform to rel / com coordinates with Moshinsky tranformation
       for (int N_ab=0; N_ab<=fab/2; ++N_ab)  // N_ab = CoM n for a,b
       {
         for (int Lam_ab=0; Lam_ab<= fab-2*N_ab; ++Lam_ab) // Lam_ab = CoM l for a,b
         {
           int Lam_cd = Lam_ab; // tcm and trel conserve lam and Lam, ie relative and com orbital angular momentum
           for (int lam_ab=(fab-2*N_ab-Lam_ab)%2; lam_ab<= (fab-2*N_ab-Lam_ab); lam_ab+=2) // lam_ab = relative l for a,b
           {
              if (Lab<abs(Lam_ab-lam_ab) or Lab>(Lam_ab+lam_ab) ) continue;
              // factor to account for antisymmetrization

              int asymm_factor = (abs(bra.op->tz2+ket.op->tz2) + abs(bra.op->tz2+ket.oq->tz2)*modelspace.phase( lam_ab + Sab ))/ 2;
              if ( asymm_factor ==0 ) continue;

              int lam_cd = lam_ab; // tcm and trel conserve lam and Lam
              int n_ab = (fab - 2*N_ab-Lam_ab-lam_ab)/2; // n_ab is determined by energy conservation

              double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);

              if (abs(mosh_ab)<1e-8) continue;

              for (int N_cd=max(0,N_ab-1); N_cd<=N_ab+1; ++N_cd) // N_cd = CoM n for c,d
              {
                int n_cd = (fcd - 2*N_cd-Lam_cd-lam_cd)/2; // n_cd is determined by energy conservation
                if (n_cd < 0) continue;
                if  (n_ab != n_cd and N_ab != N_cd) continue;

                double mosh_cd = modelspace.GetMoshinsky(N_cd,Lam_cd,n_cd,lam_cd,nc,lc,nd,ld,Lcd);
                if (abs(mosh_cd)<1e-8) continue;

                double r2cm = 0;
                double r2rel = 0;
                if (n_ab == n_cd)
                {
                  if      (N_ab == N_cd)   r2cm = (2*N_ab+Lam_ab+1.5);
                  else if (N_ab == N_cd+1) r2cm = -sqrt(N_ab*( N_ab+Lam_ab+0.5));
                  else if (N_ab == N_cd-1) r2cm = -sqrt(N_cd*( N_cd+Lam_ab+0.5));
                }
                if (N_ab == N_cd)
                {
                  if      (n_ab == n_cd)   r2rel = (2*n_ab+lam_ab+1.5);
                  else if (n_ab == n_cd+1) r2rel = -sqrt(n_ab*( n_ab+lam_ab+0.5));
                  else if (n_ab == n_cd-1) r2rel = -sqrt(n_cd*( n_cd+lam_cd+0.5));
                }
                double prefactor = njab * njcd * mosh_ab * mosh_cd * asymm_factor;
                r1r2 += (r2cm-r2rel) * prefactor;

              } // N_cd
           } // lam_ab
         } // Lam_ab
       } // N_ab

     } // Sab
   } // Lab

   double oscillator_b2 = (HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega());
   // normalize and give dimension of length
   // the 0.5 comes from the virial theorem, <V> = <T> = 1/2 (N+3/2)hw
   r1r2 *= oscillator_b2 *0.5 / sqrt((1.0+bra.delta_pq())*(1.0+ket.delta_pq()));
//   r1r2 *= 1.0 / sqrt((1.0+bra.delta_pq())*(1.0+ket.delta_pq()));
   return r1r2 ;

 }




/// Center of mass Hamiltonian
/// \f{eqnarray*}{
/// H_{CM} &= T_{CM} + \frac{1}{2} Am\omega^2 R^2 \\
///        &= T_{CM} + \frac{1}{2b^2} AR^2 \hbar\omega
/// \f}
 Operator HCM_Op(ModelSpace& modelspace)
 {
   double hw = modelspace.GetHbarOmega();
   int A = modelspace.GetTargetMass();
//   double oscillator_b2 = HBARC*HBARC/M_NUCLEON/hw;
//   Operator HcmOp = TCM_Op(modelspace) + R2CM_Op(modelspace) * (0.5*A * hw / oscillator_b2);
   Operator HcmOp = TCM_Op(modelspace) + 0.5*A*M_NUCLEON*hw*hw/HBARC/HBARC * R2CM_Op(modelspace) ;
   cout << "HcmOp: first 1b element = " << HcmOp.OneBody(0,0) << endl;
   return HcmOp;
 }



/// Returns
/// \f[ r^2 = \sum_{i} r_{i}^2 \f]
///
Operator RSquaredOp(ModelSpace& modelspace)
{
   Operator r2 = Operator(modelspace);
   r2.OneBody.zeros();
   unsigned int norbits = modelspace.GetNumberOrbits();
   double oscillator_b2 = HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega();
   for (unsigned int a=0;a<norbits;++a)
   {
      Orbit & oa = modelspace.GetOrbit(a);
      r2.OneBody(a,a) = (2*oa.n + oa.l +1.5); 
      for ( unsigned int b : r2.OneBodyChannels.at({oa.l, oa.j2, oa.tz2}) )
      {
        if ( b < a ) continue;
        Orbit & ob = modelspace.GetOrbit(b);
        {
           if (oa.n == ob.n+1)
              r2.OneBody(a,b) = -sqrt( (oa.n)*(oa.n + oa.l +0.5));
           else if (oa.n == ob.n-1)
              r2.OneBody(a,b) = -sqrt( (ob.n)*(ob.n + ob.l +0.5));
           r2.OneBody(b,a) = r2.OneBody(a,b);
        }
      }
   }
//   r2.OneBody *= (HBARC*HBARC/M_NUCLEON/hw);
   r2.OneBody *= oscillator_b2;
   return r2;
}


 Operator R2_p1_Op(ModelSpace& modelspace)
 {
   return R2_1body_Op(modelspace,"proton");
 }

/// One-body part of the proton charge radius operator.
/// Returns
/// \f[ 
/// \hat{R}^{2}_{p1} = \sum_{i} e_{i}{r}_i^2
/// \f]
 Operator R2_1body_Op(ModelSpace& modelspace,string option)
 {
   Operator r2(modelspace);
   double oscillator_b = (HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega());

   auto orbitlist = modelspace.proton_orbits;
   if (option == "neutron") orbitlist = modelspace.neutron_orbits;
   else if (option == "matter")  orbitlist.insert(orbitlist.end(),modelspace.neutron_orbits.begin(),modelspace.neutron_orbits.end());
   else if (option != "proton") cout << "!!! WARNING. BAD OPTION "  << option << " FOR imsrg_util::R2_p1_Op !!!" << endl;
 
   for (unsigned int a : orbitlist )
   {
      Orbit & oa = modelspace.GetOrbit(a);
      r2.OneBody(a,a) = (2*oa.n + oa.l +1.5); 
      for ( unsigned int b : r2.OneBodyChannels.at({oa.l, oa.j2, oa.tz2}) )
      {
        if ( b < a ) continue;
        Orbit & ob = modelspace.GetOrbit(b);
        {
           if (oa.n == ob.n+1)
              r2.OneBody(a,b) = -sqrt( (oa.n)*(oa.n + oa.l +0.5));
           else if (oa.n == ob.n-1)
              r2.OneBody(a,b) = -sqrt( (ob.n)*(ob.n + ob.l +0.5));
           r2.OneBody(b,a) = r2.OneBody(a,b);
        }
      }
   }
   r2.OneBody *= oscillator_b;
   return r2;
 }

 Operator R2_p2_Op(ModelSpace& modelspace)
 {
   return R2_2body_Op(modelspace,"proton");
 }

/// Two-body part of the proton charge radius operator.
/// Returns
/// \f[ 
/// \hat{R}^{2}_{p2} = \sum_{i\neq j} e_{i}\vec{r}_i\cdot\vec{r}_j 
/// \f]
/// evaluated in the oscillator basis.
 Operator R2_2body_Op(ModelSpace& modelspace,string option)
 {
   Operator Rp2Op(modelspace,0,0,0,2);
//   double oscillator_b = (HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega());

   int nchan = modelspace.GetNumberTwoBodyChannels();
   modelspace.PreCalculateMoshinsky();
   #pragma omp parallel for schedule(dynamic,1) 
   for (int ch=0; ch<nchan; ++ch)
   {
      TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
      int Tz = tbc.Tz;
      if (option=="proton" and Tz > 0) continue; // don't bother with nn channel
      if (option=="neutron" and Tz < 0) continue; // don't bother with pp channel
      int nkets = tbc.GetNumberKets();
      for (int ibra=0;ibra<nkets;++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         double prefactor = 1; // factor to account for double counting in pn channel.
         if (Tz==0 and (option=="proton" or option=="neutron")) prefactor = 0.5;
//         if (option=="proton" and bra.op->tz2>0) continue;
//         else if (option=="neutron" and bra.op->tz2<0) continue;
         if (option!="matter" and option!="proton" and option!="neutron") cout << "!!! WARNING. BAD OPTION "  << option << " FOR imsrg_util::R2_p2_Op !!!" << endl;
         for (int iket=ibra;iket<nkets;++iket)
         {
            Ket & ket = tbc.GetKet(iket);
//            double mat_el = Calculate_r1r2(modelspace,bra,ket,tbc.J) * oscillator_b ; 
            double mat_el = Calculate_r1r2(modelspace,bra,ket,tbc.J) * prefactor; 
            Rp2Op.TwoBody.SetTBME(ch,ibra,iket,mat_el);
            Rp2Op.TwoBody.SetTBME(ch,iket,ibra,mat_el);
         }
      }
   }
   return Rp2Op;
 }


Operator ProtonDensityAtR(ModelSpace& modelspace, double R)
{
  Operator Rho(modelspace,0,0,0,2);
  for ( auto i : modelspace.proton_orbits)
  {
    Orbit& oi = modelspace.GetOrbit(i);
    Rho.OneBody(i,i) = HO_density(oi.n,oi.l,modelspace.GetHbarOmega(),R);
  }
  return Rho;
}

Operator NeutronDensityAtR(ModelSpace& modelspace, double R)
{
  Operator Rho(modelspace,0,0,0,2);
  for ( auto i : modelspace.neutron_orbits)
  {
    Orbit& oi = modelspace.GetOrbit(i);
    Rho.OneBody(i,i) = HO_density(oi.n,oi.l,modelspace.GetHbarOmega(),R);
  }
  return Rho;
}


Operator RpSpinOrbitCorrection(ModelSpace& modelspace)
{
  Operator dr_so(modelspace,0,0,0,2);
  double M2 = M_NUCLEON*M_NUCLEON/(HBARC*HBARC);
  int norb = modelspace.GetNumberOrbits();
  for (int i=0;i<norb;i++)
  {
    Orbit& oi = modelspace.GetOrbit(i);
    double mu_i = oi.tz2<0 ? 1.79 : -1.91;
    int kappa = oi.j2 < 2*oi.l ? oi.l : -(oi.l+1);
    dr_so.OneBody(i,i) = -mu_i/M2*(kappa+1);
  }
  return dr_so;
}

// Electric monopole operator
/// Returns
/// \f[ r_{e}^2 = \sum_{i} e_{i} r_{i}^2 \f]
///
Operator E0Op(ModelSpace& modelspace)
{
   Operator e0(modelspace);
   e0.EraseZeroBody();
   e0.OneBody.zeros();
//   unsigned int norbits = modelspace.GetNumberOrbits();
   double hw = modelspace.GetHbarOmega();
   for (unsigned int a : modelspace.proton_orbits)
   {
      Orbit & oa = modelspace.GetOrbit(a);
      e0.OneBody(a,a) = (2*oa.n + oa.l +1.5); 
      for (unsigned int b : e0.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
      {
        if (b<=a) continue;
        Orbit & ob = modelspace.GetOrbit(b);
        {
           if (oa.n == ob.n+1)
              e0.OneBody(a,b) = -sqrt( (oa.n)*(oa.n + oa.l +0.5));
           else if (oa.n == ob.n-1)
              e0.OneBody(a,b) = -sqrt( (ob.n)*(ob.n + ob.l +0.5));
           e0.OneBody(b,a) = e0.OneBody(a,b);
        }
      }
   }
   e0.OneBody *= (HBARC*HBARC/M_NUCLEON/hw);
   return e0;
}

struct FBCIntegrandParameters{int n; int l; double hw;};

double FBCIntegrand(double x, void *p)
{
  struct FBCIntegrandParameters * params = (struct FBCIntegrandParameters *)p;
  return x*HO_density(params->n, params->l, params->hw, x);
}

Operator FourierBesselCoeff(ModelSpace& modelspace, int nu, double R, vector<index_t> index_list)
{
  Operator a_nu(modelspace,0,0,0,2);
  double omega = nu * M_PI / R; // coefficient of sine function, i.e. sin(omega*x)
  double L = R; // range of integration
  size_t n = 20; // number of bisections
  gsl_integration_qawo_table * table = gsl_integration_qawo_table_alloc (omega, L, GSL_INTEG_SINE, n);
  gsl_integration_workspace * workspace = gsl_integration_workspace_alloc ( n );

  const double epsabs = 1e-5; // absolute error
  const double epsrel = 1e-5; // relative error
  const size_t limit = n; // maximum number of subintervals (maybe should be different?)
  const double start = 0.0; // lower limit on integration range
  double result;
  double  abserr;
  gsl_function F;
  F.function = &FBCIntegrand;

  for (auto i : index_list )
  {
    Orbit& oi = modelspace.GetOrbit(i);
    struct FBCIntegrandParameters params = {oi.n, oi.l, modelspace.GetHbarOmega()};
    F.params = &params;
    //int status = gsl_integration_qawo (&F, start, epsabs, epsrel, limit, workspace, table, &result, &abserr);
    gsl_integration_qawo (&F, start, epsabs, epsrel, limit, workspace, table, &result, &abserr);
    a_nu.OneBody(i,i) = M_PI*M_PI/R/R/R * R/nu/M_PI*(result);
    cout << "orbit,nu = " << i << "," << nu << "  => " << a_nu.OneBody(i,i) << "  from " << result << " (" << abserr << ")" << endl;
  }
  return a_nu;
}



/// Returns the \f$ T^{2} \f$ operator
 Operator Isospin2_Op(ModelSpace& modelspace)
 {
   Operator T2 = Operator(modelspace,0,0,0,2);
   T2.OneBody.diag().fill(0.75);

   for (int ch=0; ch<T2.nChannels; ++ch)
   {
     TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
     arma::mat& TB = T2.TwoBody.GetMatrix(ch);
     // pp,nn:  2<t2.t1> = 1/(2(1+delta_ab)) along diagonal
     if (abs(tbc.Tz) == 1)
     {
        TB.diag().fill(0.5); // pp,nn TBME's
        for (int ibra=0;ibra<tbc.GetNumberKets(); ++ibra)
        {
           Ket& bra = tbc.GetKet(ibra);
           if (bra.p == bra.q)
           {
             TB(ibra,ibra) /= 2.;
           }
        }
     }
     else if (tbc.Tz == 0)
     {
        for (int ibra=0;ibra<tbc.GetNumberKets(); ++ibra)
        {
           Ket& bra = tbc.GetKet(ibra);
           Orbit& oa = modelspace.GetOrbit(bra.p);
           Orbit& ob = modelspace.GetOrbit(bra.q);
           for (int iket=ibra;iket<tbc.GetNumberKets(); ++iket)
           {
             Ket& ket = tbc.GetKet(iket);
             Orbit& oc = modelspace.GetOrbit(ket.p);
             Orbit& od = modelspace.GetOrbit(ket.q);
             if (oa.j2==oc.j2 and oa.n==oc.n and oa.l==oc.l
               and ob.j2==od.j2 and ob.n==od.n and ob.l==od.l )
             {
               // tz1 tz2 case
               if( oa.tz2 == oc.tz2 and ob.tz2==od.tz2)
               {
                 TB(ibra,iket) -= 0.5;
               }
               // t+ t- case
               if( oa.tz2 == od.tz2 and ob.tz2==oc.tz2)
               {
                 TB(ibra,iket) += 1.0;
               }
               // if a==b==c==d, we need to consider the exchange term
               if (oa.j2==ob.j2 and oa.n==ob.n and oa.l==ob.l)
               {
                  int phase = bra.Phase(tbc.J);
                  // tz1 tz2 case
                  if( oa.tz2 == oc.tz2 and ob.tz2==od.tz2)
                  {
                    TB(ibra,iket) += phase * 1.0;
                  }
                  // t+ t- case
                  if( oa.tz2 == od.tz2 and ob.tz2==oc.tz2)
                  {
                    TB(ibra,iket) -= phase * 0.5;
                  }
               }
               TB(iket,ibra) = TB(ibra,iket); // hermitian
             }

           }
        }
     }
   }
   return T2;
 }



  /// Returns a reduced electric multipole operator with units \f$ e\f$ fm\f$^{\lambda} \f$
  /// See Suhonen eq. (6.23)
  Operator ElectricMultipoleOp(ModelSpace& modelspace, int L)
  {
    Operator EL(modelspace, L,0,L%2,2);
    double bL = pow( HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega(),0.5*L); // b^L where b=sqrt(hbar/mw)
    for (int i : modelspace.proton_orbits)
    {
      Orbit& oi = modelspace.GetOrbit(i);
      double ji = 0.5*oi.j2;
      for ( int j : EL.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
      {
        if (j<i) continue;
        Orbit& oj = modelspace.GetOrbit(j);
        double jj = 0.5*oj.j2;
        double r2int = RadialIntegral(oi.n,oi.l,oj.n,oj.l,L) * bL ;
        EL.OneBody(i,j) = modelspace.phase(jj+L-0.5) * sqrt( (2*ji+1)*(2*jj+1)*(2*L+1)/4./3.1415926) * AngMom::ThreeJ(ji,jj, L, 0.5, -0.5,0) * r2int;
        EL.OneBody(j,i) = modelspace.phase((oi.j2-oj.j2)/2) * EL.OneBody(i,j);
      }
    }
    return EL;
  }

  /// Returns a reduced electric multipole operator with units \f$ e\f$ fm\f$^{\lambda} \f$
  /// See Suhonen eq. (6.23)
  Operator NeutronElectricMultipoleOp(ModelSpace& modelspace, int L)
  {
    Operator EL(modelspace, L,0,L%2,2);
    double bL = pow( HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega(),0.5*L); // b^L where b=sqrt(hbar/mw)
    for (int i : modelspace.neutron_orbits)
    {
      Orbit& oi = modelspace.GetOrbit(i);
      double ji = 0.5*oi.j2;
      for ( int j : EL.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
      {
        if (j<i) continue;
        Orbit& oj = modelspace.GetOrbit(j);
        double jj = 0.5*oj.j2;
        double r2int = RadialIntegral(oi.n,oi.l,oj.n,oj.l,L) * bL ;
        EL.OneBody(i,j) = modelspace.phase(jj+L-0.5) * sqrt( (2*ji+1)*(2*jj+1)*(2*L+1)/4./3.1415926) * AngMom::ThreeJ(ji,jj, L, 0.5, -0.5,0) * r2int;
        EL.OneBody(j,i) = modelspace.phase((oi.j2-oj.j2)/2) * EL.OneBody(i,j);
      }
    }
    return EL;
  }

  /// Returns a reduced magnetic multipole operator with units \f$ \mu_{N}\f$ fm\f$ ^{\lambda-1} \f$
  Operator MagneticMultipoleOp(ModelSpace& modelspace, int L)
  {
    return MagneticMultipoleOp_pn(modelspace,L,"both");
  }

  /// Returns a reduced magnetic multipole operator with units \f$ \mu_{N}\f$ fm\f$ ^{\lambda-1} \f$
  /// This version allows for the selection of just proton or just neutron contributions, or both.
  /// See Suhonen eq. (6.24)
  Operator MagneticMultipoleOp_pn(ModelSpace& modelspace, int L, string pn)
  {
    double bL = pow( HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega(),0.5*(L-1));
    Operator ML(modelspace, L,0,(L+1)%2,2);
    if (L<1)
    {
      cout << "A magnetic monopole operator??? Setting it to zero..." << endl;
      return ML;
    }
    int norbits = modelspace.GetNumberOrbits();
    for (int i=0; i<norbits; ++i)
    {
      Orbit& oi = modelspace.GetOrbit(i);
      if (pn=="proton" and oi.tz2>0) continue;
      if (pn=="neutron" and oi.tz2<0) continue;
      double gl = oi.tz2<0 ? 1.0 : 0.0;
      double gs = oi.tz2<0 ? 5.586 : -3.826;
      double ji = 0.5*oi.j2;
      for ( int j : ML.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
      {
        if (j<i) continue;
        Orbit& oj = modelspace.GetOrbit(j);
        double jj = 0.5*oj.j2;
        // multiply radial integral by b^L-1 = (hbar/mw)^L-1/2
        double r2int = RadialIntegral(oi.n,oi.l,oj.n,oj.l,L-1) * bL;
        double kappa =  modelspace.phase(oi.l+ji+0.5) * (ji+0.5)  +  modelspace.phase(oj.l+jj+0.5) * (jj+0.5);
        ML.OneBody(i,j) = modelspace.phase(jj+L-0.5) * sqrt( (2*ji+1)*(2*jj+1)*(2*L+1)/4./3.1415926) * AngMom::ThreeJ(ji, jj,L, 0.5,-0.5,0)
                        * (L - kappa) *(gl*(1+kappa/(L+1.))-0.5*gs )  * r2int ;
        ML.OneBody(j,i) = modelspace.phase((oi.j2-oj.j2)/2) * ML.OneBody(i,j);
      }
    }
    return ML;
  }



  /// Returns a reduced electric multipole operator with units \f$ e\f$ fm\f$^{\lambda} \f$
  /// See Suhonen eq. (6.23)
  Operator IntrinsicElectricMultipoleOp(ModelSpace& modelspace, int L)
  {
    Operator EL(modelspace, L,0,L%2,2);
    double bL = pow( HBARC*HBARC/(0.5*M_NUCLEON)/modelspace.GetHbarOmega(),0.5*L); // b^L where b=sqrt(hbar/mw)

    int emax = modelspace.GetEmax();
    int norb_rel = (2*emax+1)*(2*emax+2)/2;
    arma::mat ME_rel(norb_rel, norb_rel);
    for (int nrela=0; 2*nrela<2*emax; nrela+=1)
    {
     for (int lrela=0; lrela<2*emax-2*nrela; lrela+=1)
     {
       int ia = (2*nrela + lrela)*(2*nrela+lrela+1)/2 + lrela;
       for (int lrelb=abs(lrela-L); lrelb<=lrela+L; lrelb+=2)
       {
         for (int nrelb=0; 2*nrelb<=2*emax-lrelb; nrelb+=1)
         {
           int ib = (2*nrelb + lrelb)*(2*nrelb+lrelb+1)/2 + lrelb;
           double r2int = RadialIntegral(nrela,lrela,nrelb,lrelb,L) * bL ;
           double angint = modelspace.phase(lrela) * sqrt( (2*lrela+1)*(2*lrelb+1)*(2*L+1)/4./3.1415926) * AngMom::ThreeJ(lrela, L, lrelb, 0,0,0);
           ME_rel(ia,ib) = angint * r2int;
           printf("< %d %d || E2 || %d %d > = %f   angular: %f \n",nrela,lrela,nrelb,lrelb,ME_rel(ia,ib), angint );
         }
       }
     }
    }



    // Moshinsky transformation to lab frame bra and ket
    // lab frame: < na la ja  nb lb jb  Jab || EL || nc lc jc  nd ld jd  Jcd>
    for ( auto& itmat : EL.TwoBody.MatEl )
    {
      TwoBodyChannel& tbc_bra = modelspace.GetTwoBodyChannel(itmat.first[0]);
      TwoBodyChannel& tbc_ket = modelspace.GetTwoBodyChannel(itmat.first[1]);
      if (tbc_bra.Tz>=0) continue;
      int Jab = tbc_bra.J;
      int Jcd = tbc_ket.J;
      for (int ibra=0; ibra<tbc_bra.GetNumberKets(); ++ibra)
      {
       Ket& bra = tbc_bra.GetKet(ibra);
       int na = bra.op->n;
       int la = bra.op->l;
       double ja = bra.op->j2*0.5;
       int nb = bra.oq->n;
       int lb = bra.oq->l;
       double jb = bra.oq->j2*0.5;
       int rho_ab = 2*na+2*nb+la+lb;
       for (int iket=0; iket<tbc_ket.GetNumberKets(); ++iket)
       {
        Ket& ket = tbc_bra.GetKet(ibra);
        int nc = ket.op->n;
        int lc = ket.op->l;
        double jc = ket.op->j2*0.5;
        int nd = ket.oq->n;
        int ld = ket.oq->l;
        double jd = ket.oq->j2*0.5;
        int rho_cd = 2*nc+2*nd+lc+ld;

        double labME = 0;
        
        for (int Lab=max(abs(la-lb),abs(Jab-1)); Lab<=min(la+lb,Jab+1); Lab+=1)
        {
         for (int Lcd=max(abs(lc-ld),abs(Jcd-1))+(Lab+L)%2; Lcd<=min(lc+ld,Jcd+1); Lcd+=2)
         {
          for (int S=max(abs(Jab-Lab),abs(Jcd-Lcd)); S<=min(1,min(Jab+Lab,Jcd+Lcd)); S+=1 )
          {
            double njab = AngMom::NormNineJ( la, 0.5, ja, lb, 0.5, jb, Lab, S, Jab);
            double njcd = AngMom::NormNineJ( lc, 0.5, jc, ld, 0.5, jd, Lcd, S, Jcd);
            if (abs(njab)<1e-8 or abs(njcd)<1e-8) continue;
            for ( int nab=0; nab<=rho_ab; nab+=1)
            {
             for ( int lab=0; lab<=rho_ab-2*nab; lab+=1)
             {
              int iab = (2*nab + lab)*(2*nab+lab+1)/2 + lab;
              for (int Lam=abs(lab-Lab); Lam<=min(lab+Lab,rho_ab-2*nab-lab); Lam+=2)
              {
               int N = rho_ab - 2*nab - lab - Lam;
               double moshab = modelspace.GetMoshinsky( N, Lam, nab, lab, na,la,nb,lb,Lab);
               if (abs(moshab)<1e-8) continue;
               for ( int ncd=0; ncd<=rho_cd; ncd+=1)
               {
                for ( int lcd=abs(Lam-Lcd)+(Lcd+Lam)%2; lcd<=min(rho_cd-2*ncd,Lab+Lcd); lcd+=2)
                {
                 int icd = (2*ncd + lcd)*(2*ncd+lcd+1)/2 + lcd;
                 double moshcd = modelspace.GetMoshinsky( N, Lam, ncd, lcd, nc,lc,nd,ld,Lcd);
                 labME += njab * njcd * moshab * moshcd * ME_rel(iab,icd);
                }
               }
              }
             }
            }
          }
         }
        }
        itmat.second(ibra,iket) = labME;
       }
      }
    }
    
    cout << "done with intrinsic EL. one body = " <<  endl;
    return EL;
  }



/// Evaluate the radial integral \f[
/// \tilde{\mathcal{R}}^{\lambda}_{ab} = \int_{0}^{\infty} dx \tilde{g}_{n_a\ell_a}(x)x^{\lambda+2}\tilde{g}_{n_b\ell_b}(x)
/// \f]
/// where \f$ \tilde{g}(x) \f$ is the radial part of the harmonic oscillator wave function with unit oscillator length \f$ b=1 \f$
/// and \f$ x = r/b \f$.
/// To obtain the radial integral for some other oscillator length, multiply by \f$ b^{\lambda} \f$.
/// This implementation uses eq (6.41) from Suhonen.
/// Note this is only valid for \f$ \ell_a+\ell_b+\lambda\f$ = even.
/// If \f$ \ell_a+\ell_b+\lambda\f$ is odd, RadialIntegral_RpowK() is called.
  double RadialIntegral(int na, int la, int nb, int lb, int L)
  {
    if ((la+lb+L)%2!=0) return RadialIntegral_RpowK(na,la,nb,lb,L);
    int tau_a = max((lb-la+L)/2,0);
    int tau_b = max((la-lb+L)/2,0);
    int sigma_min = max(max(na-tau_a,nb-tau_b),0);
    int sigma_max = min(na,nb);
  
    double term1 = AngMom::phase(na+nb) * gsl_sf_fact(tau_a)*gsl_sf_fact(tau_b) * sqrt(gsl_sf_fact(na)*gsl_sf_fact(nb)
                   / (gsl_sf_gamma(na+la+1.5)*gsl_sf_gamma(nb+lb+1.5) ) );
    double term2 = 0;
    for (int sigma=sigma_min; sigma<=sigma_max; ++sigma)
    {
      term2 += gsl_sf_gamma(0.5*(la+lb+L)+sigma+1.5) / (gsl_sf_fact(sigma)*gsl_sf_fact(na-sigma)*gsl_sf_fact(nb-sigma)*gsl_sf_fact(sigma+tau_a-na)*gsl_sf_fact(sigma+tau_b-nb) );
    }
    return term1*term2;
  
  }



 double RadialIntegral_RpowK(int na, int la, int nb, int lb, int k)
 {
   double I = 0;
   int pmin = (la+lb)/2;
   int pmax = pmin + na + nb;
   for (int p=pmin;p<=pmax;++p)
   {
      I += TalmiB(na,la,nb,lb,p) * TalmiI(p,k);
   }
   return I;
 }

/// General Talmi integral for a potential r**k
/// 1/gamma(p+3/2) * 2*INT dr r**2 r**2p r**k exp(-r**2/b**2)
/// This is valid for (2p+3+k) > 0. The Gamma function diverges for non-positive integers.
 double TalmiI(int p, double k)
 {
//   return gsl_sf_gamma(p+1.5+0.5*k) / gsl_sf_gamma(p+1.5);
   return boost::math::tgamma_ratio(p+1.5+0.5*k, p+1.5);
 }

/// Calculate B coefficient for Talmi integral. Formula given in Brody and Moshinsky
/// "Tables of Transformation Brackets for Nuclear Shell-Model Calculations"
 double TalmiB(int na, int la, int nb, int lb, int p)
 {
   if ( (la+lb)%2>0 ) return 0;
   
   int q = (la+lb)/2;
//   double B1 = AngMom::phase(p-q) * gsl_sf_fact(2*p+1)/gsl_sf_fact(p)/pow(2,(na+nb))
   double B1 = AngMom::phase(p-q) * boost::math::tgamma_ratio(2*p+1, p) / pow(2,(na+nb))
              * sqrt( boost::math::tgamma_ratio(na, na+la)*boost::math::tgamma_ratio(nb, nb+lb)
                   * boost::math::factorial<double>(2*na+2*la+1) * boost::math::factorial<double>(2*nb+2*lb+1) );
   
   double B2 = 0;
   int kmin = max(0, p-q-nb);
   int kmax = min(na, p-q);
   for (int k=kmin;k<=kmax;++k)
   {
//      B2  += gsl_sf_fact(la+k) * gsl_sf_fact(p-int((la-lb)/2)-k)
//             / ( gsl_sf_fact(k) * gsl_sf_fact(2*la+2*k+1) * gsl_sf_fact(na-k) * gsl_sf_fact(2*p-la+lb-2*k+1) )
//              / ( gsl_sf_fact(nb - p + q + k) * gsl_sf_fact(p-q-k) );
      B2  += boost::math::tgamma_ratio(la+k,k) * boost::math::tgamma_ratio(p-int((la-lb)/2)-k, 2*p-la+lb-2*k+1)
             / (  boost::math::factorial<double>(2*la+2*k+1) * boost::math::factorial<double>(na-k)  
                * boost::math::factorial<double>(nb - p + q + k) * boost::math::factorial<double>(p-q-k) );
   }
   return B1 * B2;
 }

  Operator AllowedFermi_Op(ModelSpace& modelspace)
  {
    Operator Fermi(modelspace,0,1,0,2);
    Fermi.SetHermitian();
    int norbits = modelspace.GetNumberOrbits();
    for (int i=0; i<norbits; ++i)
    {
      Orbit& oi = modelspace.GetOrbit(i);
      for (int j : Fermi.OneBodyChannels[{oi.l,oi.j2,oi.tz2}] )
      {
        Orbit& oj = modelspace.GetOrbit(j);
        if (oi.n!=oj.n or oi.tz2 == oj.tz2) continue;
        Fermi.OneBody(i,j) = sqrt(oi.j2+1.0);  // Reduced matrix element
      }
    }
    return Fermi;
  }

/// Note that there is a literature convention to include the 1/sqrt(Lambda) factor
/// in the reduced matrix element rather than in the expression involving the sum
/// over one-body densities (see footnote on pg 165 of Suhonen).
/// I do not follow this convention, and instead produce the reduced matrix element
///  \f[ \langle f \| \sigma \tau_{\pm} \| i \rangle \f]
///
  Operator AllowedGamowTeller_Op(ModelSpace& modelspace)
  {
    Operator GT(modelspace,1,1,0,2);
    GT.SetHermitian();
    int norbits = modelspace.GetNumberOrbits();
    for (int i=0; i<norbits; ++i)
    {
      Orbit& oi = modelspace.GetOrbit(i);
      for (int j : GT.OneBodyChannels[{oi.l,oi.j2,oi.tz2}] )
      {
        Orbit& oj = modelspace.GetOrbit(j);
        if (oi.n!=oj.n or oi.l != oj.l or oi.tz2==oj.tz2) continue;
        double sixj = modelspace.GetSixJ(0.5,0.5,1.0,oj.j2/2.,oi.j2/2.,oi.l);
        double M_gt = 2 * modelspace.phase(oi.l+oi.j2/2.0+1.5) * sqrt((oi.j2+1)*(oj.j2+1)) * sqrt(1.5) * sixj;
        GT.OneBody(i,j) = M_gt;
      }
    }
    return GT;
  }



 /// Pauli spin operator \f[ \langle f \| \sigma \| i \rangle \f]
 Operator Sigma_Op(ModelSpace& modelspace)
 {
   return Sigma_Op_pn(modelspace,"both");
 }

 /// Pauli spin operator \f[ \langle f \| \sigma \| i \rangle \f]
 Operator Sigma_Op_pn(ModelSpace& modelspace, string pn)
 {
   Operator Sig(modelspace,1,0,0,2);
   Sig.SetHermitian();
   int norbits = modelspace.GetNumberOrbits();
   for (int i=0; i<norbits; ++i)
   {
     Orbit& oi = modelspace.GetOrbit(i);
      if (pn=="proton" and oi.tz2>0) continue;
      if (pn=="neutron" and oi.tz2<0) continue;
      for (int j : Sig.OneBodyChannels[{oi.l,oi.j2,oi.tz2}] )
      {
        Orbit& oj = modelspace.GetOrbit(j);
        if ((oi.n!=oj.n) or (oi.l != oj.l) or (oi.tz2!=oj.tz2)) continue;
        double sixj = modelspace.GetSixJ(0.5,0.5,1.0,oj.j2/2.,oi.j2/2.,oi.l);
        double M_sig = 2 * modelspace.phase(oi.l+oi.j2/2.0+1.5) * sqrt((oi.j2+1)*(oj.j2+1)) * sqrt(1.5) * sixj;
        Sig.OneBody(i,j) = M_sig;
      }
   } 
   return Sig;
 }

  void Reduce(Operator& X)
  {
    ModelSpace* modelspace = X.GetModelSpace();
    int norbits = modelspace->GetNumberOrbits();
    for (int i=0; i<norbits; ++i)
    {
      Orbit& oi = modelspace->GetOrbit(i);
      for ( int j : X.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
      {
         X.OneBody(i,j) *= sqrt(oi.j2+1.);
      }
    }

    for ( auto& itmat : X.TwoBody.MatEl )
    {
      int ch_bra = itmat.first[0];
      int J = modelspace->GetTwoBodyChannel(ch_bra).J;
      itmat.second *= sqrt(2*J+1.);
    }
  }

  void UnReduce(Operator& X)
  {
    ModelSpace* modelspace = X.GetModelSpace();
    int norbits = modelspace->GetNumberOrbits();
    for (int i=0; i<norbits; ++i)
    {
      Orbit& oi = modelspace->GetOrbit(i);
      for ( int j : X.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
      {
         X.OneBody(i,j) /= sqrt(oi.j2+1.);
      }
    }

    for ( auto& itmat : X.TwoBody.MatEl )
    {
      int ch_bra = itmat.first[0];
      int J = modelspace->GetTwoBodyChannel(ch_bra).J;
      itmat.second /= sqrt(2*J+1.);
    }

  }


  void SplitUp(Operator& OpIn, Operator& OpLow, Operator& OpHi, int ecut)
  {
    ModelSpace* modelspace = OpIn.GetModelSpace();
    OpLow = OpIn;
    int norbits = modelspace->GetNumberOrbits();
    int ncut = 0;
    for (int i=0; i<norbits; ++i)
    {
      Orbit& oi = modelspace->GetOrbit(i);
      if ( (2*oi.n + oi.l) > ecut)
      {
         ncut = i;
         break;
      }
    }

    for (int i=ncut; i<norbits; ++i)
    {
      for (int j=0; j<norbits; ++j)
      {
         OpLow.OneBody(i,j) = 0;
         OpLow.OneBody(j,i) = 0;
      }
    }

    for (auto& itmat : OpLow.TwoBody.MatEl)
    {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel( itmat.first[0] );
      int nkets = tbc.GetNumberKets();
      for (int ibra=0;ibra<nkets;++ibra)
      {
       Ket& bra = tbc.GetKet(ibra);
       for (int iket=ibra;iket<nkets;++iket)
       {
         Ket& ket = tbc.GetKet(iket);
         if ( bra.p>ncut or bra.q>ncut or ket.p>ncut or ket.q>ncut)
         {
           itmat.second(ibra,iket) = 0;
         }
       }
      }
    }
    OpHi = OpIn - OpLow;
  }


  Operator RadialOverlap(ModelSpace& modelspace)
  {
     Operator OVL(modelspace,0,1,0,1);
     index_t norb = modelspace.GetNumberOrbits();
     for (index_t i=0; i<norb; ++i)
     {
       Orbit& oi = modelspace.GetOrbit(i);
       for (index_t j=0; j<norb; ++j)
       {
         Orbit& oj = modelspace.GetOrbit(j);
         OVL.OneBody(i,j) = RadialIntegral(oi.n, oi.l, oj.n, oj.l, 0 ); // This is not quite right. Only works for li+lj=even.
       } 
     }
     return OVL;
  }


  Operator LdotS_Op(ModelSpace& modelspace)
  {
    Operator LdS(modelspace,0,0,0,2);
    index_t norbits = modelspace.GetNumberOrbits();
    for (index_t i=0;i<norbits;++i)
    {
      Orbit& oi = modelspace.GetOrbit(i);
      double ji = 0.5*oi.j2;
      int li = oi.l;
      LdS.OneBody(i,i) = 0.5* (ji*(ji+1) -li*(li+1)-0.75);
    }
    return LdS;
  }


 Operator LCM_Op(ModelSpace& modelspace)
 {

   Operator LCM(modelspace, 1,0,0,2);
/*

   unsigned int norb = modelspace.GetNumberOrbits();
   for (unsigned int i=0; i<norb; ++i)
   {
      Orbit & oi = modelspace.GetOrbit(i);
      for (auto j : LCM.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
      {
         if (j<i) continue;
         Orbit & oj = modelspace.GetOrbit(j);
         if (oi.n != oj.n or oi.l != oj.l) continue ;
         double Lij = modelspace.phase(oi.l + (oj.j2+3)/2) * sqrt(oi.j2+1)*sqrt(oj.j2+1)  * sqrt(oi.l*(oi.l+1)*(2*oi.l+1)) * modelspace.GetSixJ(oi.l,oj.l,1,oj.j2/2.,oi.j2/2.,0.5);
         LCM.OneBody(i,j) = Lij;
         LCM.OneBody(j,i) = modelspace.phase((oi.j2-oj.j2)/2) * Lij;
      }
   }

   int nchan = modelspace.GetNumberTwoBodyChannels();
   modelspace.PreCalculateMoshinsky();

*/
  return LCM;
 }




////////////////////////////// via NDBD.hh/cc (CP) //////////////////////////////


//vvDGTvv
/// This is the Double-Gamow-Teller (DGT) TBME
/// it was coded up by me, ie) Charlie Payne (CP)
/// these TBMEs are reduced by default, and I haven't added in any profiling yet
/// NOTE: it has not been tested or benchmarked as of 05/23/18 (although I'm pretty sure it should work)
  Operator DGT_Op(ModelSpace& modelspace)
  {
    Operator DGT_TBME(modelspace,0,2,0,2); // NOTE: from the constructor -- Operator::Operator(ModelSpace& ms, int Jrank, int Trank, int p, int part_rank)
    DGT_TBME.SetHermitian(); // it should be Hermitian
    // create the TBMEs of DGT
    // auto loops over the TBME channels and such
    cout<<"calculating DGT TBMEs..."<<endl;
    for (auto& itmat : DGT_TBME.TwoBody.MatEl)
    {
      int chbra = itmat.first[0]; // grab the channel count from auto
      int chket = itmat.first[1]; // " " " " " "
      TwoBodyChannel& tbc_bra = modelspace.GetTwoBodyChannel(chbra); // grab the two-body channel
      TwoBodyChannel& tbc_ket = modelspace.GetTwoBodyChannel(chket); // " " " " "
      int nbras = tbc_bra.GetNumberKets(); // get the number of bras
      int nkets = tbc_ket.GetNumberKets(); // get the number of kets
      double J = tbc_bra.J; // NOTE: by construction, J := J_ab == J_cd := J'
      double Jhat = sqrt(2*J + 1); // the hat factor of J
      #pragma omp parallel for schedule(dynamic,100)
      for (int ibra=0; ibra<nbras; ibra++)
      {
        Ket& bra = tbc_bra.GetKet(ibra); // get the final state = <ab|
        int ia = bra.p; // get the integer label a
        int ib = bra.q; // get the integer label b
        Orbit& oa = modelspace.GetOrbit(ia); // get the <a| state orbit
        Orbit& ob = modelspace.GetOrbit(ib); // get the <b| state prbit
        for (int iket=0; iket<nkets; iket++)
        {
          Ket& ket = tbc_ket.GetKet(iket); // get the initial state = |cd>
          int ic = ket.p; // get the integer label c
          int id = ket.q; // get the integer label d
          Orbit& oc = modelspace.GetOrbit(ic); // get the |c> state orbit
          Orbit& od = modelspace.GetOrbit(id); // get the |d> state orbit
          int la = oa.l; // this is just...
          int lb = ob.l;
          int lc = oc.l;
          int ld = od.l;
          double ja = oa.j2/2.0;
          double jb = ob.j2/2.0;
          double jc = oc.j2/2.0;
          double jd = od.j2/2.0; // ...for convenience
          double sumLS = 0; // for the final matrix elements
          double sumLSas = 0; // (anti-symmetric part)
          for (int S=0; S<=1; S++) // sum over total spin...
          {
            int Seval = 2*S*(S + 1) - 3; // eigenvalue of S
            double sumL = 0; // for the summation over L
            double sumLas = 0; // (anti-symmetric part)
            for (int L = abs(la-lb); L <= la+lb; L++) // ...and sum over angular momentum coupled to l_a and l_b
            {
              if ((abs(lc-ld) <= L) and (L <= lc+ld))
              {
                double tempLS = (2*L + 1)*(2*S + 1); // just for efficiency, only used in the three lines below
                double normab = sqrt(tempLS*(2*ja + 1)*(2*jb + 1)); // normalization factor for the 9j-symbol out front
                double nNJab = normab*modelspace.GetNineJ(la,lb,L,0.5,0.5,S,ja,jb,J); // the normalized 9j-symbol for a and b
                double normcd = sqrt(tempLS*(2*jc + 1)*(2*jd + 1)); // normalization factor for the second 9j-symbol
                double nNJcd = normcd*modelspace.GetNineJ(lc,ld,L,0.5,0.5,S,jc,jd,J); // the normalized 9j-symbol for c and d
                double nNJdc = normcd*modelspace.GetNineJ(ld,lc,L,0.5,0.5,S,jd,jc,J); // (anti-symmetric part)
                sumL += nNJab*nNJcd; // perform the sum over L
                sumLas += nNJab*nNJdc; // (anti-symmetric part)
              } // end of if: |lc - ld| <= L <= lc + ld
            } // end of for-loop over: L
            sumLS += Seval*sumL; // perform the sum over S
            sumLSas += Seval*sumLas; // (anti-symmetric part)
          } // end of for-loop over: S
          double Mtbme = cpNorm(ia,ib)*cpNorm(ic,id)*Jhat*(sumLS - modelspace.phase(jc + jd - J)*sumLSas); // compute the final matrix element, anti-symmetrize
          DGT_TBME.TwoBody.SetTBME(chbra,chket,ibra,iket,Mtbme); // set the two-body matrix elements (TBME) to Mtbme
        } // end of for-loop over: iket
      } // end of for-loop over: ibra
    } // end of for-loop over: auto
    cout<<"...done calculating DGT TBMEs"<<endl;
    return DGT_TBME;
  }
//^^DGT^^


//vvGT/Fvv
/// This is the M^{0\nu} TBME, for either GT or F decay type
/// it was coded up by me, ie) Charlie Payne (CP)
/// from my thesis, I employed Equations: (4.71), (4.72), (4.46), (4.40), (4.50), and (4.73)
  Operator M0nuGTF_adpt_Op(ModelSpace& modelspace, string dirname, string Type, string reduced, double Ediff, string src, string barcode)
  {
    string meth = "PSH"; // set the method for the integration over dq and dr, "PSH" or "JE"
    double t_start, t_start_tbme, t_start_omp; // profiling (v)
    t_start = omp_get_wtime(); // profiling (s)
    // run through the initial set-up routine
    // we'll call: M0nuHeader, M0nuGTF_TBME.SetHermitian, modelspace.PreCalculateMoshinsky, ndbd.PreCalcIntegrals
    if (reduced == "NR") {}
    else if (reduced == "R") {}
    else
    {
      cout<<"ERROR 7338233-1: invalid option for reduced = "<<reduced<<endl;
      exit(1);
    }
//--GT/F--
    double hw = modelspace.GetHbarOmega(); // oscillator basis frequency [MeV]
    int e2max = modelspace.GetE2max(); // 2*emax
    Operator M0nuGTF_TBME(modelspace,0,2,0,2); // NOTE: from the constructor -- Operator::Operator(ModelSpace& ms, int Jrank, int Trank, int p, int part_rank)
    NDBD ndbd(hw,e2max,Type,Ediff,src,meth); // call the Constructor NDBD::NDBD(...)
    cout<<"     reduced            =  "<<reduced<<endl;
    M0nuGTF_TBME.SetHermitian(); // it should be Hermitian
    int Anuc = modelspace.GetTargetMass(); // the mass number for the desired nucleus
    const double Rnuc = (ndbd.r0/HBARC)*pow(Anuc,1.0/3.0); // the nuclear radius [MeV^-1]
    const double prefact = 2*(2*Rnuc)/PI; // factor in-front of M0nu TBME, extra global 2 for nutbar (as confirmed by benchmarking with Ca48 NMEs) [MeV^-1]
    M0nuHeader(ndbd,dirname,reduced,Anuc,src,Rnuc,prefact,barcode); // create a header file
    modelspace.PreCalculateMoshinsky(); // pre-calculate the needed Moshinsky brackets, for efficiency
    ndbd.PreCalcIntegrals(); // pre-calculate the needed integrals over dq and dr, for efficiency
    M0nuGTF_TBME.profiler.timer["M0nuGTF_1_sur"] += omp_get_wtime() - t_start; // profiling (r)
//--GT/F--
    // create the TBMEs of M0nu
    // auto loops over the TBME channels and such
    cout<<"calculating M0nu TBMEs..."<<endl;
    t_start_tbme = omp_get_wtime(); // profiling (s)
    for (auto& itmat : M0nuGTF_TBME.TwoBody.MatEl)
    {
      int chbra = itmat.first[0]; // grab the channel count from auto
      int chket = itmat.first[1]; // " " " " " "
      TwoBodyChannel& tbc_bra = modelspace.GetTwoBodyChannel(chbra); // grab the two-body channel
      TwoBodyChannel& tbc_ket = modelspace.GetTwoBodyChannel(chket); // " " " " "
      int nbras = tbc_bra.GetNumberKets(); // get the number of bras
      int nkets = tbc_ket.GetNumberKets(); // get the number of kets
      double J = tbc_bra.J; // NOTE: by construction, J := J_ab == J_cd := J'
      double Jhat; // set below based on "reduced" variable
//--GT/F--
      if (reduced == "NR")
      {
        Jhat = 1.0; // for non-reduced elements, to compare with JE
      }
      else //if (reduced == "R")
      {
        Jhat = sqrt(2*J + 1); // the hat factor of J
      }
      t_start_omp = omp_get_wtime(); // profiling (s)
      #pragma omp parallel for schedule(dynamic,100) // need to do: PreCalculateMoshinsky() and PreCalcIntegrals() [above] and then "#pragma omp critical" [below]
      for (int ibra=0; ibra<nbras; ibra++)
      {
        Ket& bra = tbc_bra.GetKet(ibra); // get the final state = <ab|
//--GT/F--
        int ia = bra.p; // get the integer label a
        int ib = bra.q; // get the integer label b
        Orbit& oa = modelspace.GetOrbit(ia); // get the <a| state orbit
        Orbit& ob = modelspace.GetOrbit(ib); // get the <b| state prbit
        for (int iket=0; iket<nkets; iket++)
        {
          Ket& ket = tbc_ket.GetKet(iket); // get the initial state = |cd>
          int ic = ket.p; // get the integer label c
          int id = ket.q; // get the integer label d
          Orbit& oc = modelspace.GetOrbit(ic); // get the |c> state orbit
          Orbit& od = modelspace.GetOrbit(id); // get the |d> state orbit
          int na = oa.n; // this is just...
          int nb = ob.n;
          int nc = oc.n;
          int nd = od.n;
//--GT/F--
          int la = oa.l;
          int lb = ob.l;
          int lc = oc.l;
          int ld = od.l;
          double ja = oa.j2/2.0;
          double jb = ob.j2/2.0;
          double jc = oc.j2/2.0;
          double jd = od.j2/2.0; // ...for convenience
          double sumLS = 0; // for the Bessel's Matrix Elemets (BMEs)
          double sumLSas = 0; // (anti-symmetric part)
          for (int S=0; S<=1; S++) // sum over total spin...
          {
//--GT/F--
            int Seval = 1.0; // eigenvalue of S, only used for "GT" type M0nu
            if (Type == "GT")
            {
              Seval = 2*S*(S + 1) - 3;
            }
            for (int L = abs(la-lb); L <= la+lb; L++) // ...and sum over angular momentum coupled to l_a and l_b
            {
              if ((abs(lc-ld) <= L) and (L <= lc+ld))
              {
                double bulk = 0; // this is a bulk product which will only be worth calculating if the Moshinsky brackets are non-zero
                double bulkas = 0; // (anti-symmetric part)
                double sumMT = 0; // for the Moshinsky transformation
                double sumMTas = 0; // (anti-symmetric part)
                int eps_ab = 2*na + la + 2*nb + lb; // for conservation of energy in the Moshinsky brackets
                int tempmaxnr = floor((eps_ab - L)/2.0); // just for the limits below
//--GT/F--
                for (int nr = 0; nr <= tempmaxnr; nr++)
                {
                  int eps_cd = 2*nc + lc + 2*nd + ld; // for conservation of energy in the Moshinsky brackets
                  double npr = ((eps_cd - eps_ab)/2.0) + nr; // via Equation (4.73) of my thesis
                  double intnpr = npr - floor(npr); // to check if npr is an integer
                  if ((npr >= 0) and (intnpr == 0))
                  {
                    int tempmaxNcom = tempmaxnr - nr; // just for the limits below
                    for (int Ncom = 0; Ncom <= tempmaxNcom; Ncom++)
                    {
                      int tempminlr = ceil((eps_ab - L)/2.0) - (nr + Ncom); // just for the limits below
                      int tempmaxlr = floor((eps_ab + L)/2.0) - (nr + Ncom); // " " " " "
                      for (int lr = tempminlr; lr <= tempmaxlr; lr++)
                      {
//--GT/F--
                        int Lam = eps_ab - 2*(nr + Ncom) - lr; // via Equation (4.73) of my thesis
                        double tempLS = (2*L + 1)*(2*S + 1); // just for efficiency, only used in the three lines below
                        double normab = sqrt(tempLS*(2*ja + 1)*(2*jb + 1)); // normalization factor for the 9j-symbol out front
                        double nNJab = normab*modelspace.GetNineJ(la,lb,L,0.5,0.5,S,ja,jb,J); // the normalized 9j-symbol out front
                        double normcd = sqrt(tempLS*(2*jc + 1)*(2*jd + 1)); // normalization factor for the second 9j-symbol
                        double nNJcd = normcd*modelspace.GetNineJ(lc,ld,L,0.5,0.5,S,jc,jd,J); // the second normalized 9j-symbol
                        double nNJdc = normcd*modelspace.GetNineJ(ld,lc,L,0.5,0.5,S,jd,jc,J); // (anti-symmetric part)
                        bulk = Seval*nNJab*nNJcd; // bulk product of the above
                        bulkas = Seval*nNJab*nNJdc; // (anti-symmetric part)
                        double Df,Di,asDi,integral; // these *will* be set below
                        #pragma omp critical
                        {
                          Df = modelspace.GetMoshinsky(Ncom,Lam,nr,lr,na,la,nb,lb,L); // Ragnar has -- double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);
                          Di = modelspace.GetMoshinsky(Ncom,Lam,npr,lr,nc,lc,nd,ld,L); // " " " "
                          asDi = modelspace.GetMoshinsky(Ncom,Lam,npr,lr,nd,ld,nc,lc,L); // (anti-symmetric part)
//--GT/F--
                          integral = ndbd.GetIntegral(nr,lr,npr,lr); // grab the pre-calculated integral wrt dq and dr from the IntList of the NDBD class
                        }
                        sumMT += Df*Di*integral; // perform the Moshinsky transformation
                        sumMTas += Df*asDi*integral; // (anti-symmetric part)
                      } // end of for-loop over: lr
                    } // end of for-loop over: Ncom
                  } // end of if: npr \in \Nat_0
                } // end of for-loop over: nr
                sumLS += bulk*sumMT; // perform the LS-coupling sum
                sumLSas += bulkas*sumMTas; // (anti-symmetric part)
              } // end of if: |lc - ld| <= L <= lc + ld
            } // end of for-loop over: L
          } // end of for-loop over: S
//--GT/F--
          double Mtbme = cpNorm(ia,ib)*cpNorm(ic,id)*prefact*Jhat*(sumLS - modelspace.phase(jc + jd - J)*sumLSas); // compute the final matrix element, anti-symmetrize
          M0nuGTF_TBME.TwoBody.SetTBME(chbra,chket,ibra,iket,Mtbme); // set the two-body matrix elements (TBME) to Mtbme
        } // end of for-loop over: iket
      } // end of for-loop over: ibra
      M0nuGTF_TBME.profiler.timer["M0nuGTF_3_omp"] += omp_get_wtime() - t_start_omp; // profiling (r)
    } // end of for-loop over: auto
    cout<<"...done calculating M0nu TBMEs"<<endl;
    M0nuGTF_TBME.profiler.timer["M0nuGTF_2_tbme"] += omp_get_wtime() - t_start_tbme; // profiling (r)
    M0nuGTF_TBME.profiler.timer["M0nuGTF_adpt_Op"] += omp_get_wtime() - t_start; // profiling (r)
    return M0nuGTF_TBME;
  }
//^^GT/F^^


//vvTvv
/// This is the M^{0\nu} TBME, for T decay type
/// it was coded up by me, ie) Charlie Payne (CP)
/// from my thesis, I employed Equations: (4.71), (4.72), (4.31), and (4.74)
  Operator M0nuT_adpt_Op(ModelSpace& modelspace, string dirname, string reduced, double Ediff, string src, string barcode)
  {
    string meth = "PSH"; // set the method for the integration over dq and dr, "PSH" or "JE"
    string Type = "T"; // set the decay type, in this case it's for the Tensor type
    double t_start, t_start_tbme, t_start_omp; // profiling (v)
    t_start = omp_get_wtime(); // profiling (s)
    // run through the initial set-up routine
    // we'll call: M0nuHeader, M0nuT_TBME.SetHermitian, modelspace.PreCalculateMoshinsky, ndbd.PreCalcT6j, ndbd.PreCalcIntegrals
    if (reduced == "NR") {}
    else if (reduced == "R") {}
    else
    {
      cout<<"ERROR 7338233-2: invalid option for reduced = "<<reduced<<endl;
      exit(1);
    }
//--T--
    double hw = modelspace.GetHbarOmega(); // oscillator basis frequency [MeV]
    int e2max = modelspace.GetE2max(); // 2*emax
    Operator M0nuT_TBME(modelspace,0,2,0,2); // NOTE: from the constructor -- Operator::Operator(ModelSpace& ms, int Jrank, int Trank, int p, int part_rank)
    NDBD ndbd(hw,e2max,Type,Ediff,src,meth); // call the Constructor NDBD::NDBD(...)
    cout<<"     reduced            =  "<<reduced<<endl;
    M0nuT_TBME.SetHermitian(); // it should be Hermitian
    int Anuc = modelspace.GetTargetMass(); // the mass number for the desired nucleus
    const double Rnuc = (ndbd.r0/HBARC)*pow(Anuc,1.0/3.0); // the nuclear radius [MeV^-1]
    const double prefact = 2*(2*Rnuc)/PI; // factor in-front of M0nu TBME, extra global 2 for nutbar (as confirmed by benchmarking with Ca48 NMEs) [MeV^-1]
    M0nuHeader(ndbd,dirname,reduced,Anuc,src,Rnuc,prefact,barcode); // create a header file
    modelspace.PreCalculateMoshinsky(); // pre-calculate the needed Moshinsky brackets, for efficiency
    ndbd.PreCalcT6j(); // pre-calculate the expected 6j-symbols, for efficiency
    ndbd.PreCalcIntegrals(); // pre-calculate the needed integrals over dq and dr, for efficiency
//--T--
    // now we'll precalculate and cache the factors involving the Clebsch-Goradan coefficients, see "factSH" within the TMT
    double SQRT6 = sqrt(6);
    double listSH[PairFN(e2max,e2max)];
    #pragma omp parallel for schedule(dynamic,1)
    for (int lr=0; lr<=e2max; lr++)
    {
      for (int lpr=0; lpr<=e2max; lpr++)
      {
        int tempPF = PairFN(lr,lpr);
        if ((abs(lr-2) <= lpr and lpr <= lr+2) and ((lr+lpr)%2 == 0))
        {
          listSH[tempPF] = SQRT6*CG(lr,0,2,0,lpr,0)*sqrt(2*lr + 1);
        }
        else
        {
          listSH[tempPF] = 0; // just to be safe, amirite!?
        }
      }
    }
    M0nuT_TBME.profiler.timer["M0nuT_1_sur"] += omp_get_wtime() - t_start; // profiling (r)
//--T--
    // create the TBMEs of M0nu
    // auto loops over the TBME channels and such
    cout<<"calculating M0nu TBMEs..."<<endl;
    t_start_tbme = omp_get_wtime(); // profiling (s)
    double Seval = 2.0*sqrt(5.0); // eigenvalue of S for "T", where the only non-zero one is for Sf = Si = 1
    for (auto& itmat : M0nuT_TBME.TwoBody.MatEl)
    {
      int chbra = itmat.first[0]; // grab the channel count from auto
      int chket = itmat.first[1]; // " " " " " "
      TwoBodyChannel& tbc_bra = modelspace.GetTwoBodyChannel(chbra); // grab the two-body channel
      TwoBodyChannel& tbc_ket = modelspace.GetTwoBodyChannel(chket); // " " " " "
      int nbras = tbc_bra.GetNumberKets(); // get the number of bras
      int nkets = tbc_ket.GetNumberKets(); // get the number of kets
      double J = tbc_bra.J; // NOTE: by construction, J := J_ab == J_cd := J'
      double Jhat; // set below based on "reduced" variable
//--T--
      if (reduced == "NR")
      {
        Jhat = 1.0; // for non-reduced elements, to compare with JE
      }
      else //if (reduced == "R")
      {
        Jhat = sqrt(2*J + 1); // the hat factor of J
      }
      t_start_omp = omp_get_wtime(); // profiling (s)
      #pragma omp parallel for schedule(dynamic,100) // need to do: PreCalculateMoshinsky(), PreCalcT6j, and PreCalcIntegrals() [above] and then "#pragma omp critical" [below]
      for (int ibra=0; ibra<nbras; ibra++)
      {
        Ket& bra = tbc_bra.GetKet(ibra); // get the final state = <ab|
//--T--
        int ia = bra.p; // get the integer label a
        int ib = bra.q; // get the integer label b
        Orbit& oa = modelspace.GetOrbit(ia); // get the <a| state orbit
        Orbit& ob = modelspace.GetOrbit(ib); // get the <b| state prbit
        for (int iket=0; iket<nkets; iket++)
        {
          Ket& ket = tbc_ket.GetKet(iket); // get the initial state = |cd>
          int ic = ket.p; // get the integer label c
          int id = ket.q; // get the integer label d
          Orbit& oc = modelspace.GetOrbit(ic); // get the |c> state orbit
          Orbit& od = modelspace.GetOrbit(id); // get the |d> state orbit
          int na = oa.n; // this is just...
          int nb = ob.n;
          int nc = oc.n;
          int nd = od.n;
//--T--
          int la = oa.l;
          int lb = ob.l;
          int lc = oc.l;
          int ld = od.l;
          double ja = oa.j2/2.0;
          double jb = ob.j2/2.0;
          double jc = oc.j2/2.0;
          double jd = od.j2/2.0; // ...for convenience
          double sumLS = 0; // for the Bessel's Matrix Elemets (BMEs)
          double sumLSas = 0; // (anti-symmetric part)
          for (int Lf = abs(la-lb); Lf <= la+lb; Lf++) // sum over angular momentum coupled to l_a and l_b
          {
            for (int Li = abs(lc-ld); Li <= lc+ld; Li++) // sum over angular momentum coupled to l_c and l_d
            {
//--T--
              double bulk = 0; // this is a bulk product which will only be worth calculating if the Moshinsky brackets are non-zero
              double bulkas = 0; // (anti-symmetric part)
              double sumMT = 0; // for the Moshinsky transformation
              double sumMTas = 0; // (anti-symmetric part)
              if (((abs(Lf-1) <= J) and (J <= Lf+1)) and ((abs(1-Li) <= J) and (J <= 1+Li)) and ((abs(Lf-Li) <= 2) and (2 <= Lf+Li))) // so that the GetT6j in factLS/factMT below are non-zero
              {
                int eps_ab = 2*na + la + 2*nb + lb; // for conservation of energy in the Moshinsky brackets
                int tempmaxnr = floor((eps_ab - Lf)/2.0); // just for the limits below
                for (int nr = 0; nr <= tempmaxnr; nr++)
                {
                  int tempmaxNcom = tempmaxnr - nr; // just for the limits below
                  for (int Ncom = 0; Ncom <= tempmaxNcom; Ncom++)
                  {
                    int tempminlr = ceil((eps_ab - Lf)/2.0) - (nr + Ncom); // just for the limits below
                    int tempmaxlr = floor((eps_ab + Lf)/2.0) - (nr + Ncom); // " " " " "
                    for (int lr = tempminlr; lr <= tempmaxlr; lr++)
                    {
//--T--
                      int Lam = eps_ab - 2*(nr + Ncom) - lr; // via Equation (4.74) of my thesis
                      int tempminlpr = abs(Lam - Li); // just for the limits below
                      int tempmaxlpr = Lam + Li; // " " " " "
                      for (int lpr = tempminlpr; lpr <= tempmaxlpr; lpr++)
                      {
                        if ((abs(lr-2) <= lpr and lpr <= lr+2) and ((lr+lpr)%2 == 0)) // so that the CG in listSH and factMT below are non-zero
                        {
                          int eps_cd = 2*nc + lc + 2*nd + ld; // for conservation of energy in the Moshinsky brackets
                          double npr = ((eps_cd - eps_ab)/2.0) + ((lr - lpr)/2.0) + nr; // via Equation (4.74) of my thesis...
                          //double npr = ((eps_cd - lpr - Lam)/2.0) - Ncom; // ...also works
                          double intnpr = npr - floor(npr); // to check if npr is an integer
                          if ((npr >= 0) and (intnpr == 0))
                          {
//--T--
                            double tempLf = (2*Lf + 1); // just for the lines below
                            double tempLi = (2*Li + 1); // " " " " "
                            double factLS = modelspace.phase(1 + J + Li)*ndbd.GetT6j(Lf,1,J,1,Li); // factor from transforming from jj-coupling to ls-coupling
                            double normab = sqrt(tempLf*3*(2*ja + 1)*(2*jb + 1)); // normalization factor for the 9j-symbol out front
                            double nNJab = normab*modelspace.GetNineJ(la,lb,Lf,0.5,0.5,1,ja,jb,J); // the normalized 9j-symbol out front
                            double normcd = sqrt(tempLi*3*(2*jc + 1)*(2*jd + 1)); // normalization factor for the second 9j-symbol
                            double nNJcd = normcd*modelspace.GetNineJ(lc,ld,Li,0.5,0.5,1,jc,jd,J); // the second normalized 9j-symbol
                            double nNJdc = normcd*modelspace.GetNineJ(ld,lc,Li,0.5,0.5,1,jd,jc,J); // (anti-symmetric part)
                            double tempMT = sqrt(tempLf*tempLi); // the hat factor of Lf and Li
                            bulk = Seval*factLS*nNJab*nNJcd*tempMT; // bulk product of the above
                            bulkas = Seval*factLS*nNJab*nNJdc*tempMT; // (anti-symmetric part)
                            double Df,Di,asDi,integral,factMT; // these *will* be set below
                            #pragma omp critical
                            {
//--T--
                              Df = modelspace.GetMoshinsky(Ncom,Lam,nr,lr,na,la,nb,lb,Lf); // Ragnar has -- double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);
                              Di = modelspace.GetMoshinsky(Ncom,Lam,npr,lpr,nc,lc,nd,ld,Li); // " " " "
                              asDi = modelspace.GetMoshinsky(Ncom,Lam,npr,lpr,nd,ld,nc,lc,Li); // (anti-symmetric part)
                              integral = ndbd.GetIntegral(nr,lr,npr,lpr); // grab the pre-calculated integral wrt dq and dr from the IntList of the NDBD class
                              factMT = modelspace.phase(Li + lr + Lam)*ndbd.GetT6j(Lf,lr,Lam,lpr,Li); // factor from TMT
                            }
                            double tempfactMTSH = factMT*listSH[PairFN(lr,lpr)]; // just for the lines below, listSH from the spherical harmonics of order 2
                            sumMT += Df*Di*tempfactMTSH*integral; // perform the Moshinsky transformation
                            sumMTas += Df*asDi*tempfactMTSH*integral; // (anti-symmetric part)
                          } // end of if: npr \in \Nat_0
                        } // end of if: CG != 0, factMT != 0
                      } // end of for-loop over: lpr
                    } // end of for-loop over: lr
                  } // end of for-loop over: Ncom
                } // end of for-loop over: nr
//--T--
              } // end of if: factLS != 0, factMT != 0
              sumLS += bulk*sumMT; // perform the LS-coupling sum
              sumLSas += bulkas*sumMTas; // (anti-symmetric part)
            } // end of for-loop over: Li
          } // end of for-loop over: Lf
          double Mtbme = cpNorm(ia,ib)*cpNorm(ic,id)*prefact*Jhat*(sumLS - modelspace.phase(jc + jd - J)*sumLSas); // compute the final matrix element, anti-symmetrize
          M0nuT_TBME.TwoBody.SetTBME(chbra,chket,ibra,iket,Mtbme); // set the two-body matrix elements (TBME) to Mtbme
        } // end of for-loop over: iket
      } // end of for-loop over: ibra
      M0nuT_TBME.profiler.timer["M0nuT_3_omp"] += omp_get_wtime() - t_start_omp; // profiling (r)
    } // end of for-loop over: auto
    cout<<"...done calculating M0nu TBMEs"<<endl;
    M0nuT_TBME.profiler.timer["M0nuT_2_tbme"] += omp_get_wtime() - t_start_tbme; // profiling (r)
    M0nuT_TBME.profiler.timer["M0nuT_adpt_Op"] += omp_get_wtime() - t_start; // profiling (r)
    return M0nuT_TBME;
  }
//^^T^^


/// This prints the integrand(nr,lr,npr) to file as a function of q [MeV] over the interval [qmin,qmax]
/// It's been primarily used for testing/debugging/etc
  void M0nu_PrintIntegrand(ModelSpace& modelspace, string dirdat, int mode, double Ediff, string src, int nr, int lr, int npr, int lpr, double qmin, double qmax)
  {
    // create the output files which hold the integrand
    string base = dirdat;
    string Istr = "M0nu_integrand_";
    string typeGT = "GT";
    string typeF = "F";
    string typeT = "T";
    string end = ".txt";
    ofstream IGTout(base+Istr+typeGT+end);
    if (!IGTout)
    {
      cout<<"ERROR: with IGTout"<<endl;
      exit(1);
    }
    ofstream IFout(base+Istr+typeF+end);
    if (!IFout)
    {
      cout<<"ERROR: with IFout"<<endl;
      exit(1);
    }
    ofstream ITout(base+Istr+typeT+end);
    if (!ITout)
    {
      cout<<"ERROR: with ITout"<<endl;
      exit(1);
    }
    // set stuff up
    string meth = "PSH";
    double hw = modelspace.GetHbarOmega();
    int e2max = modelspace.GetE2max();
    NDBD ndbdGT(hw,e2max,typeGT,Ediff,src,meth);
    NDBD ndbdF(hw,e2max,typeF,Ediff,src,meth);
    NDBD ndbdT(hw,e2max,typeT,Ediff,src,meth);
    gsl_function FqGT,FqF,FqT;
    double fqparams[18];
    int dqswitch = 0;
    int drswitch = 0;
    ndbdGT.SetPSH(nr,lr,npr,lpr,dqswitch,drswitch); // could use either ndbdGT or ndbdF or ndbdT here
    ndbdGT.SetIparams(nr,lr,npr,lpr,drswitch,fqparams); // " " " " " " " " "
    FqGT.params = &(fqparams[0]);
    FqF.params = &(fqparams[0]);
    FqT.params = &(fqparams[0]);
    FqGT.function = &fqGT; // point to the Gamow-Teller type
    FqF.function = &fqF; // point to the Fermi type
    FqT.function = &fqT; // point to the Tensor type
    cout<<endl<<"nr = "<<nr<<", lr = "<<lr<<", npr = "<<npr<<", lpr = "<<lpr<<endl;
    cout<<"dqswitch = "<<dqswitch<<", drsqitch = "<<drswitch<<endl;
    cout<<"where: 1 = analytic"<<endl;
    cout<<"       2 = qagiu"<<endl;
    cout<<"       3 = qag"<<endl<<endl;
    if (mode == 0)
    {
      // FIRST we'll print the integrands
      cout<<"printing the integrands..."<<endl;
      int res = 40001; // plotting resolution
      double dr = (qmax - qmin)/(res - 1);
      for (int i=0; i<res; i++)
      {
        double q = i*dr + qmin;
        IGTout<<setprecision(12)<<q<<"\t"<<GSL_FN_EVAL(&FqGT,q)<<endl; // print q and integrand to the output file, GSL_FN_EVAL is found in <gsl_math.h>
        IFout<<setprecision(12)<<q<<"\t"<<GSL_FN_EVAL(&FqF,q)<<endl; // " " " " " " " " " " " " "
        ITout<<setprecision(12)<<q<<"\t"<<GSL_FN_EVAL(&FqT,q)<<endl; // " " " " " " " " " " " " "
      }
      IGTout.close();
      IFout.close();
      ITout.close();
      cout<<"...done printing the integrands"<<endl<<endl;
      // LAST we'll run the integrations and output them to screen
      cout<<"# from M0nu_PrintIntegrand"<<endl;
      cout<<"the GT integration = "<<setprecision(12)<<ndbdGT.dqdrPSH(nr,lr,npr,lpr)<<endl; // don't use GetIntgral...
      cout<<"the F integration  = "<<setprecision(12)<<ndbdF.dqdrPSH(nr,lr,npr,lpr)<<endl; // ...because ndbdGT and ndbdF share the IntList...
      cout<<"the T integration  = "<<setprecision(12)<<ndbdT.dqdrPSH(nr,lr,npr,lpr)<<endl<<endl; // ...and so does ndbdT :S
    }
    else if (mode == 1)
    {
      // FIRST we'll run the integrations and output them to screen
      cout<<"# from M0nu_PrintIntegrand"<<endl;
      cout<<"the GT integration = "<<setprecision(12)<<ndbdGT.dqdrPSH(nr,lr,npr,lpr)<<endl; // don't use GetIntgral...
      cout<<"the F integration  = "<<setprecision(12)<<ndbdF.dqdrPSH(nr,lr,npr,lpr)<<endl; // ...because ndbdGT and ndbdF share the IntList...
      cout<<"the T integration  = "<<setprecision(12)<<ndbdT.dqdrPSH(nr,lr,npr,lpr)<<endl<<endl; // ...and so does ndbdT :S
      // LAST we'll print the integrands
      cout<<"printing the integrands..."<<endl;
      int res = 40001; // plotting resolution
      double dr = (qmax - qmin)/(res - 1);
      for (int i=0; i<res; i++)
      {
        double q = i*dr + qmin;
        IGTout<<setprecision(12)<<q<<"\t"<<GSL_FN_EVAL(&FqGT,q)<<endl; // print q and integrand to the output file, GSL_FN_EVAL is found in <gsl_math.h>
        IFout<<setprecision(12)<<q<<"\t"<<GSL_FN_EVAL(&FqF,q)<<endl; // " " " " " " " " " " " " "
        ITout<<setprecision(12)<<q<<"\t"<<GSL_FN_EVAL(&FqT,q)<<endl; // " " " " " " " " " " " " "
      }
      IGTout.close();
      IFout.close();
      ITout.close();
      cout<<"...done printing the integrands"<<endl<<endl;
    }
    else
    {
      cout<<"ERROR 6633: invalid option for mode = "<<mode<<endl;
      exit(1);
    }
  }


////////////////////////// END OF: via NDBD.hh/cc (CP) //////////////////////////




 // Evaluate <bra | r1*r2 | ket>, omitting the factor (hbar * omega) /(m * omega^2)
/// Returns the normalized, anti-symmetrized, J-coupled, two-body matrix element of \f$ \frac{m\omega^2}{\hbar \omega} \vec{r}_1\cdot\vec{r}_2 \f$.
/// Calculational details are similar to Calculate_p1p2().
 double Calculate_r1xp2(ModelSpace& modelspace, Ket & bra, Ket & ket, int Jab, int Jcd)
 {
   Orbit & oa = modelspace.GetOrbit(bra.p);
   Orbit & ob = modelspace.GetOrbit(bra.q);
   Orbit & oc = modelspace.GetOrbit(ket.p);
   Orbit & od = modelspace.GetOrbit(ket.q);

   int na = oa.n;
   int nb = ob.n;
   int nc = oc.n;
   int nd = od.n;

   int la = oa.l;
   int lb = ob.l;
   int lc = oc.l;
   int ld = od.l;

   double ja = oa.j2/2.0;
   double jb = ob.j2/2.0;
   double jc = oc.j2/2.0;
   double jd = od.j2/2.0;

   int fab = 2*na + 2*nb + la + lb;
   int fcd = 2*nc + 2*nd + lc + ld;
   if (abs(fab-fcd)%2 >0) return 0; // p1*p2 only connects kets with delta N = 0,1
   if (abs(fab-fcd)>2) return 0; // p1*p2 only connects kets with delta N = 0,1


   double sa,sb,sc,sd;
   sa=sb=sc=sd=0.5;

   double r1xp2=0;

   // First, transform to LS coupling using 9j coefficients
   for (int Lab=abs(la-lb); Lab<= la+lb; ++Lab)
   {
     for (int Sab=0; Sab<=1; ++Sab)
     {
       if ( abs(Lab-Sab)>Jab or Lab+Sab<Jab) continue;

       double njab = NormNineJ(la,sa,ja, lb,sb,jb, Lab,Sab,Jab);
       if (njab == 0) continue;
       int Scd = Sab;
       int Lcd = Lab;
       double njcd = NormNineJ(lc,sc,jc, ld,sd,jd, Lcd,Scd,Jcd);
       if (njcd == 0) continue;

       // Next, transform to rel / com coordinates with Moshinsky tranformation
       for (int N_ab=0; N_ab<=fab/2; ++N_ab)  // N_ab = CoM n for a,b
       {
         for (int Lam_ab=0; Lam_ab<= fab-2*N_ab; ++Lam_ab) // Lam_ab = CoM l for a,b
         {
           int Lam_cd = Lam_ab; // Lcm and Lrel conserve lam and Lam, 
           for (int lam_ab=(fab-2*N_ab-Lam_ab)%2; lam_ab<= (fab-2*N_ab-Lam_ab); lam_ab+=2) // lam_ab = relative l for a,b
           {
              if (Lab<abs(Lam_ab-lam_ab) or Lab>(Lam_ab+lam_ab) ) continue;
              // factor to account for antisymmetrization

              int asymm_factor = (abs(bra.op->tz2+ket.op->tz2) + abs(bra.op->tz2+ket.oq->tz2)*modelspace.phase( lam_ab + Sab ))/ 2;
              if ( asymm_factor ==0 ) continue;

              int lam_cd = lam_ab; // Lcm and Lrel operators conserve lam and Lam
              int n_ab = (fab - 2*N_ab-Lam_ab-lam_ab)/2; // n_ab is determined by energy conservation

              double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);

              if (abs(mosh_ab)<1e-8) continue;

//              for (int N_cd=max(0,N_ab-1); N_cd<=N_ab+1; ++N_cd) // N_cd = CoM n for c,d
               int N_cd = N_ab;
              {
                int n_cd = (fcd - 2*N_cd-Lam_cd-lam_cd)/2; // n_cd is determined by energy conservation
                if (n_cd < 0) continue;
                if  (n_ab != n_cd ) continue;

                double mosh_cd = modelspace.GetMoshinsky(N_cd,Lam_cd,n_cd,lam_cd,nc,lc,nd,ld,Lcd);
                if (abs(mosh_cd)<1e-8) continue;

//                double lcm = 0;
//                double lrel = 0;
                double sj_rel = modelspace.GetSixJ( lam_ab, lam_cd, 1, Lcd, Lab, Lam_ab );
                double sj_cm = modelspace.GetSixJ( Lam_ab, Lam_cd, 1, Lcd, Lab, lam_ab );
                double lrel = modelspace.phase( Lab + Lcd ) * ( (2*Lab+1)*(2*Lcd+1) ) * sj_rel * sj_rel  * sqrt(lam_ab*(lam_ab+1)*(2*lam_ab+1));
                double lcm = ( (2*Lab+1)*(2*Lcd+1) ) * sj_cm * sj_cm  * sqrt(Lam_ab*(Lam_ab+1)*(2*Lam_ab+1));
                r1xp2 += (lcm-lrel) * 0.5;

              } // N_cd
           } // lam_ab
         } // Lam_ab
       } // N_ab

     } // Sab
   } // Lab

   // normalize.
   r1xp2 *= 1.0 / sqrt((1.0+bra.delta_pq())*(1.0+ket.delta_pq()));
   return r1xp2 ;

 }


 // Orbital angular momentum squared L^2 in the relative coordinate.
 // This was written with the deuteron in mind. Not sure if it will be useful for other systems.
 Operator L2rel_Op(ModelSpace& modelspace)
 {
   Operator OpOut(modelspace, 0,0,0,2);

   int nchan = modelspace.GetNumberTwoBodyChannels();
   modelspace.PreCalculateMoshinsky();
   #pragma omp parallel for schedule(dynamic,1) 
   for (int ch=0; ch<nchan; ++ch)
   {
    TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
    int J = tbc.J;
    int nkets = tbc.GetNumberKets();
    for (int ibra=0;ibra<nkets;++ibra)
    {
     Ket & bra = tbc.GetKet(ibra);
     for (int iket=ibra;iket<nkets;++iket)
     {
       Ket & ket = tbc.GetKet(iket);

       Orbit & oa = modelspace.GetOrbit(bra.p);
       Orbit & ob = modelspace.GetOrbit(bra.q);
       Orbit & oc = modelspace.GetOrbit(ket.p);
       Orbit & od = modelspace.GetOrbit(ket.q);
    
       int na = oa.n;
       int nb = ob.n;
       int nc = oc.n;
       int nd = od.n;
    
       int la = oa.l;
       int lb = ob.l;
       int lc = oc.l;
       int ld = od.l;
    
       double ja = oa.j2/2.0;
       double jb = ob.j2/2.0;
       double jc = oc.j2/2.0;
       double jd = od.j2/2.0;
    
       int fab = 2*na + 2*nb + la + lb;
       int fcd = 2*nc + 2*nd + lc + ld;
       if (fab != fcd) continue; // L2 conserves all the quantum numbers
    
       double sa,sb,sc,sd;
       sa=sb=sc=sd=0.5;
    
       double L2rel=0;
    
       for (int Lab=abs(la-lb); Lab<= la+lb; ++Lab)
       {
         for (int Sab=0; Sab<=1; ++Sab)
         {
           if ( abs(Lab-Sab)>J or Lab+Sab<J) continue;
    
           double njab = NormNineJ(la,sa,ja, lb,sb,jb, Lab,Sab,J);
           if (njab == 0) continue;
           int Scd = Sab;
           int Lcd = Lab;
           double njcd = NormNineJ(lc,sc,jc, ld,sd,jd, Lcd,Scd,J);
           if (njcd == 0) continue;
    
           // Next, transform to rel / com coordinates with Moshinsky tranformation
           for (int N_ab=0; N_ab<=fab/2; ++N_ab)  // N_ab = CoM n for a,b
           {
             for (int Lam_ab=0; Lam_ab<= fab-2*N_ab; ++Lam_ab) // Lam_ab = CoM l for a,b
             {
               int Lam_cd = Lam_ab; // tcm and trel conserve lam and Lam, ie relative and com orbital angular momentum
               for (int lam_ab=(fab-2*N_ab-Lam_ab)%2; lam_ab<= (fab-2*N_ab-Lam_ab); lam_ab+=2) // lam_ab = relative l for a,b
               {
                  if (Lab<abs(Lam_ab-lam_ab) or Lab>(Lam_ab+lam_ab) ) continue;
                  // factor to account for antisymmetrization
    
                  int asymm_factor = (abs(bra.op->tz2+ket.op->tz2) + abs(bra.op->tz2+ket.oq->tz2)*modelspace.phase( lam_ab + Sab ))/ 2;
                  if ( asymm_factor ==0 ) continue;
    
                  int lam_cd = lam_ab; // tcm and trel conserve lam and Lam
                  int n_ab = (fab - 2*N_ab-Lam_ab-lam_ab)/2; // n_ab is determined by energy conservation
    
                  double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);
    
                  if (abs(mosh_ab)<1e-8) continue;
    
                  for (int N_cd=max(0,N_ab-1); N_cd<=N_ab+1; ++N_cd) // N_cd = CoM n for c,d
                  {
                    int n_cd = (fcd - 2*N_cd-Lam_cd-lam_cd)/2; // n_cd is determined by energy conservation
                    if (n_cd < 0) continue;
                    if  (n_ab != n_cd or N_ab != N_cd) continue;
    
                    double mosh_cd = modelspace.GetMoshinsky(N_cd,Lam_cd,n_cd,lam_cd,nc,lc,nd,ld,Lcd);
                    if (abs(mosh_cd)<1e-8) continue;
    
                    double prefactor = njab * njcd * mosh_ab * mosh_cd * asymm_factor;
                    L2rel += lam_ab*(lam_ab+1) * prefactor;
    
                  } // N_cd
               } // lam_ab
             } // Lam_ab
           } // N_ab
    
         } // Sab
       } // Lab


       OpOut.TwoBody.SetTBME(ch,ibra,iket,L2rel);
       OpOut.TwoBody.SetTBME(ch,iket,ibra,L2rel);

     }
   }
  }

  return OpOut;

 }



 Operator EKKShift( Operator& Hin, int Nlower, int Nupper)
 {
   ModelSpace* modelspace = Hin.GetModelSpace();
   Operator EKKShift(*modelspace,0,0,0,2);
   double elower = 0;
   double eupper = 0;
   vector<index_t> index_lower;
   vector<index_t> index_upper;
   for (int i=0; i<modelspace->GetNumberOrbits(); ++i)
   {
     Orbit& oi = modelspace->GetOrbit(i);
     int N = 2*oi.n + oi.l;
     if (N==Nlower)
     {
       elower += Hin.OneBody(i,i);
       index_lower.push_back(i);
     }
     if (N==Nupper)
     {
       eupper += Hin.OneBody(i,i);
       index_upper.push_back(i);
     }
   }
   elower /= index_lower.size();
   eupper /= index_upper.size();
   double emid = 0.5 * (elower + eupper);
   for ( auto i : index_lower ) EKKShift.OneBody(i,i) = emid - elower;
   for ( auto i : index_upper ) EKKShift.OneBody(i,i) = emid - eupper;
   cout << "In EKKShift, elower, eupper, emid = " << elower << " , " << eupper << " , " << emid << endl;
   return EKKShift;
 }


  map<index_t,double> GetSecondOrderOccupations(Operator& H, int emax)
  {
//    ModelSpace* modelspace = H.GetModelSpace();
    map<index_t,double> hole_list;
    cout << "GetSecondOrderOccupations : Not yet implemented" << endl;
    return hole_list;
  }


  /// Embeds the one-body operator of op1 in the two-body part, using mass number A in the embedding.
  /// Note that the embedded operator is added to the two-body part, rather than overwriting.
  /// The one-body part is left as-is.
  void Embed1BodyIn2Body(Operator& op1, int A)
  {
    if (A<2)
    {
      cout << "Embed1BodyIn2Body: A = " << A << ". You clearly didn't mean to do that..." << endl;
      return;
    }
    ModelSpace* modelspace = op1.GetModelSpace();
    int Lambda = op1.GetJRank();
    for (auto& itmat : op1.TwoBody.MatEl)
    {
      index_t ch_bra = itmat.first[0];
      index_t ch_ket = itmat.first[1];
      arma::mat& TBME = itmat.second;

      TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
      index_t nbras = tbc_bra.GetNumberKets();
      index_t nkets = tbc_ket.GetNumberKets();
      int Jbra = tbc_bra.J;
      int Jket = tbc_ket.J;
      for (index_t ibra=0;ibra<nbras;++ibra)
      {
        Ket& bra = tbc_bra.GetKet(ibra);
        index_t i = bra.p;
        index_t j = bra.q;
        for (index_t iket=0;iket<nkets;++iket)
        {
          Ket& ket = tbc_ket.GetKet(iket);
          index_t k = ket.p;
          index_t l = ket.q;
          double embedded_tbme = GetEmbeddedTBME(op1,i,j,k,l,Jbra,Jket,Lambda) / (A-1.0);
          TBME(ibra,iket) += embedded_tbme;
        }
      }
    }
  }

  /// Returns a normalized TBME formed by embedding the one body part of op1 in a two body operator. Assumes A=2 in the formula.
  /// For other values of A, divide by (A-1).
  double GetEmbeddedTBME(Operator& op1, index_t i, index_t j, index_t k, index_t l, int Jbra,int Jket, int Lambda)
  {
    double embedded_tbme = 0;
    ModelSpace* modelspace = op1.GetModelSpace();
    double ji = modelspace->GetOrbit(i).j2 * 0.5;
    double jj = modelspace->GetOrbit(j).j2 * 0.5;
    double jk = modelspace->GetOrbit(k).j2 * 0.5;
    double jl = modelspace->GetOrbit(l).j2 * 0.5;
    arma::mat& OB = op1.OneBody;
    if (Lambda==0) // scalar => no sixJs, tbmes are not reduced.
    {
       if (j==l)  embedded_tbme += OB(i,k);
       if (i==k)  embedded_tbme += OB(j,l);
       if (i==l)  embedded_tbme -= OB(j,k) * modelspace->phase(ji+jj-Jbra);
       if (j==k)  embedded_tbme -= OB(i,l) * modelspace->phase(ji+jj-Jbra);
    }
    else // Tensor => slightly more complicated, tbmes are reduced.
    {
       if (j==l)  embedded_tbme += OB(i,k) * modelspace->phase(ji+jj+Jket) * SixJ(Jbra,Jket,Lambda,jk,ji,jj);
       if (i==k)  embedded_tbme += OB(j,l) * modelspace->phase(jk+jl-Jbra) * SixJ(Jbra,Jket,Lambda,jl,jj,ji);
       if (j==k)  embedded_tbme -= OB(i,l) * modelspace->phase(ji+jj+jk+jl)* SixJ(Jbra,Jket,Lambda,jl,ji,jj);
       if (i==l)  embedded_tbme -= OB(j,k) * modelspace->phase(Jket-Jbra)  * SixJ(Jbra,Jket,Lambda,jk,jj,ji);
       embedded_tbme *= sqrt((2*Jbra+1)*(2*Jket+1)) * modelspace->phase(Lambda);
    }
    if (i==j) embedded_tbme /=SQRT2;
    if (k==l) embedded_tbme /=SQRT2;
    return embedded_tbme;
  }



//  double FCcoefFunction(int n1, int l1, double hw1, int n2, int l2, double hw2, double r)
  double FCcoefIntegrand(double r, void* params)
  {
    int n1     = ((double*) params)[0];
    int l1     = ((double*) params)[1];
    double hw1 = ((double*) params)[2];
    int n2     = ((double*) params)[3];
    int l2     = ((double*) params)[4];
    double hw2 = ((double*) params)[5];
    return r*r*HO_Radial_psi( n1,l1,hw1,r)*HO_Radial_psi( n2,l2,hw2,r);
  }

  double FrequencyConversionCoeff(int n1, int l1, double hw1, int n2, int l2, double hw2)
  {
    double parameters[6] = {(double)n1, (double)l1, hw1, (double)n2, (double)l2, hw2};
    gsl_function F;
    F.function = &FCcoefIntegrand;
    F.params = &(parameters[0]);
    double r0 = 0;
    double epsabs = 1e-7;
    double epsrel = 1e-7;
    size_t limit = 1000; // is this a reasonable number?
    gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(limit);
    double result = 0;
    double abserr;

    gsl_integration_qagiu(&F, r0, epsabs, epsrel, limit, workspace, &result, &abserr );

    gsl_integration_workspace_free(workspace);
    return result;

  }


  void CommutatorTest(Operator& X, Operator& Y)
  {
    Operator Zscalar(X);
    if ( (X.IsHermitian() and Y.IsHermitian()) or (X.IsAntiHermitian() and Y.IsAntiHermitian()) ) Zscalar.SetAntiHermitian();
    if ( (X.IsHermitian() and Y.IsAntiHermitian()) or (X.IsAntiHermitian() and Y.IsHermitian()) ) Zscalar.SetHermitian();
    Zscalar.Erase();
    Operator Ztensor(Zscalar);
    Operator Yred = Y;
    Reduce(Yred);

    cout << "operator norms: " << X.Norm() << "  " << Y.Norm() << endl;
//    X.comm111ss(Y,Zscalar);
//    X.comm111st(Yred,Ztensor);
    Zscalar.comm111ss(X,Y);
    Ztensor.comm111st(X,Yred);
    Zscalar.Symmetrize();
    Ztensor.Symmetrize();
    UnReduce(Ztensor);
    cout << "comm111 norm = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << ",   " << Ztensor.OneBodyNorm() << " " << Ztensor.TwoBodyNorm() << endl;
    Zscalar -= Ztensor;
    cout << "comm111 diff = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << endl;

    Zscalar.Erase();
    Ztensor.Erase();
    cout << "121ss" << endl;
    Zscalar.comm121ss(X,Y);
    cout << "121st" << endl;
    Ztensor.comm121st(X,Yred);
    Zscalar.Symmetrize();
    Ztensor.Symmetrize();
    UnReduce(Ztensor);
    cout << "comm121 norm = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << ",   " << Ztensor.OneBodyNorm() << " " << Ztensor.TwoBodyNorm() << endl;
    Zscalar -= Ztensor;
    cout << "comm121 diff = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << endl;

    Zscalar.Erase();
    Ztensor.Erase();
    Zscalar.comm122ss(X,Y);
    Ztensor.comm122st(X,Yred);
    Zscalar.Symmetrize();
    Ztensor.Symmetrize();
    UnReduce(Ztensor);
    cout << "comm122 norm = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << ",   " << Ztensor.OneBodyNorm() << " " << Ztensor.TwoBodyNorm() << endl;
    Zscalar -= Ztensor;
    cout << "comm122 diff = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << endl;

    Zscalar.Erase();
    Ztensor.Erase();
    Zscalar.comm222_pp_hh_221ss(X,Y);
    Ztensor.comm222_pp_hh_221st(X,Yred);
    Zscalar.Symmetrize();
    Ztensor.Symmetrize();
    UnReduce(Ztensor);
    cout << "comm222_pp_hh_221 norm = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << ",   " << Ztensor.OneBodyNorm() << " " << Ztensor.TwoBodyNorm() << endl;
    Zscalar -= Ztensor;
    cout << "comm222_pp_hh_221 diff = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << endl;

    Zscalar.Erase();
    Ztensor.Erase();
    X.comm222_phss(Y,Zscalar);
//    Reduce(Y); // Not sure why I can't use Yred...
    Zscalar.comm222_phss(X,Y);
    Ztensor.comm222_phst(X,Yred);
    Zscalar.Symmetrize();
    Ztensor.Symmetrize();
    UnReduce(Ztensor);
    cout << "comm222_ph norm = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << ",   " << Ztensor.OneBodyNorm() << " " << Ztensor.TwoBodyNorm() << endl;
    Zscalar -= Ztensor;
    cout << "comm222_ph diff = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << endl;


  }



/*
  void CommutatorTest(Operator& X, Operator& Y)
  {
    Operator Zscalar(X);
    if ( (X.IsHermitian() and Y.IsHermitian()) or (X.IsAntiHermitian() and Y.IsAntiHermitian()) ) Zscalar.SetAntiHermitian();
    if ( (X.IsHermitian() and Y.IsAntiHermitian()) or (X.IsAntiHermitian() and Y.IsHermitian()) ) Zscalar.SetHermitian();
    Zscalar.Erase();
    Operator Ztensor(Zscalar);
    Operator Yred = Y;
    Reduce(Yred);

    cout << "operator norms: " << X.Norm() << "  " << Y.Norm() << endl;
//    X.comm111ss(Y,Zscalar);
//    X.comm111st(Yred,Ztensor);
    Zscalar.comm111ss(X,Y);
    Ztensor.comm111st(X,Yred);
    Zscalar.Symmetrize();
    Ztensor.Symmetrize();
    UnReduce(Ztensor);
    cout << "comm111 norm = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << ",   " << Ztensor.OneBodyNorm() << " " << Ztensor.TwoBodyNorm() << endl;
    Zscalar -= Ztensor;
    cout << "comm111 diff = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << endl;

    Zscalar.Erase();
    Ztensor.Erase();
    X.comm121ss(Y,Zscalar);
    X.comm121st(Yred,Ztensor);
    Zscalar.Symmetrize();
    Ztensor.Symmetrize();
    UnReduce(Ztensor);
    cout << "comm121 norm = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << ",   " << Ztensor.OneBodyNorm() << " " << Ztensor.TwoBodyNorm() << endl;
    Zscalar -= Ztensor;
    cout << "comm121 diff = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << endl;

    Zscalar.Erase();
    Ztensor.Erase();
    X.comm122ss(Y,Zscalar);
    X.comm122st(Yred,Ztensor);
    Zscalar.Symmetrize();
    Ztensor.Symmetrize();
    UnReduce(Ztensor);
    cout << "comm122 norm = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << ",   " << Ztensor.OneBodyNorm() << " " << Ztensor.TwoBodyNorm() << endl;
    Zscalar -= Ztensor;
    cout << "comm122 diff = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << endl;

    Zscalar.Erase();
    Ztensor.Erase();
    X.comm222_pp_hh_221ss(Y,Zscalar);
    X.comm222_pp_hh_221st(Yred,Ztensor);
    Zscalar.Symmetrize();
    Ztensor.Symmetrize();
    UnReduce(Ztensor);
    cout << "comm222_pp_hh_221 norm = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << ",   " << Ztensor.OneBodyNorm() << " " << Ztensor.TwoBodyNorm() << endl;
    Zscalar -= Ztensor;
    cout << "comm222_pp_hh_221 diff = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << endl;

    Zscalar.Erase();
    Ztensor.Erase();
    X.comm222_phss(Y,Zscalar);
    Reduce(Y); // Not sure why I can't use Yred...
    X.comm222_phst(Y,Ztensor);
    Zscalar.Symmetrize();
    Ztensor.Symmetrize();
    UnReduce(Ztensor);
    cout << "comm222_ph norm = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << ",   " << Ztensor.OneBodyNorm() << " " << Ztensor.TwoBodyNorm() << endl;
    Zscalar -= Ztensor;
    cout << "comm222_ph diff = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << endl;


  }
*/




// template <typename T>
// T VectorUnion(T& v1)
// {
//   return v1;
// }
// 
// template <typename T, typename... Args>
// T VectorUnion(T& v1, T& v2, Args... args)
// {
//   T vec(v1.size()+v2.size());
//   copy(v1.begin(),v1.end(),vec.begin());
//   copy(v2.begin(),v2.end(),vec.begin()+v1.size());
//   return VectorUnion(vec, args...);
// }
 

}// namespace imsrg_util

