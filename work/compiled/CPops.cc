/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
///                                                  ____                                         ///
///        _________________           _____________/   /\               _________________        ///
///       /____/_____/_____/|         /____/_____/ /___/  \             /____/_____/_____/|       ///
///      /____/_____/__G_ /||        /____/_____/|/   /\  /\           /____/_____/____ /||       ///
///     /____/_____/__+__/|||       /____/_____/|/ G /  \/  \         /____/_____/_____/|||       ///
///    |     |     |     ||||      |     |     |/___/   /\  /\       |     |     |     ||||       ///
///    |  I  |  M  |     ||/|      |  I  |  M  /   /\  /  \/  \      |  I  |  M  |     ||/|       ///
///    |_____|_____|_____|/||      |_____|____/ + /  \/   /\  /      |_____|_____|_____|/||       ///
///    |     |     |     ||||      |     |   /___/   /\  /  \/       |     |     |     ||||       ///
///    |  S  |  R  |     ||/|      |  S  |   \   \  /  \/   /        |  S  |  R  |  G  ||/|       ///
///    |_____|_____|_____|/||      |_____|____\ __\/   /\  /         |_____|_____|_____|/||       ///
///    |     |     |     ||||      |     |     \   \  /  \/          |     |     |     ||||       ///
///    |     |  +  |     ||/       |     |  +  |\ __\/   /           |     |  +  |  +  ||/        ///
///    |_____|_____|_____|/        |_____|_____|/\   \  /            |_____|_____|_____|/         ///
///                                               \___\/                                          ///
///                                                                                               ///
///           imsrg++ : Interface for performing standard IMSRG calculations.                     ///
///                     Usage is imsrg++  option1=value1 option2=value2 ...                       ///
///                     To get a list of options, type imsrg++ help                               ///
///                                                                                               ///
///                                                      - Ragnar Stroberg 2016                   ///
///                                                                                               ///
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <omp.h>
#include "IMSRG.hh"
#include "Parameters.hh"

using namespace imsrg_util;

int main(int argc, char** argv)
{
  // Default parameters, and everything passed by command line args.
#ifdef BUILDVERSION
  cout << "######  imsrg++ build version: " << BUILDVERSION << endl;
#endif
  Parameters parameters(argc,argv);
  if (parameters.help_mode) return 0;

  string inputtbme = parameters.s("2bme");
  string input3bme = parameters.s("3bme");
  string reference = parameters.s("reference");
  string valence_space = parameters.s("valence_space");
  string custom_valence_space = parameters.s("custom_valence_space");
  string basis = parameters.s("basis");
  string method = parameters.s("method");
  string flowfile = parameters.s("flowfile");
  string intfile = parameters.s("intfile");
  string core_generator = parameters.s("core_generator");
  string valence_generator = parameters.s("valence_generator");
  string fmt2 = parameters.s("fmt2");
  string fmt3 = parameters.s("fmt3");
  string denominator_delta_orbit = parameters.s("denominator_delta_orbit");
  string LECs = parameters.s("LECs");
  string scratch = parameters.s("scratch");
  string use_brueckner_bch = parameters.s("use_brueckner_bch");
  string valence_file_format = parameters.s("valence_file_format");
  string occ_file = parameters.s("occ_file");
  string goose_tank = parameters.s("goose_tank");

  int eMax = parameters.i("emax");
  int E3max = parameters.i("e3max");
  int lmax3 = parameters.i("lmax3");
  int targetMass = parameters.i("A");
  int nsteps = parameters.i("nsteps");
  int file2e1max = parameters.i("file2e1max");
  int file2e2max = parameters.i("file2e2max");
  int file2lmax = parameters.i("file2lmax");
  int file3e1max = parameters.i("file3e1max");
  int file3e2max = parameters.i("file3e2max");
  int file3e3max = parameters.i("file3e3max");
  int Ntb = parameters.i("NumLine"); // (CP)

  double hw = parameters.d("hw");
  double smax = parameters.d("smax");
  double ode_tolerance = parameters.d("ode_tolerance");
  double dsmax = parameters.d("dsmax");
  double ds_0 = parameters.d("ds_0");
  double domega = parameters.d("domega");
  double omega_norm_max = parameters.d("omega_norm_max");
  double denominator_delta = parameters.d("denominator_delta");
  double BetaCM = parameters.d("BetaCM");
  double hwBetaCM = parameters.d("hwBetaCM");
  double eta_criterion = parameters.d("eta_criterion");

  vector<string> opnames = parameters.v("Operators");

  vector<Operator> ops;
  vector<string> spwf = parameters.v("SPWF");

/*
  ifstream test;
  // test 2bme file
  if (inputtbme != "none")
  {
    test.open(inputtbme);
    if( not test.good() and fmt2!="oakridge")
    {
      cout << "trouble reading " << inputtbme << " exiting. " << endl;
      return 1;
    }
    test.close();
  }
  // test 3bme file
  if (input3bme != "none")
  {
    test.open(input3bme);
    if( not test.good() )
    {
      cout << "trouble reading " << input3bme << " exiting. " << endl;
      return 1;
    }
    test.close();
  }
*/


  ReadWrite rw;
  rw.SetLECs_preset(LECs);
  rw.SetScratchDir(scratch);
  rw.Set3NFormat( fmt3 );

//  ModelSpace modelspace;

  if (custom_valence_space!="") // if a custom space is defined, the input valence_space is just used as a name
  {
    if (valence_space=="") // if no name is given, then just name it "custom"
    {
      parameters.string_par["valence_space"] = "custom";
      flowfile = parameters.DefaultFlowFile();
      intfile = parameters.DefaultIntFile();
    }
    valence_space = custom_valence_space;
  }


  ModelSpace modelspace = ( reference=="default" ? ModelSpace(eMax,valence_space) : ModelSpace(eMax,reference,valence_space) );

  if (occ_file != "none" and occ_file != "" )
  {
    modelspace.Init_occ_from_file(eMax,valence_space,occ_file);
  }

  if (nsteps < 0)
    nsteps = modelspace.valence.size()>0 ? 2 : 1;


  modelspace.SetHbarOmega(hw);
  if (targetMass>0)
     modelspace.SetTargetMass(targetMass);
  modelspace.SetE3max(E3max);
  if (lmax3>0)
     modelspace.SetLmax3(lmax3);
  
/*
  cout << "Making the operator..." << endl;
  int particle_rank = input3bme=="none" ? 2 : 3;
  Operator Hbare = Operator(modelspace,0,0,0,particle_rank);
  Hbare.SetHermitian();


  if ( goose_tank == "true" or goose_tank == "True")
  {
    Hbare.SetUseGooseTank(true);
  }

  cout << "Reading interactions..." << endl;

//  int nthreads = omp_get_max_threads();
//  {
//   omp_set_num_threads(2);
//   omp_set_nested(1);
//  #pragma omp parallel sections
//  {
//    omp_set_num_threads(max(1,nthreads-3));
//    #pragma omp section
//    {
      if (inputtbme != "none")
      {
        if (fmt2 == "me2j")
          rw.ReadBareTBME_Darmstadt(inputtbme, Hbare,file2e1max,file2e2max,file2lmax);
        else if (fmt2 == "navratil" or fmt2 == "Navratil")
          rw.ReadBareTBME_Navratil(inputtbme, Hbare);
        else if (fmt2 == "oslo" )
          rw.ReadTBME_Oslo(inputtbme, Hbare);
        else if (fmt2 == "oakridge" )
        { // input format should be: singleparticle.dat,vnn.dat
          size_t comma_pos = inputtbme.find_first_of(",");
          rw.ReadTBME_OakRidge( inputtbme.substr(0,comma_pos),  inputtbme.substr( comma_pos+1 ), Hbare);
        }
        else if (fmt2 == "takayuki" )
          rw.ReadTwoBody_Takayuki( inputtbme, Hbare);
        else if (fmt2 == "nushellx" )
          rw.ReadNuShellX_int( Hbare, inputtbme );
  
        cout << "done reading 2N" << endl;
      }
      if (fmt2 != "nushellx")
        Hbare += Trel_Op(modelspace);
//    }
  
//    #pragma omp section
    if (Hbare.particle_rank >=3)
    {
//      omp_set_num_threads(1);
      rw.Read_Darmstadt_3body(input3bme, Hbare, file3e1max,file3e2max,file3e3max);
      cout << "done reading 3N" << endl;
    }
//  }
//  }
//  omp_set_num_threads(nthreads);
//  omp_set_nested(0);

//  if (fmt2 != "nushellx")
//    Hbare += Trel_Op(modelspace);

  // Add a Lawson term. If hwBetaCM is specified, use that frequency
  if (abs(BetaCM)>1e-3)
  {
    if (hwBetaCM < 0) hwBetaCM = modelspace.GetHbarOmega();
    ostringstream hcm_opname;
    hcm_opname << "HCM_" << hwBetaCM;
    Hbare += BetaCM * imsrg_util::OperatorFromString( modelspace, hcm_opname.str());
  }

  cout << "Creating HF" << endl;
  HartreeFock hf(Hbare);
  cout << "Solving" << endl;
  hf.Solve();
//  cout << "EHF = " << hf.EHF << endl;
  
  Operator HNO;
  if (basis == "HF" and method !="HF")
    HNO = hf.GetNormalOrderedH();
  else if (basis == "oscillator")
    HNO = Hbare.DoNormalOrdering();


  int n_radial_points = 40;
  double Rmax = 10.0;
  vector<index_t> spwf_indices = modelspace.String2Index(spwf);
  vector<double> R(n_radial_points);
  vector<double> PSI(n_radial_points);
  for ( index_t i=0; i< spwf.size(); ++i)
  {
    for (int rstep=0;rstep<n_radial_points;++rstep) R[rstep] = Rmax/n_radial_points * rstep;
    hf.GetRadialWF(spwf_indices[i], R, PSI);
    ofstream wf_file (intfile + "_spwf_" + spwf[i] + ".dat");
    for ( index_t rstep=0; rstep<R.size(); ++rstep)  wf_file << fixed << setw(10) << setprecision(7) << R[rstep] << "   " << setw(10) << setprecision(7) << PSI[rstep] << endl;
    cout << "About to close wf file" << endl;
//    wf_file.close();
  }
  if (spwf.size() > 0)   cout << "Done with SPWF" << endl;

  HNO -= BetaCM * 1.5*hwBetaCM;
  cout << "Hbare 0b = " << HNO.ZeroBody << endl;

  if (method != "HF")
  {
    cout << "Perturbative estimates of gs energy:" << endl;
    double EMP2 = HNO.GetMP2_Energy();
    cout << "EMP2 = " << EMP2 << endl;
    double EMP3 = HNO.GetMP3_Energy();
    cout << "EMP3 = " << EMP3 << endl;
    cout << "To 3rd order, E = " << HNO.ZeroBody+EMP2+EMP3 << endl;
  }
*/

  cout<<endl<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;
  // Calculate all the desired operators
  for (auto& opname : opnames)
  {
    ops.emplace_back( imsrg_util::OperatorFromString(modelspace,opname) );
  }

/*
  for (auto& op : ops)
  {
     if (basis == "HF") op = hf.TransformToHFBasis(op);
     op = op.DoNormalOrdering();
     if (method == "MP3")
     {
       double dop = op.MP1_Eval( HNO );
       cout << "Operator 1st order correction  " << dop << "  ->  " << op.ZeroBody + dop << endl;
     }
  }
  auto itR2p = find(opnames.begin(),opnames.end(),"Rp2");
  if (itR2p != opnames.end())
  {
    Operator& Rp2 = ops[itR2p-opnames.begin()];
    int Z = modelspace.GetTargetZ();
    int A = modelspace.GetTargetMass();
    cout << " HF point proton radius = " << sqrt( Rp2.ZeroBody ) << endl;
    cout << " HF charge radius = " << ( abs(Rp2.ZeroBody)<1e-6 ? 0.0 : sqrt( Rp2.ZeroBody + r2p + r2n*(A-Z)/Z + DF) ) << endl;
  }
  for (index_t i=0;i<ops.size();++i)
  {
    cout << opnames[i] << " = " << ops[i].ZeroBody << endl;
  }
  
  if ( method == "HF" or method == "MP3")
  {
    HNO.PrintTimes();
    return 0;
  }
  cout << "HF Single particle energies:" << endl;
  hf.PrintSPE();
  cout << endl;


  if (method == "FCI")
  {
    HNO = HNO.UndoNormalOrdering();
    rw.WriteNuShellX_int(HNO,intfile+".int");
    rw.WriteNuShellX_sps(HNO,intfile+".sp");

    for (index_t i=0;i<ops.size();++i)
    {
      ops[i] = ops[i].UndoNormalOrdering();
      if ((ops[i].GetJRank()+ops[i].GetTRank()+ops[i].GetParity())<1)
      {
        rw.WriteNuShellX_op(ops[i],intfile+opnames[i]+".int");
      }
      else
      {
        rw.WriteTensorOneBody(intfile+opnames[i]+"_1b.op",ops[i],opnames[i]);
        rw.WriteTensorTwoBody(intfile+opnames[i]+"_2b.op",ops[i],opnames[i]);
      }
    }
    HNO.PrintTimes();
    return 0;
  }

//  Operator HlowT = HNO;
//  double Temp = hw;
//  double Efermi = 0;
//  Operator Eye = HNO;
//  Eye.Eye();
//  HlowT.ScaleFermiDirac(HNO, Temp, Efermi);  // 0 is roughly fermi surface? we can do beter...
//  Eye.ScaleFermiDirac(HNO, Temp, Efermi);  // 0 is roughly fermi surface? we can do beter...
//  cout << "Initial low temp trace with T = " << Temp << " and Ef = " << Efermi << ":   " << HlowT.Trace(modelspace.GetAref(),modelspace.GetZref()) <<"  with normalization  " << Eye.Trace( modelspace.GetAref(),modelspace.GetZref() ) << endl;

  IMSRGSolver imsrgsolver(HNO);
  imsrgsolver.SetReadWrite(rw);
  imsrgsolver.SetEtaCriterion(eta_criterion);
  bool brueckner_restart = false;
  
  if (method == "NSmagnus") // "No split" magnus
  {
    omega_norm_max=500;
    method = "magnus";
  }
  if (method.find("brueckner") != string::npos)
  {
    if (method=="brueckner2") brueckner_restart=true;
    if (method=="brueckner1step")
    {
       nsteps = 1;
       core_generator = valence_generator;
    }
    use_brueckner_bch = "true";
    omega_norm_max=500;
    method = "magnus";
  }

  if (use_brueckner_bch == "true" or use_brueckner_bch == "True")
  {
//    Hbare.SetUseBruecknerBCH(true);
    HNO.SetUseBruecknerBCH(true);
    cout << "Using Brueckner flavor of BCH" << endl;
  }

  imsrgsolver.SetMethod(method);
//  imsrgsolver.SetHin(Hbare);
  imsrgsolver.SetHin(HNO);
  imsrgsolver.SetSmax(smax);
  imsrgsolver.SetFlowFile(flowfile);
  imsrgsolver.SetDs(ds_0);
  imsrgsolver.SetDsmax(dsmax);
  imsrgsolver.SetDenominatorDelta(denominator_delta);
  imsrgsolver.SetdOmega(domega);
  imsrgsolver.SetOmegaNormMax(omega_norm_max);
  imsrgsolver.SetODETolerance(ode_tolerance);
  if (denominator_delta_orbit != "none")
    imsrgsolver.SetDenominatorDeltaOrbit(denominator_delta_orbit);

  imsrgsolver.SetGenerator(core_generator);
  if (core_generator.find("imaginary")!=string::npos)
  {
   if (ds_0>1e-2)
   {
     ds_0 = 1e-4;
     dsmax = 1e-2;
     imsrgsolver.SetDs(ds_0);
     imsrgsolver.SetDsmax(dsmax);
   }
  }
  imsrgsolver.Solve();

//  HlowT = imsrgsolver.Transform(HlowT);
//  cout << "After Solve, low temp trace with T = " << Temp << " and Ef = " << Efermi << ":   " << HlowT.Trace(modelspace.GetAref(),modelspace.GetZref()) << endl;

  if (method == "magnus")
  {
//    for (size_t i=0;i<ops.size();++i)
//    {
//      Operator tmp = imsrgsolver.Transform(ops[i]);
////      rw.WriteOperatorHuman(tmp,intfile+opnames[i]+"_step1.op");
//    }
//    cout << endl;
    // increase smax in case we need to do additional steps
    smax *= 1.5;
    imsrgsolver.SetSmax(smax);
  }


  if (brueckner_restart)
  {
     arma::mat newC = hf.C * arma::expmat( -imsrgsolver.GetOmega(0).OneBody  );
//     if (input3bme != "none") Hbare.SetParticleRank(3);
     HNO = hf.GetNormalOrderedH(newC);
     imsrgsolver.SetHin(HNO);
     imsrgsolver.s = 0;
     imsrgsolver.Solve();
  }

  if (nsteps > 1) // two-step decoupling, do core first
  {
    if (method == "magnus") smax *= 2;

    imsrgsolver.SetGenerator(valence_generator);
    modelspace.ResetFirstPass();
    if (valence_generator.find("imaginary")!=string::npos)
    {
     if (ds_0>1e-2)
     {
       ds_0 = 1e-4;
       dsmax = 1e-2;
       imsrgsolver.SetDs(ds_0);
       imsrgsolver.SetDsmax(dsmax);
     }
    }
    imsrgsolver.SetSmax(smax);
    imsrgsolver.Solve();
  }



  // Transform all the operators
  if (method == "magnus")
  {
    if (ops.size()>0) cout << "transforming operators" << endl;
    for (size_t i=0;i<ops.size();++i)
    {
      cout << opnames[i] << " " << endl;
      ops[i] = imsrgsolver.Transform(ops[i]);
      cout << " (" << ops[i].ZeroBody << " ) " << endl;
//      rw.WriteOperatorHuman(ops[i],intfile+opnames[i]+"_step2.op");
    }
    cout << endl;
    // increase smax in case we need to do additional steps
    smax *= 1.5;
    imsrgsolver.SetSmax(smax);
  }


  // If we're doing targeted/ensemble normal ordering
  // we now re-normal order wrt to the core
  // and do any remaining flow.
  ModelSpace ms2(modelspace);
  bool renormal_order = false;
  if (modelspace.valence.size() > 0 )
  {
    renormal_order = modelspace.holes.size() != modelspace.core.size();
    if (not renormal_order)
    {
      for (auto c : modelspace.core)
      {
         if ( (find( modelspace.holes.begin(), modelspace.holes.end(), c) == modelspace.holes.end()) or (abs(1-modelspace.holes[c])>1e-6))
         {
           renormal_order = true;
           break;
         }
      }
    }
  }
  if ( renormal_order )
  {

    HNO = imsrgsolver.GetH_s();

    int nOmega = imsrgsolver.GetOmegaSize() + imsrgsolver.GetNOmegaWritten();
    cout << "Undoing NO wrt A=" << modelspace.GetAref() << " Z=" << modelspace.GetZref() << endl;
    HNO = HNO.UndoNormalOrdering();

    ms2.SetReference(ms2.core); // change the reference
    HNO.SetModelSpace(ms2);

    cout << "Doing NO wrt A=" << ms2.GetAref() << " Z=" << ms2.GetZref() << "  norbits = " << ms2.GetNumberOrbits() << endl;
    HNO = HNO.DoNormalOrdering();

    imsrgsolver.SetHin(HNO);
    imsrgsolver.SetEtaCriterion(1e-4);
    imsrgsolver.Solve();
    // Change operators to the new basis, then apply the rest of the transformation
    cout << "Final transformation on the operators..." << endl;
    for (auto& op : ops)
    {
//      double ZeroBody_before = op.ZeroBody;
      op = op.UndoNormalOrdering();
//      double ZeroBody_undo = op.ZeroBody;
      op.SetModelSpace(ms2);
      op = op.DoNormalOrdering();
//      double ZeroBody_mid = op.ZeroBody;
      // transform using the remaining omegas
      op = imsrgsolver.Transform_Partial(op,nOmega);
//      cout << ZeroBody_before << "   =>   " << ZeroBody_undo << "   =>   " << ZeroBody_mid<< "   =>   " << op.ZeroBody << endl;
    }
  }


  // Write the output

  // If we're doing a shell model interaction, write the
  // interaction files to disk.
  if (modelspace.valence.size() > 0)
  {
    if (valence_file_format == "antoine") // this is still being tested...
    {
      rw.WriteAntoine_int(imsrgsolver.GetH_s(),intfile+".ant");
      rw.WriteAntoine_input(imsrgsolver.GetH_s(),intfile+".inp",modelspace.GetAref(),modelspace.GetZref());
    }
    cout << "Writing files: " << intfile << endl;
    rw.WriteNuShellX_int(imsrgsolver.GetH_s(),intfile+".int");
    rw.WriteNuShellX_sps(imsrgsolver.GetH_s(),intfile+".sp");

    if (method == "magnus")
    {
       for (index_t i=0;i<ops.size();++i)
       {
          if ((ops[i].GetJRank()+ops[i].GetTRank()+ops[i].GetParity())<1)
          {
            rw.WriteNuShellX_op(ops[i],intfile+opnames[i]+".int");
          }
          else
          {
            rw.WriteTensorOneBody(intfile+opnames[i]+"_1b.op",ops[i],opnames[i]);
            rw.WriteTensorTwoBody(intfile+opnames[i]+"_2b.op",ops[i],opnames[i]);
          }
       }
    }
  }
  else // single ref. just print the zero body pieces out. (maybe check if its magnus?)
  {
    cout << "Core Energy = " << setprecision(6) << imsrgsolver.GetH_s().ZeroBody << endl;
    for (index_t i=0;i<ops.size();++i)
    {
      Operator& op = ops[i];
      cout << opnames[i] << " = " << ops[i].ZeroBody << endl;
      if ( opnames[i] == "Rp2" )
      {
         int Z = modelspace.GetTargetZ();
         int A = modelspace.GetTargetMass();
         cout << " IMSRG point proton radius = " << sqrt( op.ZeroBody ) << endl;
         cout << " IMSRG charge radius = " << sqrt( op.ZeroBody + r2p + r2n*(A-Z)/Z + DF) << endl;
      }
      if (op.GetJRank()>0) // if it's a tensor, you probably want the full operator
      {
        cout << "Writing operator to " << intfile+opnames[i]+".op" << endl;
        rw.WriteOperatorHuman(op,intfile+opnames[i]+".op");
      }
    }
  }
*/


  cout<<endl<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;
  int printCSdS = 0; // to print CharlieSdS
  int printRSdS = 0; // to print RagnarSdS
  int printBME = 0; // to print BME_rho
  int printM0nu = 0; // to print M0nu_TBME
  int printM0nuNEW = 0; // to print M0nu_qagiu
  int printM0nuJE = 0; // to print M0nu_JE
  int printM0nuAdpt = 0; // to print M0nu_adpt
  int printTime = 0;
  int printJE = 0; // to print Engel
  int printMH = 0; // to print Mihai
  for (auto& opname : opnames)
  {
    if (opname == "CharlieSdS")
    {
      printCSdS = 1;
      printTime = 1;
    }
    if (opname == "RagnarSdS")
    {
      printRSdS = 1;
      printTime = 1;
    }
    if (opname.substr(0,8) == "BME_rho_")
    {
      printBME = 1;
      printTime = 1;
    }
    if (opname.substr(0,10) == "M0nu_TBME_")
    {
      printM0nu = 1;
      printTime = 1;
    }
    if (opname.substr(0,11) == "M0nu_qagiu_")
    {
      printM0nuNEW = 1;
      printTime = 1;
    }
    if (opname.substr(0,8) == "M0nu_JE_")
    {
      printM0nuJE = 1;
      printTime = 1;
    }
    if (opname.substr(0,10) == "M0nu_adpt_")
    {
      printM0nuAdpt = 1;
      printTime = 1;
    }
    if (opname.substr(0,5) == "Engel")
    {
      printJE = 1;
      printTime = 0;
    }
    if (opname.substr(0,5) == "Mihai")
    {
      printMH = 1;
      printTime = 0;
    }
  }

  if (printCSdS == 1 && printRSdS == 1)
  {
    // print difference of operators to file
    ifstream cpsds("/home/cgpayne/cougwork/imsrg_code/ragnar_imsrg/work/scripts/output_runCPops/cpSdS.txt");
    if(!cpsds)
    {
      cout<<"Error: with opening cpsds"<<endl;
      exit(1);
    }
    cout<<endl;
    ifstream rssds("/home/cgpayne/cougwork/imsrg_code/ragnar_imsrg/work/scripts/output_runCPops/rsSdS.txt");
    if(!rssds)
    {
      cout<<"Error: with opening cssds"<<endl;
      exit(1);
    }
    double CP[Ntb];
    double RS[Ntb];
    for (int i=0; i<Ntb; i++)
    {
      cpsds>>CP[i];
      rssds>>RS[i];
    }
    ofstream diff("/home/cgpayne/cougwork/imsrg_code/ragnar_imsrg/work/scripts/output_runCPops/the_diff.txt");
    if(!diff)
    {
      cout<<"Error: with diff"<<endl;
      exit(1);
    }
    diff<<"  abs(cpSdS - rsSdS)  |  cpSdS  |  rsSdS"<<endl;
    diff<<"----------------------|---------|-------"<<endl;
    double sum = 0;
    for (int i=0; i<Ntb; i++)
    {
      double delta = abs(CP[i]-RS[i]);
      sum += delta;
      diff<<"  "<<delta<<"                   |  "<<CP[i]<<"  |  "<<RS[i]<<endl;
    }
    diff<<"~~~~~~~~~~~~~~~~~~~~~~|~~~~~~~~~|~~~~~~~"<<endl;
    diff<<"sum_i[abs(CP[i]-RS[i])] = "<<sum<<endl;
    cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
    cout<<"sum_i[abs(CP[i]-RS[i])] = "<<sum<<endl;
    cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl<<endl;
  }
  
  if (printBME == 1)
  {
    rw.WriteOperator(ops[0], intfile+opnames[0]+".txt"); // this using Ragnar's way of printing operators
  }
  
  if (printM0nu == 1)
  {
    // get the max, min, avg
    string Estr = "_emax_"+std::to_string(eMax);
    string Astr = "_A_"+std::to_string(targetMass);
    stringstream HWstream;
    HWstream<<fixed<<setprecision(2)<<modelspace.GetHbarOmega();
    string HWstr = "_hw_"+HWstream.str();
    ifstream Fin(intfile+opnames[0]+"_F"+Estr+Astr+HWstr+".txt");
    if(!Fin)
    {
      cout<<"Error: with opening Fin"<<endl;
      exit(1);
    }
    ifstream GTin(intfile+opnames[0]+"_GT"+Estr+Astr+HWstr+".txt");
    if(!GTin)
    {
      cout<<"Error: with opening GTin"<<endl;
      exit(1);
    }
    double M0F[Ntb];
    double M0GT[Ntb];
    for (int i=0; i<Ntb; i++)
    {
      Fin>>M0F[i];
      GTin>>M0GT[i];
//cout<<M0F[i]<<"\t"<<M0GT[i]<<endl;
    }
    double minF = INFINITY;
    double maxF = -INFINITY;
    double sumF = 0;
    int moticF = 0;
    double minGT = INFINITY;
    double maxGT = -INFINITY;
    double sumGT = 0;
    int moticGT = 0;
    double motol = pow(10,-7);
    for (int i=0; i<Ntb; i++)
    {
      // F
      if (M0F[i] < minF)
      {
        minF = M0F[i];
      }
      if (maxF < M0F[i])
      {
        maxF = M0F[i];
      }
      if (abs(M0F[i]) > motol)
      {
        sumF += M0F[i];
        moticF++;
      }
      // GT
      if (M0GT[i] < minGT)
      {
        minGT = M0GT[i];
      }
      if (maxGT < M0GT[i])
      {
        maxGT = M0GT[i];
      }
      if (abs(M0GT[i]) > motol)
      {
        sumGT += M0GT[i];
        moticGT++;
      }
    }
    ofstream mmaF(intfile+opnames[0]+"_F"+Estr+Astr+HWstr+"_mma.txt");
    if(!mmaF)
    {
      cout<<"Error: with mmaF"<<endl;
      exit(1);
    }
    ofstream mmaGT(intfile+opnames[0]+"_GT"+Estr+Astr+HWstr+"_mma.txt");
    if(!mmaGT)
    {
      cout<<"Error: with mmaGT"<<endl;
      exit(1);
    }
    mmaF<<setprecision(12)<<"min(M0F) = "<<minF<<endl;
    mmaF<<setprecision(12)<<"max(M0F) = "<<maxF<<endl;
    mmaF<<setprecision(12)<<"AVG(M0F) = "<<sumF/moticF<<endl;
    mmaF<<"the number of non-zero elements in M0F = "<<moticF<<endl;
    mmaGT<<setprecision(12)<<"min(M0GT) = "<<minGT<<endl;
    mmaGT<<setprecision(12)<<"max(M0GT) = "<<maxGT<<endl;
    mmaGT<<setprecision(12)<<"AVG(M0GT) = "<<sumGT/moticGT<<endl;
    mmaGT<<"the number of non-zero elements in M0GT = "<<moticGT<<endl;
    // also compare with how Ragnar outputs it for F+GT
    rw.WriteOperatorHuman(ops[0], intfile+opnames[0]+Estr+Astr+HWstr+"_ragtime.txt");
  }
  
  if (printM0nuNEW == 1)
  {
    string Estr = "_emax_"+std::to_string(eMax);
    string Astr = "_A_"+std::to_string(targetMass);
    stringstream HWstream;
    HWstream<<fixed<<setprecision(2)<<modelspace.GetHbarOmega();
    string HWstr = "_hw_"+HWstream.str();
    rw.WriteOperatorHuman(ops[0], intfile+opnames[0]+Estr+Astr+HWstr+"_ragtime.txt");
  }
  
  if (printM0nuJE == 1)
  {
    string Estr = "_emax_"+std::to_string(eMax);
    string Astr = "_A_"+std::to_string(targetMass);
    stringstream HWstream;
    HWstream<<fixed<<setprecision(2)<<modelspace.GetHbarOmega();
    string HWstr = "_hw_"+HWstream.str();
    rw.WriteOperatorHuman(ops[0], intfile+opnames[0]+Estr+Astr+HWstr+"_ragtime.txt");
  }

  if (printM0nuAdpt == 1)
  {
    rw.WriteM0nu(modelspace,ops[0],intfile,opnames[0]);
  }
  
  if (printJE == 1)
  {
    int checkJE = 0; // default to ReadTwoBodyEngel
    cout<<endl<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
    cout<<"are you using:"<<endl;
    cout<<"0 = ReadTwoBodyEngel"<<endl;
    cout<<"1 = ReadTwoBodyNewEngel"<<endl;
    cout<<"enter your choice:"<<endl;
    cin>>checkJE;
    
    Operator Engel = Operator(modelspace,0,2,0,2);
    string Estr = "_emax_"+std::to_string(eMax);
    if (checkJE == 0)
    {
      cout<<"old"<<endl<<endl;
      rw.ReadTwoBodyEngel(intfile+".op",Engel);
      rw.WriteOperatorHuman(Engel,intfile+Estr+"_ragtime.op");
    }
    else if (checkJE == 1)
    {
      cout<<"new"<<endl<<endl;
      rw.ReadTwoBodyNewEngel(intfile+".op",Engel);
      rw.WriteOperatorHuman(Engel,intfile+Estr+"_ragtime.op");
    }
    else
    {
      cout<<"JEngel conversion failed..."<<endl;
      cout<<"checkJE = "<<checkJE<<endl;
    }
    
    Engel.PrintTimes();
  }
  
  if (printMH == 1)
  {
    Operator Mihai = Operator(modelspace,0,2,0,2);
    string Estr = "_emax_"+std::to_string(eMax);
    rw.ReadTwoBodyMihai(intfile+".dat",Mihai);
    rw.WriteOperatorHuman(Mihai,intfile+Estr+"_ragtime.dat");
    
    Mihai.PrintTimes();
  }

  if (printTime == 1)
  {
    ops[0].PrintTimes();
  }
 
  return 0;
}

