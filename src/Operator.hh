
#ifndef Operator_h
#define Operator_h 1

#include "ModelSpace.hh"
#include "TwoBodyME.hh"
#include "ThreeBodyME.hh"
#include "IMSRGProfiler.hh"
#include <armadillo>
#include <string>
#include <vector>
#include <map>

using namespace std;

//typedef uint64_t orbindx3_t;

///
/// The Operator class provides a generic operator up to three-body, scalar or tensor.
/// The class contains lots of methods and overloaded operators so that the resulting
/// code that uses the operators can look as close as possible to the math that is
/// written down.
class Operator
{
 public:
  //Fields
  ModelSpace * modelspace;
  double ZeroBody; ///< The zero body piece of the operator.
  arma::mat OneBody; ///< The one body piece of the operator, stored in a single NxN armadillo matrix, where N is the number of single-particle orbits.
  TwoBodyME TwoBody; ///< The two body piece of the operator.
  ThreeBodyME ThreeBody; ///< The three body piece of the operator.

  int rank_J; ///< Spherical tensor rank of the operator
  int rank_T; ///< Isotensor rank of the operator
  int parity; ///< Parity of the operator, 0=even 1=odd
  int particle_rank; ///< Maximum particle rank. Should be 2 or 3.

  int E2max; ///< For two-body matrix elements, \f$ e_i + e_j \leq \f$ E2max
  int E3max; ///< For three-body matrix elements, \f$ e_i + e_j + e_k \leq \f$ E3max

  bool hermitian;
  bool antihermitian;
  int nChannels; ///< Number of two-body channels \f$ J,\pi,T_z \f$ associated with the model space
  static double bch_transform_threshold;
  static double bch_product_threshold;
  static map<string, double> timer; ///< For keeping timing information for various method calls
  map<array<int,3>,vector<index_t> > OneBodyChannels;
  IMSRGProfiler profiler;

  void PrintTimes();


  //Constructors
  // In the future, consider using C++11 rvalues / move constructor to avoid copies in certain cases
  ~Operator();
  Operator(); ///< Default constructor
  Operator(ModelSpace&); ///< Construct a 2-body scalar operator
  Operator(ModelSpace&, int Jrank, int Trank, int Parity, int part_rank);
  Operator( const Operator& rhs); ///< Copy constructor
  Operator( Operator&&);

  //Overloaded operators
  Operator& operator=( const Operator& rhs);
  Operator& operator+=( const Operator& rhs);
  Operator operator+( const Operator& rhs) const;
  Operator& operator-=( const Operator& rhs);
  Operator operator-( const Operator& rhs) const;
  Operator operator-( ) const;
  Operator& operator*=( const double rhs);
  Operator operator*( const double rhs) const;
  Operator& operator/=( const double rhs);
  Operator operator/( const double rhs) const;

  Operator& operator=(Operator&& rhs);

  //Methods
  void Copy(const Operator& rhs);

  // One body setter/getters
  double GetOneBody(int i,int j) {return OneBody(i,j);};
  void SetOneBody(int i, int j, double val) { OneBody(i,j) = val;};
  int GetTwoBodyDimension(int ch_bra, int ch_ket){ return TwoBody.GetMatrix(ch_bra, ch_ket).n_cols;};
  double GetTwoBody(int ch_bra, int ch_ket, int i, int j){ return TwoBody.GetMatrix(ch_bra, ch_ket)(i,j);};

  void SetE3max(int e){E3max = e;};
  int GetE3max(){return E3max;};

  // Other setter-getters
  ModelSpace * GetModelSpace();
  void SetModelSpace(ModelSpace &ms){modelspace = &ms;};

  void Erase(); ///< Set all matrix elements to zero.
  void EraseZeroBody(){ZeroBody = 0;}; // set zero-body term to zero
  void EraseOneBody(); // set all one-body terms to zero
  void EraseTwoBody(); // set all two-body terms to zero
  void EraseThreeBody(); // set all two-body terms to zero

  void SetHermitian() ;
  void SetAntiHermitian() ;
  void SetNonHermitian() ;
  bool IsHermitian()const {return hermitian;};
  bool IsAntiHermitian()const {return antihermitian;};
  bool IsNonHermitian()const {return not (hermitian or antihermitian);};
  int GetParticleRank()const {return particle_rank;};
  int GetJRank()const {return rank_J;};
  int GetTRank()const {return rank_T;};
  int GetParity()const {return parity;};
  void SetParticleRank(int pr) {particle_rank = pr;};

  void MakeReduced();
  void MakeNotReduced();

  void ScaleZeroBody(double x);
  void ScaleOneBody(double x);
  void ScaleTwoBody(double x);
  void Symmetrize(); ///< Copy the upper-half triangle to the lower-half triangle for each matrix
  void AntiSymmetrize(); ///< Copy the upper-half triangle to the lower-half triangle with a minus sign.
  void SetUpOneBodyChannels();
  int Size();

  // The actually interesting methods
  Operator DoNormalOrdering(); ///< Calls DoNormalOrdering2() or DoNormalOrdering3(), depending on the rank of the operator.
  Operator DoNormalOrdering2(); ///< Returns the normal ordered two-body operator
  Operator DoNormalOrdering3(); ///< Returns the normal ordered three-body operator
  Operator UndoNormalOrdering(); ///< Returns the operator normal-ordered wrt the vacuum

  friend Operator Commutator(const Operator& X, const Operator& Y) ; 
  friend Operator CommutatorScalarScalar( const Operator& X, const Operator& Y) ;
  friend Operator CommutatorScalarTensor( const Operator& X, const Operator& Y) ;

  Operator BCH_Product(  Operator& )  ; 
  Operator BCH_Transform( const Operator& ) ; 

  void CalculateKineticEnergy();
  void Eye(); ///< set to identity operator

  double GetMP2_Energy();


  double Norm() const;
  double OneBodyNorm() const;
  double TwoBodyNorm() const;


  void PrintOneBody() const {OneBody.print();};
  void PrintTwoBody(int ch) const {TwoBody.PrintMatrix(ch,ch);};


  //Methods

  static void Set_BCH_Transform_Threshold(double x){bch_transform_threshold=x;};
  static void Set_BCH_Product_Threshold(double x){bch_product_threshold=x;};

  void DoPandyaTransformation(vector<arma::mat>&, vector<arma::mat>&) const ;
  void AddInversePandyaTransformation(vector<arma::mat>&);


  void comm110ss( const Operator& X, const Operator& Y) ; 
  void comm220ss( const Operator& X, const Operator& Y) ;
  void comm111ss( const Operator& X, const Operator& Y) ;
  void comm121ss( const Operator& X, const Operator& Y) ;
  void comm221ss( const Operator& X, const Operator& Y) ;
  void comm122ss( const Operator& X, const Operator& Y) ;
  void comm222_pp_hhss( const Operator& X, const Operator& Y) ;
  void comm222_phss( const Operator& X, const Operator& Y) ;
  void comm222_pp_hh_221ss( const Operator& X, const Operator& Y) ;

// scalar-tensor commutators

  void DoTensorPandyaTransformation(map<array<int,2>,arma::mat>&, map<array<int,2>,arma::mat>&) const;
  void AddInverseTensorPandyaTransformation(map<array<int,2>,arma::mat>&);

  void comm111st( const Operator& X, const Operator& Y) ;
  void comm121st( const Operator& X, const Operator& Y) ;
  void comm122st( const Operator& X, const Operator& Y) ;
  void comm222_pp_hh_221st( const Operator& X, const Operator& Y) ;
  void comm222_phst( const Operator& X, const Operator& Y) ;

};

/// Non member function, multiply by scalar from left side
Operator operator*(const double lhs, const Operator& rhs);
Operator operator*(const double lhs, const Operator&& rhs);



#endif

