//
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                          HEP Random
//                       --- MixMaxRng ---
//                       class header file
// -----------------------------------------------------------------------
//
// This file interfaces the MixMax PseudoRandom Number Generator 
// proposed by:
//
// G.K.Savvidy and N.G.Ter-Arutyunian,
//   On the Monte Carlo simulation of physical systems,
//   J.Comput.Phys. 97, 566 (1991);
//   Preprint EPI-865-16-86, Yerevan, Jan. 1986
//   http://dx.doi.org/10.1016/0021-9991(91)90015-D
//
// K.Savvidy
//   "The MIXMAX random number generator"
//   Comp. Phys. Commun. (2015)
//   http://dx.doi.org/10.1016/j.cpc.2015.06.003
//
// K.Savvidy and G.Savvidy
//   "Spectrum and Entropy of C-systems. MIXMAX random number generator"
//   Chaos, Solitons & Fractals, Volume 91, (2016) pp. 33-38
//   http://dx.doi.org/10.1016/j.chaos.2016.05.003
//
// =======================================================================
// Implementation by Konstantin Savvidy - Copyright 2004-2023
// July 2023 - Updated class structure upon suggestions from Marco Barbone
// =======================================================================

#ifndef MixMaxRng_h
#define MixMaxRng_h 1

#include <array>
#include <cstdint>

#include "CLHEP/Random/RandomEngine.h"

namespace CLHEP {
    
/**
  * @author  K.Savvidy
  * @ingroup random
  */

using myID_t = std::uint32_t;
using myuint_t = std::uint64_t;

class alignas(128) MixMaxRng: public HepRandomEngine {

  static const int N = 17;

public:

  MixMaxRng(std::istream& is);
  MixMaxRng();
  MixMaxRng(long seed);
  ~MixMaxRng();
  // Constructors and destructor.

  MixMaxRng(const MixMaxRng& rng);
  MixMaxRng& operator=(const MixMaxRng& rng);
  // Copy constructor and assignment operator.

  inline double flat()
  {
    if (counter >= N) iterate();  
#if defined(__x86_64__)
      return convert1double(V[counter++]);
#else
      return INV_M61*static_cast<double>(V[counter++]);
#endif
  }
  // Returns a pseudo random number between 0 and 1
  // excluding the zero: in (0,1] 
  // smallest number which it will give is approximately 10^-19

  void flatArray (const int size, double* vect);
  // Fills the array "vect" of specified size with flat random values.

  inline void setSeed(long longSeed, int = 0 /* extraSeed */)
  {
    seed_spbox(theSeed = longSeed);
  }
  // Sets the state of the algorithm according to seed.

  void setSeeds(const long * seeds, int seedNum=0);
  // Sets the initial state of the engine according to the array of between one and four 32-bit seeds.
  // If the size of long is greater on the platform, only the lower 32-bits are used.
  // Streams created from seeds differing by at least one bit somewhere are guaranteed absolutely
  // to be independent and non-colliding for at least the next 10^100 random numbers

  void saveStatus( const char filename[] = "MixMaxRngState.conf" ) const;
  // Saves the the current engine state in the file given, by default MixMaxRngState.conf

  void restoreStatus( const char filename[] = "MixMaxRngState.conf" );
  // Reads a valid engine state from a given file, by default MixMaxRngState.conf
  // and restores it.

  void showStatus() const;
  // Dumps the engine status on the screen.

  operator double(){ return flat(); }
  // Returns same as flat()
  operator float();
  // less precise flat, faster if possible
  operator unsigned int();
  // 32-bit flat

  virtual std::ostream & put (std::ostream & os) const;
  virtual std::istream & get (std::istream & is);
  static  std::string beginTag ( );
  virtual std::istream & getState ( std::istream & is );

  std::string name() const { return "MixMaxRng"; }
  static std::string engineName();

  std::vector<unsigned long> put () const;
  bool get (const std::vector<unsigned long> & v);
  bool getState (const std::vector<unsigned long> & v);

private:

  static constexpr long long int SPECIAL   = ((N==17)? 0 : ((N==240)? 487013230256099140ULL:0) ); // etc...
  static constexpr long long int SPECIALMUL= ((N==17)? 36: ((N==240)? 51                   :53) ); // etc...
  // Note the potential for confusion...
  static constexpr int BITS=61;
  static constexpr myuint_t M61=2305843009213693951ULL;
  static constexpr double INV_M61=0.43368086899420177360298E-18;
  static constexpr unsigned int VECTOR_STATE_SIZE = 2*N+4; // 2N+4 for MIXMAX

  inline myuint_t MIXMAX_MOD_MERSENNE(myuint_t k)
  {
    return ((((k)) & M61) + (((k)) >> BITS) );
  }

  static constexpr int rng_get_N();
  static constexpr long long int rng_get_SPECIAL();
  static constexpr int rng_get_SPECIALMUL();
  void seed_uniquestream( myID_t clusterID, myID_t machineID, myID_t runID, myID_t  streamID );
  void seed_spbox(myuint_t seed);
  void print_state() const;
  myuint_t precalc();
  myuint_t get_next();
  
  MixMaxRng Branch();
  void BranchInplace(int id);

  MixMaxRng(myID_t clusterID, myID_t machineID, myID_t runID, myID_t  streamID ); // Constructor with four 32-bit seeds
  inline void seed64(myuint_t seedval)
  {
    seed_uniquestream( 0, 0, (myID_t)(seedval>>32), (myID_t)seedval );
    // seed with one 64-bit seed
  }

  inline void iterate()
  {
    myuint_t  tempP, tempV;
    V[0] = ( tempV = sumtot);
    myuint_t insumtot=V[0], ovflow = 0; // will keep a running sum of all new elements
    tempP = 0;              // will keep a partial sum of all old elements
    myuint_t tempPO;
    tempPO = MULWU(tempP); tempP = modadd(tempP, V[1] ); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); V[1]  = tempV; insumtot += tempV; if (insumtot < tempV) {ovflow++;};
    tempPO = MULWU(tempP); tempP = modadd(tempP, V[2] ); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); V[2]  = tempV; insumtot += tempV; if (insumtot < tempV) {ovflow++;};
    tempPO = MULWU(tempP); tempP = modadd(tempP, V[3] ); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); V[3]  = tempV; insumtot += tempV; if (insumtot < tempV) {ovflow++;};
    tempPO = MULWU(tempP); tempP = modadd(tempP, V[4] ); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); V[4]  = tempV; insumtot += tempV; if (insumtot < tempV) {ovflow++;};
    tempPO = MULWU(tempP); tempP = modadd(tempP, V[5] ); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); V[5]  = tempV; insumtot += tempV; if (insumtot < tempV) {ovflow++;};
    tempPO = MULWU(tempP); tempP = modadd(tempP, V[6] ); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); V[6]  = tempV; insumtot += tempV; if (insumtot < tempV) {ovflow++;};
    tempPO = MULWU(tempP); tempP = modadd(tempP, V[7] ); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); V[7]  = tempV; insumtot += tempV; if (insumtot < tempV) {ovflow++;};
    tempPO = MULWU(tempP); tempP = modadd(tempP, V[8] ); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); V[8]  = tempV; insumtot += tempV; if (insumtot < tempV) {ovflow++;};
    tempPO = MULWU(tempP); tempP = modadd(tempP, V[9] ); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); V[9]  = tempV; insumtot += tempV; if (insumtot < tempV) {ovflow++;};
    tempPO = MULWU(tempP); tempP = modadd(tempP, V[10]); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); V[10] = tempV; insumtot += tempV; if (insumtot < tempV) {ovflow++;};
    tempPO = MULWU(tempP); tempP = modadd(tempP, V[11]); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); V[11] = tempV; insumtot += tempV; if (insumtot < tempV) {ovflow++;};
    tempPO = MULWU(tempP); tempP = modadd(tempP, V[12]); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); V[12] = tempV; insumtot += tempV; if (insumtot < tempV) {ovflow++;};
    tempPO = MULWU(tempP); tempP = modadd(tempP, V[13]); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); V[13] = tempV; insumtot += tempV; if (insumtot < tempV) {ovflow++;};
    tempPO = MULWU(tempP); tempP = modadd(tempP, V[14]); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); V[14] = tempV; insumtot += tempV; if (insumtot < tempV) {ovflow++;};
    tempPO = MULWU(tempP); tempP = modadd(tempP, V[15]); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); V[15] = tempV; insumtot += tempV; if (insumtot < tempV) {ovflow++;};
    tempPO = MULWU(tempP); tempP = modadd(tempP, V[16]); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); V[16] = tempV; insumtot += tempV; if (insumtot < tempV) {ovflow++;};
    sumtot = MIXMAX_MOD_MERSENNE(MIXMAX_MOD_MERSENNE(insumtot) + (ovflow <<3 ));

    counter=1;
  }

#if defined __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#pragma GCC diagnostic ignored "-Wuninitialized"
#endif
  inline double convert1double(myuint_t u)
  {
    const double one = 1;
    const myuint_t onemask = *(myuint_t*)&one;
    myuint_t tmp = (u>>9) | onemask; // bits between 52 and 62 dont affect the result!
    double d = *(double*)&tmp;
    return d-1.0;
  }
#if defined __GNUC__
#pragma GCC diagnostic pop
#endif
  myuint_t MOD_MULSPEC(myuint_t k);
  inline myuint_t MULWU(myuint_t k)
  {
    return (( (k)<<(SPECIALMUL) & M61) ^ ( (k) >> (BITS-SPECIALMUL))  );
  }
  myuint_t iterate_raw_vec(myuint_t* Y, myuint_t sumtotOld);
  myuint_t apply_bigskip(myuint_t* Vout, myuint_t* Vin, myID_t clusterID, myID_t machineID, myID_t runID, myID_t  streamID );

  inline myuint_t modadd(myuint_t xfoo, myuint_t xbar)
  {
#if defined(__x86_64__) && defined(__GNUC__) && (!defined(__ICC))
    myuint_t out;
    /* Assembler trick suggested by Andrzej Goerlich     */
    __asm__ ("addq %2, %0; "
             "btrq $61, %0; "
             "adcq $0, %0; "
             :"=r"(out)
             :"0"(xfoo), "r"(xbar)
            );
    return out;
#else
    return MIXMAX_MOD_MERSENNE(xfoo+xbar);
#endif
  }

#if defined(__x86_64__)
  myuint_t mod128(__uint128_t s);
  myuint_t fmodmulM61(myuint_t cum, myuint_t a, myuint_t b);
#else // on all other platforms, including 32-bit linux, PPC and PPC64, ARM and all Windows
  myuint_t fmodmulM61(myuint_t cum, myuint_t s, myuint_t a);
#endif

private:

  myuint_t V[N];
  myuint_t sumtot;
  int counter;
};

}  // namespace CLHEP

#endif
