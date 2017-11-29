
#ifndef __testing_h__
#define __testing_h__

// Order of header inclusion matters here: ShadowNumbers need to
// precede DualNumbers which need to precede *Vector

const unsigned int NDIM = 2;

#define TESTPRINT(a) do { std::cout << #a " = " << a << std::endl; } while (0)

#ifdef USE_SPARSE
  #ifdef USE_SHADOW
    #if defined(USE_STRUCT)
      #include "metaphysicl/dualshadowsparsestruct.h"
    #elif defined(USE_DYNAMIC)
      #include "metaphysicl/dualshadowdynamicsparsevector.h"
    #else // !USE_STRUCT && !USE_DYNAMIC
      #include "metaphysicl/dualshadowsparsevector.h"
    #endif // USE_STRUCT, USE_DYNAMIC
    typedef MetaPhysicL::ShadowNumber<double, long double> RawScalar;
  #else // USE_SHADOW
    #if defined(USE_STRUCT)
      #include "metaphysicl/dualsparsenumberstruct.h"
    #elif defined(USE_DYNAMIC)
      #include "metaphysicl/dualdynamicsparsenumbervector.h"
    #else // !USE_STRUCT && !USE_DYNAMIC
      #include "metaphysicl/dualsparsenumbervector.h"
    #endif // USE_STRUCT, USE_DYNAMIC
    typedef double RawScalar;
  #endif // USE_SHADOW
  #if defined(USE_STRUCT)
    typedef
    MetaPhysicL::SetConstructor<MetaPhysicL::UnsignedIntType<0,RawScalar>, MetaPhysicL::UnsignedIntType<1,RawScalar> >::type IndexSet;
    typedef MetaPhysicL::SparseNumberStruct<IndexSet> RawVector;
    #define VectorUnitVector SparseNumberStructUnitVector
    #define VectorOf SparseNumberStructOf
  #elif defined(USE_DYNAMIC)
    typedef MetaPhysicL::DynamicSparseNumberVector<RawScalar, unsigned int> RawVector;
    #define VectorUnitVector DynamicSparseNumberVectorUnitVector
    #define VectorOf DynamicSparseNumberVectorOf
  #else // !USE_STRUCT && !USE_DYNAMIC
    typedef MetaPhysicL::SetConstructor<MetaPhysicL::UnsignedIntType<0>, MetaPhysicL::UnsignedIntType<1> >::type IndexSet;
    typedef MetaPhysicL::SparseNumberVector<RawScalar, IndexSet> RawVector;
    #define VectorUnitVector SparseNumberVectorUnitVector
    #define VectorOf SparseNumberVectorOf
  #endif // USE_STRUCT, USE_DYNAMIC
#else // USE_SPARSE
  #ifdef USE_SHADOW
    #include "metaphysicl/dualshadowvector.h"
    typedef MetaPhysicL::ShadowNumber<double, long double> RawScalar;
  #else // USE_SHADOW
    #include "metaphysicl/dualnumbervector.h"
    typedef double RawScalar;
  #endif // USE_SHADOW
  typedef MetaPhysicL::NumberVector<NDIM, RawScalar> RawVector;
  #define VectorUnitVector NumberVectorUnitVector
  #define VectorOf NumberVectorOf
#endif // USE_SPARSE

#endif // __testing_h__
