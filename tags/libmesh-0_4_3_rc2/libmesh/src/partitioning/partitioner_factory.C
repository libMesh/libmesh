// $Id: partitioner_factory.C,v 1.7 2004-04-25 05:43:33 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
  
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// C++ Includes   -----------------------------------

// // Local Includes -----------------------------------
// #include "libmesh_config.h"
// #include "partitioner.h"
// #include "centroid_partitioner.h"
// #include "metis_partitioner.h"
// #include "parmetis_partitioner.h"
// #include "linear_partitioner.h"
// #include "sfc_partitioner.h"
// #include "hilbert_sfc_partitioner.h"
// #include "morton_sfc_partitioner.h"
// #include "factory.h"



// // ------------------------------------------------------------
// // Explicit instantiation of a Partitioner factory
// // (this is becoming complicated due to compiler sematics)
// #if  defined(__IBMCPP__) || defined(__HP_aCC)

//   template class Factory<Partitioner;
//   template <class Partitioner> std::map<std::string, FactoryBase* > Factory<Partitioner>::factory_map;

//   // GCC - different syntax possibly required 
// #elif defined(__GNUC__) && defined (__GNUC_MINOR__)
// # if (__GNUC__ >= 3) && (__GNUC_MINOR__ >= 4)

//   template <class Partitioner>
//   std::map<std::string, FactoryBase* > Factory<Partitioner>::factory_map;

// # else
        
//   std::map<std::string, FactoryBase* > Factory<Partitioner>::factory_map;

// # endif

//   // All others
// #else
      
//   std::map<std::string, FactoryBase* > Factory<Partitioner>::factory_map;
  
// #endif



// // ------------------------------------------------------------
// // Register Partitioning classes with the factory.  These will never
// // be called from user code, they just need to get instantiated.  Hide
// // them in an anonymous namespace to prevent name clashes
// namespace {
  
// #ifdef HAVE_METIS
//   FactoryImp<MetisPartitioner,      Partitioner> metis    ("Metis");
// #endif
  
// #ifdef HAVE_PARMETIS
//   FactoryImp<ParmetisPartitioner,   Partitioner> parmetis ("Parmetis");
// #endif

// #ifdef HAVE_SFCURVES
//   FactoryImp<SFCPartitioner,        Partitioner> sfc      ("SFCurves");
//   FactoryImp<HilbertSFCPartitioner, Partitioner> hilbert  ("Hilbert");
//   FactoryImp<MortonSFCPartitioner,  Partitioner> morton   ("Morton");
// #endif
  
//   FactoryImp<LinearPartitioner,     Partitioner> linear   ("Linear");
//   FactoryImp<CentroidPartitioner,   Partitioner> centroid ("Centroid");
  
// }
