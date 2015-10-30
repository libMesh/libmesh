#include <hilbert.h>

// user-provided functions for performing MPI reductions on 
// HilbertIndices datatypes.  Note that these can actually 
// be implemented without any knowledge of MPI.
void __hilbert_max_op (Hilbert::HilbertIndices *in, Hilbert::HilbertIndices *inout, int *len, void *)
{
  Hilbert::BitVecType 
    a(3*sizeof(double)*sizeof(Hilbert::inttype)),
    b(3*sizeof(double)*sizeof(Hilbert::inttype));
  
  assert (a.rackCount() == 3);
  assert (b.rackCount() == 3);
  
  for (int i=0; i<*len; i++, in++, inout++)
    {
      a = *in; /**/ b = *inout;
      
      // if (a < b), then inout already contains max(a,b)
      if (b < a) *inout = *in;
    }
}



void __hilbert_min_op (Hilbert::HilbertIndices *in, Hilbert::HilbertIndices *inout, int *len, void *)
{
  Hilbert::BitVecType 
    a(3*sizeof(double)*sizeof(Hilbert::inttype)),
    b(3*sizeof(double)*sizeof(Hilbert::inttype));
  
  assert (a.rackCount() == 3);
  assert (b.rackCount() == 3);
  
  for (int i=0; i<*len; i++, in++, inout++)
    {
      a = *in; /**/ b = *inout;
      
      // if (b < a), then inout already contains min(a,b)
      if (a < b)  *inout = *in;
    }
}
