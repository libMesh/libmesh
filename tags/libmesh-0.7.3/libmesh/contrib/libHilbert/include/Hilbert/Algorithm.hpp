/*
 * Copyright (C) 2006-2007 Chris Hamilton <chamilton@cs.dal.ca>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef _ALGORITHM_HPP_
#define _ALGORITHM_HPP_


#include "Hilbert/Common.hpp"
#include "Hilbert/BigBitVec.hpp"
#include "Hilbert/GetLocation.hpp"
#include "Hilbert/SetLocation.hpp"
#include "Hilbert/GetBits.hpp"
#include "Hilbert/SetBits.hpp"
#include "Hilbert/GrayCodeRank.hpp"
#include <string.h>


// Templated Hilbert functions.
// P - is the class used to represent each dimension
//     of a multidimensional point.
// H - is the class used to represent the point as a Hilbert
//     index.
// I - is the class used to represent interim variables.
//     Needs to have as many bits as there are dimensions.
//
// In general, each of P,H,I should be a FixBitVec if they
// fit in that fixed precision.  Otherwise, use a BigBitVec.
// Whatever you use, they must all be of consistent underlying
// storage types.


// The dimension across which the Hilbert curve travels
// principally.
// D0 is d + 1, where d is the actual dimension you want to
// walk across.
// MUST HAVE 0 <= D0 < n
#define D0	1


namespace Hilbert
{
	// 'Transforms' a point.
	template<class I>
	H_INLINE
	void
	transform(
		const I &e,
		int d,
		int n,
		I &a
		)
	{
		a ^= e;
		a.rotr( d, n );//#D d+1, n );
		return;
	}

	// Inverse 'transforms' a point.
	template<class I>
	H_INLINE
	void
	transformInv(
		const I &e,
		int d,
		int n,
		I &a
		)
	{
		a.rotl( d, n );//#D d+1, n );
		a ^= e;
		return;
	}

	// Update for method 1 (GrayCodeInv in the loop)
	template<class I>
	H_INLINE
	void
	update1(
		const I &l,
		const I &t,
		const I &w,
		int n,
		I &e,
		int &d
		)
	{
		assert( 0 <= d && d < n );
		e = l;
		e.toggleBit( d ); //#D d == n-1 ? 0 : d+1 );
	
		// Update direction
		d += 1 + t.fsb();
		if ( d >= n ) d -= n;
		if ( d >= n ) d -= n;
		assert( 0 <= d && d < n );

		if ( ! (w.rack() & 1) )
			e.toggleBit( d == 0 ? n-1 : d-1 ); //#D d );

		return;
	}

	// Update for method 2 (GrayCodeInv out of loop)
	template<class I>
	H_INLINE
	void
	update2(
		const I &l,
		const I &t,
		const I & /* w */,
		int n,
		I &e,
		int &d
		)
	{
		assert( 0 <= d && d < n );
		e = l;
		e.toggleBit( d );//#D d == n-1 ? 0 : d+1 );
	
		// Update direction
		d += 1 + t.fsb();
		if ( d >= n ) d -= n;
		if ( d >= n ) d -= n;
		assert( 0 <= d && d < n );

		return;
	}

	template <class P,class H,class I>
	H_INLINE
	void
	_coordsToIndex(
		const P *p,
		int m,
		int n,
		H &h,
		int *ds = NULL // #HACK
		)
	{
		I e(n), l(n), t(n), w(n);
		int d, i;
		int ho = m*n;

		// Initialize
		e.zero();
		d = D0;
		l.zero();
		h.zero();

		int r, b;
		BBV_MODSPLIT(r,b,n-1);
		FBV_UINT bm = (1<<b);

		// Work from MSB to LSB
		for ( i = m-1; i >= 0; i-- )
		{
			// #HACK
			if ( ds ) ds[i] = d;
			
			// Get corner of sub-hypercube where point lies.
			getLocation<P,I>(p,n,i,l);

			// Mirror and reflect the location.
			// t = T_{(e,d)}(l)
			t = l;
			transform<I>(e,d,n,t);

			w = t;
			if ( i < m-1 )
				w.racks()[r] ^= bm;

			// Concatenate to the index.
			ho -= n;
			setBits<H,I>(h,n,ho,w);

			// Update the entry point and direction.
			update2<I>(l,t,w,n,e,d);
		}

		h.grayCodeInv();

		return;
	}

	// This is wrapper to the basic Hilbert curve index
	// calculation function.  It will support fixed or
	// arbitrary precision, templated.  Depending on the
	// number of dimensions, it will use the most efficient
	// representation for interim variables.
	// Assumes h is big enough for the output (n*m bits!)
	template<class P,class H>
	H_INLINE
	void
	coordsToIndex(
		const P *p,	// [in ] point
		int m,		// [in ] precision of each dimension in bits
		int n,		// [in ] number of dimensions
		H &h		// [out] Hilbert index
		)
	{
		// Intermediate variables will fit in fixed width?
		if ( n <= FBV_BITS )
			_coordsToIndex<P,H,CFixBitVec>(p,m,n,h);
		// Otherwise, they must be BigBitVecs.
		else
			_coordsToIndex<P,H,CBigBitVec>(p,m,n,h);

		return;
	}


	template <class P,class H,class I>
	H_INLINE
	void
	_indexToCoords(
		P *p,
		int m,
		int n,
		const H &h
		)
	{
		I e(n), l(n), t(n), w(n);
		int d, i, j, ho;

		// Initialize
		e.zero();
		d = D0;
		l.zero();
		for ( j = 0; j < n; j++ )
			p[j].zero();

		ho = m*n;
		
		// Work from MSB to LSB
		for ( i = m-1; i >= 0; i-- )
		{
			// Get the Hilbert index bits
			ho -= n;
			getBits<H,I>(h,n,ho,w);

			// t = GrayCode(w)
			t = w;
			t.grayCode();

			// Reverse the transform
			// l = T^{-1}_{(e,d)}(t)
			l = t;
			transformInv<I>(e,d,n,l);

			// Distribute these bits
			// to the coordinates.
			setLocation<P,I>(p,n,i,l);

			// Update the entry point and direction.
			update1<I>(l,t,w,n,e,d);
		}

		return;
	}

	// This is wrapper to the basic Hilbert curve inverse
	// index function.  It will support fixed or
	// arbitrary precision, templated.  Depending on the
	// number of dimensions, it will use the most efficient
	// representation for interim variables.
	// Assumes each entry of p is big enough to hold the
	// appropriate variable.
	template<class P,class H>
	H_INLINE
	void
	indexToCoords(
		P *p,				// [out] point
		int m,			// [in ] precision of each dimension in bits
		int n,			// [in ] number of dimensions
		const H &h	// [out] Hilbert index
		)
	{
		// Intermediate variables will fit in fixed width?
		if ( n <= FBV_BITS )
			_indexToCoords<P,H,CFixBitVec>(p,m,n,h);
		// Otherwise, they must be BigBitVecs.
		else
			_indexToCoords<P,H,CBigBitVec>(p,m,n,h);

		return;
	}

	template <class P,class HC,class I>
	H_INLINE
	void
	_coordsToCompactIndex(
		const P *p,
		const int *ms,
		int n,
		HC &hc,
		int M = 0,
		int m = 0
		)
	{
		int i, mn;
		int *ds;
		
		// Get total precision and max precision
		// if not supplied
		if ( M == 0 || m == 0 )
		{
			M = m = 0;
			for ( i = 0; i < n; i++ )
			{
				if ( ms[i] > m ) m = ms[i];
				M += ms[i];
			}
		}

		mn = m*n;

		// If we could avoid allocation altogether (ie: have a
		// fixed buffer allocated on the stack) then this increases
		// speed by a bit (4% when n=4, m=20)
		ds = new int [ m ];

		if ( mn > FBV_BITS )
		{
			CBigBitVec h(mn);
			_coordsToIndex<P,CBigBitVec,I>(p,m,n,h,ds);
			compactIndex<CBigBitVec,HC>(ms,ds,n,m,h,hc);
		}
		else
		{
			CFixBitVec h;
			_coordsToIndex<P,CFixBitVec,I>(p,m,n,h,ds);
			compactIndex<CFixBitVec,HC>(ms,ds,n,m,h,hc);
		}

		delete [] ds;

		return;
	}

	// This is wrapper to the basic Hilbert curve index
	// calculation function.  It will support fixed or
	// arbitrary precision, templated.  Depending on the
	// number of dimensions, it will use the most efficient
	// representation for interim variables.
	// Assumes h is big enough for the output (n*m bits!)
	template<class P,class HC>
	H_INLINE
	void
	coordsToCompactIndex(
		const P *p,		// [in ] point
		const int *ms,// [in ] precision of each dimension in bits
		int n,				// [in ] number of dimensions
		HC &hc,					// [out] Hilbert index
		int M = 0,
		int m = 0
		)
	{
		// Intermediate variables will fit in fixed width?
		if ( n <= FBV_BITS )
			_coordsToCompactIndex<P,HC,CFixBitVec>(p,ms,n,hc,M,m);
		// Otherwise, they must be BigBitVecs.
		else
			_coordsToCompactIndex<P,HC,CBigBitVec>(p,ms,n,hc,M,m);

		return;
	}

	template <class P,class HC,class I>
	H_INLINE
	void
	_compactIndexToCoords(
		P *p,
		const int *ms,
		int n,
		const HC &hc,
		int M = 0,
		int m = 0
		)
	{
		I e(n), l(n), t(n), w(n), r(n), mask(n), ptrn(n);
		int d, i, j, b;

		// Get total precision and max precision
		// if not supplied
		if ( M == 0 || m == 0 )
		{
			M = m = 0;
			for ( i = 0; i < n; i++ )
			{
				if ( ms[i] > m ) m = ms[i];
				M += ms[i];
			}
		}
		
		// Initialize
		e.zero();
		d = D0;
		l.zero();
		for ( j = 0; j < n; j++ )
			p[j].zero();
		
		// Work from MSB to LSB
		for ( i = m-1; i >= 0; i-- )
		{
			// Get the mask and ptrn
			extractMask<I>(ms,n,d,i,mask,b);
			ptrn = e;
			ptrn.rotr(d,n);//#D ptrn.Rotr(d+1,n);

			// Get the Hilbert index bits
			M -= b;
			r.zero(); // GetBits doesn't do this
			getBits<HC,I>(hc,b,M,r);

			// w = GrayCodeRankInv(r)
			// t = GrayCode(w)
			grayCodeRankInv<I>(mask,ptrn,r,n,b,t,w);

			// Reverse the transform
			// l = T^{-1}_{(e,d)}(t)
			l = t;
			transformInv<I>(e,d,n,l);

			// Distribute these bits
			// to the coordinates.
			setLocation<P,I>(p,n,i,l);

			// Update the entry point and direction.
			update1<I>(l,t,w,n,e,d);
		}

		return;
	}

	// This is wrapper to the basic Hilbert curve inverse
	// index function.  It will support fixed or
	// arbitrary precision, templated.  Depending on the
	// number of dimensions, it will use the most efficient
	// representation for interim variables.
	// Assumes each entry of p is big enough to hold the
	// appropriate variable.
	template<class P,class HC>
	H_INLINE
	void
	compactIndexToCoords(
		P *p,					// [out] point
		const int *ms,// [in ] precision of each dimension in bits
		int n,				// [in ] number of dimensions
		const HC &hc,		// [out] Hilbert index
		int M = 0,
		int m = 0
		)
	{
		// Intermediate variables will fit in fixed width?
		if ( n <= FBV_BITS )
			_compactIndexToCoords<P,HC,CFixBitVec>(p,ms,n,hc,M,m);
		// Otherwise, they must be BigBitVecs.
		else
			_compactIndexToCoords<P,HC,CBigBitVec>(p,ms,n,hc,M,m);

		return;
	}

}


#endif
