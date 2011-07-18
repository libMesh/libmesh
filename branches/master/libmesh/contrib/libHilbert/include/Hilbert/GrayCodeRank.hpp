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

#ifndef _GRAYCODERANK_HPP_
#define _GRAYCODERANK_HPP_


#include "Hilbert/Common.hpp"
#include "Hilbert/BigBitVec.hpp"


namespace Hilbert
{
	// This is the bulk of the cost in calculating
	// a Compact Hilbert Index.  It compresses a previously
	// calculated index when provided the rotation
	// at each level of precision.
	template<class H,class HC>
	H_INLINE
	void
	compactIndex(
		const int *ms,
		const int *ds,
		int n,
		int m,
		H &h,
		HC &hc
		)
	{
		int i, j, hr, hcr;
		FBV_UINT hm, hcm;

		hc.zero();

		hr = hcr = 0;
		hm = hcm = 1;

		// Run through the levels of precision
		for ( i = 0; i < m; i++ )
		{
			// #D
			// int k = ds[i] + 1;
			// if ( k == n ) k = 0;
			// j = k;

/*#define COMP	{ \
	if(ms[j]>i) { \
		if(h.racks()[hr]&hm)hc.racks()[hcr]|=hcm; \
		hcm<<=1;if(hcm==0){hcm=1;++hcr;} \
		} \
	if((++j)==n)j=0; \
	hm<<=1;if(hm==0){hm=1;++hr;} \
	}
#define COMP1(a)	case a: COMP;
#define COMP2(a)	COMP1(a+1); \
		COMP1(a);
#define COMP4(a)	COMP2(a+2); \
		COMP2(a);
#define COMP8(a)	COMP4(a+4); \
		COMP4(a);
#define COMP16(a)	COMP8(a+8); \
		COMP8(a);
#define COMP32(a)	COMP16(a+16); \
		COMP16(a);

			// This complicated mess only buys a marginal performance increase.
			int k, kd;
			k = n;
			j = ds[i];
			while ( k )
			{
				kd = k;
				switch ( kd )
				{
					default: kd = 32;
					COMP32(1);
				}
				k -= kd;
			}*/
			
			// Run through the dimensions
			j = ds[i];
			do
			{
				// This dimension contributes a bit?
				if ( ms[j] > i )
				{
					if ( h.racks()[hr] & hm )
						hc.racks()[hcr] |= hcm;
					hcm<<=1; if (hcm==0) {hcm=1;++hcr;}
				}

				if ( ++j == n ) j = 0;
				hm<<=1; if (hm==0) {hm=1;++hr;}
			}
			while ( j != ds[i] );
		}

		return;
	}

	template<class I>
	H_INLINE
	void
	grayCodeRank(
		const I &mask,
		const I &gi,
		int n,
		I &r
		)
	{
		r.zero();
		int i, ir, jr;
		FBV_UINT im, jm;

		jr = 0; jm = 1;
		ir = 0; im = 1;
		for ( i = 0; i < n; ++i )
		{
			if ( mask.racks()[ir] & im )
			{
				if ( gi.racks()[ir] & im )
					r.racks()[jr] |= jm;
				jm<<=1; if ( jm == 0 ) { jm = 1; ++jr; }
			}
			
			im<<=1; if ( im == 0 ) { im = 1; ++ir; }
		}

		return;
	}


	template<class I>
	H_INLINE
	void
	grayCodeRankInv(
		const I &mask,
		const I &ptrn,
		const I &r,
		int n,
		int b,
		I &g,
		I &gi
		)
	{
		g.zero();
		gi.zero();

		int i, j, ir, jr;
		FBV_UINT im, jm;

		i = n-1;
		BBV_MODSPLIT(ir,im,i);
		im = (FBV1<<im);

		j = b-1;
		BBV_MODSPLIT(jr,jm,j);
		jm = (FBV1<<jm);

		FBV_UINT gi0, gi1, g0;
		gi1 = gi0 = g0 = 0;

		for ( ; i >= 0; --i )
		{
			// Unconstrained bit?
			if ( mask.racks()[ir] & im )
			{
				gi1 = gi0;
				gi0 = (r.racks()[jr] & jm) > 0;
				g0 = gi0 ^ gi1;
				if ( gi0 ) gi.racks()[ir] |= im;
				if ( g0 ) g.racks()[ir] |= im;
				jm>>=1; if ( jm == 0 ) { jm=((FBV_UINT)1)<<(FBV_BITS-1); --jr; }
			}
			else
			{
				g0 = (ptrn.racks()[ir] & im) > 0;
				gi1 = gi0;
				gi0 = g0 ^ gi1;
				if ( gi0 ) gi.racks()[ir] |= im;
				if ( g0 ) g.racks()[ir] |= im;
			}

			im>>=1; if ( im == 0 ) { im=((FBV_UINT)1)<<(FBV_BITS-1); --ir; }
		}

		return;
	}

	
	template <class I>
	H_INLINE
	void
	extractMask(
		const int *ms,
		int n,
		int d,
		int i,
		I &mask,
		int &b
		)
	{
		assert( 0 <= d && d < n );
		int j, jr;
		FBV_UINT jm;

		mask.zero();
		b = 0;
		
		jm = 1; jr = 0;
		j = d; // #D j = (d==n-1) ? 0 : d+1;
		do
		{
			if ( ms[j] > i )
			{
				mask.racks()[jr] |= jm;
				++b;
			}
		
			jm <<= 1; if ( jm == 0 ) { jm=1; ++jr; }
			if ( ++j == n ) j = 0;
		}
		while ( j != d );

		return;
	}
}


#endif
