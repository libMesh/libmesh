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

#ifndef _SETLOCATION_HPP_
#define _SETLOCATION_HPP_


#include <Hilbert/Common.hpp>
#include <Hilbert/BigBitVec.hpp>


namespace Hilbert
{
	template <class P,class I>
	H_INLINE
	void
	setLocation(
		P *p,
		int n,
		int i,
		const I &l
		)
	{
		// Easy to understand implementation
		/*int j;
		for ( j = 0; j < n; j++ )
			p[j].setBit(i,l.getBit(j));
		return;*/

		// Much faster loop-unrolled implementation.
		int ir, ib;
		FBV_UINT im;
		BBV_MODSPLIT(ir,ib,i);
		im = (FBV1<<ib);

#define SETBIT p->racks()[ir]|=im; if ((*lr&jm)==0) p->racks()[ir]^=im; jm<<=1; ++p;
#define SETBIT1(a) \
		case (a+1): SETBIT;
#define SETBIT2(a) \
		SETBIT1(a+1); \
		SETBIT1(a);
#define SETBIT4(a) \
		SETBIT2(a+2); \
		SETBIT2(a);
#define SETBIT8(a) \
		SETBIT4(a+4); \
		SETBIT4(a);
#define SETBIT16(a) \
		SETBIT8(a+8); \
		SETBIT8(a);
#define SETBIT32(a) \
		SETBIT16(a+16); \
		SETBIT16(a);

		int j = 0;
		FBV_UINT jm = 1;
		const FBV_UINT *lr = l.racks();
		while ( j < n )
		{
			if ( jm == 0 )
			{
				jm = 1;
				++lr;
			}
			
			int dj = n - j;
			switch ( n - j )
			{
				default: dj = 32;
				SETBIT32(0);
			}
			j += dj;
		}
		
		return;
	}

}


#endif
