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

#ifndef _SETBITS_HPP_
#define _SETBITS_HPP_


#include <Hilbert/Common.hpp>
#include <Hilbert/BigBitVec.hpp>


namespace Hilbert
{

	template<class H,class I>
	H_INLINE
	void
	setBits(
		H &h,				// destination
		int n,			// number of bits
		int i,			// bit position
		const I &w	// bits to place
		)
	{
		// This is terribly inefficient.
		int j;
		for ( j = 0; j < n; j++ )
			h.setBit(i+j,w.getBit(j));
	}


	// <CBigBitVec,CBigBitVec>
	// #TODO


	// <CBigBitVec,CFixBitVec>
	template<>
	H_INLINE
	void
	setBits(
		CBigBitVec &h,
		int n,
		int i,
		const CFixBitVec &w
		)
	{
		int ir, ib, t;
		BBV_MODSPLIT(ir,ib,i);
		h.racks()[ir] |= w.rack() << ib;
		t = ib+n;
		if ( t > FBV_BITS )
		{
			t -= FBV_BITS;
			h.racks()[ir+1] |= w.rack() >> t;
		}
		return;
	}
	
	
	// <CFixBitVec,CFixBitVec>
	template<>
	H_INLINE
	void
	setBits(
		CFixBitVec &h,
		int /* n */,
		int i,
		const CFixBitVec &w
		)
	{
		h |= w << i;
		return;
	}

}


#endif
