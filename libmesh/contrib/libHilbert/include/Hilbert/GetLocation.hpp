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

#ifndef _GETLOCATION_HPP_
#define _GETLOCATION_HPP_


#include <Hilbert/Common.hpp>
#include <Hilbert/BigBitVec.hpp>


namespace Hilbert
{
	
	template<class P,class I>
	H_INLINE
	void
	_getLocation(
		const P *p,
		int jo,
		int jn,
		int ir,
		FBV_UINT im,
		FBV_UINT &l
		)
	{
		l = 0;
    switch ( jn )
    {
#define GETLOC_CASE(i)		case ((i)+1): if (p[jo+(i)].racks()[ir]&im) l|=(FBV1<<(i))
#define GETLOC_CASE2(i) \
			GETLOC_CASE((i)+1); \
			GETLOC_CASE(i)
#define GETLOC_CASE4(i) \
			GETLOC_CASE2((i)+2); \
			GETLOC_CASE2(i)
#define GETLOC_CASE8(i) \
			GETLOC_CASE4((i)+4); \
			GETLOC_CASE4(i)
#define GETLOC_CASE16(i) \
			GETLOC_CASE8((i)+8); \
			GETLOC_CASE8(i)
#define GETLOC_CASE32(i) \
			GETLOC_CASE16((i)+16); \
			GETLOC_CASE16(i)
#if FBV_BITS == 64
			GETLOC_CASE32(32);
#endif
			GETLOC_CASE32(0);
		}
		return;
	}
	
	
	template<class P,class I>
	H_INLINE
	void
	getLocation(
		const P *p,
		int n,
		int i,
		I	&l
		)
	{
		/*int j;
		for ( j = n-1; j >= 0; --j )
			l.setBit(j,p[j].getBit(i));
		return;*/
		
		int j, jo, ir;
		FBV_UINT im;
		
		if ( P::type() == eBig )
		{
			ir = i / FBV_BITS;
			im = FBV1 << (i-ir*FBV_BITS);
		}
		else
		{
			ir = 0;
			im = FBV1 << i;
		}
		
		j = 0;
		jo = 0;
		if ( I::type() == eBig )
		{
			for ( ; j < l.rackCount()-1; ++j, jo += FBV_BITS )
				_getLocation<P,I>(p,jo,FBV_BITS,ir,im,l.racks()[j]);
		}
		_getLocation<P,I>(p,jo,n-jo,ir,im,l.racks()[j]);
		
		return;
	}
	
}


#endif
