// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// Local includes
#include "libmesh/quadrature_nodal.h"
#include "libmesh/quadrature_trap.h"
#include "libmesh/quadrature_simpson.h"
#include "libmesh/enum_to_string.h"

namespace libMesh
{

void QNodal::init_3D(const ElemType, unsigned int)
{
#if LIBMESH_DIM == 3

  switch (_type)
    {
    case TET4:
    case PRISM6:
    case HEX8:
    case PYRAMID5:
      {
        // Construct QTrap rule that matches our own nodal pyramid quadrature permissions
        QTrap rule(/*dim=*/3, /*ignored*/_order);
        rule.allow_nodal_pyramid_quadrature = this->allow_nodal_pyramid_quadrature;
        rule.init(_type, /*ignored*/_p_level);
        _points.swap (rule.get_points());
        _weights.swap(rule.get_weights());
        return;
      }

    case PRISM15:
      {
        // A rule with 15 points which is exact for linears, and
        // naturally produces a lumped approximation to the mass
        // matrix. The quadrature points are numbered the same way as
        // the reference element nodes.
        _points =
          {
            Point(0.,0.,-1), Point(+1,0.,-1), Point(0.,+1,-1),
            Point(0.,0.,+1), Point(+1,0.,+1), Point(0.,+1,+1),
            Point(.5,0.,-1), Point(.5,.5,-1), Point(0.,.5,-1),
            Point(0.,0.,0.), Point(+1,0.,0.), Point(0.,+1,0.),
            Point(.5,0.,+1), Point(.5,.5,+1), Point(0.,.5,+1),
          };

        // vertex (wv), tri edge (wt), and quad edge (wq) weights are
        // obtained using the same approach that was used for the Quad8,
        // see quadrature_nodal_2D.C for details.
        Real wv = Real(1) / 34;
        Real wt = Real(4) / 51;
        Real wq = Real(2) / 17;

        _weights = {wv, wv, wv, wv, wv, wv,
                    wt, wt, wt,
                    wq, wq, wq,
                    wt, wt, wt};

        return;
      }

    case HEX20:
      {
        // A rule with 20 points which is exact for linears, and
        // naturally produces a lumped approximation to the mass
        // matrix. The quadrature points are numbered the same way as
        // the reference element nodes.
        _points =
          {
            Point(-1,-1,-1), Point(+1,-1,-1), Point(+1,+1,-1), Point(-1,+1,-1),
            Point(-1,-1,+1), Point(+1,-1,+1), Point(+1,+1,+1), Point(-1,+1,+1),
            Point(0.,-1,-1), Point(+1,0.,-1), Point(0.,+1,-1), Point(-1,0.,-1),
            Point(-1,-1,0.), Point(+1,-1,0.), Point(+1,+1,0.), Point(-1,+1,0.),
            Point(0.,-1,+1), Point(+1,0.,+1), Point(0.,+1,+1), Point(-1,0.,+1)
          };

        // vertex (wv), and edge (we) weights are obtained using the
        // same approach that was used for the Quad8, see
        // quadrature_nodal_2D.C for details.
        Real wv = Real(7) / 31;
        Real we = Real(16) / 31;

        _weights = {wv, wv, wv, wv, wv, wv, wv, wv,
                    we, we, we, we, we, we, we, we, we, we, we, we};

        return;
      }

    case TET10:
    case PRISM18:
    case HEX27:
    case PYRAMID13:
    case PYRAMID14:
      {
        // Construct QSimpson rule that matches our own nodal pyramid quadrature permissions
        QSimpson rule(/*dim=*/3, /*ignored*/_order);
        rule.allow_nodal_pyramid_quadrature = this->allow_nodal_pyramid_quadrature;
        rule.init(_type, /*ignored*/_p_level);
        _points.swap (rule.get_points());
        _weights.swap(rule.get_weights());

        // We can't do a proper Simpson rule for pyramids regardless
        if (_type == PYRAMID13)
          {
            _points.resize(13);
            _weights.resize(13);
          }

        return;
      }

    case PRISM20:
      {
        _points =
          {
            Point(0.,0.,-1), Point(+1,0.,-1), Point(0.,+1,-1),
            Point(0.,0.,+1), Point(+1,0.,+1), Point(0.,+1,+1),
            Point(.5,0.,-1), Point(.5,.5,-1), Point(0.,.5,-1),
            Point(0.,0.,0.), Point(+1,0.,0.), Point(0.,+1,0.),
            Point(.5,0.,+1), Point(.5,.5,+1), Point(0.,.5,+1),
            Point(.5,0.,0.), Point(.5,.5,0.), Point(0.,.5,0.),
            Point(1/Real(3),1/Real(3),-1), Point(1/Real(3),1/Real(3),+1)
          };

        // Symmetry gives us identical weights on vertices, triangle
        // edges, square edges, triangle faces, square faces; then we
        // have the midnode.  Solving for weights which exactly
        // integrate cubics and xi^2*zeta^2 would give a unique answer
        // ... with wse=0 and wte<0.  See Quad8 in
        // quadrature_nodal_2D.C for discussion of this problem.
        //
        // Dropping the xi^2 zeta^2 constraint gives us a rank-1 null
        // space to play with ... but adding anything from that null
        // space to push wte closer to positive just pushes wse
        // negative.
        //
        // Dropping exact integration of xi^3 type terms gives us a
        // rank-2 null space.  I just found the max(min(weight)) from
        // that solution space using Octave.  Someone who cares more
        // than me might want to repeat this exercise with better than
        // double precision...

        constexpr Real wv = 8.546754839782711e-03;
        constexpr Real wte = 8.548403936750599e-03;
        constexpr Real wse = wv; // here's our min(weight) constraint...
        constexpr Real wtf = 1.153811903370667e-01;
        constexpr Real wsf = 2.136754673824396e-01;

        _weights = {wv, wv, wv, wv, wv, wv,
                    wte, wte, wte, wse, wse, wse, wte, wte, wte,
                    wsf, wsf, wsf, wtf, wtf};

        return;
      }

    case PRISM21:
      {
        _points =
          {
            Point(0.,0.,-1), Point(+1,0.,-1), Point(0.,+1,-1),
            Point(0.,0.,+1), Point(+1,0.,+1), Point(0.,+1,+1),
            Point(.5,0.,-1), Point(.5,.5,-1), Point(0.,.5,-1),
            Point(0.,0.,0.), Point(+1,0.,0.), Point(0.,+1,0.),
            Point(.5,0.,+1), Point(.5,.5,+1), Point(0.,.5,+1),
            Point(.5,0.,0.), Point(.5,.5,0.), Point(0.,.5,0.),
            Point(1/Real(3),1/Real(3),-1), Point(1/Real(3),1/Real(3),+1),
            Point(1/Real(3),1/Real(3),0.)
          };

        // Symmetry gives us identical weights on vertices, triangle
        // edges, square edges, triangle faces, square faces; then we
        // have the midnode.  Solving for weights which exactly
        // integrate cubics, xi^2*zeta^2, and xi^3*zeta^2, gives a
        // unique answer:
        constexpr Real wv = Real(1)/120;
        constexpr Real wte = Real(1)/45;
        constexpr Real wse = Real(1)/30;
        constexpr Real wtf = Real(3)/40;
        constexpr Real wsf = Real(4)/45;

        _weights = {wv, wv, wv, wv, wv, wv,
                    wte, wte, wte, wse, wse, wse, wte, wte, wte,
                    wsf, wsf, wsf, wtf, wtf, Real(3)/10};

        return;
      }

    case PYRAMID18:
      {
        libmesh_error_msg_if(!allow_nodal_pyramid_quadrature,
                             "Nodal quadrature on Pyramid elements is not allowed by default since\n"
                             "the Jacobian of the inverse element map is not well-defined at the Pyramid apex.\n"
                             "Set the QBase::allow_nodal_pyramid_quadrature flag to true to ignore skip this check.");

        _points =
          {
            Point(-1,-1,0.), Point(+1,-1,0.), Point(+1,+1,0.),
            Point(-1,+1,0.), Point(0.,0.,+1), Point(0.,-1,0.),
            Point(+1,0.,0.), Point(0.,+1,0.), Point(-1,0.,0.),
            Point(-.5,-.5,.5), Point(-.5,.5,.5), Point(.5,.5,.5),
            Point(-.5,.5,.5), Point(0.,0.,0.), Point(0.,-2/Real(3),1/Real(3)),
            Point(2/Real(3),0.,1/Real(3)), Point(0.,2/Real(3),1/Real(3)),
            Point(-2/Real(3),0.,1/Real(3))
          };

        // Even with triangle faces to play with, I can't seem to get
        // exact integrals of any higher order functions without
        // negative weights.  So I punt and just use QTrap weights.
        _weights = {1/Real(4), 1/Real(4), 1/Real(4), 1/Real(4), 1/Real(3),
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

        return;
      }

    case TET14:
      {
        _points.resize(14);
        _weights.resize(14);

        _points[0](0) = 0.;   _points[5](0) = .5;
        _points[0](1) = 0.;   _points[5](1) = .5;
        _points[0](2) = 0.;   _points[5](2) = 0.;

        _points[1](0) = 1.;   _points[6](0) = 0.;
        _points[1](1) = 0.;   _points[6](1) = .5;
        _points[1](2) = 0.;   _points[6](2) = 0.;

        _points[2](0) = 0.;   _points[7](0) = 0.;
        _points[2](1) = 1.;   _points[7](1) = 0.;
        _points[2](2) = 0.;   _points[7](2) = .5;

        _points[3](0) = 0.;   _points[8](0) = .5;
        _points[3](1) = 0.;   _points[8](1) = 0.;
        _points[3](2) = 1.;   _points[8](2) = .5;

        _points[4](0) = .5;   _points[9](0) = 0.;
        _points[4](1) = 0.;   _points[9](1) = .5;
        _points[4](2) = 0.;   _points[9](2) = .5;


        _points[10](0) = 1/Real(3); _points[11](0) = 1/Real(3);
        _points[10](1) = 1/Real(3); _points[11](1) = 0.;
        _points[10](2) = 0.;        _points[11](2) = 1/Real(3);

        _points[12](0) = 1/Real(3); _points[13](0) = 0.;
        _points[12](1) = 1/Real(3); _points[13](1) = 1/Real(3);
        _points[12](2) = 1/Real(3); _points[13](2) = 1/Real(3);

        // RHS:
        //
        // The ``optimal'' nodal quadrature here, the one that
        // integrates every Lagrange polynomial on these nodes
        // exactly, produces an indefinite mass matrix ...
        //
        // const Real wv = 1/Real(240);
        // const Real we = 0;
        // const Real wf = 3/Real(80);

        // We could average that with our Tet10 rule and get:
        //
        // const Real wv = (1/Real(240)+1/Real(192))/2;
        // const Real we = Real(14)/576/2;
        // const Real wf = 3/Real(80)/2;
        //
        // Except our Tet10 rule won't actually exactly integrate
        // quadratics! (exactly integrating quadratics wouldn't even
        // have given a positive semidefinite mass matrix there...)
        //
        // John derived equations for wv and we based on symmetry and
        // the requirement to exactly integrate quadratics; within
        // those constraints we might pick the wf that maximizes the
        // minimum nodal weight:
        // const Real wf = Real(15)/440;
        // const Real wv = -1/Real(120) + wf/3;
        // const Real we = 1/Real(30) - wf*(Real(8)/9);
        //
        // But John also did the same clever optimization trick that
        // quadrature_nodal_2D.C discusses in the context of Quad8 and
        // Hex20 outputs, and for Tet14 that gives us:
        const Real wv = Real(87)/47120;
        const Real we = Real(164)/26505;
        const Real wf = Real(1439)/47120;

        _weights = {wv, wv, wv, wv, we, we, we, we, we, we, wf, wf, wf, wf};
        return;
      }

    default:
      libmesh_error_msg("ERROR: Unsupported type: " << Utility::enum_to_string(_type));
    }
#endif
}

} // namespace libMesh
