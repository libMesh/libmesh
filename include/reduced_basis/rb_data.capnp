@0xffc2b6c54145b6c4;
# Schema's unique identifier

# Import C++-specific functionalities (e.g. namespace generation)
using Cxx = import "/capnp/c++.capnp";

# Declare a namespace for the symbols defined in the schema
$Cxx.namespace("RBData");

using Integer = UInt32;
using Real = Float64;

struct Complex @0xda65bfa81ea2ce0b {
  real @0 :Real;
  imag @1 :Real;
}

using Number = <NUMBER>;

struct RBParameter @0xb73e071e0a405648 {
  name  @0 :Text;
  value @1 :Real;
}

struct ParameterRanges @0xeccae95c5c74616e {
  name     @0 :List(Text);
  minValue @1 :List(Real);
  maxValue @2 :List(Real);
}

struct DiscreteParameterList @0xc7c5f4dfa33dfb27 {
  name   @0 :List(Text);
  values @1 :List(List(Real));
}

struct Output @0xeb6a19346114773a {
  outputDualInnerprods @0 :List(Number);
  outputVectors        @1 :List(List(Number));
}

struct RBEvaluation @0xa459b0816a4ad3e3 {
  nBfs                 @0 :Integer;
  parameterRanges      @1 :ParameterRanges;
  discreteParameters   @2 :DiscreteParameterList;
  aqAqInnerprods       @3 :List(Number);
  fqAqInnerprods       @4 :List(Number);
  fqInnerprods         @5 :List(Number);
  rbFqVectors          @6 :List(List(Number));
  rbAqMatrices         @7 :List(List(Number));
  outputs              @8 :List(Output);
  rbInnerProductMatrix @9 :List(Number);
}

struct TransientRBEvaluation @0xf0ee71757fa42963 {
  rbEvaluation      @0  :RBEvaluation;
  deltaT            @1  :Real;
  eulerTheta        @2  :Real;
  nTimeSteps        @3  :Integer;
  timeStep          @4  :Integer;
  rbL2Matrix        @5  :List(Number);
  rbMqMatrices      @6  :List(List(Number));
  initialL2Errors   @7  :List(Real);
  initialConditions @8  :List(List(Number));
  fqMqInnerprods    @9  :List(Number);
  mqMqInnerprods    @10 :List(Number);
  aqMqInnerprods    @11 :List(Number);
}

struct Point3D @0xde7dfcf8ecccaa90 {
  x @0 :Real;
  y @1 :Real;
  z @2 :Real;
}

struct MeshElem @0xcd01b7bd6045605d {
  type         @0 :Text;
  points       @1 :List(Point3D);
  subdomainId  @2 :Integer;
}

struct RBEIMEvaluation @0xf8121d2237427a80 {
  rbEvaluation                @0 :RBEvaluation;
  interpolationElemIds        @1 :List(Integer);
  interpolationMatrix         @2 :List(Number);
  interpolationPointsElems    @3 :List(MeshElem);
  interpolationPointsVar      @4 :List(Integer);
  interpolationPoints         @5 :List(Point3D);
  extraInterpolationMatrixRow @6 :List(Number);
  extraInterpolationPointElem @7 :MeshElem;
  extraInterpolationPointVar  @8 :Integer;
  extraInterpolationPoint     @9 :Point3D;
}

struct RBSCMEvaluation @0xb8dd038628a64b16 {
  parameterRanges    @0 :ParameterRanges;
  discreteParameters @1 :DiscreteParameterList;
  bMin               @2 :List(Real);
  bMax               @3 :List(Real);
  cJStabilityVector  @4 :List(Real);
  cJ                 @5 :List(List(RBParameter));
  scmUbVectors       @6 :List(Real);
}
