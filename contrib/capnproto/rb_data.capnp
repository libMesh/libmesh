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

struct RBParameter @0xb73e071e0a405648 {
  name  @0 :Text;
  value @1 :Real;
}

struct ParameterRanges @0xeccae95c5c74616e {
  names     @0 :List(Text);
  minValues @1 :List(Real);
  maxValues @2 :List(Real);
}

struct DiscreteParameterList @0xc7c5f4dfa33dfb27 {
  names  @0 :List(Text);
  values @1 :List(List(Real));
}

struct Point3D @0xde7dfcf8ecccaa90 {
  x @0 :Real;
  y @1 :Real;
  z @2 :Real;
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

struct RBEvaluationReal @0xa459b0816a4ad3e3 {
  nBfs                 @0  :Integer;
  parameterRanges      @1  :ParameterRanges;
  discreteParameters   @2  :DiscreteParameterList;
  aqAqInnerprods       @3  :List(Real);
  fqAqInnerprods       @4  :List(Real);
  fqInnerprods         @5  :List(Real);
  rbFqVectors          @6  :List(List(Real));
  rbAqMatrices         @7  :List(List(Real));
  rbInnerProductMatrix @8  :List(Real);
  outputDualInnerprods @9  :List(List(Real));
  outputVectors        @10 :List(List(List(Real)));
}
struct RBEvaluationComplex @0x9a82baf67396076d {
  nBfs                 @0  :Integer;
  parameterRanges      @1  :ParameterRanges;
  discreteParameters   @2  :DiscreteParameterList;
  aqAqInnerprods       @3  :List(Complex);
  fqAqInnerprods       @4  :List(Complex);
  fqInnerprods         @5  :List(Complex);
  rbFqVectors          @6  :List(List(Complex));
  rbAqMatrices         @7  :List(List(Complex));
  rbInnerProductMatrix @8  :List(Complex);
  outputDualInnerprods @9  :List(List(Complex));
  outputVectors        @10 :List(List(List(Complex)));
}

struct TransientRBEvaluationReal @0xf0ee71757fa42963 {
  rbEvaluation      @0  :RBEvaluationReal;
  deltaT            @1  :Real;
  eulerTheta        @2  :Real;
  nTimeSteps        @3  :Integer;
  timeStep          @4  :Integer;
  rbL2Matrix        @5  :List(Real);
  rbMqMatrices      @6  :List(List(Real));
  initialL2Errors   @7  :List(Real);
  initialConditions @8  :List(List(Real));
  fqMqInnerprods    @9  :List(Real);
  mqMqInnerprods    @10 :List(Real);
  aqMqInnerprods    @11 :List(Real);
}
struct TransientRBEvaluationComplex @0x856e1656058c2e19 {
  rbEvaluation      @0  :RBEvaluationComplex;
  deltaT            @1  :Real;
  eulerTheta        @2  :Real;
  nTimeSteps        @3  :Integer;
  timeStep          @4  :Integer;
  rbL2Matrix        @5  :List(Complex);
  rbMqMatrices      @6  :List(List(Complex));
  initialL2Errors   @7  :List(Real);
  initialConditions @8  :List(List(Complex));
  fqMqInnerprods    @9  :List(Complex);
  mqMqInnerprods    @10 :List(Complex);
  aqMqInnerprods    @11 :List(Complex);
}

struct RBEIMEvaluationReal @0xf8121d2237427a80 {
  nBfs                     @0 :Integer;
  parameterRanges          @1 :ParameterRanges;
  discreteParameters       @2 :DiscreteParameterList;
  interpolationXyz         @3 :List(Point3D);
  interpolationComp        @4 :List(Integer);
  interpolationSubdomainId @5 :List(Integer);
  interpolationElemId      @6 :List(Integer);
  interpolationQp          @7 :List(Integer);
  interpolationXyzPerturb  @8 :List(List(Point3D));
  interpolationMatrix      @9 :List(Real);
}
struct RBEIMEvaluationComplex @0xc35a5eb004965455 {
  nBfs                     @0 :Integer;
  parameterRanges          @1 :ParameterRanges;
  discreteParameters       @2 :DiscreteParameterList;
  interpolationXyz         @3 :List(Point3D);
  interpolationComp        @4 :List(Integer);
  interpolationSubdomainId @5 :List(Integer);
  interpolationElemId      @6 :List(Integer);
  interpolationQp          @7 :List(Integer);
  interpolationXyzPerturb  @8 :List(List(Point3D));
  interpolationMatrix      @9 :List(Complex);
}
