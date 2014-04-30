/***************************************************************************\
|* Function Parser for C++ v4.5.1                                          *|
|*-------------------------------------------------------------------------*|
|* Function optimizer                                                      *|
|*-------------------------------------------------------------------------*|
|* Copyright: Joel Yliluoma                                                *|
|*                                                                         *|
|* This library is distributed under the terms of the                      *|
|* GNU Lesser General Public License version 3.                            *|
|* (See lgpl.txt and gpl.txt for the license text.)                        *|
\***************************************************************************/

/* NOTE:
 This file contains generated code (from the optimizer sources) and is
 not intended to be modified by hand. If you want to modify the optimizer,
 download the development version of the library.
*/

#include "fpconfig.hh"
#ifdef FP_SUPPORT_OPTIMIZER
#include "fparser.hh"
#include "extrasrc/fptypes.hh"
#include "extrasrc/fpaux.hh"
#define lZ4 );lC c3
#define lY4 )n4}if(
#define lX4 b.Opcode
#define lW4 ;tmp.nJ
#define lV4 &&p.yG2
#define lU4 iD|tL1)
#define lT4 Immeds
#define lS4 xK 2)t72
#define lR4 FP_GetOpcodeName
#define lQ4 ){case
#define lP4 tB info
#define lO4 b.Value)
#define lN4 );}xO1
#define lM4 .end()
#define lL4 iI l8 a),
#define lK4 if(n5 iK3
#define lJ4 :start_at()
#define lI4 ParamHolder
#define lH4 NumConstant:
#define lG4 info.SaveMatchedParamIndex(
#define lF4 ||tree.GetOpcode
#define lE4 },{{1,0
#define lD4 ),Value
#define lC4 (xW n93
#define lB4 Finite
#define lA4 tI1 yD1
#define l94 };enum
#define l84 xD lT 2,
#define l74 ;}void
#define l64 )l74
#define l54 :if lI
#define l44 "Found "
#define l34 stackpos
#define l24 sim.nX 1,
#define l14 "%d, cost"
#define l04 "dup(%u) "
#define iZ3 ;unsigned
#define iY3 (exponent
#define iX3 eQ{assert
#define iW3 "PUSH " cL3
#define iV3 "immed "<<
#define iU3 mFuncParsers
#define iT3 param.data
#define iS3 stderr
#define iR3 sep2=" "
#define iQ3 FPHASH_CONST
#define iP3 cache_needed[
#define iO3 fprintf
#define iN3 ::cout<<"Applying "
#define iM3 FUNCTIONPARSER_INSTANTIATE_OPTIMIZE
#define iL3 FUNCTIONPARSER_INSTANTIATE_EMPTY_OPTIMIZE
#define iK3 HANDLE_UNARY_CONST_FUNC
#define iJ3 .yB t02
#define iI3 {data->
#define iH3 AddFrom(
#define iG3 (param.
#define iF3 ,i2 2,
#define iE3 {if(n21
#define iD3 within,
#define iC3 )t2 tA2
#define iB3 c_count
#define iA3 s_count
#define i93 MaxOp
#define i83 2)lS 2*
#define i73 tmp2.nJ
#define i63 else nM
#define i53 eM1 val
#define i43 sim.x81
#define i33 },0,0x4
#define i23 yD2 iJ
#define i13 .empty()
#define i03 ].swap(
#define tZ3 =synth.
#define tY3 codes[b
#define tX3 whydump
#define tW3 middle2
#define tV3 ::string
#define tU3 {switch
#define tT3 nparams
#define tS3 ){std::
#define tR3 l4 4,1,
#define tQ3 l3 2,2,
#define tP3 l2 2,2,
#define tO3 l4 16,1,
#define tN3 );o<<"\n";
#define tM3 l4 0,1,
#define tL3 ,l1 0x12 xI2
#define tK3 iT 1,0,
#define tJ3 nQ 0,
#define tI3 cAbs nQ
#define tH3 ;case
#define tG3 iH tH3
#define tF3 fp_pow(
#define tE3 typename
#define tD3 false;}
#define tC3 cAbsIf)
#define tB3 size_t>
#define tA3 l71++b)
#define t93 .second
#define t83 ]t93
#define t73 ].first
#define t63 cH r);}
#define t53 Ne_Mask
#define t43 Gt_Mask
#define t33 Lt_Mask
#define t23 opcode,
#define t13 FindPos
#define t03 )i01 xH=
#define eZ3 tree l8 2
#define eY3 xC cond nC
#define eX3 tree xC
#define eW3 ,xD1);lC
#define eV3 xC cMul);
#define eU3 resize(
#define eT3 public:
#define eS3 pclone
#define eR3 for(eP3<
#define eQ3 a=0 eP3
#define eP3 ;a
#define eO3 :if(&*tG1){
#define eN3 :{n61 r=
#define eM3 --c81.
#define eL3 eS n32
#define eK3 xJ(tZ2
#define eJ3 newpow
#define eI3 change
#define eH3 (i53
#define eG3 (count
#define eF3 133,2,
#define eE3 Params
#define eD3 Needs
#define eC3 byteCode
#define eB3 t2 cV1);
#define eA3 lY1 nC==
#define e93 cLog2by
#define e83 iP2 n11
#define e73 factor_t
#define e63 value1
#define e53 fp_mod(
#define e43 else{if(
#define e33 xJ());nD
#define e23 nN1);}if
#define e13 Ge0Lt1
#define e03 Others
#define cZ3 tG1=r.specs;if(r.found){
#define cY3 info=(*x5)[
#define cX3 ,tG1,info
#define cW3 ,tree,info
#define cV3 ,cLess eS
#define cU3 ;size_t
#define cT3 &&IsLogicalValue(
#define cS3 ,leaf1 l8
#define cR3 ;}static yV1
#define cQ3 );if(
#define cP3 iI))cJ
#define cO3 cG ifp2
#define cN3 cAbsNot
#define cM3 IsImmed(
#define cL3 ;DumpTree
#define cK3 switch(tJ
#define cJ3 }switch
#define cI3 stackptr
#define cH3 cLog);xK
#define cG3 l8 0));
#define cF3 opcodes
#define cE3 did_muli
#define cD3 Value){
#define cC3 yG const
#define cB3 used[b]
#define cA3 );exponent
#define c93 if(a>0){
#define c83 ,iE,1,lA1+1);
#define c73 break;}
#define c63 ;break;
#define c53 (constraints&
#define c43 c63 cJ3(
#define c33 param=*(
#define c23 sizeof(
#define c13 cAbsIf,
#define c03 cNotNot,
#define yZ3 l2 18,2,
#define yY3 cLess,cF
#define yX3 cTan iE1
#define yW3 450998,
#define yV3 cLog iE1
#define yU3 cExp2 nQ
#define yT3 ,1,562 nP
#define yS3 ,cEqual eS
#define yR3 7168,
#define yQ3 lJ 2},0,
#define yP3 Gt0Le1
#define yO3 IsLogicalValue iI l8
#define yN3 p1 cH p1)
#define yM3 return p.yW
#define yL3 ))return
#define yK3 ;if(half
#define yJ3 xH l64
#define yI3 default:
#define yH3 range nY3
#define yG3 ==cOr)?0:1))
#define yF3 range<yY2
#define yE3 range xJ
#define yD3 cAdd lS2
#define yC3 if(op1==
#define yB3 tI1 nE1
#define yA3 iterator
#define y93 begin();
#define y83 TreeSet
#define y73 parent
#define y63 insert(i
#define y53 newrel
#define y43 IsNever
#define y33 void set
#define y23 b_needed
#define y13 cachepos
#define y03 half=
#define xZ3 131,4,1,
#define xY3 131,8,1,
#define xX3 4,1,2,1,
#define xW3 src_pos
#define xV3 reserve(
#define xU3 n01 void
#define xT3 treeptr
#define xS3 tV1 void
#define xR3 ImmedTag
#define xQ3 a,const
#define xP3 RefCount
#define xO3 Birth();
#define xN3 mulgroup
#define xM3 cost_t
#define xL3 iftree
#define xK3 fpdata
#define xJ3 middle
#define xI3 nG1 1,
#define xH3 sqrt_cost
#define xG3 const int
#define xF3 mul_count
#define xE3 if iI eQ1
#define xD3 maxValue1
#define xC3 minValue1
#define xB3 maxValue0
#define xA3 minValue0
#define x93 ValueType
#define x83 PlusInf
#define x73 1),l03(1))
#define x63 result
#define x53 x63;}case
#define x43 x63 x42
#define x33 ;range.xS2
#define x23 x63 cQ
#define x13 .value;else
#define x03 x63 cT=false;if(
#define nZ3 x63 tZ1
#define nY3 xJ x63
#define nX3 max.known
#define nW3 abs_mul
#define nV3 pos_set
#define nU3 y3 cAdd);
#define nT3 subtree
#define nS3 invtree
#define nR3 MakeHash(
#define nQ3 parampair
#define nP3 rulenumit
#define nO3 cAnd l3
#define nN3 MakeEqual
#define nM3 MakeTrue,
#define nL3 newbase
#define nK3 fp_equal(
#define nJ3 branch1op
#define nI3 branch2op
#define nH3 )val=l33
#define nG3 ),cH1 y8
#define nF3 ,cAnd,l6
#define nE3 ::vector<unsigned>&t51
#define nD3 ::ByteCodeSynth xJ
#define nC3 std::vector<bool>
#define nB3 tree l8 a)
#define nA3 ;for nD2 a=
#define n93 l8 a)xE
#define n83 overlap
#define n73 truth_b
#define n63 truth_a
#define n53 found_dup
#define n43 )lS 3*3*
#define n33 continue;
#define n23 lH3 size(
#define n13 ContainsOtherCandidates
#define n03 eS1 n41
#define lZ3 rangeutil
#define lY3 synth lP2
#define lX3 Plan_Has(
#define lW3 StackMax)
#define lV3 iH true;}
#define lU3 x63 xP
#define lT3 const xM2
#define lS3 namespace
#define lR3 DelParam(
#define lQ3 ByteCode[
#define lP3 inverted
#define lO3 )iH lH
#define lN3 y43:
#define lM3 .known&&
#define lL3 depcodes
#define lK3 explicit
#define lJ3 cCosh nQ
#define lI3 VarBegin
#define lH3 eE3.
#define lG3 )){data x7
#define lF3 ].data);
#define lE3 (nK3
#define lD3 default_function_handling
#define lC3 sim.Eat(
#define lB3 .what nV1
#define lA3 ){lC3
#define l93 begin(),
#define l83 cond_add
#define l73 cond_mul
#define l63 cond_and
#define l53 IsAlways
#define l43 lQ4 cImmed:
#define l33 func lM1
#define l23 const eN
#define l13 .size();++
#define l03 Value_t
#define iZ2 ,cPow,
#define iY2 ,cGreater
#define iX2 t93);
#define iW2 l03&
#define iV2 bool eB1
#define iU2 Optimize()
#define iT2 costree
#define iS2 sintree
#define iR2 leaf_count
#define iQ2 sub_params
#define iP2 x63.
#define iO2 printf(
#define iN2 swap(tmp);
#define iM2 cbrt_count
#define iL2 sqrt_count
#define iK2 );nD lC
#define iJ2 eM1 n3 0),
#define iI2 (x63
#define iH2 min.n3 0),
#define iG2 fp_abs(yG2))
#define iF2 fp_abs(yW)
#define iE2 ,std::cout
#define iD2 p1 cG ifp1
#define iC2 pcall_tree
#define iB2 after_powi
#define iA2 (half&63)-1;
#define i92 GetHash().
#define i82 TreeCountItem
#define i72 if(op==
#define i62 info.lQ[b].
#define i52 info=info;
#define i42 yK tD3
#define i32 grammar
#define i22 cEqual,
#define i12 cLog nQ
#define i02 cNeg,lT 1,
#define tZ2 ),0},{
#define tY2 tree nC
#define tX2 MakeNEqual
#define tW2 l7 op1 tI1
#define tV2 Dump(std::
#define tU2 isInteger(
#define tT2 Comparison
#define tS2 ,bool abs)
#define tR2 needs_flip
#define tQ2 (l03
#define tP2 iQ1 apos==
#define tO2 value]
#define tN2 ~size_t(0)
#define tM2 xJ1 xR+1);
#define tL2 ,lA 0x12 nH
#define tK2 {std::cout<<
#define tJ2 const std::eO
#define tI2 Rule&rule,
#define tH2 ,const e0&
#define tG2 >::res,b8<
#define tF2 n33}
#define tE2 cG tree);
#define tD2 mul_item
#define tC2 innersub
#define tB2 xN3.
#define tA2 xN3)
#define t92 cbrt_cost
#define t82 best_cost
#define t72 y3 cPow)
#define t62 condition
#define t52 per_item
#define t42 item_type
#define t32 }data;data.
#define t22 (eR)));nZ
#define t12 eZ2 x63(
#define t02 l03(2)));
#define eZ2 ){l03
#define eY2 nR eZ2 tmp=
#define eX2 ,tree nW
#define eW2 first2
#define eV2 l4 18,1,
#define eU2 cIf,iT 3,
#define eT2 ,cMul x4
#define eS2 lJ 1},0,
#define eR2 tI 1},0,
#define eQ2 Decision
#define eP2 not_tree
#define eO2 group_by
#define eN2 ;eW iT1
#define eM2 (list.first
#define eL2 );tB2
#define eK2 );t3=!t3;}
#define eJ2 exponent=
#define eI2 (leaf1 l8 1)xF1
#define eH2 ->second
#define eG2 targetpos
#define eF2 ParamSpec
#define eE2 rhs.hash2;}
#define eD2 rhs.hash1
#define eC2 struct
#define eB2 Forget()
#define eA2 source_tree
#define e92 GetParam tM
#define e82 :xK l03(
#define e72 ;sim.Push(
#define e62 .cM3
#define e52 1)e62
#define e42 nC==cLog2&&
#define e32 (xN3
#define e22 <tW,xM3>
#define e12 CodeTree lW
#define e02 p1_evenness
#define cZ2 isNegative(
#define cY2 (lR)iK2
#define cX2 neg_set
#define cW2 );sim.nX 2,
#define cV2 nC==cNot||
#define cU2 if iI nC==
#define cT2 cNop,cNop}}
#define cS2 cTanh,cNop,
#define cR2 NewHash
#define cQ2 >eC2 c7<
#define cP2 matches
#define cO2 else{x5=new
#define cN2 .match_tree
#define cM2 eS1 void*)&
#define cL2 cGreater,cF
#define cK2 cTan nQ
#define cJ2 cCos nQ
#define cI2 (std::move(
#define cH2 if iI l8 0)
#define cG2 cM3)||
#define cF2 +=1 iH nC1;
#define cE2 negated
#define cD2 Specializer
#define cC2 params
#define cB2 coshtree
#define cA2 sinhtree
#define c92 best_score
#define c82 mulvalue
#define c72 pow_item
#define c62 subgroup
#define c52 nC==cPow&&tL
#define c42 PowiResult
#define c32 maxValue
#define c22 minValue
#define c12 fp_min(yU,
#define c02 div_tree
#define yZ2 pow_tree
#define yY2 l03 nU
#define yX2 preserve
#define yW2 PullResult()
#define yV2 dup_or_fetch
#define yU2 nominator]
#define yT2 ;TriTruthValue
#define yS2 )c63}
#define yR2 ;pow xC cPow);pow
#define yQ2 Rehash(false
#define yP2 test_order
#define yO2 nQ3,
#define yN2 .param_count
#define yM2 shift(index)
#define yL2 rulenumber
#define yK2 cTanh nQ
#define yJ2 cSinh nQ
#define yI2 cInv,lT 1,
#define yH2 constraints=
#define yG2 max.val
#define yF2 factor_immed
#define yE2 changes
#define yD2 .Rehash()
#define yC2 yT1 c73}
#define yB2 cG leaf2 l8
#define yA2 cG leaf1 l8
#define y92 cG cond l8
#define y82 for(tE3
#define y72 exp_diff
#define y62 ExponentInfo
#define y52 lower_bound(
#define y42 factor
#define y32 is_logical
#define y22 newrel_and
#define y12 .IsDefined(
#define y02 ;iE.Remember(
#define xZ2 iH Unknown;}
#define xY2 res_stackpos
#define xX2 half_pos
#define xW2 ifdata.ofs
#define xV2 >>1)):(
#define xU2 CodeTreeData
#define xT2 exponent)
#define xS2 multiply(
#define xR2 (IfData&ifdata
#define xQ2 iB push_back(
#define xP2 iB size()
#define xO2 var_trees
#define xN2 const n41
#define xM2 CodeTree&
#define xL2 Value(Value::
#define xK2 tree l8 1)
#define xJ2 yD2;
#define xI2 },{{3,
#define xH2 parent_opcode
#define xG2 cM3)){if
#define xF2 log2_exponent
#define xE2 lR,eR)iK2
#define xD2 dup_fetch_pos
#define xC2 {e0 start_at;
#define xB2 y43 cR lC
#define xA2 cSin nQ
#define x92 Value_EvenInt
#define x82 AddCollection
#define x72 ConditionType
#define x62 SpecialOpcode
#define x52 =i eH2.
#define x42 cQ=false
#define x32 fp_max(yU);
#define x22 assimilated
#define x12 denominator
#define x02 fraction
#define nZ2 synth.Find
#define nY2 .GetDepth()
#define nX2 );p2 cH p2)xX1
#define nW2 nB OPCODE
#define nV2 ,eP,synth);
#define nU2 DUP_BOTH();
#define nT2 template lL
#define nS2 0x80000000u
#define nR2 IsDescendantOf
#define nQ2 found_log2
#define nP2 :lD m.min.set(
#define nO2 TopLevel)
#define nN2 a;if(nF2){x5=(
#define nM2 );nZ l5::
#define nL2 ,l03(-1)))xB
#define nK2 synth.PushImmed(
#define nJ2 cL3 cP
#define nI2 :lC3 1,
#define nH2 ,l03(1))){
#define nG2 bool t3=false;
#define nF2 &*start_at
#define nE2 ,yA yZ3
#define nD2 (size_t
#define nC2 nD2 a y0
#define nB2 SetOpcode(
#define nA2 div_params
#define n92 immed_sum
#define n82 lQ3++IP]
#define n72 OPCODE(opcode)
#define n62 FactorStack xJ
#define n52 l53 cR lC
#define n42 cLessOrEq,
#define n32 282870 nP
#define n22 cNotNot nQ
#define n12 cNot nQ
#define n02 DumpHashesFrom
#define lZ2 replacing_slot
#define lY2 RefParams
#define lX2 if_always[
#define lW2 WhatDoWhenCase
#define lV2 exponent_immed
#define lU2 new_base_immed
#define lT2 base_immed
#define lS2 ||op1==
#define lR2 data[a t83
#define lQ2 if(newrel_or==
#define lP2 .AddOperation(
#define lO2 DUP_ONE(apos);
#define lN2 flipped
#define lM2 .UseGetNeeded(
#define lL2 eA 2,131,
#define lK2 [xR-1-offset].
#define lJ2 lQ3 a
#define lI2 y6 Immed.size());
#define lH2 1 y6 xP2
#define lG2 OptimizedUsing
#define lF2 Var_or_Funcno
#define lE2 lF2;
#define lD2 GetParams(
#define lC2 crc32_t
#define lB2 signed_chain
#define lA2 MinusInf
#define l92 n_immeds
#define l82 stack.size()
#define l72 FindClone(xH
#define l62 );xH.SetParamsMove(
#define l52 lQ3 IP]
#define l42 GetOpcode())
#define l32 needs_rehash
#define l22 AnyWhere_Rec
#define l12 minimum_need
#define l02 ~unsigned(0)
#define iZ1 41,42,43,44,
#define iY1 yC xA,l5::i7
#define iX1 p1_logical_b
#define iW1 p0_logical_b
#define iV1 p1_logical_a
#define iU1 p0_logical_a
#define iT1 xC tY2);
#define iS1 2*2)lS 3*
#define iR1 ,PowiCache&iE,
#define iQ1 else if(
#define iP1 synth.DoDup(
#define iO1 cache_needed
#define iN1 eA 2,1,eA 2,
#define iM1 treelist
#define iL1 has_bad_balance
#define iK1 e73 y42
#define iJ1 set(fp_ceil);cX
#define iI1 std::cout<<"POP "
#define iH1 stack[l82-
#define iG1 stack.push_back
#define iF1 lQ4 l53:
#define iE1 ,l4 2,1,
#define iD1 );synth.StackTopIs(
#define iC1 divgroup
#define iB1 iI2))break;x63*=
#define iA1 {if(GetOpcode()
#define i91 cNEqual
#define i81 ));nB1 cG y5 l8
#define i71 tI 2},0,0x0},{{
#define i61 Oneness_NotOne|
#define i51 Value_IsInteger
#define i41 Constness_Const
#define i31 lG2(
#define i21 reltype
#define i11 SequenceOpcodes
#define i01 {CodeTree xJ
#define tZ1 .yG2)
#define tY1 goto fail;}
#define tX1 l1 0x4 nH
#define tW1 template<
#define tV1 lY2);
#define tU1 template lX
#define tT1 TreeCountType xJ
#define tS1 );lR3 a);}}
#define tR1 >tQ2(1),
#define tQ1 l03(0.0)){xN
#define tP1 n72);
#define tO1 0.5))t72;lC
#define tN1 MaxChildDepth
#define tM1 repl_param_list,
#define tL1 (unsigned
#define tK1 ;for tL1
#define tJ1 tI1 SetParam(
#define tI1 );tree.
#define tH1 synth.FindAndDup iI);
#define tG1 (*x5 t01
#define tF1 std::pair<It,It>
#define tE1 cPow,lA
#define tD1 Sign_Negative
#define tC1 Value_Logical
#define tB1 yC MakeFalse,{l5
#define tA1 new_factor_immed
#define t91 if(remaining[a])
#define t81 tU1 void
#define t71 e9(),std::vector<
#define t61 ,l0 2,
#define t51 ByteCode,size_t&IP,size_t limit,size_t y1
#define t41 yL3 true;
#define t31 synth.xL 1
#define t21 ,iE nV2
#define t11 :public e9,public std::vector<
#define t01 )[a].start_at
#define eZ1 for nD2 b=0;b<
#define eY1 .GetImmed()
#define eX1 CollectionSet xJ
#define eW1 .back().thenbranch
#define eV1 (!tree e62)cR
#define eU1 eT1.erase(cs_it);
#define eT1 TreeCounts
#define eS1 (const
#define eR1 eS1 iW2 v
#define eQ1 l8 e52)&&
#define eP1 )nQ3 t93
#define eO1 ,cMul l3 16,
#define eN1 .min.set(fp_floor);
#define eM1 m.max.
#define eL1 <<tree.i92
#define eK1 tree l8 0)
#define eJ1 .Become iI l8
#define eI1 occurance_pos
#define eH1 exponent_hash
#define eG1 exponent_list
#define eF1 CollectMulGroup(
#define eE1 source_set
#define eD1 exponent,y83
#define eC1 *const func)(
#define eB1 operator
#define eA1 retry_anyparams_3
#define e91 retry_anyparams_2
#define e81 needlist_cached_t
#define e71 grammar_rules[*r]
#define e61 yQ3 0x4},{{1,
#define e51 ,cPow iE1
#define e41 eS2 0x4},{{1,
#define e31 CodeTreeImmed xJ(
#define e21 by_float_exponent
#define e11 nK3 exponent
#define e01 new_exp
#define cZ1 =comp.AddItem(atree
#define cY1 return BecomeZero;
#define cX1 return BecomeOne;
#define cW1 if(lQ.size()<=n1)
#define cV1 addgroup
#define cU1 found_log2by
#define cT1 nD2 a=c5 a-->0;)
#define cS1 }c73 case
#define cR1 nC==cN3)
#define cQ1 if(keep_powi
#define cP1 ParsePowiMuli(
#define cO1 lF2)
#define cN1 eS 529654 nP
#define cM1 branch1_backup
#define cL1 l71 a-->0;)if(
#define cK1 branch2_backup
#define cJ1 exponent_map
#define cI1 plain_set
#define cH1 rangehalf
#define cG1 LightWeight(
#define cF1 ;iQ1!x63
#define cE1 lS3 FPoptimizer_Optimize
#define cD1 if(value
#define cC1 tU1 yP
#define cB1 tU1 static
#define cA1 ::MakeTrue
#define c91 (cond yN&&cond eH))
#define c81 NeedList
#define c71 lM4&&i->first==
#define c61 should_regenerate=true;
#define c51 should_regenerate,
#define c41 Collection
#define c31 CodeTree xJ r;r xC
#define c21 RelationshipResult
#define c11 Subdivide_Combine(
#define c01 long value
#define yZ1 )const yK
#define yY1 rhs yZ1 hash1
#define yX1 best_sep_factor
#define yW1 needlist_cached
#define yV1 inline unsigned
#define yU1 t23 bool pad
#define yT1 .lR3 a);
#define yS1 .GetParamCount()
#define yR1 tK1 eQ3<xU;++a)
#define yQ1 changed=true;
#define yP1 MakesInteger(
#define yO1 const iW2 value
#define yN1 best_sep_cost
#define yM1 MultiplicationRange
#define yL1 pihalf_limits
#define yK1 n_stacked
#define yJ1 cR2.hash1
#define yI1 AnyParams_Rec
#define yH1 Become(value l8 0))
#define yG1 PositionalParams,0}
#define yF1 ByteCodeSynth xJ&synth)
#define yE1 .sep_list[n]
#define yD1 AddParamMove(
#define yC1 always_sincostan
#define yB1 Recheck_RefCount_Div
#define yA1 Recheck_RefCount_Mul
#define y91 xN3;xN3 xC
#define y81 MultiplyAndMakeLong(
#define y71 l03(0)
#define y61 covers_plus1
#define y51 max.template set_if<
#define y41 l8 a)e62))
#define y31 tI1 lR3
#define y21 xJ())y3 cMul);lC
#define y11 if(synth.FindAndDup(
#define y01 SynthesizeParam(
#define xZ1 ;std::cout<<
#define xY1 grammar_func
#define xX1 ;eX3 xL3 nC);cJ}
#define xW1 cOr l3 16,1,
#define xV1 252421 nP 24830,
#define xU1 l2 0,2,165888 nP
#define xT1 cCos iE1
#define xS1 l1 0x12 nH
#define xR1 Modulo_Radians},
#define xQ1 PositionType
#define xP1 CollectionResult
#define xO1 tU1 bool
#define xN1 const_offset
#define xM1 inline TriTruthValue
#define xL1 stacktop_desired
#define xK1 int mStackPtr=0;
#define xJ1 SetStackTop(
#define xI1 FPoptimizer_ByteCode
#define xH1 1)?(poly^(
#define xG1 y71)
#define xF1 xE leaf2 l8
#define xE1 (p1 eY1
#define xD1 cond_type
#define xC1 fphash_value_t
#define xB1 Recheck_RefCount_RDiv
#define xA1 tmp.AddParamMove iI);
#define x91 cMul)lW4 0));tmp.
#define x81 SwapLastTwoInStack();
#define x71 ParamSpec_Extract xJ(
#define x61 fPExponentIsTooLarge(
#define x51 CollectMulGroup_Item(
#define x41 pair<l03,y83>
#define x31 nL xJ1 xR-1);
#define x21 covers_full_cycle
#define x11 AssembleSequence(
#define x01 252180 nP 281854,
#define nZ1 {DataP slot_holder(xX[
#define nY1 <<std::dec<<")";}
#define nX1 yC MakeNotP1,l5::
#define nW1 yC MakeNotP0,l5::
#define nV1 !=xA)if(TestCase(
#define nU1 }inline
#define nT1 std::pair<T1,T2>&
#define nS1 tW1 tE3
#define nR1 has_good_balance_found
#define nQ1 n_occurrences
#define nP1 found_log2_on_exponent
#define nO1 covers_minus1
#define nN1 tree.lR3 a
#define nM1 iP1 found[data.
#define nL1 nB3 eY1;
#define nK1 needs_resynth
#define nJ1 immed_product
#define nI1 ,2,1 xM if(found[data.
#define nH1 c73 switch(bitmask&
#define nG1 cMul l3 0,
#define nF1 Sign_Positive
#define nE1 SetParamMove(
#define nD1 CodeTreeImmed tQ2(
#define nC1 Suboptimal
#define nB1 changed_if
#define nA1 n_as_tanh_param
#define n91 cA3 xJ2
#define n81 opposite=
#define n71 xC1(
#define n61 MatchResultType
#define n51 (long double)
#define n41 CodeTree xJ&
#define n31 yK CodeTree xJ(
#define n21 needs_sincos
#define n11 resulting_exponent
#define n01 ;tU1
#define lZ1 Unknown:yI3;}
#define lY1 GetParam(a)
#define lX1 inverse_nominator]
#define lW1 cSin iE1
#define lV1 cM3)eZ2
#define lU1 void FunctionParserBase
#define lT1 ,cIf,l0 3,
#define lS1 xC cLog);tree eV3
#define lR1 SetParams(lD2));
#define lQ1 GetPositivityInfo iI)!=
#define lP1 o<<"("<<std::hex<<data.
#define lO1 ;cR2.hash2+=
#define lN1 for nD2 eQ3<c5++a)
#define lM1 (val);else*this=model;}
#define lL1 IfBalanceGood(
#define lK1 n_as_tan_param
#define lJ1 changed_exponent
#define lI1 inverse_denominator
#define lH1 retry_positionalparams_2
#define lG1 unsigned index
#define lF1 situation_flags&
#define lE1 518 nP 400412,
#define lD1 7168 nP 401798
#define lC1 data.subfunc_opcode
#define lB1 CopyOnWrite();
#define lA1 recursioncount
#define l91 PlanNtimesCache(
#define l81 FPoptimizer_Grammar
#define l71 GetParamCount();
#define l61 AddOperation(cInv,1,1 xM}
#define l51 ParamSpec_SubFunctionData
#define l41 lV4<l03(
#define l31 ;lB1 lR3 a);}
#define l21 );xN2
#define l11 PositionalParams_Rec
#define l01 yC MakeNotNotP1,l5::
#define iZ yC MakeNotNotP0,l5::
#define iY xJ2 eX3 op1 tI1 DelParams()
#define iX lQ3 xW2+
#define iW DumpTreeWithIndent(*this);
#define iV tW1 unsigned Compare>
#define iU >=xG1
#define iT lA 0x4},{{
#define iS yQ3 0x0},{{1,
#define iR edited_powgroup
#define iQ has_unknown_max
#define iP has_unknown_min
#define iO static const yE3
#define iN p0 cT&&p0.yW>=l03(0.0))
#define iM synthed_tree
#define iL SelectedParams,0},0,0x0},{{
#define iK tU3(type lQ4 cond_or:
#define iJ ;CodeTree xJ
#define iI (tree
#define iH ;return
#define iG )lV3
#define iF collections
#define iE cache
#define iD ;xQ2 nS2
#define iC )iD);
#define iB ByteCode.
#define iA goto ReplaceTreeWithOne tH3
#define i9 ]);lY3
#define i8 )xZ1 std::endl;DumpHashes(
#define i7 Never yC xA,l5::Never}
#define i6 !=xA)return lX2
#define i5 e21.data
#define i4 lK3 xU2(
#define i3 needs_sinhcosh
#define i2 cAdd l3 0,
#define i1 tZ2 l03(
#define i0 tU1 nA
#define tZ MakeFalse,l5::
#define tY ].relationship
#define tX }},{ProduceNewTree,2,1,
#define tW int_exponent_t
#define tV ,cAdd x4
#define tU CalculateResultBoundaries iI l8
#define tT p0=CalculateResultBoundaries(
#define tS ,ByteCode,IP,limit,y1,stack);
#define tR 408964 nP 24963,
#define tQ 528503 nP 24713,
#define tP matched_params
#define tO [n1 t73=true;lQ[n1 t83
#define tN l81::Grammar*
#define tM (a);bool needs_cow=GetRefCount()>1;
#define tL powgroup l8
#define tK i01 tmp;tmp xC
#define tJ GetLogicalValue iI l8
#define tI xD AnyParams,
#define tH relationships
#define tG tN2&&found[data.
#define tF nD1(
#define tE has_mulgroups_remaining
#define tD by_exponent
#define tC const l51
#define tB MatchInfo xJ&
#define tA Rehash();iQ2.push_back(
#define t9 best_factor
#define t8 RootPowerTable xJ::RootPowers[
#define t7 MatchPositionSpec_AnyParams xJ
#define t6 lS3 FPoptimizer_CodeTree
#define t5 n_as_sinh_param
#define t4 n_as_cosh_param
#define t3 is_signed
#define t2 ;yD1
#define t1 lD2));xN3 xJ2
#define t0 result_positivity
#define eZ biggest_minimum
#define eY 122999 nP 139399,
#define eX 142455 nP 141449,
#define eW cond_tree
#define eV else_tree
#define eU then_tree
#define eT nM l03(-x73;
#define eS ,yZ3
#define eR xK2 eY1
#define eQ n41 tree)
#define eP sequencing
#define eO string lR4(
#define eN std::vector<CodeTree xJ>
#define eM if_stack
#define eL n_as_sin_param
#define eK n_as_cos_param
#define eJ PowiResolver::
#define eI cIf,tM3
#define eH .BalanceGood
#define eG yD1 c62
#define eF valueType
#define eE back().endif_location
#define eD xC1 key
#define eC xJ2 eX3 leaf1 nC tI1
#define eB yD1 mul);
#define eA 130,1,
#define e9 MatchPositionSpecBase
#define e8 lK3 CodeTree(
#define e7 smallest_maximum
#define e6 }PACKED_GRAMMAR_ATTRIBUTE;
#define e5 ))i01
#define e4 :{AdoptChildrenWithSameOpcode iI);
#define e3 :goto ReplaceTreeWithZero tH3
#define e2 ReplaceTreeWithParam0;
#define e1 factor_needs_rehashing
#define e0 MatchPositionSpecBaseP
#define cZ tE3 tT1::yA3
#define cY x71 nN.param_list,
#define cX return m;}case
#define cW 243,244,245,246,249,250,251,253,255,256,257,258,259}};}
#define cV ];};extern"C"{
#define cU l03(1.5)*fp_const_pi xJ()
#define cT .min.known
#define cS cQ&&p0.yG2<=fp_const_negativezero xJ())
#define cR )return false;
#define cQ .nX3
#define cP iI)xZ1"\n";
#define cO 79,122,123,160,161,163,164,165,166,167,168,169,178,179,180,200,204,212,216,224,236,237,239,240,
#define cN 27,28,29,30,31,32,33,35,36,
#define cM const ParamSpec_SubFunction
#define cL const ParamSpec_ParamHolder
#define cK otherhalf
#define cJ goto redo;
#define cI StackState
#define cH .Rehash(lA4
#define cG .AddParam(
#define cF l2 16,2,
#define cE const SequenceOpCode xJ
#define cD paramholder_matches[x0]
#define cC nE1 0,xT2;lR3 1);
#define cB MatchPositionSpec_PositionalParams xJ
#define cA xN2 tree,std::ostream&o
#define c9 paramholder_matches.
#define c8 CalculatePowiFactorCost(
#define c7 ImmedHashGenerator
#define c6 );void c1 unsigned t23 cD2<
#define c5 tree.l71
#define c4 ::map<fphash_t,std::set<std tV3> >
#define c3 ComparisonSetBase::
#define c2 yD1 comp.cI1[a].value);
#define c1 AddFunctionOpcode(
#define c0 T1,tE3 T2>inline iV2()(
#define yZ has_nonlogical_values
#define yY from_logical_context)
#define yX CodeTree xJ tmp,tmp2;tmp2 xC
#define yW min.val
#define yV +=fp_const_twopi xJ();
#define yU fp_sin(min),fp_sin(max))
#define yT fp_const_twopi xJ()cQ3
#define yS AnyParams,0}},{ProduceNewTree,
#define yR for nD2 a=xW.l71 a-->0;)
#define yQ POWI_CACHE_SIZE
#define yP static inline CodeTree xJ
#define yO ++IP;tF2 if(l52==cF3.
#define yN .FoundChild
#define yM BalanceResultType
#define yL xP3(0),Opcode(
#define yK {return
#define yJ const yK data->
#define yI lP2 GetOpcode(),
#define yH for nD2 eQ3<l71++a){if(
#define yG static void nR3 nB fphash_t&cR2,
#define yF MatchPositionSpec_AnyWhere
#define yE if iG3 data.match_type==
#define yD void OutFloatHex(std::ostream&o,
#define yC },{l5::
#define yB AddParam(CodeTreeImmed(
#define yA cGreaterOrEq,
#define y9 ,tE3 CodeTree xJ::
#define y8 xJ model=cH1 xJ()){if(known
#define y7 AssembleSequence_Subdivide(
#define y6 ]=nS2|unsigned(
#define y5 branch2
#define y4 !=tN2){nM1
#define y3 ;lC3 2,
#define y2 unsigned c iZ3 short l[
#define y1 factor_stack_base
#define y0 =0 eP3<y73.l71++a)if(
#define xZ ,tE1 0x4 nH
#define xY (nP3 r=range.first;r!=range t93;++r){
#define xX data->eE3
#define xW branch1
#define xV );tF2 if eM2 eY1==l03(
#define xU nN yN2
#define xT =fp_cosh(m.yW);i53=fp_cosh eH3);
#define xS {eT1.erase(i);tF2
#define xR StackTop
#define xQ FPOPT_autoptr
#define xP +=x63 iH x63;}tU1 inline l03
#define xO int_exponent
#define xN tree.ReplaceWithImmed(
#define xM iD1*this)iH;}
#define xL GetStackTop()-
#define xK sim.AddConst(
#define xJ <l03>
#define xI .SetParam(0,xL3 l8 0))iJ p1;p1 xC
#define xH newnode
#define xG has_highlevel_opcodes
#define xF eS2 0x0},{{
#define xE .IsIdenticalTo(
#define xD ,cAdd,
#define xC .nB2
#define xB {if(needs_cow){lB1 goto
#define xA Unchanged
#define x9 cA=std::cout
#define x8 best_selected_sep
#define x7 ->Recalculate_Hash_NoRecursion();}
#define x6 l71++a)if(ApplyGrammar(i32,nB3,
#define x5 position
#define x4 l3 2,1,
#define x3 std::vector<CodeTree>
#define x2 ;lN1{yE3
#define x1 TestImmedConstraints iG3 constraints,tree)cR
#define x0 paramholder_index
#define nZ return true tH3
#define nY occurance_counts
#define nX x81 lC3
#define nW )){tree.FixIncompleteHashes();}
#define nV );cQ1 lA3 1,cInv yS2 xK-1)t72;lC
#define nU >p=tU a)cQ3 p.
#define nT ;i73 0));tmp xC cInv);tmp.yD1 tmp2)iH
#define nS -->0;){xN2 powgroup=lY1;if(powgroup
#define nR cH2 e62)
#define nQ ,l0 1,
#define nP ,{2,
#define nO const FPoptimizer_CodeTree::n41 tree
#define nN model_tree
#define nM return yE3(
#define nL ){using lS3 FUNCTIONPARSERTYPES;
#define nK eN&lY2
#define nJ AddParam iI l8
#define nI ConstantFolding_LogicCommon iI,c3
#define nH },{{2,
#define nG nS1 Ref>inline void xQ<Ref>::
#define nF AnyParams,1},0,0x0},{{
#define nE ):data(new xU2 xJ(
#define nD goto do_return;}
#define nC .GetOpcode()
#define nB FUNCTIONPARSERTYPES::
#define nA xU2 xJ::xU2(
#define n9 b;}};tW1>eC2 Comp<nB
#define n8 lF2(),eE3(),Hash(),Depth(1),i31 0){}
#define n7 SynthesizeByteCode(synth);
#define n6 while(ApplyGrammar(cM2
#define n5 GetIntegerInfo iI l8 0))==l53)goto e2
#define n4 ;tree.yD1 nB1 iG
#define n3 template set_if<cGreater>tQ2(
#define n2 DumpParams xJ iG3 data.param_list,iT3 yN2,o);
#define n1 restholder_index
#define n0 CodeTree xJ exponent;exponent xC cMul cA3 cG
#define lZ lR cQ3 fp_nequal(tmp,xG1){xN l03(1)/tmp);nD}lC
#define lY :if(ParamComparer xJ()(eE3[1],eE3[0])tS3 swap(eE3[0],eE3[1]);Opcode=
#define lX <tE3 l03>
#define lW xJ tmp;tmp xC cPow)lW4 0));tmp.yB l03(
#define lV i41,0x0},
#define lU yD1 pow l8 1));pow.lR3 1);pow.Rehash(yB3 0,pow);goto NowWeAreMulGroup;}
#define lT GroupFunction,0},lV{{
#define lS ,l03(1)/l03(
#define lR eK1 eY1
#define lQ restholder_matches
#define lP yJ1|=key;xC1 crc=(key>>10)|(key<<(64-10))lO1((~n71 crc))*3)^1234567;}};
#define lO nB1;nB1 iT1 nB1.AddParamMove iI cG3 nB1 cG xW l8
#define lN tU1 CodeTree xJ::CodeTree(
#define lM tree.SetParam(0,eK1 l8 0)tJ1 1,CodeTreeImmed(
#define lL lX void ByteCodeSynth xJ::c1 unsigned t23 cD2<
#define lK cMul,lT 2,
#define lJ cMul,AnyParams,
#define lI iI l8 0)e62)&&xK2 e62)){xN
#define lH CalculateResultBoundaries(tmp);}case
#define lG :eI3=comp.AddRelationship(atree l8 0),atree l8 1),c3
#define lF cPow,l0 2
#define lE tE3 l03>inline iV2()eS1 iW2 xQ3 iW2 b)yK a
#define lD {yE3 m=tU 0));
#define lC break tH3
#define lB t81 CodeTree xJ::
#define lA yG1,0,
#define l9 l1 0x0 nH
#define l8 .GetParam(
#define l7 iJ nB1;nB1 iT1 nB1.SetParamsMove iI.lD2));nB1 xJ2 eX3
#define l6 SelectedParams,0},0,0x0 nH
#define l5 RangeComparisonData
#define l4 yG1},{ProduceNewTree,
#define l3 ,AnyParams,0}},{ReplaceParams,
#define l2 yG1},{ReplaceParams,
#define l1 cMul,SelectedParams,0},0,
#define l0 lA 0x0},{{
#ifdef _MSC_VER
typedef
unsigned
int
lC2;
#else
#include <stdint.h>
typedef
uint_least32_t
lC2;
#endif
lS3
crc32{enum{startvalue=0xFFFFFFFFUL,poly=0xEDB88320UL}
;tW1
lC2
crc>eC2
b8{enum{b1=(crc&xH1
crc
xV2
crc>>1),b2=(b1&xH1
b1
xV2
b1>>1),b3=(b2&xH1
b2
xV2
b2>>1),b4=(b3&xH1
b3
xV2
b3>>1),b5=(b4&xH1
b4
xV2
b4>>1),b6=(b5&xH1
b5
xV2
b5>>1),b7=(b6&xH1
b6
xV2
b6>>1),res=(b7&xH1
b7
xV2
b7>>1)}
;}
;inline
lC2
update(lC2
crc,unsigned
b){
#define B4(n) b8<n tG2 n+1 tG2 n+2 tG2 n+3>::res
#define R(n) B4(n),B4(n+4),B4(n+8),B4(n+12)
static
const
lC2
table[256]={R(0x00),R(0x10),R(0x20),R(0x30),R(0x40),R(0x50),R(0x60),R(0x70),R(0x80),R(0x90),R(0xA0),R(0xB0),R(0xC0),R(0xD0),R(0xE0),R(0xF0)}
;
#undef R
#undef B4
return((crc>>8))^table[(crc^b)&0xFF];nU1
lC2
calc_upd(lC2
c,const
unsigned
char*buf,size_t
size){lC2
value=c;for
nD2
p=0;p<size;++p)value=update(value,buf[p])iH
value;nU1
lC2
calc
eS1
unsigned
char*buf,size_t
size)yK
calc_upd(startvalue,buf,size);}
}
#ifndef FPOptimizerAutoPtrHH
#define FPOptimizerAutoPtrHH
nS1
Ref>class
xQ{eT3
xQ():p(0){}
xQ(Ref*b):p(b){xO3}
xQ
eS1
xQ&b):p(b.p){xO3
nU1
Ref&eB1*(yZ1*p;nU1
Ref*eB1->(yZ1
p;}
xQ&eB1=(Ref*b){Set(b)iH*this;}
xQ&eB1=eS1
xQ&b){Set(b.p)iH*this;}
#ifdef __GXX_EXPERIMENTAL_CXX0X__
xQ(xQ&&b):p(b.p){b.p=0;}
xQ&eB1=(xQ&&b){if(p!=b.p){eB2;p=b.p;b.p=0;}
return*this;}
#endif
~xQ(){eB2
l74
UnsafeSetP(Ref*newp){p=newp
l74
swap(xQ<Ref>&b){Ref*tmp=p;p=b.p;b.p=tmp;}
private:inline
static
void
Have(Ref*p2);inline
void
eB2;inline
void
xO3
inline
void
Set(Ref*p2);private:Ref*p;}
;nG
eB2{if(!p)return;p->xP3-=1;if(!p->xP3)delete
p;}
nG
Have(Ref*p2){if(p2)++(p2->xP3);}
nG
Birth(){Have(p);}
nG
Set(Ref*p2){Have(p2);eB2;p=p2;}
#endif
#include <utility>
eC2
Compare2ndRev{nS1
T>inline
iV2()eS1
T&xQ3
T&b
yZ1
a
t93>b
t93;}
}
;eC2
Compare1st{nS1
c0
const
nT1
xQ3
nT1
b
yZ1
a.first<b.first;}
nS1
c0
const
nT1
a,T1
b
yZ1
a.first<b;}
nS1
c0
T1
xQ3
nT1
b
yZ1
a<b.first;}
}
;
#ifndef FPoptimizerHashHH
#define FPoptimizerHashHH
#ifdef _MSC_VER
typedef
unsigned
long
long
xC1;
#define FPHASH_CONST(x) x##ULL
#else
#include <stdint.h>
typedef
uint_fast64_t
xC1;
#define FPHASH_CONST(x) x##ULL
#endif
lS3
FUNCTIONPARSERTYPES{eC2
fphash_t{xC1
hash1,hash2;fphash_t():hash1(0),hash2(0){}
fphash_t
eS1
xC1&xQ3
xC1&b):hash1(a),hash2(b){}
iV2==eS1
fphash_t&yY1==eD2&&hash2==eE2
iV2!=eS1
fphash_t&yY1!=eD2||hash2!=eE2
iV2<eS1
fphash_t&yY1!=eD2?hash1<eD2:hash2<eE2}
;}
#endif
#ifndef FPOptimizer_CodeTreeHH
#define FPOptimizer_CodeTreeHH
#ifdef FP_SUPPORT_OPTIMIZER
#include <vector>
#include <utility>
lS3
l81{eC2
Grammar;}
lS3
xI1{tU1
class
ByteCodeSynth;}
t6{tU1
class
CodeTree
n01
eC2
xU2
n01
class
CodeTree{typedef
xQ<xU2
xJ>DataP;DataP
data;eT3
CodeTree();~CodeTree();eC2
OpcodeTag{}
;e8
nW2
o,OpcodeTag);eC2
FuncOpcodeTag{}
;e8
nW2
o,unsigned
f,FuncOpcodeTag);eC2
xR3{}
;e8
const
iW2
v,xR3);
#ifdef __GXX_EXPERIMENTAL_CXX0X__
e8
l03&&v,xR3);
#endif
eC2
VarTag{}
;e8
unsigned
varno,VarTag);eC2
CloneTag{}
;e8
lT3
b,CloneTag);void
GenerateFrom
eS1
tE3
FunctionParserBase
xJ::Data&data,bool
keep_powi=false);void
GenerateFrom
eS1
tE3
FunctionParserBase
xJ::Data&data,const
x3&xO2,bool
keep_powi=false);void
SynthesizeByteCode(std::vector<unsigned>&eC3,std::vector
xJ&immed,size_t&stacktop_max);void
SynthesizeByteCode(xI1
nD3&synth,bool
MustPopTemps=true)const
cU3
SynthCommonSubExpressions(xI1::yF1
const;void
SetParams
eS1
x3&xS3
SetParamsMove(x3&tV1
CodeTree
GetUniqueRef();
#ifdef __GXX_EXPERIMENTAL_CXX0X__
void
SetParams(x3&&tV1
#endif
void
SetParam
nD2
which,lT3
b);void
nE1
size_t
which,xM2
b);void
AddParam
eS1
xM2
param);void
yD1
xM2
param);void
AddParams
eS1
x3&xS3
AddParamsMove(x3&xS3
AddParamsMove(x3&lY2,size_t
lZ2);void
DelParam
nD2
index);void
DelParams();void
Become
eS1
xM2
b);inline
size_t
GetParamCount(yZ1
lD2).size();nU1
xM2
GetParam
nD2
n)yK
lD2)[n];nU1
lT3
GetParam
nD2
n
yZ1
lD2)[n];nU1
void
nB2
nW2
o)iI3
Opcode=o;nU1
nW2
GetOpcode()yJ
Opcode;nU1
nB
fphash_t
GetHash()yJ
Hash;nU1
const
x3&lD2
yZ1
xX;nU1
x3&lD2)yK
xX;nU1
size_t
GetDepth()yJ
Depth;nU1
const
iW2
GetImmed()yJ
Value;nU1
unsigned
GetVar()yJ
lE2
nU1
unsigned
GetFuncNo()yJ
lE2
nU1
bool
IsDefined(yZ1
GetOpcode()!=nB
cNop;nU1
bool
cM3
yZ1
GetOpcode()==nB
cImmed;nU1
bool
IsVar(yZ1
GetOpcode()==nB
lI3;nU1
unsigned
GetRefCount()yJ
xP3
l74
ReplaceWithImmed
eS1
iW2
i);void
Rehash(bool
constantfolding=true);void
Sort();inline
void
Mark_Incompletely_Hashed()iI3
Depth=0;nU1
bool
Is_Incompletely_Hashed()yJ
Depth==0;nU1
const
tN
GetOptimizedUsing()yJ
lG2;nU1
void
SetOptimizedUsing
eS1
tN
g)iI3
lG2=g;}
bool
RecreateInversionsAndNegations(bool
prefer_base2=false);void
FixIncompleteHashes();void
swap(xM2
b){data.swap(b.data);}
bool
IsIdenticalTo
eS1
xM2
b)const;void
lB1}
n01
eC2
xU2{int
xP3;nW2
Opcode;l03
Value
iZ3
lE2
eN
eE3;nB
fphash_t
Hash
cU3
Depth;const
tN
lG2;xU2();xU2
eS1
xU2&b);i4
nW2
o);i4
nW2
o,unsigned
f);i4
const
iW2
i);
#ifdef __GXX_EXPERIMENTAL_CXX0X__
i4
l03&&i);xU2(xU2&&b);
#endif
bool
IsIdenticalTo
eS1
xU2&b)const;void
Sort();void
Recalculate_Hash_NoRecursion();private:void
eB1=eS1
xU2&b);}
n01
yP
CodeTreeImmed
eS1
iW2
i)n31
i
y9
xR3());}
#ifdef __GXX_EXPERIMENTAL_CXX0X__
cC1
CodeTreeImmed
tQ2&&i)n31
std::move(i)y9
xR3());}
#endif
cC1
CodeTreeOp(nW2
opcode)n31
opcode
y9
OpcodeTag());}
cC1
CodeTreeFuncOp(nW2
t23
unsigned
f)n31
t23
f
y9
FuncOpcodeTag());}
cC1
CodeTreeVar
tL1
varno)n31
varno
y9
VarTag());}
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
t81
DumpHashes(x9)xU3
DumpTree(x9)xU3
DumpTreeWithIndent(x9,const
std
tV3&indent="\\"
);
#endif
}
#endif
#endif
#ifndef FPOptimizer_GrammarHH
#define FPOptimizer_GrammarHH
#include <iostream>
t6{tU1
class
CodeTree;}
lS3
l81{enum
ImmedConstraint_Value{ValueMask=0x07,Value_AnyNum=0x0,x92=0x1,Value_OddInt=0x2,i51=0x3,Value_NonInteger=0x4,tC1=0x5
l94
ImmedConstraint_Sign{SignMask=0x18,Sign_AnySign=0x00,nF1=0x08,tD1=0x10,Sign_NoIdea=0x18
l94
ImmedConstraint_Oneness{OnenessMask=0x60,Oneness_Any=0x00,Oneness_One=0x20,Oneness_NotOne=0x40
l94
ImmedConstraint_Constness{ConstnessMask=0x180,Constness_Any=0x00,i41=0x80,Constness_NotConst=0x100
l94
Modulo_Mode{Modulo_None=0,Modulo_Radians=1
l94
Situation_Flags{LogicalContextOnly=0x01,NotForIntegers=0x02,OnlyForIntegers=0x04,OnlyForComplex=0x08,NotForComplex=0x10
l94
x62{NumConstant,lI4,SubFunction
l94
ParamMatchingType{PositionalParams,SelectedParams,AnyParams,GroupFunction
l94
RuleType{ProduceNewTree,ReplaceParams}
;
#ifdef __GNUC__
# define PACKED_GRAMMAR_ATTRIBUTE __attribute__((packed))
#else
# define PACKED_GRAMMAR_ATTRIBUTE
#endif
typedef
std::pair<x62,const
void*>eF2
n01
eF2
ParamSpec_Extract
tL1
paramlist,lG1)n01
bool
ParamSpec_Compare
eS1
void*xQ3
void*b,x62
type)iZ3
ParamSpec_GetDepCode
eS1
eF2&b);eC2
ParamSpec_ParamHolder{lG1:8
iZ3
constraints:9
iZ3
depcode:15;e6
tU1
eC2
ParamSpec_NumConstant{l03
constvalue
iZ3
modulo;}
;eC2
l51{unsigned
param_count:2
iZ3
param_list:30;nW2
subfunc_opcode:8;ParamMatchingType
match_type:3
iZ3
n1:5;e6
eC2
ParamSpec_SubFunction{l51
data
iZ3
constraints:9
iZ3
depcode:7;e6
eC2
Rule{RuleType
ruletype:2
iZ3
situation_flags:5
iZ3
repl_param_count:2+9
iZ3
repl_param_list:30;l51
match_tree;e6
eC2
Grammar{unsigned
rule_count
iZ3
short
rule_list[999
cV
extern
const
Rule
grammar_rules[];}
t81
DumpParam
eS1
eF2&p,std::ostream&o=std::cout)xU3
DumpParams
tL1
paramlist,unsigned
count,std::ostream&o=std::cout);}
#endif
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif
#define CONSTANT_POS_INF HUGE_VAL
#define CONSTANT_NEG_INF (-HUGE_VAL)
lS3
FUNCTIONPARSERTYPES{tU1
inline
l03
fp_const_pihalf()yK
fp_const_pi
xJ()*l03(0.5);}
tU1
inline
l03
fp_const_twopi(t12
fp_const_pi
xJ());lU3
fp_const_twoe(t12
fp_const_e
xJ());lU3
fp_const_twoeinv(t12
fp_const_einv
xJ());lU3
fp_const_negativezero()yK-Epsilon
xJ::value;}
}
#ifdef FP_SUPPORT_OPTIMIZER
#include <vector>
#include <utility>
#include <iostream>
cE1{using
lS3
l81;using
t6;using
lS3
FUNCTIONPARSERTYPES
n01
class
MatchInfo{eT3
std::vector<std::pair<bool,eN> >lQ;eN
paramholder_matches;std::vector<unsigned>tP;eT3
MatchInfo():lQ(),paramholder_matches(),tP(){}
eT3
bool
SaveOrTestRestHolder
tL1
n1,l23&iM1){cW1{lQ.eU3
n1+1);lQ
tO=iM1
lV3
if(lQ[n1
t73==false){lQ
tO=iM1
lV3
l23&found=lQ[n1
t83;if(iM1.size()!=found.size()cR
for
nD2
eQ3<iM1
l13
a)if(!iM1[a]xE
found[a])cR
return
true
l74
SaveRestHolder
tL1
n1,eN&iM1){cW1
lQ.eU3
n1+1);lQ
tO.swap(iM1);}
bool
SaveOrTestParamHolder
tL1
x0,xN2
xT3){if(c9
size()<=x0){c9
xV3
x0+1);c9
eU3
x0);c9
push_back(xT3
iG
if(!cD
y12)){cD=xT3
lV3
return
xT3
xE
cD
l64
SaveMatchedParamIndex(lG1){tP.push_back(index);}
xN2
GetParamHolderValueIfFound
tL1
x0)const{static
const
CodeTree
xJ
dummytree;if(c9
size()<=x0)return
dummytree
iH
cD;}
xN2
GetParamHolderValue
tL1
x0
yZ1
cD;}
bool
HasRestHolder
tL1
n1
yZ1
lQ.size()>n1&&lQ[n1
t73==true;}
l23&GetRestHolderValues
tL1
n1)const{static
l23
empty_result;cW1
return
empty_result
iH
lQ[n1
t83;}
const
std::vector<unsigned>&GetMatchedParamIndexes(yZ1
tP
l74
swap(tB
b){lQ.swap(b.lQ);c9
swap(b.paramholder_matches);tP.swap(b.tP);}
tB
eB1=eS1
tB
b){lQ=b.lQ;paramholder_matches=b.paramholder_matches;tP=b.tP
iH*this;}
}
;class
e9;typedef
xQ<e9>e0;class
e9{eT3
int
xP3;eT3
e9():xP3(0){}
virtual~e9(){}
}
;eC2
n61{bool
found;e0
specs;n61(bool
f):found(f),specs(){}
n61(bool
f
tH2
s):found(f),specs(s){}
}
xU3
SynthesizeRule
eS1
tI2
n41
tree,lP4)n01
n61
TestParam
eS1
eF2&yO2
xN2
tree
tH2
start_at,lP4)n01
n61
TestParams(tC&nN,xN2
tree
tH2
start_at,lP4,bool
nO2
n01
bool
ApplyGrammar
eS1
Grammar&i32,FPoptimizer_CodeTree::n41
tree,bool
from_logical_context=false)xU3
ApplyGrammars(FPoptimizer_CodeTree::eQ
n01
bool
IsLogisticallyPlausibleParamsMatch(tC&cC2,const
eQ;}
lS3
l81{t81
DumpMatch
eS1
tI2
nO,const
FPoptimizer_Optimize::lP4,bool
DidMatch,std::ostream&o=std::cout)xU3
DumpMatch
eS1
tI2
nO,const
FPoptimizer_Optimize::lP4,const
char*tX3,std::ostream&o=std::cout);}
#endif
#include <string>
tJ2
l81::x62
yU1=false);tJ2
nW2
yU1=false);
#include <string>
#include <sstream>
#include <assert.h>
#include <iostream>
using
lS3
l81;using
lS3
FUNCTIONPARSERTYPES;tJ2
l81::x62
yU1){
#if 1
const
char*p=0;switch(opcode
lQ4
lH4
p="NumConstant"
;lC
lI4:p="ParamHolder"
;lC
SubFunction:p="SubFunction"
c63}
std::ostringstream
tmp;assert(p);tmp<<p;if(pad)while(tmp.str().size()<12)tmp<<' '
iH
tmp.str();
#else
std::ostringstream
tmp;tmp<<opcode;if(pad)while(tmp.str().size()<5)tmp<<' '
iH
tmp.str();
#endif
}
tJ2
nW2
yU1){
#if 1
const
char*p=0;switch(opcode
lQ4
cAbs:p="cAbs"
;lC
cAcos:p="cAcos"
;lC
cAcosh:p="cAcosh"
;lC
cArg:p="cArg"
;lC
cAsin:p="cAsin"
;lC
cAsinh:p="cAsinh"
;lC
cAtan:p="cAtan"
;lC
cAtan2:p="cAtan2"
;lC
cAtanh:p="cAtanh"
;lC
cCbrt:p="cCbrt"
;lC
cCeil:p="cCeil"
;lC
cConj:p="cConj"
;lC
cCos:p="cCos"
;lC
cCosh:p="cCosh"
;lC
cCot:p="cCot"
;lC
cCsc:p="cCsc"
;lC
cExp:p="cExp"
;lC
cExp2:p="cExp2"
;lC
cFloor:p="cFloor"
;lC
cHypot:p="cHypot"
;lC
cIf:p="cIf"
;lC
cImag:p="cImag"
;lC
cInt:p="cInt"
;lC
cLog:p="cLog"
;lC
cLog2:p="cLog2"
;lC
cLog10:p="cLog10"
;lC
cMax:p="cMax"
;lC
cMin:p="cMin"
;lC
cPolar:p="cPolar"
;lC
cPow:p="cPow"
;lC
cReal:p="cReal"
;lC
cSec:p="cSec"
;lC
cSin:p="cSin"
;lC
cSinh:p="cSinh"
;lC
cSqrt:p="cSqrt"
;lC
cTan:p="cTan"
;lC
cTanh:p="cTanh"
;lC
cTrunc:p="cTrunc"
;lC
cImmed:p="cImmed"
;lC
cJump:p="cJump"
;lC
cNeg:p="cNeg"
;lC
cAdd:p="cAdd"
;lC
cSub:p="cSub"
;lC
cMul:p="cMul"
;lC
cDiv:p="cDiv"
;lC
cMod:p="cMod"
;lC
cEqual:p="cEqual"
;lC
i91:p="cNEqual"
;lC
cLess:p="cLess"
;lC
cLessOrEq:p="cLessOrEq"
;lC
cGreater:p="cGreater"
;lC
cGreaterOrEq:p="cGreaterOrEq"
;lC
cNot:p="cNot"
;lC
cAnd:p="cAnd"
;lC
cOr:p="cOr"
;lC
cDeg:p="cDeg"
;lC
cRad:p="cRad"
;lC
cFCall:p="cFCall"
;lC
cPCall:p="cPCall"
c63
#ifdef FP_SUPPORT_OPTIMIZER
case
cFetch:p="cFetch"
;lC
cPopNMov:p="cPopNMov"
;lC
e93:p="cLog2by"
;lC
cNop:p="cNop"
c63
#endif
case
cSinCos:p="cSinCos"
;lC
cSinhCosh:p="cSinhCosh"
;lC
cN3:p="cAbsNot"
;lC
cAbsNotNot:p="cAbsNotNot"
;lC
cAbsAnd:p="cAbsAnd"
;lC
cAbsOr:p="cAbsOr"
;lC
cAbsIf:p="cAbsIf"
;lC
cDup:p="cDup"
;lC
cInv:p="cInv"
;lC
cSqr:p="cSqr"
;lC
cRDiv:p="cRDiv"
;lC
cRSub:p="cRSub"
;lC
cNotNot:p="cNotNot"
;lC
cRSqrt:p="cRSqrt"
;lC
lI3:p="VarBegin"
c63}
std::ostringstream
tmp;assert(p);tmp<<p;if(pad)while(tmp.str().size()<12)tmp<<' '
iH
tmp.str();
#else
std::ostringstream
tmp;tmp<<opcode;if(pad)while(tmp.str().size()<5)tmp<<' '
iH
tmp.str();
#endif
}
#ifdef FP_SUPPORT_OPTIMIZER
#include <vector>
#include <utility>
#ifndef FP_GENERATING_POWI_TABLE
enum{MAX_POWI_BYTECODE_LENGTH=20}
;
#else
enum{MAX_POWI_BYTECODE_LENGTH=999}
;
#endif
enum{MAX_MULI_BYTECODE_LENGTH=3}
;lS3
xI1{tU1
class
ByteCodeSynth{eT3
ByteCodeSynth():ByteCode(),Immed(),cI(),xR(0),StackMax(0){iB
xV3
64);Immed.xV3
8);cI.xV3
16
l64
Pull(std::vector<unsigned>&bc,std::vector
xJ&imm,size_t&StackTop_max){for
tL1
eQ3<xP2;++a){lJ2]&=~nS2;}
iB
swap(bc);Immed.swap(imm);StackTop_max=StackMax;}
size_t
GetByteCodeSize(yZ1
xP2;}
size_t
GetStackTop(yZ1
xR
l74
PushVar
tL1
varno){xQ2
varno);tM2}
void
PushImmed
tQ2
immed
nL
xQ2
cImmed);Immed.push_back(immed);tM2}
void
StackTopIs(nO,int
offset=0){if((int)xR>offset){cI
lK2
first=true;cI
lK2
second=tree;}
}
bool
IsStackTop(nO,int
offset=0
yZ1(int)xR>offset&&cI
lK2
first&&cI
lK2
second
xE
tree);nU1
void
EatNParams
tL1
eat_count){xR-=eat_count
l74
ProducedNParams
tL1
produce_count){xJ1
xR+produce_count
l64
DoPopNMov
nD2
eG2,size_t
srcpos
nL
xQ2
cPopNMov)lU4
eG2)lU4
srcpos);xJ1
srcpos+1);cI[eG2]=cI[srcpos];xJ1
eG2+1
l64
DoDup
nD2
xW3
nL
if(xW3==xR-1){xQ2
cDup);}
else{xQ2
cFetch)lU4
xW3);}
tM2
cI[xR-1]=cI[xW3];}
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
tW1
int>void
Dump(tS3
ostream&o=std::cout;o<<"Stack state now("
<<xR<<"):\n"
nA3
0
eP3<xR;++a){o<<a<<": "
;if(cI[a
t73){nO=cI[a
t83;o<<'['<<std::hex<<(void*)(&tree.lD2))<<std::dec<<','<<tree.GetRefCount()<<']'
cL3
iI,o);}
else
o<<"?"
;o<<"\n"
;}
o<<std::flush;}
#endif
size_t
t13(nO)const{for
nD2
a=xR
eP3-->0;)if(cI[a
t73&&cI[a
t83
xE
tree
yL3
a
iH
tN2;}
bool
Find(nO
yZ1
t13
iI)!=tN2;}
bool
FindAndDup(nO){size_t
pos=t13
iI
cQ3
pos!=tN2){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<l44"duplicate at ["
<<pos<<"]: "
cL3
iI)xZ1" -- issuing cDup or cFetch\n"
;
#endif
DoDup(pos
iG
return
tD3
eC2
IfData{size_t
ofs;}
;void
SynthIfStep1
xR2,nW2
op
x31
xW2=xP2;xQ2
op
iC
xQ2
nS2
l64
SynthIfStep2
xR2
x31
iX
lH2+2);iX
2
lI2
xW2=xP2;xQ2
cJump
iC
xQ2
nS2
l64
SynthIfStep3
xR2
x31
iB
back()|=nS2;iX
lH2-1);iX
2
lI2
tM2
for
nD2
eQ3<xW2;++a){if(lJ2]==cJump&&lJ2+1]==(nS2|(xW2-1))){lJ2+lH2-1);lJ2+2
lI2
cJ3(lJ2]lQ4
cAbsIf:case
cIf:case
cJump:case
cPopNMov:a+=2;lC
cFCall:case
cPCall:case
cFetch:a+=1
c63
yI3
c73}
}
protected:void
xJ1
size_t
value){xR=value;if(xR>lW3{StackMax=xR;cI.eU3
lW3;}
}
protected:std::vector<unsigned>ByteCode;std::vector
xJ
Immed;std::vector<std::pair<bool,FPoptimizer_CodeTree::CodeTree
xJ> >cI
cU3
xR
cU3
StackMax;private:void
incStackPtr(){if(xR+2>lW3
cI.eU3
StackMax=xR+2);}
tW1
bool
IsIntType,bool
IsComplexType>eC2
cD2{}
;eT3
void
AddOperation
tL1
t23
unsigned
eat_count,unsigned
produce_count=1){EatNParams(eat_count);c1
opcode);ProducedNParams(produce_count
l64
c1
unsigned
t23
cD2<false,false>c6
false,true>c6
true,false>c6
true,true>);inline
void
c1
unsigned
opcode){c1
t23
cD2<bool(nB
IsIntType
xJ::x63),bool(nB
IsComplexType
xJ::x63)>());}
}
n01
eC2
SequenceOpCode
n01
eC2
i11{static
cE
AddSequence;static
cE
MulSequence;}
xU3
x11
long
count,cE&eP,yF1;}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
lS3
FUNCTIONPARSERTYPES;lS3
xI1{tU1
eC2
SequenceOpCode{l03
basevalue
iZ3
op_flip
iZ3
op_normal,op_normal_flip
iZ3
op_inverse,op_inverse_flip;}
n01
cE
i11
xJ::AddSequence={y71,cNeg
xD
cAdd,cSub,cRSub}
n01
cE
i11
xJ::MulSequence={l03(1),cInv,cMul,cMul,cDiv,cRDiv}
;
#define findName(a,b,c) "var"
#define TryCompilePowi(o) false
#define mData this
#define mByteCode ByteCode
#define mImmed Immed
nT2
false,false>){xK1
# define FP_FLOAT_VERSION 1
# define FP_COMPLEX_VERSION 0
# include "extrasrc/fp_opcode_add.inc"
# undef FP_COMPLEX_VERSION
# undef FP_FLOAT_VERSION
}
nT2
true,false>){xK1
# define FP_FLOAT_VERSION 0
# define FP_COMPLEX_VERSION 0
# include "extrasrc/fp_opcode_add.inc"
# undef FP_COMPLEX_VERSION
# undef FP_FLOAT_VERSION
}
#ifdef FP_SUPPORT_COMPLEX_NUMBERS
nT2
false,true>){xK1
# define FP_FLOAT_VERSION 1
# define FP_COMPLEX_VERSION 1
# include "extrasrc/fp_opcode_add.inc"
# undef FP_COMPLEX_VERSION
# undef FP_FLOAT_VERSION
}
nT2
true,true>){xK1
# define FP_FLOAT_VERSION 0
# define FP_COMPLEX_VERSION 1
# include "extrasrc/fp_opcode_add.inc"
# undef FP_COMPLEX_VERSION
# undef FP_FLOAT_VERSION
}
#endif
#undef findName
#undef mImmed
#undef mByteCode
#undef mData
#undef TryCompilePowi
}
using
lS3
xI1;
#define POWI_TABLE_SIZE 256
#define POWI_WINDOW_SIZE 3
lS3
xI1{
#ifndef FP_GENERATING_POWI_TABLE
extern
const
unsigned
char
powi_table[POWI_TABLE_SIZE];const
#endif
unsigned
char
powi_table[POWI_TABLE_SIZE]={0,1,1,1,2,1,2,1,xX3
4,1,2,xY3
2,1,xX3
8,eF3
xZ3
15,1,16,1,2,1,4,1,2,xY3
2,1,4,eF3
1,16,1,25,xZ3
27,5,8,3,2,1,30,1,31,3,32,1,2,1,xX3
8,1,2,xZ3
39,1,16,137,2,1,4,eF3
xY3
45,135,4,31,2,5,32,1,2,131,50,1,51,1,8,3,2,1,54,1,55,3,16,1,57,133,4,137,2,135,60,1,61,3,62,133,63,1,iN1
131,iN1
139,lL2
eA
30,1,130,137,2,31,lL2
eA
eA
130,eF3
1,eA
eA
2,1,130,133,iN1
61,130,133,62,139,130,137,eA
lL2
eA
eA
iN1
131,eA
eA
130,131,2,133,lL2
130,141,eA
130,eF3
1,eA
5,135,eA
lL2
eA
lL2
130,133,130,141,130,131,eA
eA
2,131}
;}
static
xG3
yQ=256;
#define FPO(x)
lS3{class
PowiCache{private:int
iE[yQ];int
iO1[yQ];eT3
PowiCache():iE(),iO1(){iE[1]=1;}
bool
Plan_Add(c01,int
count){cD1>=yQ
cR
iO1[tO2+=count
iH
iE[tO2!=0
l74
lX3
c01){cD1<yQ)iE[tO2=1
l74
Start
nD2
value1_pos){for(int
n=2;n<yQ;++n)iE[n]=-1;Remember(1,value1_pos);DumpContents();}
int
Find(c01)const{cD1<yQ){if(iE[tO2>=0){FPO(iO3(iS3,"* I found %ld from cache (%u,%d)\n",value,(unsigned)cache[value],iP3 value]))iH
iE[tO2;}
}
return-1
l74
Remember(c01,size_t
l34){cD1>=yQ)return;FPO(iO3(iS3,"* Remembering that %ld can be found at %u (%d uses remain)\n",value,(unsigned)l34,iP3 value]));iE[tO2=(int)l34
l74
DumpContents()const{FPO(for(int a=1;a<POWI_CACHE_SIZE;++a)if(cache[a]>=0||iP3 a]>0){iO3(iS3,"== cache: sp=%d, val=%d, needs=%d\n",cache[a],a,iP3 a]);})}
int
UseGetNeeded(c01){cD1>=0&&value<yQ)return--iO1[tO2
iH
0;}
}
n01
size_t
y7
long
count
iR1
cE&eP,yF1
xU3
c11
size_t
apos,long
aval,size_t
bpos,long
bval
iR1
unsigned
cumulation_opcode,unsigned
cimulation_opcode_flip,yF1;void
l91
c01
iR1
int
need_count,int
lA1=0){cD1<1)return;
#ifdef FP_GENERATING_POWI_TABLE
if(lA1>32)throw
false;
#endif
if(iE.Plan_Add(value,need_count
yL3;long
y03
1;cD1<POWI_TABLE_SIZE){y03
powi_table[tO2
yK3&128){half&=127
yK3&64)y03-iA2
FPO(iO3(iS3,"value=%ld, half=%ld, otherhalf=%ld\n",value,half,value/half));l91
half
c83
iE.lX3
half)iH;}
iQ1
half&64){y03-iA2}
}
else
cD1&1)y03
value&((1<<POWI_WINDOW_SIZE)-1);else
y03
value/2;long
cK=value-half
yK3>cK||half<0)std::swap(half,cK);FPO(iO3(iS3,"value=%ld, half=%ld, otherhalf=%ld\n",value,half,otherhalf))yK3==cK){l91
half,iE,2,lA1+1);}
else{l91
half
c83
l91
cK>0?cK:-cK
c83}
iE.lX3
value);}
tU1
size_t
y7
c01
iR1
cE&eP,yF1{int
y13=iE.Find(value
cQ3
y13>=0)yK
y13;}
long
y03
1;cD1<POWI_TABLE_SIZE){y03
powi_table[tO2
yK3&128){half&=127
yK3&64)y03-iA2
FPO(iO3(iS3,"* I want %ld, my plan is %ld * %ld\n",value,half,value/half))cU3
xX2=y7
half
t21
if(iE
lM2
half)>0||xX2!=t31){iP1
xX2)y02
half,t31);}
x11
value/half
nV2
size_t
l34=t31
y02
value,l34);iE.DumpContents()iH
l34;}
iQ1
half&64){y03-iA2}
}
else
cD1&1)y03
value&((1<<POWI_WINDOW_SIZE)-1);else
y03
value/2;long
cK=value-half
yK3>cK||half<0)std::swap(half,cK);FPO(iO3(iS3,"* I want %ld, my plan is %ld + %ld\n",value,half,value-half))yK3==cK){size_t
xX2=y7
half
t21
c11
xX2,half,xX2,half,iE,eP.op_normal,eP.op_normal_flip,synth);}
else{long
part1=half;long
part2=cK>0?cK:-cK
cU3
part1_pos=y7
part1
t21
size_t
part2_pos=y7
part2
t21
FPO(iO3(iS3,"Subdivide(%ld: %ld, %ld)\n",value,half,otherhalf));c11
part1_pos,part1,part2_pos,part2,iE,cK>0?eP.op_normal:eP.op_inverse,cK>0?eP.op_normal_flip:eP.op_inverse_flip,synth);}
size_t
l34=t31
y02
value,l34);iE.DumpContents()iH
l34;}
t81
c11
size_t
apos,long
aval,size_t
bpos,long
bval
iR1
unsigned
cumulation_opcode,unsigned
cumulation_opcode_flip,yF1{int
a_needed=iE
lM2
aval);int
y23=iE
lM2
bval);bool
lN2=false;
#define DUP_BOTH() do{if(apos<bpos){size_t tmp=apos;apos=bpos;bpos=tmp;lN2=!lN2;}FPO(iO3(iS3,"-> " l04 l04"op\n",(unsigned)apos,(unsigned)bpos));iP1 apos);iP1 apos==bpos?t31:bpos);}while(0)
#define DUP_ONE(p) do{FPO(iO3(iS3,"-> " l04"op\n",(unsigned)p));iP1 p);}while(0)
if(a_needed>0){if(y23>0){nU2}
e43
bpos!=t31)nU2
else{lO2
lN2=!lN2;}
}
}
iQ1
y23>0){if(apos!=t31)nU2
else
DUP_ONE(bpos);}
e43
apos==bpos&&apos==t31)lO2
tP2
t31&&bpos==synth.xL
2){FPO(iO3(iS3,"-> op\n"));lN2=!lN2;}
tP2
synth.xL
2&&bpos==t31)FPO(iO3(iS3,"-> op\n"));tP2
t31)DUP_ONE(bpos);iQ1
bpos==t31){lO2
lN2=!lN2;}
else
nU2}
lY3
lN2?cumulation_opcode_flip:cumulation_opcode,2);}
t81
cG1
long
count,cE&eP,yF1{while
eG3<256){int
y03
xI1::powi_table[count]yK3&128){half&=127;cG1
half
nV2
count/=half;}
else
c73
if
eG3==1)return;if(!eG3&1)){lY3
cSqr,1);cG1
count/2
nV2}
else{iP1
t31);cG1
count-1
nV2
lY3
cMul,2);}
}
}
lS3
xI1{t81
x11
long
count,cE&eP,yF1{if
eG3==0)nK2
eP.basevalue);else{bool
tR2=false;if
eG3<0){tR2=true;count=-count;}
if(false)cG1
count
nV2
iQ1
count>1){PowiCache
iE;l91
count,iE,1)cU3
xL1
tZ3
GetStackTop();iE.Start(t31);FPO(iO3(iS3,"Calculating result for %ld...\n",count))cU3
xY2=y7
count
t21
size_t
n_excess
tZ3
xL
xL1;if(n_excess>0||xY2!=xL1-1){synth.DoPopNMov(xL1-1,xY2);}
}
if(tR2)lY3
eP.op_flip,1);}
}
}
#endif
#ifndef FPOptimizer_ValueRangeHH
#define FPOptimizer_ValueRangeHH
t6{lS3
lZ3{iV
eC2
Comp{}
;tW1>eC2
Comp<nB
cLess>{tW1
lE<n9
cLessOrEq>{tW1
lE<=n9
cGreater>{tW1
lE>n9
cGreaterOrEq>{tW1
lE>=n9
cEqual>{tW1
lE==n9
i91>{tW1
lE!=b;}
}
;}
tU1
eC2
cH1{l03
val;bool
known;cH1():val(),known(false){}
cH1
eR1):val(v),known(true){nU1
y33
eR1){known=true;val=v;}
y33
tQ2(eC1
l03
nG3
nH3
y33
tQ2(eC1
const
iW2
nG3
nH3
iV
void
set_if
tQ2
v,l03(eC1
l03
nG3&&lZ3::Comp<Compare>()(val,v)nH3
iV
void
set_if
eR1,l03(eC1
const
iW2
nG3&&lZ3::Comp<Compare>()(val,v)nH3}
n01
eC2
range{cH1
xJ
min,max;range():min(),max(){}
range
tQ2
mi,l03
ma):min(mi),max(ma){}
range(bool,l03
ma):min(),max(ma){}
range
tQ2
mi,bool):min(mi),max(){}
void
set_abs();void
set_neg();}
n01
bool
IsLogicalTrueValue
eS1
yE3&p
tS2
n01
bool
IsLogicalFalseValue
eS1
yE3&p
tS2;}
#endif
#ifndef FPOptimizer_RangeEstimationHH
#define FPOptimizer_RangeEstimationHH
t6{enum
TriTruthValue{l53,y43,Unknown}
n01
yE3
CalculateResultBoundaries
eS1
eQ
n01
bool
IsLogicalValue
eS1
eQ
n01
TriTruthValue
GetIntegerInfo
eS1
eQ
n01
xM1
GetEvennessInfo
eS1
eQ{if(!tree
e62
yL3
Unknown;yO1=tree
eY1;if(nB
isEvenInteger(value
yL3
l53;if(nB
isOddInteger(value
yL3
y43
xZ2
tU1
xM1
GetPositivityInfo
eS1
eQ{yE3
p=CalculateResultBoundaries
iI
cQ3
p
cT&&p.yW>=l03(yL3
l53;if(p
cQ
l41
yL3
y43
xZ2
tU1
xM1
GetLogicalValue
n03
tree
tS2{yE3
p=CalculateResultBoundaries
iI
cQ3
IsLogicalTrueValue(p,abs
yL3
l53;if(IsLogicalFalseValue(p,abs
yL3
y43
xZ2}
#endif
#ifndef FPOptimizer_ConstantFoldingHH
#define FPOptimizer_ConstantFoldingHH
t6{t81
ConstantFolding(eQ;}
#endif
lS3{using
lS3
FUNCTIONPARSERTYPES;using
t6;eC2
ComparisonSetBase{enum{t33=0x1,Eq_Mask=0x2,Le_Mask=0x3,t43=0x4,t53=0x5,Ge_Mask=0x6}
;static
int
Swap_Mask(int
m)yK(m&Eq_Mask)|((m&t33)?t43:0)|((m&t43)?t33:0);}
enum
c21{Ok,BecomeZero,BecomeOne,nC1
l94
x72{cond_or,l63,l73,l83}
;}
n01
eC2
ComparisonSet:public
ComparisonSetBase{eC2
tT2
i01
a
iJ
b;int
relationship;tT2():a(),b(),relationship(){}
}
;std::vector<tT2>tH;eC2
Item
i01
value;bool
cE2;Item():value(),cE2(false){}
}
;std::vector<Item>cI1;int
xN1;ComparisonSet():tH(),cI1(),xN1(0){}
c21
AddItem
n03
a,bool
cE2,x72
type){for
nD2
c=0;c<cI1
l13
c)if(cI1[c].value
xE
a)){if(cE2!=cI1[c].cE2)iK
cX1
case
l83:cI1.erase(cI1.begin()+c);xN1
cF2
case
l63:case
l73:cY1}
}
return
nC1;}
Item
pole;pole.value=a;pole.cE2=cE2;cI1.push_back(pole)iH
Ok;}
c21
AddRelationship(CodeTree
xJ
a,CodeTree
xJ
b,int
i21,x72
type)iK
if(i21==7)cX1
lC
l83:if(i21==7){xN1
cF2}
lC
l63:case
l73:if(i21==0)cY1
c73
if(!(a.GetHash()<b.GetHash())){a.swap(b);i21=Swap_Mask(i21);}
for
nD2
c=0;c<tH
l13
c){if(tH[c].a
xE
a)&&tH[c].b
xE
b))iK{int
y53=tH[c
tY|i21;if(y53==7)cX1
tH[c
tY=y53
c63}
case
l63:case
l73:{int
y53=tH[c
tY&i21;if(y53==0)cY1
tH[c
tY=y53
c63}
case
l83:{int
newrel_or=tH[c
tY|i21;int
y22=tH[c
tY&i21;lQ2
5&&y22==0){tH[c
tY=t53
iH
nC1;}
lQ2
7&&y22==0){xN1+=1;tH.erase(tH.begin()+c)iH
nC1;}
lQ2
7&&y22==Eq_Mask){tH[c
tY=Eq_Mask;xN1
cF2}
tF2}
return
nC1;}
}
tT2
comp;comp.a=a;comp.b=b;comp.relationship=i21;tH.push_back(comp)iH
Ok;}
}
;nS1
l03,tE3
CondType>bool
ConstantFolding_LogicCommon(n41
tree,CondType
xD1,bool
y32){bool
should_regenerate=false;ComparisonSet
xJ
comp;lN1{tE3
c3
c21
eI3=c3
Ok;xN2
atree=nB3;switch(atree
nC
lQ4
cEqual
lG
Eq_Mask
eW3
i91
lG
t53
eW3
cLess
lG
t33
eW3
cLessOrEq
lG
Le_Mask
eW3
cGreater
lG
t43
eW3
cGreaterOrEq
lG
Ge_Mask
eW3
cNot:eI3
cZ1
l8
0),true
eW3
cNotNot:eI3
cZ1
l8
0),false,xD1)c63
yI3
if(y32||IsLogicalValue(atree))eI3
cZ1,false,xD1);cJ3(eI3){ReplaceTreeWithZero:xN
0)iH
true;ReplaceTreeWithOne:xN
1);nZ
c3
Ok:lC
c3
BecomeZero
e3
c3
BecomeOne:iA
c3
nC1:c61
c73}
if(should_regenerate){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before ConstantFolding_LogicCommon: "
nJ2
#endif
if(y32){tree.DelParams();}
else{for
cT1{xN2
atree=nB3;if(IsLogicalValue(atree))nN1);}
}
for
nD2
eQ3<comp.cI1
l13
a){if(comp.cI1[a].cE2){c31
cNot);r.c2
r
t63
iQ1!y32){c31
cNotNot);r.c2
r
t63
else
tree.c2}
for
nD2
eQ3<comp.tH
l13
a){c31
cNop);switch(comp.tH[a
tY
lQ4
c3
t33:r
xC
cLess
lZ4
Eq_Mask:r
xC
cEqual
lZ4
t43:r
xC
cGreater
lZ4
Le_Mask:r
xC
cLessOrEq
lZ4
t53:r
xC
i91
lZ4
Ge_Mask:r
xC
cGreaterOrEq
yS2
r.yD1
comp.tH[a].a);r.yD1
comp.tH[a].b);r
t63
if(comp.xN1!=0)tree.yB
l03(comp.xN1)));
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After ConstantFolding_LogicCommon: "
nJ2
#endif
return
true;}
return
tD3
xO1
ConstantFolding_AndLogic(iX3(tree.GetOpcode()==cAnd
lF4()==cAbsAnd)iH
nI
l63,true
lN4
ConstantFolding_OrLogic(iX3(tree.GetOpcode()==cOr
lF4()==cAbsOr)iH
nI
cond_or,true
lN4
ConstantFolding_AddLogicItems(iX3(tree.GetOpcode()==cAdd)iH
nI
l83,false
lN4
ConstantFolding_MulLogicItems(iX3(tree.GetOpcode()==cMul)iH
nI
l73,false);}
}
#include <vector>
#include <map>
#include <algorithm>
lS3{using
lS3
FUNCTIONPARSERTYPES;using
t6;eC2
CollectionSetBase{enum
xP1{Ok,nC1}
;}
n01
eC2
CollectionSet:public
CollectionSetBase{eC2
c41
i01
value
iJ
y42;bool
e1;c41():value(),y42(),e1(false){}
c41
n03
v,xN2
f):value(v),y42(f),e1(false){}
}
;std::multimap<fphash_t,c41>iF;typedef
tE3
std::multimap<fphash_t,c41>::yA3
xQ1;CollectionSet():iF(){}
xQ1
FindIdenticalValueTo
n03
value){fphash_t
hash=value.GetHash();for(xQ1
i=iF.y52
hash);i!=iF
c71
hash;++i){cD1
xE
i
eH2.value
yL3
i;}
return
iF
lM4;}
bool
Found
eS1
xQ1&b)yK
b!=iF
lM4;}
xP1
AddCollectionTo
n03
y42,const
xQ1&into_which){c41&c=into_which
eH2;if(c.e1)c.y42
cG
y42);else
i01
add;add
xC
cAdd);add.yD1
c.y42);add
cG
y42);c.y42.swap(add);c.e1=true;}
return
nC1;}
xP1
x82
n03
value,xN2
y42){const
fphash_t
hash=value.GetHash();xQ1
i=iF.y52
hash);for(;i!=iF
c71
hash;++i){if(i
eH2.value
xE
value
yL3
AddCollectionTo(y42,i);}
iF.y63,std::make_pair(hash,c41(value,y42)))iH
Ok;}
xP1
x82
n03
a)yK
x82(a,nD1
1)));}
}
n01
eC2
ConstantExponentCollection{typedef
eN
y83;typedef
std::x41
y62;std::vector<y62>data;ConstantExponentCollection():data(){}
void
MoveToSet_Unique
eS1
iW2
eD1&eE1){data.push_back(std::x41(eD1()));data.back()t93.swap(eE1
l64
MoveToSet_NonUnique
eS1
iW2
eD1&eE1){tE3
std::vector<y62>::yA3
i=std::y52
data.l93
data
lM4,exponent,Compare1st()cQ3
i!=data
c71
xT2{i
eH2.y63
eH2
lM4,eE1.l93
eE1
lM4);}
else{data.y63,std::x41
iY3,eE1));}
}
bool
iU2{bool
changed=false;std::sort(data.l93
data
lM4,Compare1st());redo:for
nD2
eQ3<data
l13
a
eZ2
exp_a=data[a
t73;if
lE3
exp_a,l03(1)))n33
for
nD2
b=a+1;b<data
l13
b
eZ2
exp_b=data[b
t73;l03
y72=exp_b-exp_a;if(y72>=fp_abs(exp_a))break;l03
exp_diff_still_probable_integer=y72*l03(16
cQ3
tU2
exp_diff_still_probable_integer)&&!(tU2
exp_b)&&!tU2
y72))){y83&a_set=lR2;y83&b_set=data[b
t83;
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before ConstantExponentCollection iteration:\n"
;tV2
cout);
#endif
if(isEvenInteger(exp_b)&&!isEvenInteger(y72+exp_a
e5
tmp2;tmp2
eV3
tmp2.SetParamsMove(b_set);tmp2
i23
tmp;tmp
xC
cAbs);tmp.yD1
tmp2);tmp
xJ2
b_set.eU3
1);b_set[0].iN2}
a_set.insert(a_set
lM4,b_set.l93
b_set
lM4);y83
b_copy=b_set;data.erase(data.begin()+b);MoveToSet_NonUnique(y72,b_copy);yQ1
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After ConstantExponentCollection iteration:\n"
;tV2
cout);
#endif
cJ}
}
}
return
changed;}
#ifdef DEBUG_SUBSTITUTIONS
void
tV2
ostream&out){for
nD2
eQ3<data
l13
a){out.precision(12);out<<data[a
t73<<": "
;eZ1
lR2
l13
b){if(b>0)out<<'*'
cL3(lR2[b],out);}
out<<std::endl;}
}
#endif
}
n01
static
CodeTree
xJ
x51
n41
value,bool&xG)tU3(value
nC
lQ4
cPow:i01
eJ2
value
l8
1);value.yH1
iH
exponent;}
case
cRSqrt:value.yH1;xG=true
iH
nD1-0.5))tH3
cInv:value.yH1;xG=true
iH
nD1-1));yI3
c73
return
nD1
1));}
cB1
void
eF1
eX1&mul,xN2
tree,xN2
y42,bool&c51
bool&xG){lN1
i01
value
iI
l8
a))iJ
exponent(x51
value,xG)cQ3!y42.cG2
y42
eY1!=l03(1.0
e5
e01;e01
eV3
e01
cG
xT2;e01
cG
y42);e01
xJ2
exponent.swap(e01);}
#if 0 /* FIXME: This does not work */
cD1
nC==cMul){if(1){bool
exponent_is_even=exponent
e62)&&isEvenInteger
iY3
eY1);eZ1
value.tA3{bool
tmp=false
iJ
val(value
l8
b))iJ
exp(x51
val,tmp)cQ3
exponent_is_even||(exp
e62)&&isEvenInteger(exp
eY1)e5
e01;e01
eV3
e01
cG
xT2;e01.yD1
exp);e01.ConstantFolding(cQ3!e01.cG2!isEvenInteger(e01
eY1)){goto
cannot_adopt_mul;}
}
}
}
eF1
mul,value,exponent,c51
xG);}
else
cannot_adopt_mul:
#endif
{if(mul.x82(value,xT2==CollectionSetBase::nC1)c61}
}
}
xO1
ConstantFolding_MulGrouping(eQ{bool
xG=false;bool
should_regenerate=false;eX1
mul;eF1
mul,tree,nD1
1)),c51
xG);typedef
std::pair<CodeTree
xJ,eN>eG1;typedef
std::multimap<fphash_t,eG1>cJ1;cJ1
tD;y82
eX1::xQ1
j=mul.iF.y93
j!=mul.iF
lM4;++j){n41
value=j
eH2.value
iJ&eJ2
j
eH2.y42;if(j
eH2.e1)exponent
xJ2
const
fphash_t
eH1=exponent.GetHash();tE3
cJ1::yA3
i=tD.y52
eH1);for(;i!=tD
c71
eH1;++i)if(i
eH2.first
xE
xT2){if(!exponent.cG2!e11
eY1,l03(1)))c61
i
eH2
t93.push_back(value);goto
skip_b;}
tD.y63,std::make_pair(eH1,std::make_pair
iY3,eN
nD2(1),value))));skip_b:;}
#ifdef FP_MUL_COMBINE_EXPONENTS
ConstantExponentCollection
xJ
e21;y82
cJ1::yA3
j,i=tD.y93
i!=tD
lM4;i=j){j=i;++j;eG1&list=i
eH2;if
eM2.lV1
eJ2
list.first
eY1;if(!iY3==xG1)e21.MoveToSet_Unique
iY3,list
iX2
tD.erase(i);}
}
if(e21.iU2)c61
#endif
if(should_regenerate)i01
before=tree;before.lB1
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before ConstantFolding_MulGrouping: "
cL3(before)xZ1"\n"
;
#endif
tree.DelParams();y82
cJ1::yA3
i=tD.y93
i!=tD
lM4;++i){eG1&list=i
eH2;
#ifndef FP_MUL_COMBINE_EXPONENTS
if
eM2.lV1
eJ2
list.first
eY1;if
iY3==xG1
n33
if(e11
nH2
tree.AddParamsMove(list
iX2
tF2}
#endif
CodeTree
xJ
mul;mul
eV3
mul.SetParamsMove(list
iX2
mul
xJ2
if(xG&&list.first.xG2
eM2
eY1==l03(1)/l03(3
e5
cbrt;cbrt
xC
cCbrt);cbrt.eB
cbrt
cH
cbrt
xV
0.5
e5
sqrt;sqrt
xC
cSqrt);sqrt.eB
sqrt
cH
sqrt
xV-0.5
e5
rsqrt;rsqrt
xC
cRSqrt);rsqrt.eB
rsqrt
cH
rsqrt
xV-1
e5
inv;inv
xC
cInv);inv.eB
inv
cH
inv);tF2}
CodeTree
xJ
pow
yR2.eB
pow.yD1
list.first);pow
cH
pow);}
#ifdef FP_MUL_COMBINE_EXPONENTS
tD.clear()nA3
0
eP3<i5
l13
a
eZ2
eJ2
i5[a
t73;if(e11
nH2
tree.AddParamsMove(i5[a]iX2
tF2
CodeTree
xJ
mul;mul
eV3
mul.SetParamsMove(i5[a]iX2
mul
i23
pow
yR2.eB
pow.yB
xT2);pow
cH
pow);}
#endif
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After ConstantFolding_MulGrouping: "
nJ2
#endif
return!tree
xE
before);}
return
tD3
xO1
ConstantFolding_AddGrouping(eQ{bool
should_regenerate=false;eX1
add;lN1{if
iI
l8
a)nC==cMul)n33
if(add.x82
iI
l8
a))==CollectionSetBase::nC1)c61}
nC3
remaining
iI
yS1)cU3
tE=0;lN1{xN2
xN3=nB3;if
e32
nC==cMul){eZ1
tB2
tA3{if
e32
l8
b)e62))n33
tE3
eX1::xQ1
c=add.FindIdenticalValueTo
e32
l8
b)cQ3
add.Found(c
e5
tmp
e32
y9
CloneTag());tmp.lR3
b);tmp
xJ2
add.AddCollectionTo(tmp,c);c61
goto
done_a;}
}
remaining[a]=true;tE+=1;done_a:;}
}
if(tE>0){if(tE>1
tS3
vector<std::pair<CodeTree
xJ,tB3>nY;std::multimap<fphash_t,tB3
eI1;bool
n53=false;lN1
t91{eZ1
nB3.tA3{xN2
p=nB3
l8
b);const
fphash_t
p_hash=p.GetHash();for(std::multimap<fphash_t,tB3::const_iterator
i=eI1.y52
p_hash);i!=eI1
c71
p_hash;++i){if(nY[i
eH2
t73
xE
p)){nY[i
eH2
t83+=1;n53=true;goto
found_mulgroup_item_dup;}
}
nY.push_back(std::make_pair(p,size_t(1)));eI1.insert(std::make_pair(p_hash,nY.size()-1));found_mulgroup_item_dup:;}
}
if(n53)i01
eO2;{size_t
max=0;for
nD2
p=0;p<nY
l13
p)if(nY[p
t83<=1)nY[p
t83=0;else{nY[p
t83*=nY[p
t73
nY2;if(nY[p
t83>max){eO2=nY[p
t73;max=nY[p
t83;}
}
}
CodeTree
xJ
group_add;group_add
xC
cAdd);
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Duplicate across some trees: "
cL3(eO2)xZ1" in "
nJ2
#endif
lN1
t91
eZ1
nB3.tA3
if(eO2
xE
nB3
l8
b)e5
tmp
iI
l8
a)y9
CloneTag());tmp.lR3
b);tmp
xJ2
group_add.yD1
tmp);remaining[a]=false
c63}
group_add
i23
group;group
eV3
group.yD1
eO2);group.yD1
group_add);group
xJ2
add.x82(group);c61}
}
lN1
t91{if(add.x82
iI
l8
a))==CollectionSetBase::nC1)c61}
}
if(should_regenerate){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before ConstantFolding_AddGrouping: "
nJ2
#endif
tree.DelParams();y82
eX1::xQ1
j=add.iF.y93
j!=add.iF
lM4;++j){n41
value=j
eH2.value
iJ&coeff=j
eH2.y42;if(j
eH2.e1)coeff
xJ2
if(coeff.xG2
lE3
coeff
eY1,xG1)n33
if
lE3
coeff
eY1
nH2
tree.yD1
value);tF2}
CodeTree
xJ
mul;mul
eV3
mul.yD1
value);mul.yD1
coeff);mul
cH
mul);}
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After ConstantFolding_AddGrouping: "
nJ2
#endif
return
true;}
return
tD3}
lS3{using
lS3
FUNCTIONPARSERTYPES;using
t6
n01
bool
ConstantFolding_IfOperations(iX3(tree.GetOpcode()==cIf
lF4()==cAbsIf);for(;;){cH2
nC==cNot){eX3
cIf);eK1
eJ1
0)cG3
xK2.swap
iI
l8
2));}
else
cH2
cR1{eX3
tC3;eK1
eJ1
0)cG3
xK2.swap
iI
l8
2));}
else
c73
cK3
0),tY2==tC3
iF1
tree
eJ1
1));nZ
lN3
tree
eJ1
2));nZ
lZ1
cH2
nC==cIf||eK1
nC==tC3
i01
cond=eK1
iJ
n63;n63
eY3==cIf?cNotNot:cAbsNotNot);n63
y92
1));ConstantFolding(n63)iJ
n73;n73
eY3==cIf?cNotNot:cAbsNotNot);n73
y92
2));ConstantFolding(n73
cQ3
n63.cG2
n73
e62
e5
eU;eU
eY3);eU
y92
1));eU.nJ
1));eU.nJ
2));eU
i23
eV;eV
eY3);eV
y92
2));eV.nJ
1));eV.nJ
2));eV
xJ2
eX3
cond
nC
tJ1
0,cond
l8
0)yB3
1,eU
yB3
2,eV
iG}
if
iI
l8
1)nC==eZ3)nC&&iI
l8
1)nC==cIf||xK2
nC==cAbsIf
e5&leaf1=xK2
iJ&leaf2=eZ3
cQ3
leaf1
l8
0)xF1
0))&&eI2
1))||leaf1
l8
2)xF1
2))e5
eU;eU
iT1
eU.nJ
0));eU
yA2
1));eU
yB2
1));eU
i23
eV;eV
iT1
eV.nJ
0));eV
yA2
2));eV
yB2
2));eV
eC
SetParam(0
cS3
0)yB3
1,eU
yB3
2,eV
iG
if
eI2
1))&&leaf1
l8
2)xF1
2)e5
eW
eN2
eW.AddParamMove
iI
cG3
eW
yA2
0));eW
yB2
0));eW
eC
nE1
0,eW
tJ1
2
cS3
2)tJ1
1
cS3
1)iG
if
eI2
2))&&leaf1
l8
2)xF1
1)e5
eP2;eP2
xC
leaf2
nC==cIf?cNot:cN3);eP2
yB2
0));eP2
i23
eW
eN2
eW.AddParamMove
iI
cG3
eW
yA2
0));eW.yD1
eP2);eW
eC
nE1
0,eW
tJ1
2
cS3
2)tJ1
1
cS3
1)iG}
n41
xW=xK2
iJ&y5=eZ3
cQ3
xW
xE
y5)){tree
eJ1
1)iG
const
OPCODE
op1=xW
nC;const
OPCODE
op2=y5
nC;yC3
op2){if(xW
yS1==1)i01
lO
0
i81
0));nB1
iY
n4
if(xW
yS1==2&&y5
yS1==2){if(xW
l8
0)xE
y5
l8
0)e5
param0=xW
l8
0)iJ
lO
1
i81
1));nB1
iY;tree.yD1
param0)n4
if(xW
l8
1)xE
y5
l8
1)e5
param1=xW
l8
1)iJ
lO
0
i81
0));nB1
iY;tree.yD1
nB1
lA4
param1
iG}
yC3
yD3
cMul
lS2
cAnd
lS2
cOr
lS2
cAbsAnd
lS2
cAbsOr
lS2
cMin
lS2
cMax){eN
n83;yR{for
nD2
b=y5.l71
b-->0;){if
lC4
y5
l8
b))){if(n83
i13){xW.lB1
y5.lB1}
n83.push_back(xW
l8
a));y5.lR3
b);xW
yC2}
if(!n83
i13){xW
xJ2
y5
yD2
tW2
SetParamsMove(n83)n4}
}
yC3
yD3
cMul||(op1==cAnd
cT3
y5))||(op1==cOr
cT3
y5))){yR
if
lC4
y5)){xW.lB1
xW
yT1
xW
i23
cK1=y5;y5=tF
op1==yD3
cOr)?0:1))tW2
yD1
cK1
lY4(op1==cAnd
lS2
cOr)&&op2==cNotNot){n41
nI3=y5
l8
0);yR
if
lC4
nI3)){xW.lB1
xW
yT1
xW
i23
cK1=nI3;y5=tF
op1
yG3
tW2
yD1
cK1
lY4
op2==cAdd||op2==cMul||(op2==cAnd
cT3
xW))||(op2==cOr
cT3
xW))){for
nD2
a=y5.cL1
y5
n93
xW)){y5.lB1
y5
yT1
y5
i23
cM1=xW;xW=tF
op2==cAdd||op2
yG3
l7
op2
lA4
cM1
lY4(op2==cAnd||op2==cOr)&&op1==cNotNot){n41
nJ3=xW
l8
0)nA3
y5.cL1
y5
n93
nJ3)){y5.lB1
y5
yT1
y5
i23
cM1=nJ3;xW=tF
op2
yG3
l7
op2
lA4
cM1)n4}
return
tD3}
#include <limits>
lS3{using
lS3
FUNCTIONPARSERTYPES;using
t6
n01
int
maxFPExponent()yK
std::numeric_limits
xJ::max_exponent;}
xO1
x61
l03
base,l03
xT2{if(base<xG1
return
true;if
lE3
base,xG1||nK3
base,l03(1))cR
return
exponent>=l03(maxFPExponent
xJ())/fp_log2(base
lN4
ConstantFolding_PowOperations(iX3(tree.GetOpcode()==cPow);nR&&xK2.lV1
const_value=tF3
lR,eR);xN
const_value)iH
tD3
xE3
nK3
eR
nH2
tree
eJ1
0)iG
nR&&nK3
lR
nH2
xN
1)iH
tD3
nR&&xK2
nC==cMul){bool
yE2=false;l03
lT2=lR
iJ
xN3=xK2
nA3
tB2
cL1
xN3
l8
a).lV1
imm=xN3
l8
a)eY1;{if(x61
lT2,imm))break;l03
lU2=tF3
lT2,imm);if
lE3
lU2,xG1)break;if(!yE2){yE2=true;tB2
lB1}
lT2=lU2;xN3
yC2
if(yE2){xN3
xJ2
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before pow-mul change: "
nJ2
#endif
eK1.Become(e31
lT2));xK2.Become
e32);
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After pow-mul change: "
nJ2
#endif
}
}
xE3
eK1
nC==cMul
eZ2
lV2=eR;l03
yF2=1.0;bool
yE2=false
iJ&xN3=eK1
nA3
tB2
cL1
xN3
l8
a).lV1
imm=xN3
l8
a)eY1;{if(x61
imm,lV2))break;l03
tA1=tF3
imm,lV2);if
lE3
tA1,xG1)break;if(!yE2){yE2=true;tB2
lB1}
yF2*=tA1;xN3
yC2
if(yE2){xN3
i23
eJ3;eJ3
xC
cPow);eJ3.SetParamsMove
iI.lD2));eJ3.yQ2);eX3
cMul
lA4
eJ3);tree
cG
e31
yF2)iG}
cH2
nC==cPow&&tree
eQ1
eK1
l8
1).lV1
a=eK1
l8
1)eY1;l03
b=eR;l03
c=a*b;if(isEvenInteger(a)&&!isEvenInteger(c
e5
nL3;nL3
xC
cAbs);nL3.nJ
0)cG3
nL3.Rehash(yB3
0,nL3);}
else
tree.SetParam(0,eK1
l8
0)tJ1
1,e31
c));}
return
tD3}
lS3{using
lS3
FUNCTIONPARSERTYPES;using
t6;eC2
l5{enum
eQ2{MakeFalse=0,MakeTrue=1,tX2=2,nN3=3,MakeNotNotP0=4,MakeNotNotP1=5,MakeNotP0=6,MakeNotP1=7,xA=8
l94
lW2{Never=0,Eq0=1,Eq1=2,yP3=3,e13=4}
;eQ2
if_identical;eQ2
lX2
4];eC2{eQ2
what:4;lW2
when:4;}
iU1,iV1,iW1,iX1
n01
eQ2
Analyze
n03
a,xN2
b)const{if(a
xE
b
yL3
if_identical;yE3
tT
a);yE3
p1=CalculateResultBoundaries(b
cQ3
p0
cQ&&p1
cT){if(p0.yG2<p1.yW&&lX2
0]i6
0];if(p0.yG2<=p1.yW&&lX2
1]i6
1];}
if(p0
cT&&p1
cQ){if(p0.yW>p1.yG2&&lX2
2]i6
2];if(p0.yW>=p1.yG2&&lX2
3]i6
3];}
if(IsLogicalValue(a)){if(iU1
lB3
iU1.when,p1
yL3
iU1.what;if(iW1
lB3
iW1.when,p1
yL3
iW1.what;}
if(IsLogicalValue(b)){if(iV1
lB3
iV1.when,p0
yL3
iV1.what;if(iX1
lB3
iX1.when,p0
yL3
iX1.what;}
return
xA;}
cB1
bool
TestCase(lW2
when,const
yE3&p){if(!p
cT||!p
cQ
cR
switch(when
lQ4
Eq0:yM3==l03(0.0)lV4==p.yW
tH3
Eq1:yM3==l03(1.0)lV4==p.yG2
tH3
yP3:yM3>y71
lV4<=l03(1)tH3
e13:yM3>=y71
l41
1);yI3;}
return
tD3}
;lS3
RangeComparisonsData{static
const
l5
Data[6]={{l5
cA1,{l5::tZ
xA,l5::tZ
xA
iZ
Eq1
l01
Eq1
nW1
Eq0
nX1
Eq0}
tB1
cA1,l5::xA,l5
cA1,l5::xA
iZ
Eq0
l01
Eq0
nW1
Eq1
nX1
Eq1}
tB1
cA1,l5::tX2,l5::tZ
MakeFalse
nW1
yP3
l01
e13
iY1
yC
nM3{l5::xA,l5
cA1,l5::tZ
nN3
nW1
e13
l01
yP3
iY1
tB1::tZ
tZ
nM3
l5::tX2
iZ
e13
nX1
yP3
iY1
yC
nM3{l5::tZ
nN3,l5::xA,l5
cA1
iZ
yP3
nX1
e13
iY1}
}
;}
xO1
ConstantFolding_Comparison(eQ{using
lS3
RangeComparisonsData;assert(tree.GetOpcode()>=cEqual&&tree.GetOpcode()<=cGreaterOrEq);switch(Data[tY2-cEqual].Analyze
iI
l8
0),xK2)lQ4
l5::MakeFalse:xN
0);nZ
l5
cA1:xN
1
nM2
nN3:eX3
cEqual
nM2
tX2:eX3
i91
nM2
MakeNotNotP0:eX3
cNotNot
y31
1
nM2
MakeNotNotP1:eX3
cNotNot
y31
0
nM2
MakeNotP0:eX3
cNot
y31
1
nM2
MakeNotP1:eX3
cNot
y31
0
nM2
xA:;}
if
iI
l8
e52))switch
iI
l8
0)nC
lQ4
cAsin:lM
fp_sin
t22
cAcos:lM
fp_cos(eR)));eX3
tY2==cLess?cGreater:tY2==cLessOrEq?cGreaterOrEq:tY2==cGreater?cLess:tY2==cGreaterOrEq?cLessOrEq:tY2);nZ
cAtan:lM
fp_tan
t22
cLog:lM
fp_exp
t22
cSinh:lM
fp_asinh
t22
cTanh:if(fp_less(fp_abs(eR)nH2
lM
fp_atanh(eR))iG
break;yI3
c73
return
tD3}
#include <list>
#include <algorithm>
#ifdef FP_SUPPORT_OPTIMIZER
using
lS3
FUNCTIONPARSERTYPES;lS3{
#ifdef DEBUG_SUBSTITUTIONS
yD
double
d){union{double
d;uint_least64_t
h;t32
d=d;lP1
h
nY1
#ifdef FP_SUPPORT_FLOAT_TYPE
yD
float
f){union{float
f;uint_least32_t
h;t32
f=f;lP1
h
nY1
#endif
#ifdef FP_SUPPORT_LONG_DOUBLE_TYPE
yD
long
double
ld){union{long
double
ld;eC2{uint_least64_t
a
iZ3
short
b;}
s;t32
ld=ld;lP1
s.b<<data.s.a
nY1
#endif
#ifdef FP_SUPPORT_LONG_INT_TYPE
yD
long
ld){o<<"("
<<std::hex<<ld
nY1
#endif
#endif
}
t6{lN
nE)){}
lN
const
iW2
i
y9
xR3
nE
i
lG3
#ifdef __GXX_EXPERIMENTAL_CXX0X__
lN
l03&&i
y9
xR3
nE
std::move(i)lG3
#endif
lN
unsigned
v
y9
VarTag
nE
lI3,v
lG3
lN
nW2
o
y9
OpcodeTag
nE
o
lG3
lN
nW2
o,unsigned
f
y9
FuncOpcodeTag
nE
o,f
lG3
lN
xN2
b
y9
CloneTag
nE*b.data)){}
tU1
CodeTree
xJ::~CodeTree(){}
lB
ReplaceWithImmed
eS1
iW2
i){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Replacing "
cL3(*this
cQ3
cM3))OutFloatHex(std::cout,GetImmed())xZ1" with const value "
<<i;OutFloatHex(std::cout,i)xZ1"\n"
;
#endif
data=new
xU2
xJ(i);}
tU1
eC2
ParamComparer{iV2()n03
a,xN2
b)const{if(a
nY2!=b
nY2)return
a
nY2<b
nY2
iH
a.GetHash()<b.GetHash();}
}
xU3
xU2
xJ::Sort()tU3(Opcode
lQ4
cAdd:case
cMul:case
cMin:case
cMax:case
cAnd:case
cAbsAnd:case
cOr:case
cAbsOr:case
cHypot:case
cEqual:case
i91:std::sort(lH3
l93
lH3
end(),ParamComparer
xJ());lC
cLess
lY
cGreater;}
lC
cLessOrEq
lY
cGreaterOrEq;}
lC
cGreater
lY
cLess;}
lC
cGreaterOrEq
lY
cLessOrEq;}
break;yI3
c73}
lB
AddParam
n03
param){xX.push_back(param);}
lB
yD1
n41
param){xX.push_back(CodeTree
xJ());xX.back().swap(param);}
lB
SetParam
nD2
which,xN2
b)nZ1
which
lF3
xX[which]=b;}
lB
nE1
size_t
which,n41
b)nZ1
which
lF3
xX[which
i03
b);}
lB
AddParams
eS1
nK){xX.insert(xX
lM4,lY2.l93
lY2
lM4);}
lB
AddParamsMove(nK){size_t
endpos=xX.size(),added=lY2.size();xX.eU3
endpos+added,CodeTree
xJ());for
nD2
p=0;p<added;++p)xX[endpos+p
i03
lY2[p]);}
lB
AddParamsMove(nK,size_t
lZ2)nZ1
lZ2
lF3
lR3
lZ2);AddParamsMove(tV1}
lB
SetParams
eS1
nK){eN
tmp(tV1
xX.iN2}
lB
SetParamsMove(nK){xX.swap(tV1
lY2.clear();}
#ifdef __GXX_EXPERIMENTAL_CXX0X__
lB
SetParams(eN&&lY2){SetParamsMove(tV1}
#endif
lB
DelParam
nD2
index){eN&eE3=xX;
#ifdef __GXX_EXPERIMENTAL_CXX0X__
lH3
erase(lH3
begin()+index);
#else
eE3[index].data=0;for
nD2
p=index;p+1<eE3
l13
p)eE3[p].data.UnsafeSetP(&*eE3[p+1
lF3
eE3[n23)-1].data.UnsafeSetP(0);lH3
eU3
n23)-1);
#endif
}
lB
DelParams(){xX.clear(lN4
CodeTree
xJ::IsIdenticalTo
n03
b)const{if(&*data==&*b.data)return
true
iH
data->IsIdenticalTo(*b.data
lN4
xU2
xJ::IsIdenticalTo
eS1
xU2
xJ&b)const{if(Hash!=b.Hash
cR
if(Opcode!=lX4
cR
switch(Opcode
l43
return
nK3
Value,lO4
tH3
lI3:return
lF2==b.lE2
case
cFCall:case
cPCall:if(lF2!=b.lF2
cR
break;yI3
c73
if(n23)!=b.n23)cR
for
nD2
eQ3<eE3
l13
a){if(!eE3[a]xE
b.eE3[a])cR}
return
true;}
lB
Become
n03
b){if(&b!=this&&&*data!=&*b.data){DataP
tmp=b.data;lB1
data.iN2}
}
lB
CopyOnWrite(){if(GetRefCount()>1)data=new
xU2
xJ(*data);}
tU1
CodeTree
xJ
CodeTree
xJ::GetUniqueRef(){if(GetRefCount()>1)return
CodeTree
xJ(*this,CloneTag())iH*this;}
i0):yL
cNop
lD4(),n8
i0
const
xU2&b):yL
lX4
lD4(lO4,lF2(b.cO1,eE3(b.eE3),Hash(b.Hash),Depth(b.Depth),i31
b.lG2){}
i0
const
iW2
i):yL
cImmed
lD4(i),n8
#ifdef __GXX_EXPERIMENTAL_CXX0X__
i0
xU2
xJ&&b):yL
lX4
lD4
cI2
lO4),lF2(b.cO1,eE3
cI2
b.eE3)),Hash(b.Hash),Depth(b.Depth),i31
b.lG2){}
i0
l03&&i):yL
cImmed
lD4
cI2
i)),n8
#endif
i0
nW2
o):yL
o
lD4(),n8
i0
nW2
o,unsigned
f):yL
o
lD4(),lF2(f),eE3(),Hash(),Depth(1),i31
0){}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <iostream>
using
lS3
FUNCTIONPARSERTYPES;
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
lS3{t81
n02(nO,std
c4&done,std::ostream&o){lN1
n02
lL4
done,o);std::ostringstream
buf
cL3
iI,buf);done[tree.GetHash()].insert(buf.str());}
}
#endif
t6{
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
t81
DumpHashes(cA){std
c4
done;n02
iI,done,o);for(std
c4::const_iterator
i=done.y93
i!=done
lM4;++i){const
std::set<std
tV3>&flist=i
eH2;if(flist.size()!=1)o<<"ERROR - HASH COLLISION?\n"
;for(std::set<std
tV3>::const_iterator
j=flist.y93
j!=flist
lM4;++j){o<<'['<<std::hex<<i->first.hash1<<','<<i->first.hash2<<']'<<std::dec;o<<": "
<<*j<<"\n"
;}
}
}
t81
DumpTree(cA){const
char*iR3;switch
iI
nC
l43
o<<tree
eY1
tG3
lI3:o<<"Var"
<<iI.GetVar()-lI3)tG3
cAdd:iR3"+"
;lC
cMul:iR3"*"
;lC
cAnd:iR3"&"
;lC
cOr:iR3"|"
;lC
cPow:iR3"^"
c63
yI3
iR3;o<<lR4
iI
nC);cU2
cFCall||tY2==cPCall)o<<':'<<tree.GetFuncNo();}
o<<'(';if
iI
yS1<=1&&sep2[1])o<<(sep2+1)<<' ';lN1{if(a>0)o<<' '
cL3
lL4
o
cQ3
a+1<tree
yS1)o<<sep2;}
o<<')';}
t81
DumpTreeWithIndent(cA,const
std
tV3&indent){o<<'['<<std::hex<<(void*)(&tree.lD2))<<std::dec<<','<<tree.GetRefCount()<<']';o<<indent<<'_';switch
iI
nC
l43
o<<"cImmed "
<<tree
eY1;o<<'\n'
tG3
lI3:o<<"VarBegin "
<<iI.GetVar()-lI3);o<<'\n'
iH;yI3
o<<lR4
iI
nC);cU2
cFCall||tY2==cPCall)o<<':'<<tree.GetFuncNo();o<<'\n';}
lN1{std
tV3
ind=indent;for
nD2
p=0;p<ind.size();p+=2)if(ind[p]=='\\')ind[p]=' ';ind+=(a+1<tree
yS1)?" |"
:" \\"
;DumpTreeWithIndent
lL4
o,ind);}
o<<std::flush;}
#endif
}
#endif
using
lS3
l81;using
lS3
FUNCTIONPARSERTYPES;
#include <cctype>
lS3
l81{unsigned
ParamSpec_GetDepCode
eS1
eF2&b)tU3(b.first
lQ4
lI4:{cL*s=(cL*)b
t93
iH
s->depcode;}
case
SubFunction:{cM*s=(cM*)b
t93
iH
s->depcode;}
yI3
c73
return
0;}
t81
DumpParam
eS1
eF2&yO2
std::ostream&o){static
const
char
ParamHolderNames[][2]={"%"
,"&"
,"x"
,"y"
,"z"
,"a"
,"b"
,"c"
}
iZ3
yH2
0;switch(nQ3.first
lQ4
lH4{const
ParamSpec_NumConstant
xJ&param=*eS1
ParamSpec_NumConstant
xJ*eP1;using
lS3
FUNCTIONPARSERTYPES;o.precision(12);o<<param.constvalue
c63}
case
lI4:{cL&c33
cL*eP1;o<<ParamHolderNames[param.index];yH2
param.constraints
c63}
case
SubFunction:{cM&c33
cM*eP1;yH2
param.constraints;yE
GroupFunction){if
iG3
lC1==cNeg){o<<"-"
;n2}
iQ1
param.lC1==cInv){o<<"/"
;n2}
else{std
tV3
opcode=lR4((nW2)param.lC1).substr(1)nA3
0
eP3<opcode
l13
a)opcode[a]=(char)std::toupper(opcode[a]);o<<opcode<<"( "
;n2
o<<" )"
;}
}
else{o<<'('<<lR4((nW2)param.lC1)<<' ';yE
PositionalParams)o<<'[';yE
SelectedParams)o<<'{';n2
if
iG3
data.n1!=0)o<<" <"
<<iT3.n1<<'>';yE
PositionalParams)o<<"]"
;yE
SelectedParams)o<<"}"
;o<<')';}
c73
cJ3(ImmedConstraint_Value
c53
ValueMask)lQ4
ValueMask:lC
Value_AnyNum:lC
x92:o<<"@E"
;lC
Value_OddInt:o<<"@O"
;lC
i51:o<<"@I"
;lC
Value_NonInteger:o<<"@F"
;lC
tC1:o<<"@L"
c43
ImmedConstraint_Sign
c53
SignMask)lQ4
SignMask:lC
Sign_AnySign:lC
nF1:o<<"@P"
;lC
tD1:o<<"@N"
c43
ImmedConstraint_Oneness
c53
OnenessMask)lQ4
OnenessMask:lC
Oneness_Any:lC
Oneness_One:o<<"@1"
;lC
Oneness_NotOne:o<<"@M"
c43
ImmedConstraint_Constness
c53
ConstnessMask)lQ4
ConstnessMask:lC
i41:if(nQ3.first==lI4){cL&c33
cL*eP1;if
iG3
index<2)c73
o<<"@C"
;lC
Constness_NotConst:o<<"@V"
;lC
Oneness_Any:c73}
t81
DumpParams
tL1
paramlist,unsigned
count,std::ostream&o){for
tL1
eQ3<count;++a){if(a>0)o<<' ';const
eF2&param=x71
paramlist,a);DumpParam
xJ(param,o)iZ3
depcode=ParamSpec_GetDepCode(param
cQ3
depcode!=0)o<<"@D"
<<depcode;}
}
}
#include <algorithm>
using
lS3
l81;using
lS3
FUNCTIONPARSERTYPES;lS3{cL
plist_p[37]={{2,0,0x0}
nP
0,0x4}
nP
nF1,0x0}
nP
tD1|Constness_NotConst,0x0}
nP
Sign_NoIdea,0x0}
nP
tC1,0x0}
,{3,Sign_NoIdea,0x0}
,{3,0,0x0}
,{3,tC1,0x0}
,{3,0,0x8}
,{3,Value_OddInt,0x0}
,{3,Value_NonInteger,0x0}
,{3,x92,0x0}
,{3,nF1,0x0}
,{0,tD1|lV{0,lV{0,nF1|lV{0,x92|lV{0,i41,0x1}
,{0,i51|nF1|lV{0,i61
i41,0x1}
,{0,i61
lV{0,Oneness_One|lV{0,tC1|lV{1,lV{1,x92|lV{1,i61
lV{1,i51|lV{1,nF1|lV{1,tD1|lV{6,0,0x0}
,{4,0,0x0}
,{4,i51,0x0}
,{4,lV{4,0,0x16}
,{5,0,0x0}
,{5,lV}
n01
eC2
plist_n_container{static
const
ParamSpec_NumConstant
xJ
plist_n[20];}
n01
const
ParamSpec_NumConstant
xJ
plist_n_container
xJ::plist_n[20]={{l03(-2
i1-1
i1-0.5
i1-0.25
i1
0
tZ2
fp_const_deg_to_rad
eK3
fp_const_einv
eK3
fp_const_log10inv
xJ(i1
0.5
tZ2
fp_const_log2
xJ(i1
1
tZ2
fp_const_log2inv
xJ(i1
2
tZ2
fp_const_log10
eK3
fp_const_e
eK3
fp_const_rad_to_deg
eK3-fp_const_pihalf
xJ(),xR1{y71,xR1{fp_const_pihalf
xJ(),xR1{fp_const_pi
xJ(),xR1}
;cM
plist_s[517]={{{1,15,i02
398,i02
477,i02
15,cNeg,GroupFunction,0}
,i41,0x1}
,{{1,15,yI2
24,yI2
465,yI2
466,yI2
498,cInv,lT
2,327995
xD
l0
2,48276
xD
l6
260151
xD
l6
470171
xD
l6
169126
xD
l6
48418
xD
l6
1328
xD
l6
283962
xD
l6
169275
xD
l6
39202
xD
l6
283964
xD
l6
283973
xD
l6
476619
xD
l6
296998
xD
l6
47
xD
SelectedParams,0
i33
nH
161839
xD
l6
25036
xD
l6
35847
xD
l6
60440
xD
l6
30751
xD
l6
270599
xD
l6
60431
xD
l6
259119
xD
l6
183474
xD
l6
332066
xD
l6
7168
xD
l6
197632
xD
l6
291840
xD
l6
283648
xD
l6
238866
xD
l6
239902
xD
l6
31751
xD
l6
244743
xD
l6
384022
xD
SelectedParams,0
i33
nH
385262
xD
l6
386086
xD
l6
393254
xD
SelectedParams,0}
,0,0x5
nH
393254
xD
l6
386095
xD
l6
387312
xD
l6
18662
xD
l6
61670
xD
l6
387397
xD
l6
247855
xD
SelectedParams,0}
,0,0x1
nH
342063
xD
l6
297007
xD
l6
15820
xD
l6
393263
xD
l6
393263
xD
SelectedParams,0}
,0,0x5
nH
161847
xD
l6
258103
xD
l6
249073
xD
l6
249076
xD
iL
0,0
xD
nF
0,0
i71
1,45
xD
nF
1,53
xD
nF
1,54
xD
nF
1,55
xD
nF
1,56
xD
nF
1,26
xD
nF
1,259
eR2
0x16}
,{{1,253
xD
nF
1,272
i71
1,323
eR2
0x16
lE4
xD
nF
1,21
xD
nF
1,447
eR2
0x4}
,{{1,449
eR2
0x4
lE4
eR2
0x4
lE4
tI
2
i33}
,{{1,15
xD
nF
1,24
tI
2}
,0,0x0
nH
58392
i71
0,0
tI
1}
,nF1,0x0
nH
24591
l84
33807
l84
48143
l84
285720
l84
290840
l84
305152,l9
312400,l9
39202,l9
121894,l9
421926,l9
429094,l9
443430,l9
317834,l9
329098,l9
7633,l9
7706,l9
7730,l9
38,l9
50587,l9
406528,l9
24583,l9
31751,l9
405511,l9
321551,xS1
327713,l9
322596,l9
88361,l9
335174,l9
327050,l9
493606,l9
496678,l9
503846,l9
516134,l9
7217,l9
333875,l9
336896,l9
524326,l9
509952,l9
286727,l9
90127,l9
131087,l9
296976,tX1
324623,l1
0x14
nH
332815,l1
0x10
xI2
7340056,tX1
289092,l9
92176,xS1
337935,l1
0x0
xI2
7340060
tL3
7340176,l9
338959,l1
0x0
xI2
7340061,xS1
7206,l9
yR3
l9
357414,l9
368678,l9
370745,l1
0x7
xI2
7340177,l9
39277,tX1
426398
tL3
40272286,xS1
490910
tL3
40336798,xS1
50600,l9
426462,xS1
490974,xS1
370726,l1
0x6
nH
371750,l1
0x6
nH
428070,l1
0x0
xI2
40336862,xS1
38378,l9
50671,l1
0x0
xI2
47662080,l9
477184,l9
568320,l9
371727,l1
0x7
xI2
15779306,l9
370703,l1
0x7
nH
39277,l9
39279,l1
0x4
xI2
15779238,l9
39338,tX1
436262,l9
508966,l9
39409,tX1
296998,tX1
35847,l9
15,tX1
377894,l9
386063,l1
0x1
nH
15,l9
7192,l9
122904,l9
121880,l9
30751,l9
57,l9
7456,l9
15674,l1
0x0
xI2
67579935,l9
39237,l9
58768,l9
62924,l9
121856,l9
15760,l1
0x0
xI2
64009216,l1
0x0}
,{{0,0,xF
0,0,iS
2,e41
2,e61
3,e41
3,e61
38,xF
1,38,iS
14,xF
1,57,xF
1,16,eS2
0x0
nH
471103,eS2
0x1}
,{{1,303,xF
1,323,yQ3
0x0
nH
471363,eS2
0x16}
,{{1,293,e41
294,e61
295,xF
1,296,iS
400,xF
1,0,xF
1,460,xF
1,465,xF
1,16,eS2
0x1}
,{{1,57,yQ3
0x1
lE4,iS
21,xF
1,15,eS2
0x0
nH
24591,xF
1,24,iS
517,yQ3
0x0
nH
46095,lK
46104,lK
15397,lK
287789,lK
66584,lK
404763,lK
62504,lK
15409,lK
39951,lK
24591,lK
33807,lK
50200,lK
62509,lK
50176,lF,178176,tE1
0x12
nH
283648,lF,19456,lF,27648,lF,89088,lF,86016,lF,488448,lF,14342,lF,58375,lF,46147
xZ
46151,lF,284679,lF,7183,lF,46159
xZ
38993
xZ
50265,lF,50249,lF,283808,lF,284835,lF,24822,lF,10240,lF,11264,lF,7170,lF,yR3
lF,17408,lF,164864,lF,237568,lF,242688,tE1
0x14
nH
476160,lF,25607,lF,121871,lF,50252,lF,39374,lF,50183,lF,7192,lF,121887,lF,252979,lF,46155,lF,38919,lF,50267,lF,50268,lF,50253,lF,46190,lF,50295,lF,7563,tE1
0x10
nH
416811,lF,416819,lF,40046,lF,46191
xZ
415795,lF,40047
xZ
415787,lF,39015,tE1
0x5
nH
39326
xZ
39326,lF,39332,tE1
0x5
nH
39333,tE1
0x1
nH
50590
xZ
50590,lF,39338
xZ
39338,lF,39335,tE1
0x5
nH
15786
xZ
146858,lF,39372,lF,39379,lF,39380,lF,39390
xZ
50654
xZ
50654,lF,24,tE1
0x6
nH
62,lF,24,lF,62,tE1
0x6
nH
43,lF,43
xZ
51,lF,51
xZ
50269,lF,50176
xZ
50270,lF,39159,lF,39183
xZ
7168
xZ
31744,lF,99328,lF,31746,lF,100376,lF,39409
xZ
39411
xZ
39411,lF,39420,lF,39420
xZ
15,lF,39025,tE1
0x5
nH
39422,lF,16384,lF,62853,lF,15360,lF,15,tE1
0x1
nH
16,lF,7183,tE1
0x1
nH
7172
iZ2
yG1,nF1,0x0
nH
24591
iZ2
lT
2,50200
iZ2
lT
2,63521
iZ2
lT
2,62500
iZ2
lT
2,50453
iZ2
lT
2,62488
iZ2
lT
1,0,tI3
7,tI3
194,tI3
0,cAcos
tJ3
cAcosh
tJ3
cAsin
tJ3
cAsinh
nQ
119,cAsinh
tJ3
cAtan
t61
306176,cAtan2
t61
yR3
cAtan2
tJ3
cAtanh
nQ
246,cCeil
tJ3
cCeil,tK3
cJ2
0,cCos,iT
1,7,cJ2
91,cJ2
92,cJ2
119,cJ2
236,cJ2
255,cJ2
214,lJ3
236,lJ3
464,lJ3
0,cCosh,tK3
lJ3
0,cExp
nQ
7,cExp
nQ
91,cExp
tJ3
yU3
7,yU3
91,yU3
246,cFloor
tJ3
cFloor,lA
0x4
nH
309540,cHypot
t61
316708,cHypot
t61
316724,cHypot,l0
3,32513024,eU2
34627584
lT1
31493120,eU2
89213952
lT1
149042176
lT1
246647808
lT1
301234176
lT1
494360576
lT1
498558976
lT1
62933520
lT1
62933520,eU2
62933526
lT1
62933526,eU2
24670208
lT1
579378176
lT1
573578240
lT1
32513024
lT1
566254592
lT1
7900160
lT1
588822528,cIf
nQ
119,cInt
nQ
246,i12
0,i12
7,i12
31,i12
194,i12
363,i12
15,cLog,lT
1,24,cLog,lT
1,0,cLog10
tJ3
cLog2
t61
yR3
cMax
t61
35847,cMax
t61
30751,cMax
tJ3
cMax,AnyParams,1
i33
nH
yR3
cMin
t61
35847,cMin
t61
30751,cMin
tJ3
cMin,AnyParams,1
i33
nH
24591,cMin,lT
1,0,xA2
7,xA2
91,xA2
92,xA2
119,xA2
149,xA2
231,cSin,lA
0x5}
,{{1,246,xA2
255,xA2
254,xA2
0,cSin,iT
1,273,cSin,lA
0x1}
,{{1,214,yJ2
231,cSinh,lA
0x5}
,{{1,246,yJ2
254,yJ2
255,yJ2
464,yJ2
0,cSinh,tK3
yJ2
15,cSqrt,lT
1,0,cK2
0,cTan,iT
1,115,cTan,iT
1,116,cK2
231,cK2
246,cK2
273,cK2
254,cK2
255,cK2
0,yK2
0,cTanh,iT
1,213,yK2
231,yK2
246,yK2
254,yK2
255,yK2
0,cTrunc
t61
15384,cSub,lT
2,15384,cDiv,lT
2,476626,cDiv,lT
2,121913,cEqual
t61
yR3
cEqual
tL2
yR3
cEqual
t61
31744,i22
lA
0x20
nH
31751,i22
lA
0x24
nH
31751,cEqual
t61
121913,i91
t61
yR3
cLess
tL2
41984,cLess,lA
0x4
nH
41984,cLess
t61
7,cLess
t61
yR3
cLessOrEq
t61
296182,cLessOrEq
t61
7168
iY2
tL2
41984
iY2,lA
0x4
nH
41984
iY2
t61
7
iY2
t61
yR3
yA
l0
2,296182,cGreaterOrEq
tJ3
n12
245,n12
7,n12
550,n12
553,n12
554,n12
556,n12
31,n12
559,n12
15,n12
560,cNot
t61
7706
nF3
7168
nF3
35847
nF3
30751
nF3
463903
nF3
466975,cAnd,iL
0,0,cAnd,nF
2,yR3
cOr,l6
7706,cOr,l6
35847,cOr,l6
463903,cOr,l6
466975,cOr,l6
30751,cOr,iL
1,0,n22
91,n22
131,n22
245,n22
215,n22
246,cDeg
nQ
246,cRad
t61
yR3
cAbsAnd,l6
yR3
cAbsOr,iL
1,0,cN3
tJ3
cAbsNotNot,l0
3,32513024,c13
lA
0x0}
,}
;}
lS3
l81{const
Rule
grammar_rules[262]={{ProduceNewTree,17,1,0,{1,0,cAbs,eV2
409,{1,146,cAtan,eV2
403
nP
1324,cAtan2,eV2
405
nP
307201,cAtan2
eS
253174
nP
255224,cAtan2
eS
259324
nP
257274,cAtan2,eV2
152,{1,252,cCeil
iE1
486,{1,68,xT1
482,{1,122,xT1
483,{1,124,xT1
151,{1,125,xT1
419,{1,123,xT1
0,{1,403,cCos,l2
2,1,246,{1,252,cCos,l2
18,1,0,{1,400,xT1
301,{1,406,cCosh,l2
2,1,246,{1,252,cCosh,l2
18,1,0,{1,400,cCosh
iE1
458,{1,121,cFloor,eV2
150,{1,252,cFloor,tM3
156,{3,7382016,eI
549,{3,8430592,eI
556,{3,8436736,eI
157,{3,42998784,eI
550,{3,42999808,eI
562,{3,43039744,eI
557,{3,49291264,eI
538,{3,49325056,eI
469,{3,1058318,eI
473,{3,1058324,eI
473,{3,9438734,eI
469,{3,9438740,cIf,l2
0,3,32542225,{3,36732434,cIf,l2
0,3,32542231,{3,36732440,cIf,tO3
573,{3,32513026,cIf,tO3
515,{3,455505423,cIf,tO3
515,{3,433506837,cIf
iE1
78,{1,256,yV3
69,{1,258,yV3
404,{1,72,yV3
159,{1,147,cLog,l2
0,1,0
nP
487425,cMax
l3
16,1,445
nP
yW3
cMax
l3
0,1,0
nP
483329,cMin
l3
16,1,446
nP
yW3
cMin,yS
0,1,153
nP
24832
iZ2
tM3
153
nP
25854
iZ2
tM3
154
nP
129039
iZ2
xU1
32055
iZ2
xU1
32056
iZ2
xU1
32057
iZ2
l2
0,2,166288
nP
32137
iZ2
xU1
33082
iZ2
l2
0,2,7168
nP
12688
iZ2
l2
0,2,7434
nP
12553
e51
435
nP
46146
e51
436
nP
46154
e51
437
nP
46150
e51
169
nP
83983
e51
168
nP
130082
e51
175
nP
133154
iZ2
tP3
476160
nP
471055
iZ2
tP3
274432
nP
273423
iZ2
tP3
251904
nP
266274
iZ2
tP3
251904
nP
263186
e51
171,{1,252,lW1
421,{1,68,lW1
151,{1,122,lW1
419,{1,124,lW1
170,{1,125,lW1
482,{1,123,lW1
0,{1,405,lW1
172,{1,252,cSinh
iE1
328,{1,404,cSinh
iE1
173,{1,252,yX3
0,{1,408,yX3
176,{1,410,yX3
177,{1,252,cTanh,l2
0,1,442
nP
449551,i2
1,441
nP
yW3
i2
1,167
nP
268549,i2
1,181
nP
276749,i2
1,180
nP
276500
iF3
190770
nP
189622
iF3
194748
nP
193723
iF3
202943
nP
196795
iF3
59699
nP
298148
iF3
59714
nP
325815
iF3
59724
nP
343224
xD
yS
2,1,337,{1,333
tI
1
tX
336,{1,338
tI
1}
}
,{ReplaceParams,2,1,340
nP
1363
tV
342
nP
1365
tV
463
nP
472524
tV
47
nP
356711
tV
349
nP
200751
tV
360
nP
199727
tV
480
nP
207053
tV
481
nP
208077
tV
417
nP
211144
tV
209
nP
211145
tV
418
nP
215240
tV
212
nP
212329
tV
204
nP
373097
tV
211
nP
372944
tV
217
nP
201944
tV
221
nP
223448
tV
367
nP
508329
tV
219
nP
508126
tV
224
nP
225705
tV
223
nP
225776
tV
365
nP
230825
tV
426
nP
377057
tV
497
nP
377054
tV
497
nP
204201
tV
426
nP
375280
tV
224
nP
375006,cAdd
tQ3
407781
nP
233698,cAdd
tQ3
59763
nP
233842,i2
1,372
nP
1397,xI3
95
nP
24705,xI3
96
nP
24708,xI3
444
nP
449551,xI3
443
nP
yW3
xI3
100
nP
101750,xI3
108
nP
106821,xI3
105
nP
103749,nG1
2,110607
nP
108869,nG1
2,107535
nP
109893,lJ
0
tX
112
nP
111634,cMul,SelectedParams,0
tX
567,{1,52,lJ
1
tX
568,{1,42,lJ
1}
}
,{ReplaceParams,2,1,467
nP
45516
eT2
356
nP
51555
eT2
468
nP
49612
eT2
357
nP
47459
eT2
429
nP
438699
eT2
432
nP
441774
eT2
486
nP
498726
eT2
494
nP
504870
eT2
382
nP
435579
eT2
497
nP
435709
eT2
426
nP
508287
eT2
414
nP
500092
eT2
499
nP
352744
eT2
345
nP
367092
eT2
381
nP
425318
eT2
478
nP
425460
eT2
47
nP
512501
eT2
505
nP
355817
eT2
47
nP
516598
eT2
507
nP
518182
eT2
508
nP
358896
eT2
351
nP
388605
eT2
511
nP
360939
eT2
503
nP
354788
eT2
514
nP
525350
eT2
510
nP
394342
eT2
386
nP
351346,cMul
tQ3
363004
nP
361968
eO1
1,117
nP
1157
eO1
1,118
nP
1158
eO1
1,402
nP
411024
eO1
2,58768
nP
1472
eO1
2,15760
nP
1474,cMul
l3
17,1,0,{1,400,cMul
l3
17,1,57,{1,14,lJ
0}
}
,{ProduceNewTree,4,1,538
nP
41,i22
tR3
0
nP
5167,i22
cF
41984
nP
409641,i22
cF
tQ
i22
cF
eX
i22
cF
eY
cEqual
cN1
24849
yS3
tR
cEqual
eL3
281873
yS3
lD1
yS3
lE1
i22
l4
4
yT3
41,i91,tR3
538
nP
5167,i91,cF
41984
nP
409641,i91,cF
tQ
i91,cF
eX
i91,cF
eY
i91
cN1
24849,i91
eS
tR
i91
eL3
281873,i91
eS
lD1,i91
eS
lE1
i91,cF
tQ
yY3
eX
yY3
eY
cLess,eV2
571
nP
46080,cLess
cN1
24832
cV3
xV1
cLess
eS
tR
cLess
eL3
281856
cV3
x01
cLess
eS
lD1
cV3
lE1
cLess,l4
20
yT3
409641,yY3
tQ
n42
cF
eX
n42
cF
eY
n42
eV2
565
nP
409615,cLessOrEq
cN1
24832,cLessOrEq
eS
xV1
cLessOrEq
eS
tR
cLessOrEq
eL3
281856,cLessOrEq
eS
x01
cLessOrEq
eS
lD1,cLessOrEq
eS
lE1
n42
l4
20
yT3
409647,n42
cF
tQ
cL2
eX
cL2
eY
cGreater,eV2
539
nP
409615
iY2
cN1
24832
iY2
eS
xV1
cGreater
eS
tR
cGreater
eL3
281856
iY2
eS
x01
cGreater
eS
lD1
iY2
eS
lE1
cGreater,l4
20,1,538
nP
409647,cL2
tQ
yA
cF
eX
yA
cF
eY
yA
eV2
572
nP
46080
nE2
529654
nP
24832
nE2
xV1
yA
yZ3
tR
yA
yZ3
n32
281856
nE2
x01
yA
yZ3
lD1
nE2
lE1
yA
l4
20,1,538
nP
409641,yA
tR3
519,{1,137,cNot,tO3
571,{1,2,cNot,l2
0,1,452
nP
yW3
nO3
0,2,537097,{3,547892744,cAnd,yS
16,1,566,{1,5,cAnd,AnyParams,1}
}
,{ReplaceParams,16,1,569
nP
13314,nO3
16,1,544
nP
553498,nO3
16,1,546
nP
462369,nO3
16,1,548
nP
466465,nO3
0,1,457
nP
yW3
xW1
570
nP
13314,xW1
563
nP
8197,xW1
541
nP
553498,xW1
542
nP
462369,xW1
543
nP
466465,xW1
564
nP
143365,cOr,yS
4,1,525,{1,137,c03
tO3
572,{1,2,c03
l4
17,1,0,{1,0,c03
eV2
537,{1,256,cAbsNotNot,yS
18,1,531,{1,254,cAbsNotNot,yS
0,1,572,{3,43039744,c13
tM3
571,{3,49325056,c13
tO3
454,{3,32513586,c13
l2
16,3,32542225,{3,36732434,c13
yG1}
,}
;eC2
grammar_optimize_abslogical_type{y2
9
cV
grammar_optimize_abslogical_type
grammar_optimize_abslogical={9,{34,192,228,238,242,247,254,260,261}
}
;}
eC2
grammar_optimize_ignore_if_sideeffects_type{y2
59
cV
grammar_optimize_ignore_if_sideeffects_type
grammar_optimize_ignore_if_sideeffects={59,{0,20,21,22,23,24,25,26,cN
iZ1
78,cO
cW
eC2
grammar_optimize_nonshortcut_logical_evaluation_type{y2
56
cV
grammar_optimize_nonshortcut_logical_evaluation_type
grammar_optimize_nonshortcut_logical_evaluation={56,{0,25,cN
iZ1
78,cO
241,243,244,245,246,248,249,250,251,252,253,255,256,257,258,259}
}
;}
eC2
grammar_optimize_recreate_type{y2
22
cV
grammar_optimize_recreate_type
grammar_optimize_recreate={22,{18,55,56,57,80,81,82,83,84,85,117,118,120,121,130,131,132,133,134,135,136,137}
}
;}
eC2
grammar_optimize_round1_type{y2
125
cV
grammar_optimize_round1_type
grammar_optimize_round1={125,{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,19,25,cN
37,38,iZ1
45,46,47,48,49,50,51,52,53,54,58,59,60,61,62,63,64,65,66,67,68,69,70,71,78,79,80,81,82,83,84,85,86,87,88,93,94,95,96,97,98,99,100,101,117,118,119,120,121,122,123,124,125,126,127,128,129,138,160,161,162,163,164,165,166,167,168,169,178,179,180,200,204,212,216,224,236,237,239,240,cW
eC2
grammar_optimize_round2_type{y2
103
cV
grammar_optimize_round2_type
grammar_optimize_round2={103,{0,15,16,17,25,cN
39,40,iZ1
45,46,47,48,49,50,51,52,53,54,59,60,72,73,78,79,86,87,88,89,90,91,92,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,119,122,123,124,125,126,127,128,139,159,160,161,162,163,164,165,166,167,168,169,178,179,180,200,204,212,216,224,236,237,239,240,cW
eC2
grammar_optimize_round3_type{y2
79
cV
grammar_optimize_round3_type
grammar_optimize_round3={79,{74,75,76,77,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,170,171,172,173,174,175,176,177,181,182,183,184,185,186,187,188,189,190,191,193,194,195,196,197,198,199,201,202,203,205,206,207,208,209,210,211,213,214,215,217,218,219,220,221,222,223,225,226,227,229,230,231,232,233,234,235}
}
;}
eC2
grammar_optimize_round4_type{y2
12
cV
grammar_optimize_round4_type
grammar_optimize_round4={12,{18,55,56,57,130,131,132,133,134,135,136,137}
}
;}
eC2
grammar_optimize_shortcut_logical_evaluation_type{y2
53
cV
grammar_optimize_shortcut_logical_evaluation_type
grammar_optimize_shortcut_logical_evaluation={53,{0,25,cN
iZ1
78,cO
cW}
lS3
l81{tU1
eF2
ParamSpec_Extract
tL1
paramlist,lG1){index=(paramlist>>(index*10))&1023;if(index>=57)return
eF2(SubFunction,cM2
plist_s[index-57]cQ3
index>=37)return
eF2(NumConstant,cM2
plist_n_container
xJ::plist_n[index-37])iH
eF2(lI4,cM2
plist_p[index]);}
}
#ifdef FP_SUPPORT_OPTIMIZER
#include <stdio.h>
#include <algorithm>
#include <map>
#include <sstream>
using
lS3
FUNCTIONPARSERTYPES;using
lS3
l81;using
t6;using
cE1;lS3{nS1
It,tE3
T,tE3
Comp>tF1
MyEqualRange(It
first,It
last,const
T&val,Comp
comp){size_t
len=last-first;while(len>0){size_t
y03
len/2;It
xJ3(first);xJ3+=half;if(comp(*xJ3,val)){first=xJ3;++first;len=len-half-1;}
iQ1
comp(val,*xJ3)){len=half;}
else{It
left(first);{It&eW2=left;It
last2(xJ3)cU3
len2=last2-eW2;while(len2>0){size_t
half2=len2/2;It
tW3(eW2);tW3+=half2;if(comp(*tW3,val)){eW2=tW3;++eW2;len2=len2-half2-1;}
else
len2=half2;}
}
first+=len;It
right(++xJ3);{It&eW2=right;It&last2=first
cU3
len2=last2-eW2;while(len2>0){size_t
half2=len2/2;It
tW3(eW2);tW3+=half2;if(comp(val,*tW3))len2=half2;else{eW2=tW3;++eW2;len2=len2-half2-1;}
}
}
return
tF1(left,right);}
}
return
tF1(first,first);}
tU1
eC2
OpcodeRuleCompare{iV2()n03
tree,unsigned
yL2)const{const
Rule&rule=grammar_rules[yL2]iH
tY2<rule
cN2.subfunc_opcode;}
iV2()tL1
yL2,const
eQ
const{const
Rule&rule=grammar_rules[yL2]iH
rule
cN2.subfunc_opcode<tY2;}
}
n01
bool
TestRuleAndApplyIfMatch
eS1
tI2
n41
tree,bool
yY{MatchInfo
xJ
info;n61
found(false,e0()cQ3(rule.lF1
LogicalContextOnly)&&!yY{tY1
if(nB
IsIntType
xJ::x63){if(rule.lF1
NotForIntegers)tY1
e43
rule.lF1
OnlyForIntegers)tY1
if(nB
IsComplexType
xJ::x63){if(rule.lF1
NotForComplex)tY1
e43
rule.lF1
OnlyForComplex)tY1
for(;;){
#ifdef DEBUG_SUBSTITUTIONS
#endif
found=TestParams(rule
cN2,tree,found.specs,info,true
cQ3
found.found)break;if(!&*found.specs){fail:;
#ifdef DEBUG_SUBSTITUTIONS
DumpMatch(rule
cW3,false);
#endif
return
tD3}
#ifdef DEBUG_SUBSTITUTIONS
DumpMatch(rule
cW3,true);
#endif
SynthesizeRule(rule
cW3
iG}
cE1{xO1
ApplyGrammar
eS1
Grammar&i32,n41
tree,bool
yY{if
iI.GetOptimizedUsing()==&i32){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Already optimized:  "
cL3
iI)xZ1"\n"
<<std::flush;
#endif
return
tD3
if(true){bool
changed=false;switch
iI
nC
lQ4
cNot:case
cNotNot:case
cAnd:case
cOr:for
nD2
eQ3<tree.x6
true))yQ1
lC
cIf:case
cAbsIf:if(ApplyGrammar(i32,eK1,tY2==cIf))yQ1
for
nD2
a=1
eP3<tree.x6
yY)yQ1
break;yI3
for
nD2
eQ3<tree.x6
false))yQ1}
if(changed){tree.Mark_Incompletely_Hashed(iG}
typedef
const
unsigned
short*nP3;std::pair<nP3,nP3>range=MyEqualRange(i32.rule_list,i32.rule_list+i32.rule_count,tree,OpcodeRuleCompare
xJ());std::vector<unsigned
short>rules;rules.xV3
range
t93-range.first);for
xY
if(IsLogisticallyPlausibleParamsMatch(e71
cN2,tree))rules.push_back(*r);}
range.first=!rules
i13?&rules[0]:0;range
t93=!rules
i13?&rules[rules.size()-1]+1:0;if(range.first!=range
t93){
#ifdef DEBUG_SUBSTITUTIONS
if(range.first!=range
t93)tK2"Input ("
<<lR4
iI
nC)<<")["
<<tree
yS1<<"]"
;if(yY
std::cout<<"(Logical)"
iZ3
first=l02,prev=l02;const
char*sep=", rules "
;for
xY
if(first==l02)first=prev=*r;iQ1*r==prev+1)prev=*r;else
tK2
sep<<first;sep=","
;if(prev!=first)std::cout<<'-'<<prev;first=prev=*r;}
}
if(first!=l02)tK2
sep<<first;if(prev!=first)std::cout<<'-'<<prev;}
std::cout<<": "
cL3
iI)xZ1"\n"
<<std::flush;}
#endif
bool
changed=false;for
xY
#ifndef DEBUG_SUBSTITUTIONS
if(!IsLogisticallyPlausibleParamsMatch(e71
cN2,tree))n33
#endif
if(TestRuleAndApplyIfMatch(e71,tree,yY){yQ1
c73}
if(changed){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Changed."
<<std::endl
xZ1"Output: "
cL3
iI)xZ1"\n"
<<std::flush;
#endif
tree.Mark_Incompletely_Hashed(iG}
tree.SetOptimizedUsing(&i32)iH
tD3
xO1
ApplyGrammar
eS1
void*p,FPoptimizer_CodeTree::eQ
yK
ApplyGrammar(*eS1
Grammar*)p,tree);}
t81
ApplyGrammars(FPoptimizer_CodeTree::eQ{
#ifdef DEBUG_SUBSTITUTIONS
std
iN3"grammar_optimize_round1\n"
;
#endif
n6
grammar_optimize_round1
eX2
#ifdef DEBUG_SUBSTITUTIONS
std
iN3"grammar_optimize_round2\n"
;
#endif
n6
grammar_optimize_round2
eX2
#ifdef DEBUG_SUBSTITUTIONS
std
iN3"grammar_optimize_round3\n"
;
#endif
n6
grammar_optimize_round3
eX2
#ifndef FP_ENABLE_SHORTCUT_LOGICAL_EVALUATION
#ifdef DEBUG_SUBSTITUTIONS
std
iN3"grammar_optimize_nonshortcut_logical_evaluation\n"
;
#endif
n6
grammar_optimize_nonshortcut_logical_evaluation
eX2
#endif
#ifdef DEBUG_SUBSTITUTIONS
std
iN3"grammar_optimize_round4\n"
;
#endif
n6
grammar_optimize_round4
eX2
#ifdef FP_ENABLE_SHORTCUT_LOGICAL_EVALUATION
#ifdef DEBUG_SUBSTITUTIONS
std
iN3"grammar_optimize_shortcut_logical_evaluation\n"
;
#endif
n6
grammar_optimize_shortcut_logical_evaluation
eX2
#endif
#ifdef FP_ENABLE_IGNORE_IF_SIDEEFFECTS
#ifdef DEBUG_SUBSTITUTIONS
std
iN3"grammar_optimize_ignore_if_sideeffects\n"
;
#endif
n6
grammar_optimize_ignore_if_sideeffects
eX2
#endif
#ifdef DEBUG_SUBSTITUTIONS
std
iN3"grammar_optimize_abslogical\n"
;
#endif
n6
grammar_optimize_abslogical
eX2
#undef C
}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
#include <algorithm>
#include <assert.h>
#include <cstring>
#include <cmath>
#include <memory> /* for auto_ptr */
using
lS3
FUNCTIONPARSERTYPES;using
lS3
l81;using
t6;using
cE1;lS3{xO1
TestImmedConstraints
tL1
bitmask,const
eQ
tU3(bitmask&ValueMask
lQ4
Value_AnyNum:case
ValueMask:lC
x92:if(GetEvennessInfo
iI)!=n52
Value_OddInt:if(GetEvennessInfo
iI)!=xB2
i51:if(GetIntegerInfo
iI)!=n52
Value_NonInteger:if(GetIntegerInfo
iI)!=xB2
tC1:if(!IsLogicalValue
iI)cR
nH1
SignMask
lQ4
Sign_AnySign:lC
nF1:if(lQ1
n52
tD1:if(lQ1
xB2
Sign_NoIdea:if(lQ1
Unknown
cR
nH1
OnenessMask
lQ4
Oneness_Any:case
OnenessMask:lC
Oneness_One:if
eV1
if(!nK3
fp_abs
iI
eY1),l03(1))cR
lC
Oneness_NotOne:if
eV1
if
lE3
fp_abs
iI
eY1),l03(1))cR
nH1
ConstnessMask
lQ4
Constness_Any:lC
i41:if
eV1
lC
Constness_NotConst:if
iI
e62)cR
c73
return
true;}
tW1
unsigned
extent,unsigned
nbits,tE3
t42=unsigned
int>eC2
nbitmap{private:static
const
unsigned
bits_in_char=8;static
const
unsigned
t52=(c23
t42)*bits_in_char)/nbits;t42
data[(extent+t52-1)/t52];eT3
void
inc(lG1,int
by=1){data[pos(index)]+=by*t42(1<<yM2);nU1
void
dec(lG1){inc(index,-1);}
int
get(lG1
yZ1(data[pos(index)]>>yM2)&mask()cR3
pos(lG1)yK
index/t52
cR3
shift(lG1)yK
nbits*(index%t52)cR3
mask()yK(1<<nbits)-1
cR3
mask(lG1)yK
mask()<<yM2;}
}
;eC2
eD3{int
SubTrees:8;int
e03:8;int
l12:8;int
lT4:8;nbitmap<lI3,2>SubTreesDetail;eD3(tS3
memset(this,0,c23*this));}
eD3
eS1
eD3&b
tS3
memcpy(this,&b,c23
b));}
eD3&eB1=eS1
eD3&b
tS3
memcpy(this,&b,c23
b))iH*this;}
}
n01
eD3
CreateNeedList_uncached(tC&cC2){eD3
c81
tK1
eQ3<cC2
yN2;++a){const
eF2&nQ3=x71
cC2.param_list,a);switch(nQ3.first
lQ4
SubFunction:{cM&c33
cM*eP1;yE
GroupFunction)++c81.lT4;else{++c81.SubTrees;assert(iT3.subfunc_opcode<VarBegin);c81.SubTreesDetail.inc
iG3
lC1);}
++c81.l12
c63}
case
lH4
case
lI4:++c81.e03;++c81.l12
c63}
}
return
c81;}
tU1
eD3&CreateNeedList(tC&cC2){typedef
std::map<tC*,eD3>e81;static
e81
yW1;e81::yA3
i=yW1.y52&cC2
cQ3
i!=yW1
c71&cC2)return
i
eH2
iH
yW1.y63,std::make_pair(&cC2,CreateNeedList_uncached
xJ(cC2)))eH2;}
tU1
CodeTree
xJ
CalculateGroupFunction
eS1
eF2&yO2
const
lP4)tU3(nQ3.first
lQ4
lH4{const
ParamSpec_NumConstant
xJ&param=*eS1
ParamSpec_NumConstant
xJ*eP1
iH
CodeTreeImmed
iG3
constvalue);}
case
lI4:{cL&c33
cL*eP1
iH
info.GetParamHolderValueIfFound
iG3
index);}
case
SubFunction:{cM&c33
cM*eP1
iJ
x63;x63
xC
param.lC1);iP2
lD2).xV3
iT3
yN2)tK1
eQ3<iT3
yN2;++a)i01
tmp(CalculateGroupFunction(x71
iT3.param_list,a),info));iP2
yD1
tmp);}
x63
yD2
iH
x63;}
}
return
CodeTree
xJ();}
}
cE1{xO1
IsLogisticallyPlausibleParamsMatch(tC&cC2,const
eQ{eD3
c81(CreateNeedList
xJ(cC2))cU3
tT3=c5
if(tT3<size_t(c81.l12))i42
for
nD2
eQ3<tT3;++a){unsigned
opcode=nB3
nC;switch(opcode
l43
if(c81.lT4>0)eM3
lT4;else
eM3
e03;lC
lI3:case
cFCall:case
cPCall:eM3
e03
c63
yI3
assert(opcode<VarBegin);if(c81.SubTrees>0&&c81.SubTreesDetail.get(opcode)>0){eM3
SubTrees;c81.SubTreesDetail.dec(opcode);}
else
eM3
e03;}
}
if(c81.lT4>0||c81.SubTrees>0||c81.e03>0)i42
if(cC2.match_type!=AnyParams){if(0||c81.SubTrees<0||c81.e03<0)i42}
return
true;}
tU1
n61
TestParam
eS1
eF2&yO2
xN2
tree
tH2
start_at,lP4)tU3(nQ3.first
lQ4
lH4{const
ParamSpec_NumConstant
xJ&param=*eS1
ParamSpec_NumConstant
xJ*eP1;if
eV1
l03
imm=tree
eY1;switch
iG3
modulo
lQ4
Modulo_None:lC
Modulo_Radians:imm=e53
imm,yT
imm<xG1
imm
yV
if(imm>fp_const_pi
xJ())imm-=fp_const_twopi
xJ(yS2
return
nK3
imm,param.constvalue);}
case
lI4:{cL&c33
cL*eP1;if(!x1
return
info.SaveOrTestParamHolder
iG3
index,tree);}
case
SubFunction:{cM&c33
cM*eP1;yE
GroupFunction){if(!x1
CodeTree
xJ
xY1=CalculateGroupFunction(yO2
info);
#ifdef DEBUG_SUBSTITUTIONS
DumpHashes(xY1)xZ1*eS1
void**)&xY1
eY1
xZ1"\n"
xZ1*eS1
void**)&tree
eY1
xZ1"\n"
;DumpHashes
iI)xZ1"Comparing "
cL3(xY1)xZ1" and "
cL3
iI)xZ1": "
xZ1(xY1
xE
tree)?"true"
:"false"
)xZ1"\n"
;
#endif
return
xY1
xE
tree);}
e43!nF2){if(!x1
if
iI
nC!=param.lC1
cR}
return
TestParams
iG3
data,tree,start_at,info,false);}
}
}
return
tD3
tU1
eC2
l11
xC2
MatchInfo
xJ
info;l11()lJ4,info(){}
}
n01
class
MatchPositionSpec_PositionalParams
t11
l11
xJ>{eT3
lK3
MatchPositionSpec_PositionalParams
nD2
n):t71
l11
xJ>(n){}
}
;eC2
l22
xC2
l22()lJ4{}
}
;class
yF
t11
l22>{eT3
unsigned
trypos;lK3
yF
nD2
n):t71
l22>(n),trypos(0){}
}
n01
n61
TestParam_AnyWhere
eS1
eF2&yO2
xN2
tree
tH2
start_at,lP4,nC3&used,bool
nO2{xQ<yF>x5
iZ3
nN2
yF*)nF2
eP3=x5->trypos;goto
retry_anywhere_2;}
cO2
yF
iI
yS1)eP3=0;}
eR3
c5++a){if(used[a])n33
retry_anywhere
eN3
TestParam(yO2
nB3
cX3);cZ3
used[a]=true;if(nO2
lG4
a);x5->trypos=a
iH
n61(true,&*x5);}
}
retry_anywhere_2
eO3
goto
retry_anywhere;}
}
return
tD3
tU1
eC2
yI1
xC2
MatchInfo
xJ
info;nC3
used;lK3
yI1
nD2
tT3)lJ4,info(),used(tT3){}
}
n01
class
MatchPositionSpec_AnyParams
t11
yI1
xJ>{eT3
lK3
MatchPositionSpec_AnyParams
nD2
n,size_t
m):t71
yI1
xJ>(n,yI1
xJ(m)){}
}
n01
n61
TestParams(tC&nN,xN2
tree
tH2
start_at,lP4,bool
nO2{if(nN.match_type!=AnyParams){if(xU!=tree
yS1
cR}
if(!IsLogisticallyPlausibleParamsMatch(nN,tree))i42
switch(nN.match_type
lQ4
PositionalParams:{xQ<cB>x5
iZ3
nN2
cB*)nF2
eP3=xU-1;goto
lH1;}
cO2
cB(xU)eP3=0;}
eR3
xU;++a){(*x5)[a].i52
retry_positionalparams
eN3
TestParam(cY
a),nB3
cX3);cZ3
tF2}
lH1
eO3
cY3
a].info;goto
retry_positionalparams;}
c93--a;goto
lH1;}
cY3
0].info
iH
tD3
if(nO2
for
tL1
eQ3<xU;++a)lG4
a)iH
n61(true,&*x5);}
case
SelectedParams:case
AnyParams:{xQ<t7>x5;nC3
used
iI
yS1);std::vector<unsigned>lL3(xU);std::vector<unsigned>yP2(xU)yR1{const
eF2
nQ3=cY
a);lL3[a]=ParamSpec_GetDepCode(nQ3);}
{unsigned
b=0
yR1
if(lL3[a]!=0)yP2[b++]=a
yR1
if(lL3[a]==0)yP2[b++]=a;}
unsigned
nN2
t7*)nF2;if(xU==0){a=0;goto
retry_anyparams_4;}
a=xU-1;goto
e91;}
cO2
t7(xU,tree
yS1)eP3=0;if(xU!=0){(*x5)[0].i52(*x5)[0].used=used;}
}
eR3
xU;++a){c93(*x5)[a].i52(*x5)[a].used=used;}
retry_anyparams
eN3
TestParam_AnyWhere
xJ(cY
yP2[a]),tree
cX3,used,nO2;cZ3
tF2}
e91
eO3
cY3
a].info;used=(*x5)[a].used;goto
retry_anyparams;}
eA1:c93--a;goto
e91;}
cY3
0].info
iH
tD3
retry_anyparams_4:if(nN.n1!=0){if(!TopLevel||!info.HasRestHolder(nN.n1)){eN
cP2;cP2.reserve
iI
yS1)tK1
b=0;b<c5++b){if(cB3)n33
cP2.push_back
iI
l8
b));cB3=true;if(nO2
lG4
b);}
if(!info.SaveOrTestRestHolder(nN.n1,cP2)){goto
eA1;}
}
else{l23&cP2=info.GetRestHolderValues(nN.n1)nA3
0
eP3<cP2
l13
a){bool
found=false
tK1
b=0;b<c5++b){if(cB3)n33
if(cP2[a]xE
tree
l8
b))){cB3=true;if(nO2
lG4
b);found=true
c63}
}
if(!found){goto
eA1;}
}
}
}
return
n61(true,xU?&*x5:0);}
case
GroupFunction:c73
return
tD3}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
#include <algorithm>
#include <assert.h>
using
t6;using
cE1;lS3{tU1
CodeTree
xJ
y01
const
eF2&yO2
lP4,bool
inner=true)tU3(nQ3.first
lQ4
lH4{const
ParamSpec_NumConstant
xJ&param=*eS1
ParamSpec_NumConstant
xJ*eP1
iH
CodeTreeImmed
iG3
constvalue);}
case
lI4:{cL&c33
cL*eP1
iH
info.GetParamHolderValue
iG3
index);}
case
SubFunction:{cM&c33
cM*eP1
iJ
tree;eX3
param.lC1)tK1
eQ3<iT3
yN2;++a)i01
nparam=y01
x71
iT3.param_list,a),info,true
lA4
nparam);}
if
iG3
data.n1!=0){eN
trees(info.GetRestHolderValues
iG3
data.n1)tI1
AddParamsMove(trees);if
iI
yS1==1){assert(tree.GetOpcode()==cAdd lF4()==cMul lF4()==cMin lF4()==cMax lF4()==cAnd lF4()==cOr lF4()==cAbsAnd lF4()==cAbsOr);tree
eJ1
0));}
else
if
iI
yS1==0)tU3
iI
nC
lQ4
cAdd:case
cOr:tree=nD1
0));lC
cMul:case
cAnd:tree=nD1
1));yI3
c73}
}
if(inner)tree
yD2
iH
tree;}
}
return
CodeTree
xJ();}
}
cE1{t81
SynthesizeRule
eS1
tI2
n41
tree,lP4)tU3(rule.ruletype
lQ4
ProduceNewTree:{tree.Become(y01
x71
rule.tM1
0),info,false)yS2
case
ReplaceParams:yI3{std::vector<unsigned>list=info.GetMatchedParamIndexes();std::sort(list.l93
list
lM4)nA3
list.size()eP3-->0;)tree.lR3
list[a])tK1
eQ3<rule.repl_param_count;++a)i01
nparam=y01
x71
rule.tM1
a),info,true
lA4
nparam);}
c73}
}
}
#endif
#ifdef DEBUG_SUBSTITUTIONS
#include <sstream>
#include <cstring>
using
lS3
FUNCTIONPARSERTYPES;using
lS3
l81;using
t6;using
cE1;lS3
l81{t81
DumpMatch
eS1
tI2
xN2
tree,const
lP4,bool
DidMatch,std::ostream&o){DumpMatch(rule
cW3,DidMatch?l44"match"
:l44"mismatch"
,o);}
t81
DumpMatch
eS1
tI2
xN2
tree,const
lP4,const
char*tX3,std::ostream&o){static
const
char
ParamHolderNames[][2]={"%"
,"&"
,"x"
,"y"
,"z"
,"a"
,"b"
,"c"
}
;o<<tX3<<" (rule "
<<(&rule-grammar_rules)<<")"
<<":\n  Pattern    : "
;{eF2
tmp;tmp.first=SubFunction;ParamSpec_SubFunction
tmp2;tmp2.data=rule
cN2;tmp
t93=cM2
tmp2;DumpParam
xJ(tmp,o);}
o<<"\n  Replacement: "
;DumpParams
xJ(rule.tM1
rule.repl_param_count,o
tN3
o<<"  Tree       : "
cL3
iI,o
tN3
if(!std::strcmp(tX3,l44"match"
))DumpHashes
iI,o)nA3
0
eP3<info.c9
size();++a){if(!info.paramholder_matches[a]y12))n33
o<<"           "
<<ParamHolderNames[a]<<" = "
cL3(info.paramholder_matches[a],o
tN3}
eZ1
info.lQ
l13
b){if(!i62
first)n33
for
nD2
eQ3<i62
second
l13
a){o<<"         <"
<<b<<"> = "
cL3(i62
second[a],o);o<<std::endl;}
}
o<<std::flush;}
}
#endif
#include <list>
#include <algorithm>
#ifdef FP_SUPPORT_OPTIMIZER
using
lS3
FUNCTIONPARSERTYPES;lS3{xO1
MarkIncompletes(FPoptimizer_CodeTree::eQ{if
iI.Is_Incompletely_Hashed(t41
bool
l32=false;lN1
l32|=MarkIncompletes
iI
l8
a)cQ3
l32)tree.Mark_Incompletely_Hashed()iH
l32;}
t81
FixIncompletes(FPoptimizer_CodeTree::eQ{if
iI.Is_Incompletely_Hashed()){lN1
FixIncompletes
iI
l8
a)tI1
Rehash();}
}
}
t6{lB
Sort()iI3
Sort();}
lB
Rehash(bool
constantfolding){if(constantfolding)ConstantFolding(*this);else
Sort();data
x7
tU1
eC2
c7{cC3
iW2
cD3
yJ1=0;
#if 0
long
double
value=Value;eD=crc32::calc(eS1
unsigned
char*)&value,c23
value));key^=(key<<24);
#elif 0
union{eC2{unsigned
char
filler1[16];l03
v
iZ3
char
filler2[16];}
buf2;eC2{unsigned
char
filler3[sizeof
tQ2)+16-c23
xC1)];eD;}
buf1;}
data;memset(&data,0,c23
data));data.buf2.v=Value;eD=data.buf1.key;
#else
int
exponent;l03
x02=std::frexp(Value,&xT2;eD=tL1
iY3+0x8000)&0xFFFF
cQ3
x02<0){x02=-x02;key=key^0xFFFF;}
else
key+=0x10000;x02-=l03(0.5);key<<=39;key|=n71(x02+x02)*l03(1u<<31))<<8;
#endif
lP
#ifdef FP_SUPPORT_COMPLEX_NUMBERS
nS1
T
cQ2
std::complex<T> >{cC3
std::complex<T>&cD3
c7<T>::nR3
cR2,Value.real());nB
fphash_t
temp;c7<T>::nR3
temp,Value.imag());yJ1^=temp.hash2;cR2.hash2^=temp.hash1;}
}
;
#endif
#ifdef FP_SUPPORT_LONG_INT_TYPE
tW1
cQ2
long>{yG
long
cD3
eD=Value;lP
#endif
#ifdef FP_SUPPORT_GMP_INT_TYPE
tW1
cQ2
GmpInt>{cC3
GmpInt&cD3
eD=Value.toInt();lP
#endif
t81
xU2
xJ::Recalculate_Hash_NoRecursion(){fphash_t
cR2(n71
Opcode)<<56,Opcode*iQ3(0x1131462E270012B));Depth=1;switch(Opcode
l43{c7
xJ::nR3
cR2,Value
yS2
case
lI3:{yJ1|=n71
cO1<<48
lO1((n71
cO1)*11)^iQ3(0x3A83A83A83A83A0)c63}
case
cFCall:case
cPCall:{yJ1|=n71
cO1<<48
lO1((~n71
cO1)*7)^3456789;}
yI3{size_t
tN1=0
nA3
0
eP3<eE3
l13
a){if(eE3[a]nY2>tN1)tN1=eE3[a]nY2;yJ1+=((eE3[a].i92
hash1*(a+1))>>12)lO1
eE3[a].i92
hash1
lO1(3)*iQ3(0x9ABCD801357);cR2.hash2*=iQ3(0xECADB912345)lO1(~eE3[a].i92
hash2)^4567890;}
Depth+=tN1;}
}
if(Hash!=cR2){Hash=cR2;lG2=0;}
}
lB
FixIncompleteHashes(){MarkIncompletes(*this);FixIncompletes(*this);}
}
#endif
#include <cmath>
#include <list>
#include <cassert>
#ifdef FP_SUPPORT_OPTIMIZER
using
lS3
FUNCTIONPARSERTYPES;lS3{using
t6
n01
bool
x11
xN2
tree,long
count,const
xI1::SequenceOpCode
xJ&eP,xI1
nD3&synth,size_t
max_bytecode_grow_length);static
const
eC2
SinCosTanDataType{OPCODE
whichopcode;OPCODE
inverse_opcode;enum{nominator,x12,inverse_nominator,lI1}
;OPCODE
codes[4];}
SinCosTanData[12]={{cTan,cCot,{cSin,cCos,cCsc,cSec}
}
,{cCot,cCot,{cCos,cSin,cSec,cCsc}
}
,{cCos,cSec,{cSin,cTan,cCsc,cCot}
}
,{cSec,cCos,{cTan,cSin,cCot,cCsc}
}
,{cSin,cCsc,{cCos,cCot,cSec,cTan}
}
,{cCsc,cSin,{cCot,cCos,cTan,cSec}
}
,{cS2{cSinh,cCosh,cT2,{cSinh,cNop,{cS2
cNop,cCosh}
}
,{cCosh,cNop,{cSinh,cS2
cNop}
}
,{cNop,cTanh,{cCosh,cSinh,cT2,{cNop,cSinh,{cNop,cTanh,cCosh,cNop}
}
,{cNop,cCosh,{cTanh,cSinh,cT2}
;}
t6{lB
SynthesizeByteCode(std::vector<unsigned>&ByteCode,std::vector
xJ&Immed,size_t&stacktop_max){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Making bytecode for:\n"
;iW
#endif
while(RecreateInversionsAndNegations()){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"One change issued, produced:\n"
;iW
#endif
FixIncompleteHashes();using
cE1;using
lS3
l81;const
void*g=cM2
grammar_optimize_recreate;while(ApplyGrammar(*eS1
Grammar*)g,*this)){FixIncompleteHashes();}
}
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Actually synthesizing, after recreating inv/neg:\n"
;iW
#endif
xI1
nD3
synth;SynthesizeByteCode(synth,false);synth.Pull(ByteCode,Immed,stacktop_max);}
lB
SynthesizeByteCode(xI1
nD3&synth,bool
MustPopTemps)const{y11*this))yK;}
for
nD2
eQ3<12;++a){const
SinCosTanDataType&data=SinCosTanData[a];if(data.whichopcode!=cNop)iA1!=data.whichopcode)n33
CodeTree
nS3;nS3.lR1
nS3
xC
data.inverse_opcode);nS3.yQ2);y11
nS3)){synth.l61
else
iA1!=cInv)n33
if(GetParam(0)nC!=data.inverse_opcode)n33
y11
GetParam(0))){synth.l61
size_t
found[4];eZ1
4;++b){CodeTree
tree;if(data.tY3]==cNop){eX3
cInv);CodeTree
nT3;nT3.lR1
nT3
xC
data.tY3^2]);nT3.yQ2
lA4
nT3);}
else{tree.lR1
eX3
data.tY3]);}
tree.yQ2);found[b]tZ3
t13
iI);}
if(found[data.yU2!=tG
x12]y4
yU2);nM1
x12
i9
cDiv
nI1
yU2!=tG
lI1]y4
yU2);nM1
lI1
i9
cMul
nI1
lX1!=tG
lI1]y4
lX1);nM1
lI1
i9
cRDiv
nI1
lX1!=tG
x12]y4
lX1);nM1
x12
i9
cMul,2,1);synth.l61
size_t
n_subexpressions_synthesized=SynthCommonSubExpressions(synth);switch(l42{case
lI3:synth.PushVar(GetVar());lC
cImmed:nK2
GetImmed());lC
cAdd:case
cMul:case
cMin:case
cMax:case
cAnd:case
cOr:case
cAbsAnd:case
cAbsOr:iA1==cMul){bool
cE3=false;yH
lY1
e62)&&isLongInteger(lY1
eY1)){c01=makeLongInteger(lY1
eY1);CodeTree
tmp(*this,tE3
CodeTree::CloneTag());tmp
yT1
tmp
xJ2
if(x11
tmp,value,xI1::i11
xJ::AddSequence,synth,MAX_MULI_BYTECODE_LENGTH)){cE3=true
c63}
}
}
if(cE3)c73
int
yK1=0;nC3
done(GetParamCount(),false);CodeTree
iM;iM
xC
l42;for(;;){bool
found=false;yH
done[a])n33
if(synth.IsStackTop(lY1)){found=true;done[a]=true;lY1.n7
iM
cG
lY1
cQ3++yK1>1){synth
yI
2);iM.yQ2
iD1
iM);yK1=yK1-2+1;}
}
}
if(!found)c73
yH
done[a])n33
lY1.n7
iM
cG
lY1
cQ3++yK1>1){synth
yI
2);iM.yQ2
iD1
iM);yK1=yK1-2+1;}
}
if(yK1==0)tU3(l42{case
cAdd:case
cOr:case
cAbsOr:nK2
0);lC
cMul:case
cAnd:case
cAbsAnd:nK2
1);lC
cMin:case
cMax:nK2
0)c63
yI3
c73++yK1;}
assert(n_stacked==1)c63}
case
cPow:{lT3
p0=GetParam(0);lT3
p1=GetParam(1
cQ3!p1.cG2!isLongInteger
xE1)||!x11
p0,makeLongInteger
xE1),xI1::i11
xJ::MulSequence,synth,MAX_POWI_BYTECODE_LENGTH)){p0.n7
p1.n7
synth
yI
2);cS1
cIf:case
cAbsIf:{tE3
xI1
nD3::IfData
ifdata;GetParam(0).n7
synth.SynthIfStep1(ifdata,l42;GetParam(1).n7
synth.SynthIfStep2(ifdata);GetParam(2).n7
synth.SynthIfStep3(ifdata
yS2
case
cFCall:case
cPCall:{for
nD2
eQ3<l71++a)lY1.n7
synth
yI
tL1)GetParamCount());lY3
nS2|GetFuncNo(),0,0
yS2
yI3{for
nD2
eQ3<l71++a)lY1.n7
synth
yI
tL1)GetParamCount()yS2}
synth.StackTopIs(*this
cQ3
MustPopTemps&&n_subexpressions_synthesized>0){size_t
top
tZ3
GetStackTop();synth.DoPopNMov(top-1-n_subexpressions_synthesized,top-1);}
}
}
lS3{xO1
x11
xN2
tree,long
count,const
xI1::SequenceOpCode
xJ&eP,xI1
nD3&synth,size_t
max_bytecode_grow_length){if
eG3!=0){xI1
nD3
backup=synth;tree.n7
size_t
bytecodesize_backup
tZ3
GetByteCodeSize();xI1::x11
count
nV2
size_t
bytecode_grow_amount
tZ3
GetByteCodeSize()-bytecodesize_backup;if(bytecode_grow_amount>max_bytecode_grow_length){synth=backup
iH
tD3
return
true;}
else{xI1::x11
count,eP,synth
iG}
}
#endif
#include <cmath>
#include <cassert>
#ifdef FP_SUPPORT_OPTIMIZER
using
lS3
FUNCTIONPARSERTYPES;lS3{using
t6;
#define FactorStack std::vector
const
eC2
PowiMuliType{unsigned
opcode_square
iZ3
opcode_cumulate
iZ3
opcode_invert
iZ3
opcode_half
iZ3
opcode_invhalf;}
iseq_powi={cSqr,cMul,cInv,cSqrt,cRSqrt}
,iseq_muli={l02
xD
cNeg,l02,l02}
n01
l03
cP1
const
PowiMuliType&cF3,const
std
nE3,n62&stack
t12
1);while(IP<limit){if(l52==cF3.opcode_square){if(!isInteger
iB1
2;yO
opcode_invert){x63=-x63;yO
opcode_half){if
iI2>y71&&isEvenInteger
iB1
l03(0.5);yO
opcode_invhalf){if
iI2>y71&&isEvenInteger
iB1
l03(-0.5);++IP;tF2
size_t
xD2=IP;l03
lhs(1
cQ3
l52==cFetch){lG1=n82;if(index<y1||size_t(index-y1)>=l82){IP=xD2
c63}
lhs=stack[index-y1];goto
yV2;}
if(l52==cDup){lhs=x63;goto
yV2;yV2:iG1
iI2);++IP;l03
subexponent=cP1
cF3
tS
if(IP>=limit||l52!=cF3.opcode_cumulate){IP=xD2
c63}
++IP;stack.pop_back();x63+=lhs*subexponent;tF2
c73
return
x63;}
tU1
l03
ParsePowiSequence
eS1
std
nE3){n62
stack;iG1
tQ2(1))iH
cP1
iseq_powi
tS}
tU1
l03
ParseMuliSequence
eS1
std
nE3){n62
stack;iG1
tQ2(1))iH
cP1
iseq_muli
tS}
tU1
class
CodeTreeParserData{eT3
lK3
CodeTreeParserData(bool
k_powi):stack(),clones(),keep_powi(k_powi){}
void
Eat
nD2
tT3,OPCODE
opcode)i01
xH;xH
xC
opcode);eN
cC2=Pop(tT3
l62
cC2
cQ3!keep_powi)switch(opcode
lQ4
cTanh:i01
sinh,cosh;sinh
xC
cSinh);sinh
cG
xH
cG3
sinh
xJ2
cosh
xC
cCosh);cosh.yD1
xH
cG3
cosh
i23
pow
yR2.yD1
cosh);pow.yB
l03(-1)));pow
xJ2
xH
eV3
xH.nE1
0,sinh);xH.yD1
pow
yS2
case
cTan:i01
sin,cos;sin
xC
cSin);sin
cG
xH
cG3
sin
xJ2
cos
xC
cCos);cos.yD1
xH
cG3
cos
i23
pow
yR2.yD1
cos);pow.yB
l03(-1)));pow
xJ2
xH
eV3
xH.nE1
0,sin);xH.yD1
pow
yS2
case
cPow:{xN2
p0=xH
l8
0
l21
p1=xH
l8
1
cQ3
p1
nC==cAdd){eN
xN3(p1
yS1)nA3
0
eP3<p1.l71++a)i01
pow
yR2
cG
p0);pow
cG
p1
l8
a));pow
xJ2
xN3[a
i03
pow);}
xH
xC
cMul
l62
tA2;}
c73
yI3
c73
xH.Rehash(!keep_powi);l72,false);
#ifdef DEBUG_SUBSTITUTIONS
iI1<<tT3<<", "
<<lR4(opcode)<<"->"
<<lR4(xH
nC)<<": "
iW3(xH
i8
xH);
#endif
iG1(yJ3
EatFunc
nD2
tT3,OPCODE
t23
unsigned
funcno
t03
CodeTreeFuncOp
xJ(t23
funcno);eN
cC2=Pop(tT3
l62
cC2);xH.yQ2);
#ifdef DEBUG_SUBSTITUTIONS
iI1<<tT3<<", "
iW3(xH
i8
xH);
#endif
l72);iG1(yJ3
AddConst(yO1
t03
CodeTreeImmed(value);l72);Push(yJ3
AddVar
tL1
varno
t03
CodeTreeVar
xJ(varno);l72);Push(yJ3
SwapLastTwoInStack(){iH1
1
i03
iH1
2]l64
Dup(){Fetch(l82-1
l64
Fetch
nD2
which){Push(stack[which]);}
nS1
T>void
Push(T
tree){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<iW3
iI
i8
tree);
#endif
iG1
iI
l64
PopNMov
nD2
target,size_t
source){stack[target]=stack[source];stack.eU3
target+1);}
CodeTree
xJ
yW2{clones.clear()iJ
x63(stack.back());stack.eU3
l82-1)iH
x63;}
eN
Pop
nD2
n_pop){eN
x63(n_pop)tK1
n=0;n<n_pop;++n)x63[n
i03
iH1
n_pop+n]);
#ifdef DEBUG_SUBSTITUTIONS
for
nD2
n=n_pop;n-->0;){iI1
cL3
iI2[n]i8
x63[n]);}
#endif
stack.eU3
l82-n_pop)iH
x63;}
size_t
GetStackTop(yZ1
l82;}
private:void
FindClone(n41,bool=true)yK;}
private:eN
stack;std::multimap<fphash_t,CodeTree
xJ>clones;bool
keep_powi;private:CodeTreeParserData
eS1
CodeTreeParserData&);CodeTreeParserData&eB1=eS1
CodeTreeParserData&);}
n01
eC2
IfInfo
i01
t62
iJ
thenbranch
cU3
endif_location;IfInfo():t62(),thenbranch(),endif_location(){}
}
;}
t6{lB
GenerateFrom
eS1
tE3
FunctionParserBase
xJ::Data&xK3,bool
keep_powi){eN
xO2;xO2.xV3
xK3.mVariablesAmount)tK1
n=0;n<xK3.mVariablesAmount;++n){xO2.push_back(CodeTreeVar
xJ(n+lI3));}
GenerateFrom(xK3,xO2,keep_powi);}
lB
GenerateFrom
eS1
tE3
FunctionParserBase
xJ::Data&xK3,const
x3&xO2,bool
keep_powi){const
std::vector<unsigned>&ByteCode=xK3.mByteCode;const
std::vector
xJ&Immed=xK3.mImmed;
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"ENTERS GenerateFrom()\n"
;
#endif
CodeTreeParserData
xJ
sim(keep_powi);std::vector<IfInfo
xJ>eM;for
nD2
IP=0,DP=0;;++IP){iB2:while(!eM
i13&&(eM.eE==IP||(IP<xP2&&l52==cJump&&eM
eW1
y12)))){CodeTree
elsebranch=sim.yW2
e72
eM.back().t62)e72
eM
eW1)e72
elsebranch);lC3
3,cIf);eM.pop_back();}
if(IP>=xP2)break
iZ3
opcode=l52;if((opcode==cSqr||opcode==cDup||(opcode==cInv&&!IsIntType
xJ::x63)||opcode==cNeg||opcode==cSqrt||opcode==cRSqrt||opcode==cFetch)){size_t
was_ip=IP;l03
eJ2
ParsePowiSequence
xJ(ByteCode,IP,eM
i13?xP2:eM.eE,sim.xL
1
cQ3
exponent!=l03(1.0)){xK
xT2
t72;goto
iB2;}
if(opcode==cDup||opcode==cFetch||opcode==cNeg
eZ2
y42=ParseMuliSequence
xJ(ByteCode,IP,eM
i13?xP2:eM.eE,sim.xL
1
cQ3
y42!=l03(1.0)){xK
y42)y3
cMul);goto
iB2;}
}
IP=was_ip;}
if(n72>=lI3){lG1=opcode-lI3
e72
xO2[index]);}
else
tU3(n72
lQ4
cIf:case
cAbsIf:{eM.eU3
eM.size()+1);CodeTree
res(sim.yW2);eM.back().t62.swap(res);eM.eE=xP2;IP+=2;tF2
case
cJump:{CodeTree
res(sim.yW2);eM
eW1.swap(res);eM.eE=lQ3
IP+1]+1;IP+=2;tF2
case
cImmed:xK
Immed[DP++]);lC
cDup:sim.Dup();lC
cNop:lC
cFCall:{unsigned
funcno=n82;assert(funcno<fpdata.mFuncPtrs.size())iZ3
cC2=xK3.mFuncPtrs[funcno].mParams;sim.EatFunc(cC2,n72,funcno
yS2
case
cPCall:{unsigned
funcno=n82;assert(funcno<fpdata.iU3.size());const
FunctionParserBase
xJ&p=*xK3.iU3[funcno].mParserPtr
iZ3
cC2=xK3.iU3[funcno].mParams;x3
paramlist=sim.Pop(cC2);CodeTree
iC2;iC2.GenerateFrom(*p.mData,paramlist)e72
iC2
yS2
case
cInv:xK
1
cW2
cDiv);lC
cNeg
nI2
cNeg)c63
xK
0
cW2
cSub);lC
cSqr:lS4;lC
cSqrt
e82
tO1
cRSqrt
e82-tO1
cCbrt
e82
1)/l03(3))t72;lC
cDeg:xK
fp_const_rad_to_deg
y21
cRad:xK
fp_const_deg_to_rad
y21
cExp:cQ1)goto
lD3;xK
fp_const_e
xJ()cW2
cPow);lC
cExp2:cQ1)goto
lD3;xK
2.0
cW2
cPow);lC
cCot
nI2
cTan
nV
cCsc
nI2
cSin
nV
cSec
nI2
cCos
nV
cInt:
#ifndef __x86_64
cQ1
lA3
1,cInt
yS2
#endif
xK
l03(0.5))nU3
lC3
1,cFloor);lC
cLog10
nI2
cH3
fp_const_log10inv
y21
cLog2
nI2
cH3
fp_const_log2inv
y21
e93:l24
cH3
fp_const_log2inv
xJ());lC3
3,cMul);lC
cHypot:lS4;i43
lS4
nU3
xK
l03(tO1
cSinCos:sim.Dup();lC3
1,cSin);l24
cCos);lC
cSinhCosh:sim.Dup();lC3
1,cSinh);l24
cCosh);lC
cRSub:i43
case
cSub:cQ1
lA3
2,cSub
yS2
xK-1)y3
cMul)nU3
lC
cRDiv:i43
case
cDiv:cQ1||IsIntType
xJ::x63
lA3
2,cDiv
yS2
xK-1)t72
y3
cMul);lC
cAdd:case
cMul:case
cMod:case
cPow:case
cEqual:case
cLess:case
cGreater:case
i91:case
cLessOrEq:case
cGreaterOrEq:case
cAnd:case
cOr:case
cAbsAnd:case
cAbsOr:lC3
2,tP1
lC
cNot:case
cNotNot:case
cN3:case
cAbsNotNot
nI2
tP1
lC
cFetch:sim.Fetch(n82);lC
cPopNMov:{unsigned
stackOffs_target=n82
iZ3
stackOffs_source=n82;sim.PopNMov(stackOffs_target,stackOffs_source
yS2
yI3
lD3:iZ3
funcno=opcode-cAbs;assert(funcno<FUNC_AMOUNT);const
FuncDefinition&func=Functions[funcno];lC3
func.cC2,tP1
c73}
}
Become(sim.yW2);
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Produced tree:\n"
;iW
#endif
}
}
#endif
#include <algorithm>
#ifdef FP_SUPPORT_OPTIMIZER
#include <assert.h>
#define FP_MUL_COMBINE_EXPONENTS
lS3{using
lS3
FUNCTIONPARSERTYPES;using
t6
n01
static
void
AdoptChildrenWithSameOpcode(eQ{
#ifdef DEBUG_SUBSTITUTIONS
bool
x22=false;
#endif
for
cT1
if
iI
l8
a)nC==tY2){
#ifdef DEBUG_SUBSTITUTIONS
if(!x22)tK2"Before assimilation: "
nJ2
x22=true;}
#endif
tree.AddParamsMove
iI
l8
a).GetUniqueRef().lD2),a);}
#ifdef DEBUG_SUBSTITUTIONS
if(x22)tK2"After assimilation:   "
nJ2}
#endif
}
}
t6{t81
ConstantFolding(eQ{tree.Sort();
#ifdef DEBUG_SUBSTITUTIONS
void*cI3=0
xZ1"["
<<(&cI3)<<"]Runs ConstantFolding for: "
nJ2
DumpHashes
iI)xZ1
std::flush;
#endif
if(false){redo:;tree.Sort();
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"["
<<(&cI3)<<"]Re-runs ConstantFolding: "
nJ2
DumpHashes
iI);
#endif
}
if
iI
nC!=cImmed){yE3
p=CalculateResultBoundaries
iI
cQ3
p
cT&&p
cQ&&p.yW==p
tZ1{xN
p.yW);nD}
if(false){ReplaceTreeWithOne:xN
l03(1));goto
do_return;ReplaceTreeWithZero:xN
xG1;goto
do_return;ReplaceTreeWithParam0:
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before replace: "
xZ1
std::hex<<'['
eL1
hash1<<','
eL1
hash2<<']'<<std::dec
nJ2
#endif
tree
eJ1
0));
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After replace: "
xZ1
std::hex<<'['
eL1
hash1<<','
eL1
hash2<<']'<<std::dec
nJ2
#endif
cJ
cJ3
iI
nC
l43
lC
lI3:lC
cAnd:case
cAbsAnd
e4
bool
yZ=false;for
cT1{if(!yO3
a)))yZ=true;cK3
a),tY2==cAbsAnd)lQ4
y43
e3
l53:nN1);lC
lZ1
cJ3
iI
yS1
lQ4
0:iA
1:eX3
tY2==cAnd?cNotNot:cAbsNotNot);cJ
yI3
cU2
cAnd||!yZ)if(ConstantFolding_AndLogic
cP3
cS1
cOr:case
cAbsOr
e4
bool
yZ=false;for
cT1{if(!yO3
a)))yZ=true;cK3
a),tY2==cAbsOr)iF1
iA
lN3
nN1);lC
lZ1
cJ3
iI
yS1
lQ4
0
e3
1:eX3
tY2==cOr?cNotNot:cAbsNotNot);cJ
yI3
cU2
cOr||!yZ)if(ConstantFolding_OrLogic
cP3
cS1
cNot:case
cN3:{unsigned
n81
0;switch
iI
l8
0)nC
lQ4
cEqual:n81
i91;lC
i91:n81
cEqual;lC
cLess:n81
cGreaterOrEq;lC
cGreater:n81
cLessOrEq;lC
cLessOrEq:n81
cGreater;lC
cGreaterOrEq:n81
cLess;lC
cNotNot:n81
cNot;lC
cNot:n81
cNotNot;lC
cN3:n81
cAbsNotNot;lC
cAbsNotNot:n81
cN3
c63
yI3
c73
if(opposite){eX3
OPCODE(opposite)tI1
SetParamsMove
iI
l8
0).GetUniqueRef().lD2));cJ
cJ3(tJ
0),tree
cR1
lQ4
l53
e3
lN3
iA
lZ1
cU2
cNot&&GetPositivityInfo
iI
l8
0))==l53)eX3
cN3);cH2
nC==cIf||eK1
nC==tC3
i01
xL3=tree
l8
0
l21
ifp1=xL3
l8
1
l21
ifp2=xL3
l8
2
cQ3
ifp1
cV2
ifp1
cR1{tree
xI
ifp1
nC==cNot?cNotNot:cAbsNotNot);iD2
cG3
yN3
iJ
p2;p2
iT1
p2
cO3
nX2
if(ifp2
cV2
ifp2
cR1{tree
xI
tY2);iD2);yN3
iJ
p2;p2
xC
ifp2
nC==cNot?cNotNot:cAbsNotNot);p2
cO3
l8
0)nX2
cS1
cNotNot:case
cAbsNotNot:{if(yO3
0)))goto
e2
cK3
0),tY2==cAbsNotNot)lQ4
y43
e3
l53:iA
lZ1
cU2
cNotNot&&GetPositivityInfo
iI
l8
0))==l53)eX3
cAbsNotNot);cH2
nC==cIf||eK1
nC==tC3
i01
xL3=tree
l8
0
l21
ifp1=xL3
l8
1
l21
ifp2=xL3
l8
2
cQ3
ifp1
cV2
ifp1
cR1{tree.SetParam(0,xL3
cG3
tree
cG
ifp1)iJ
p2;p2
iT1
p2
cO3
nX2
if(ifp2
cV2
ifp2
cR1{tree
xI
tY2);iD2);yN3;tree
cO3)xX1
cS1
cIf:case
cAbsIf:{if(ConstantFolding_IfOperations
cP3
c73
case
cMul:{NowWeAreMulGroup:;AdoptChildrenWithSameOpcode
iI);l03
nJ1=l03(1)cU3
l92=0;bool
nK1=false;lN1{if(!tree
y41
n33
l03
immed=nL1
if(immed==xG1
goto
ReplaceTreeWithZero;nJ1*=immed;++l92;}
if(l92>1||(l92==1&&nK3
nJ1,l03(1))))nK1=true;if(nK1){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"cMul: Will add new "
iV3
nJ1<<"\n"
;
#endif
for
cT1
if
iI
y41{
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<" - For that, deleting "
iV3
nL1
std::cout<<"\n"
;
#endif
e23(!nK3
nJ1,l03(1)))tree
cG
e31
nJ1));cJ3
iI
yS1
lQ4
0:iA
1:goto
e2
yI3
if(ConstantFolding_MulGrouping
cP3
if(ConstantFolding_MulLogicItems
cP3
cS1
cAdd
e4
l03
n92=0.0
cU3
l92=0;bool
nK1=false;lN1{if(!tree
y41
n33
l03
immed=nL1
n92+=immed;++l92;}
if(l92>1||(l92==1&&n92==xG1)nK1=true;if(nK1){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"cAdd: Will add new "
iV3
n92<<"\n"
xZ1"In: "
nJ2
#endif
for
cT1
if
iI
y41{
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<" - For that, deleting "
iV3
nL1
std::cout<<"\n"
;
#endif
e23(!(n92==l03(0.0)))tree
cG
e31
n92));cJ3
iI
yS1
lQ4
0
e3
1:goto
e2
yI3
if(ConstantFolding_AddGrouping
cP3
if(ConstantFolding_AddLogicItems
cP3
cS1
cMin
e4
size_t
yX2=0;yE3
e7;lN1{while(a+1<tree
yS1&&tree
n93
tree
l8
a+1)))nN1+1);yF3
max
lM3(!e7
cQ||(p
tZ1<e7
tZ1){e7.yG2=p.yG2;e7
cQ=true;yX2=a;}
}
if(e7
cQ)for
cT1{yF3
min
lM3
a!=yX2&&p.yW>=e7
tZ1
e23
iI
yS1==1){goto
e2
cS1
cMax
e4
size_t
yX2=0;yE3
eZ;lN1{while(a+1<tree
yS1&&tree
n93
tree
l8
a+1)))nN1+1);yF3
min
lM3(!eZ
cT||p.yW>eZ.yW)){eZ.yW=p.yW;eZ
cT=true;yX2=a;}
}
if(eZ
cT){for
cT1{yF3
max
lM3
a!=yX2&&(p
tZ1<eZ.yW){nN1);}
}
}
if
iI
yS1==1){goto
e2
cS1
cEqual:case
i91:case
cLess:case
cGreater:case
cLessOrEq:case
cGreaterOrEq:if(ConstantFolding_Comparison
cP3
lC
cAbs:{yE3
tT
eK1
cQ3
iN
goto
e2
if(p0
cS{eX3
cMul
tI1
yB
l03(1)));goto
NowWeAreMulGroup;}
cH2
nC==cMul){xN2
p=eK1;eN
nV3;eN
cX2
nA3
0
eP3<p.l71++a){tT
p
l8
a)cQ3
iN{nV3.push_back(p
l8
a));}
if(p0
cS{cX2.push_back(p
l8
a));}
}
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Abs: mul group has "
<<nV3.size()<<" pos, "
<<cX2.size()<<"neg\n"
;
#endif
if(!nV3
i13||!cX2
i13){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"AbsReplace-Before: "
cL3
iI)xZ1"\n"
<<std::flush;DumpHashes
iI
iE2);
#endif
CodeTree
xJ
eS3;eS3
xC
cMul)nA3
0
eP3<p.l71++a){tT
p
l8
a)cQ3(iN||(p0
cS){}
else
eS3
cG
p
l8
a));}
eS3
i23
nW3;nW3
xC
cAbs);nW3.yD1
eS3);nW3
i23
y91
cMul
eL2
yD1
nW3
eL2
AddParamsMove(nV3
cQ3!cX2
i13){if(cX2.size()%2)tB2
yB
l03(-1))eL2
AddParamsMove(cX2);}
tree.Become
e32);
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"AbsReplace-After: "
cL3
iI
iE2)xZ1"\n"
<<std::flush;DumpHashes
iI
iE2);
#endif
goto
NowWeAreMulGroup;}
}
c73
#define HANDLE_UNARY_CONST_FUNC(funcname) nR){xN funcname(lR));nD
case
cLog:iK3(fp_log);cH2
nC==cPow)i01
pow=eK1;if(GetPositivityInfo(pow
l8
0))==l53){pow.lB1
pow
lS1
tree.lU
if(GetEvennessInfo(pow
l8
1))==l53){pow.CopyOnWrite()iJ
abs;abs
xC
cAbs);abs.yD1
pow
cG3
abs
xJ2
pow
lS1
pow.nE1
0,abs
tI1
lU}
else
cH2
nC==cAbs)i01
pow=eK1
l8
0
cQ3
pow
nC==cPow){pow.CopyOnWrite()iJ
abs;abs
xC
cAbs);abs.yD1
pow
cG3
abs
xJ2
pow
lS1
pow.nE1
0,abs
tI1
lU}
lC
cAcosh:iK3(fp_acosh);lC
cAsinh:iK3(fp_asinh);lC
cAtanh:iK3(fp_atanh);lC
cAcos:iK3(fp_acos);lC
cAsin:iK3(fp_asin);lC
cAtan:iK3(fp_atan);lC
cCosh:iK3(fp_cosh);lC
cSinh:iK3(fp_sinh);lC
cTanh:iK3(fp_tanh);lC
cSin:iK3(fp_sin);lC
cCos:iK3(fp_cos);lC
cTan:iK3(fp_tan);lC
cCeil:lK4(fp_ceil);lC
cTrunc:lK4(fp_trunc);lC
cFloor:lK4(fp_floor);lC
cInt:lK4(fp_int);lC
cCbrt:iK3(fp_cbrt);lC
cSqrt:iK3(fp_sqrt);lC
cExp:iK3(fp_exp);lC
cLog2:iK3(fp_log2);lC
cLog10:iK3(fp_log10);lC
e93
l54
fp_log2(lR)*eR
iK2
cArg:iK3(fp_arg);lC
cConj:iK3(fp_conj);lC
cImag:iK3(fp_imag);lC
cReal:iK3(fp_real);lC
cPolar
l54
fp_polar(xE2
cMod
l54
e53
xE2
cAtan2:{yE3
tT
eK1);yE3
p1=tU
1));nR&&nK3
lR,xG1){if(p1
cQ&&(p1
tZ1<xG1{xN
fp_const_pi
e33
if(p1
cT&&p1.yW>=tQ1
xG1;nD}
xE3
nK3
eR,xG1){if(p0
cQ&&(p0
tZ1<xG1{xN-fp_const_pihalf
e33
if(p0
cT&&p0.yW>xG1{xN
fp_const_pihalf
e33}
if
lI
fp_atan2(lR,eR));nD
if((p1
cT&&p1.yW>xG1||(p1
cQ&&(p1
tZ1<fp_const_negativezero
xJ()e5
yZ2;yZ2
xC
cPow);yZ2.AddParamMove
iI
l8
1));yZ2.yB
l03(-1)));yZ2
i23
c02;c02
eV3
c02.AddParamMove
iI
cG3
c02.yD1
yZ2);c02
xJ2
eX3
cAtan
yB3
0,c02
y31
1);cS1
cPow:{if(ConstantFolding_PowOperations
cP3
c73
case
cDiv:nR&&tree
eQ1
eR!=tQ1
lR/eR
iK2
cInv:nR&&lR!=tQ1
l03(1)/lR
iK2
cSub
l54
lR-eR
iK2
cNeg:nR){xN-lR
iK2
cRad:nR){xN
RadiansToDegrees
cY2
cDeg:nR){xN
DegreesToRadians
cY2
cSqr:nR){xN
lR*lR
iK2
cExp2:iK3(fp_exp2);lC
cRSqrt:nR){xN
l03(1)/fp_sqrt
cY2
cCot:eY2
fp_tan(lZ
cSec:eY2
fp_cos(lZ
cCsc:eY2
fp_sin(lZ
cHypot
l54
fp_hypot(xE2
cRDiv:case
cRSub:case
cDup:case
cFetch:case
cPopNMov:case
cSinCos:case
cSinhCosh:case
cNop:case
cJump:lC
cPCall:case
cFCall:c73
do_return:;
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"["
<<(&cI3)<<"]Done ConstantFolding, result: "
nJ2
DumpHashes
iI);
#endif
}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
t6{t81
yE3::set_abs(nL
bool
has_negative=!min.known||yW<l03();bool
has_positive=!nX3||yG2>l03();bool
crosses_axis=has_negative&&has_positive;cH1
xJ
newmax;if(min
lM3
nX3)newmax.set(fp_max(iF2,iG2
cQ3
crosses_axis)min.set
tQ2());e43
min
lM3
nX3)min.set(fp_min(iF2,iG2);iQ1
min.known)min.set(iF2);else
min.set(iG2;}
max=newmax;}
t81
yE3::set_neg(tS3
swap(min,max);yW=-yW;yG2=-yG2;}
xO1
IsLogicalTrueValue
eS1
yE3&p
tS2{if(nB
IsIntType
xJ::x63){if(p
cT&&p.yW>=l03(1
t41
if(!abs&&p
cQ
lV4<=l03(-1
t41}
e43
p
cT&&p.yW>=l03(0.5
t41
if(!abs&&p
cQ
lV4<=l03(-0.5
t41}
return
tD3
xO1
IsLogicalFalseValue
eS1
yE3&p
tS2{if(nB
IsIntType
xJ::x63){if(abs)return
p
cQ
l41
1);else
return
p
cT&&p
cQ&&p.yW>l03(-1)l41
1);}
e43
abs)return
p
cQ
l41
0.5);else
return
p
cT&&p
cQ&&p.yW>l03(-0.5)l41
0.5);}
}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
lS3
FUNCTIONPARSERTYPES;using
t6;t6{tU1
yE3
CalculateResultBoundaries
eS1
eQ
#ifdef DEBUG_SUBSTITUTIONS_extra_verbose
{using
lS3
FUNCTIONPARSERTYPES;yE3
tmp=CalculateResultBoundaries_do
iI)xZ1"Estimated boundaries: "
;if(tmp
cT)std::cout<<tmp.yW;else
std::cout<<"-inf"
xZ1" .. "
;if(tmp
cQ)std::cout<<tmp.yG2;else
std::cout<<"+inf"
xZ1": "
cL3
iI)xZ1
std::endl
iH
tmp;}
tU1
yE3
CodeTree
xJ::CalculateResultBoundaries_do
eS1
eQ
#endif
{iO
yL1(-fp_const_pihalf
xJ(),fp_const_pihalf
xJ());iO
pi_limits(-fp_const_pi
xJ(),fp_const_pi
xJ());iO
abs_pi_limits(y71,fp_const_pi
xJ());iO
plusminus1_limits
tQ2(-x73;using
lS3
std;switch
iI
nC
l43
nM
tree
eY1,tree
eY1)tH3
cAnd:case
cAbsAnd:case
cOr:case
cAbsOr:case
cNot:case
cN3:case
cNotNot:case
cAbsNotNot:case
cEqual:case
i91:case
cLess:case
cLessOrEq:case
cGreater:case
cGreaterOrEq:{nM
y71,l03(1));}
case
cAbs:lD
m.set_abs();cX
cLog:lD
m.iH2
fp_log);iJ2
fp_log);cX
cLog2:lD
m.iH2
fp_log2);iJ2
fp_log2);cX
cLog10:lD
m.iH2
fp_log10);iJ2
fp_log10);cX
cAcosh:lD
m.min.template
set_if<cGreaterOrEq
tR1
fp_acosh);m.y51
cGreaterOrEq
tR1
fp_acosh);cX
cAsinh
nP2
fp_asinh);eM1
set(fp_asinh);cX
cAtanh:lD
m.min.n3-1),fp_atanh);m.y51
cLess
tR1
fp_atanh);cX
cAcos:lD
nM(m
cQ&&eH3)<l03(1))?fp_acos
eH3):y71,(m
cT&&(m.yW)>=l03(-1))?fp_acos(m.yW):fp_const_pi
xJ());}
case
cAsin:lD
m.min.n3-1),fp_asin,yL1.yW);m.y51
cLess
tR1
fp_asin,yL1
tZ1;cX
cAtan
nP2
fp_atan,yL1.yW);eM1
set(fp_atan,yL1
tZ1;cX
cAtan2:{nR&&nK3
lR,xG1)yK
abs_pi_limits;}
xE3
nK3
eR,xG1)yK
yL1;}
return
pi_limits;}
case
cSin:lD
bool
x21=!m
cT||!m
cQ||eH3-m.yW)>=(yT
x21)eT
l03
min=e53
m.yW,yT
min<xG1
min
yV
l03
max=e53
i53,yT
max<xG1
max
yV
if(max<min)max
yV
bool
y61=(min<=fp_const_pihalf
xJ()&&max>=fp_const_pihalf
xJ());bool
nO1=(min<=cU&&max>=cU
cQ3
y61&&nO1)eT
if(nO1)nM
l03(-1),x32
if(y61)nM
c12
l03(1));nM
c12
x32}
case
cCos:lD
if(m
cT)m.yW+=fp_const_pihalf
xJ(cQ3
m
cQ)i53+=fp_const_pihalf
xJ();bool
x21=!m
cT||!m
cQ||eH3-m.yW)>=(yT
x21)eT
l03
min=e53
m.yW,yT
min<xG1
min
yV
l03
max=e53
i53,yT
max<xG1
max
yV
if(max<min)max
yV
bool
y61=(min<=fp_const_pihalf
xJ()&&max>=fp_const_pihalf
xJ());bool
nO1=(min<=cU&&max>=cU
cQ3
y61&&nO1)eT
if(nO1)nM
l03(-1),x32
if(y61)nM
c12
l03(1));nM
c12
x32}
case
cTan:{nM);}
case
cCeil:lD
eM1
iJ1
cFloor:lD
m
eN1
cX
cTrunc:lD
m
eN1
eM1
iJ1
cInt:lD
m
eN1
eM1
iJ1
cSinh
nP2
fp_sinh);eM1
set(fp_sinh);cX
cTanh
nP2
fp_tanh,plusminus1_limits.min);eM1
set(fp_tanh,plusminus1_limits.max);cX
cCosh:lD
if(m
cT){if(m
cQ){if(m.yW>=y71&&i53
iU{m.yW
xT}
iQ1(m.yW)<y71&&i53
iU{l03
tmp
xT
if(tmp>i53)i53=tmp;m.yW=l03(1);}
else{m.yW
xT
std::swap(m.yW,i53);}
}
e43
m.yW
iU{m
x42;m.yW=fp_cosh(m.yW);}
else{m
x42;m.yW=l03(1);}
}
}
else{m
cT=true;m.yW=l03(1
cQ3
m
cQ){m.yW=fp_cosh
eH3);m
x42;}
else
m
x42;}
cX
cIf:case
cAbsIf:{yE3
res1=tU
1));yE3
res2=tU
2)cQ3!res2
cT)res1
cT=false;iQ1
res1
cT&&(res2.yW)<res1.yW)res1.yW=res2.yW;iQ1
std::isnan(res2.yW))res1.yW=res2.yW;if(!res2
cQ)res1
x42;iQ1
res1
cQ&&(res2
tZ1>res1
tZ1
res1.yG2=res2.yG2;iQ1
std::isnan(res2
tZ1)res1.yG2=res2.yG2
iH
res1;}
case
cMin:{bool
iP=false;bool
iQ=false;yH3
x2
m=tU
a)cQ3!m
cT)iP=true
cF1
cT||(m.yW)<iP2
yW)iP2
yW=m.yW;if(!m
cQ)iQ=true
cF1
cQ||eH3)<nZ3
iP2
yG2=i53;}
if(iP)x03
iQ)x43
iH
x53
cMax:{bool
iP=false;bool
iQ=false;yH3
x2
m=tU
a)cQ3!m
cT)iP=true
cF1
cT||m.yW>iP2
yW)iP2
yW=m.yW;if(!m
cQ)iQ=true
cF1
cQ||i53>nZ3
iP2
yG2=i53;}
if(iP)x03
iQ)x43
iH
x53
cAdd:{yH3(y71,xG1
x2
item=tU
a)cQ3
item
cT)iP2
yW+=item.yW;else
x03
item
cQ)iP2
yG2+=item.yG2;else
x43;if(!x63
cT&&!x23)c73
if
iI2
cT&&x23&&iP2
yW>nZ3
std::swap
iI2.yW,nZ3
iH
x53
cMul:{eC2
Value{enum
x93{lB4,lA2,x83}
;x93
eF;l03
value;Value(x93
t):eF(t),value(0){}
Value
tQ2
v):eF(lB4),value(v){}
bool
cZ2
yZ1
eF==lA2||(eF==lB4&&value<xG1
l74
eB1*=eS1
Value&rhs){if(eF==lB4&&rhs.eF==lB4)value*=rhs
x13
eF=(cZ2)!=rhs.cZ2)?lA2:x83);}
iV2<eS1
Value&rhs
yZ1(eF==lA2&&rhs.eF!=lA2)||(eF==lB4&&(rhs.eF==x83||(rhs.eF==lB4&&value<rhs.value)));}
}
;eC2
yM1{Value
c22,c32;yM1():c22(Value::x83),c32(Value::lA2){}
void
xS2
Value
e63,const
Value&value2){e63*=value2;if(e63<c22)c22=e63;if(c32<e63)c32=e63;}
}
;yH3
tQ2(x73
x2
item=tU
a)cQ3!item
cT&&!item
cQ)nM);Value
xA3=x63
cT?Value
iI2.yW):xL2
lA2);Value
xB3=x23?Value
iI2
tZ1:xL2
x83);Value
xC3=item
cT?Value(item.yW):xL2
lA2);Value
xD3=item
cQ?Value(item
tZ1:xL2
x83);yM1
range
x33
xA3,xC3)x33
xA3,xD3)x33
xB3,xC3)x33
xB3,xD3
cQ3
range.c22.eF==Value::lB4)iP2
yW=range.c22
x13
x03
range.c32.eF==Value::lB4)iP2
yG2=range.c32
x13
x43;if(!x63
cT&&!x23)c73
if
iI2
cT&&x23&&iP2
yW>nZ3
std::swap
iI2.yW,nZ3
iH
x53
cMod:{yE3
x=tU
0));yE3
y=tU
1)cQ3
y
cQ){if(y.yG2
iU{if(!x
cT||(x.yW)<xG1
nM-y.yG2,y
tZ1;i63
y71,y
tZ1;}
e43!x
cQ||(x
tZ1
iU
nM
y.yG2,-y
tZ1;i63
y.yG2,fp_const_negativezero
xJ());}
}
i63);}
case
cPow:{xE3
eR==xG1{nM
l03(x73;}
nR&&lR==xG1{nM
y71,xG1;}
nR&&nK3
lR
nH2
nM
l03(x73;}
xE3
eR>y71&&GetEvennessInfo
iI
l8
1))==l53
eZ2
eJ2
eR;yE3
tmp=tU
0));yH3;x63
cT=true;iP2
yW=0;if(tmp
cT&&tmp.yW
iU
iP2
yW=tF3
tmp.yW,xT2;iQ1
tmp
cQ&&tmp.yG2<=xG1
iP2
yW=tF3
tmp.yG2,xT2;x43;if(tmp
cT&&tmp
cQ){x23=true;iP2
yG2=fp_max(fp_abs(tmp.yW),fp_abs(tmp
tZ1);iP2
yG2=fp_pow
iI2.yG2,xT2;}
return
x63;}
yE3
tT
eK1);yE3
p1=tU
1))yT2
p0_positivity=(p0
cT&&(p0.yW)iU?l53:(p0
cQ&&(p0
tZ1<y71?lN3
Unknown)yT2
e02=GetEvennessInfo
iI
l8
1))yT2
t0=Unknown;switch(p0_positivity
iF1
t0=l53;lC
lN3{t0=e02
c63}
yI3
switch(e02
iF1
t0=l53;lC
lN3
lC
Unknown:{xE3!tU2
eR)&&eR
iU{t0=l53;}
c73}
cJ3(t0
iF1{l03
min=y71;if(p0
cT&&p1
cT){min=tF3
p0.yW,p1.yW
cQ3
p0.yW<y71&&(!p1
cQ||p1.yG2
iU&&min
iU
min=y71;}
if(p0
cT&&p0.yW>=y71&&p0
cQ&&p1
cQ
eZ2
max=tF3
p0.yG2,p1
tZ1;if(min>max)std::swap(min,max);nM
min,max);}
nM
min,false);}
case
lN3{nM
false,fp_const_negativezero
xJ());}
yI3{c73
cS1
cNeg:lD
m.set_neg();cX
cSub:{yX
cNeg);i73
1));tmp
xC
cAdd)lW4
0));tmp.yD1
tmp2
lO3
cInv:{e12-1))lO3
cDiv:{yX
cInv);i73
1));tmp
xC
x91
yD1
tmp2
lO3
cRad:tK
x91
yB
fp_const_rad_to_deg
xJ())lO3
cDeg:tK
x91
yB
fp_const_deg_to_rad
xJ())lO3
cSqr:{e12
2))lO3
cExp:tK
cPow);tmp.yB
fp_const_e
xJ()))lW4
0)lO3
cExp2:tK
cPow);tmp
iJ3
tmp.nJ
0)lO3
cCbrt
nP2
fp_cbrt);eM1
set(fp_cbrt);cX
cSqrt:lD
if(m
cT)m.yW=(m.yW)<y71?0:fp_sqrt(m.yW
cQ3
m
cQ)i53=eH3)<y71?0:fp_sqrt
eH3);cX
cRSqrt:{e12-0.5))lO3
cHypot:i01
xsqr,ysqr,add,sqrt;xsqr.nJ
0));xsqr
iJ3
ysqr.nJ
1));ysqr
iJ3
xsqr
xC
cPow);ysqr
xC
cPow);add.yD1
xsqr);add.yD1
ysqr);add
xC
cAdd);sqrt.yD1
add);sqrt
xC
cSqrt)iH
CalculateResultBoundaries(sqrt);}
case
e93:{yX
cLog2);i73
0));tmp
eV3
tmp.yD1
tmp2)lW4
1)lO3
cCot:{yX
cTan)nT
lH
cSec:{yX
cCos)nT
lH
cCsc:{yX
cSin)nT
CalculateResultBoundaries(tmp);}
lC
cRDiv:case
cRSub:case
cDup:case
cFetch:case
cPopNMov:case
cSinCos:case
cSinhCosh:case
cNop:case
cJump:case
lI3:lC
cArg:case
cConj:case
cImag:case
cReal:case
cPolar:lC
cPCall:lC
cFCall:c73
nM);}
tU1
TriTruthValue
GetIntegerInfo
eS1
eQ
tU3
iI
nC
l43
return
isInteger
iI
eY1)?l53:y43
tH3
cFloor:case
cCeil:case
cTrunc:case
cInt:return
l53
tH3
cAnd:case
cOr:case
cNot:case
cNotNot:case
cEqual:case
i91:case
cLess:case
cLessOrEq:case
cGreater:case
cGreaterOrEq:return
l53
tH3
cIf:{TriTruthValue
a=GetIntegerInfo
iI
l8
1))yT2
b=GetIntegerInfo
iI
l8
2)cQ3
a==b)return
a
xZ2
case
cAdd:case
cMul:{for
cT1
if(GetIntegerInfo
iI
l8
a))!=l53)return
Unknown
iH
l53;}
yI3
c73
return
Unknown;}
xO1
IsLogicalValue
eS1
eQ
tU3
iI
nC
l43
return
fp_equal
iI
eY1,xG1||fp_equal
iI
eY1,l03(1))tH3
cAnd:case
cOr:case
cNot:case
cNotNot:case
cAbsAnd:case
cAbsOr:case
cN3:case
cAbsNotNot:case
cEqual:case
i91:case
cLess:case
cLessOrEq:case
cGreater:case
cGreaterOrEq:nZ
cMul:{for
cT1
if(!yO3
a))cR
return
true;}
case
cIf:case
cAbsIf:yK
yO3
1))&&yO3
2));}
yI3
c73
return
tD3}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
lS3
FUNCTIONPARSERTYPES;
#if defined(__x86_64) || !defined(FP_SUPPORT_CPLUSPLUS11_MATH_FUNCS)
# define CBRT_IS_SLOW
#endif
#if defined(DEBUG_POWI) || defined(DEBUG_SUBSTITUTIONS)
#include <cstdio>
#endif
lS3
xI1{extern
const
unsigned
char
powi_table[256];}
lS3{using
t6
n01
bool
IsOptimizableUsingPowi(long
immed,long
penalty=0){xI1
nD3
synth;synth.PushVar(lI3)cU3
bytecodesize_backup
tZ3
GetByteCodeSize();xI1::x11
immed,xI1::i11
xJ::MulSequence,synth)cU3
bytecode_grow_amount
tZ3
GetByteCodeSize()-bytecodesize_backup
iH
bytecode_grow_amount<size_t(MAX_POWI_BYTECODE_LENGTH-penalty);}
t81
ChangeIntoRootChain(n41
tree,bool
lP3,long
iL2,long
iM2){while(iM2>0)tK
cCbrt);xA1
tmp.Rehash(tI1
iN2--iM2;}
while(iL2>0)tK
cSqrt
cQ3
lP3){tmp
xC
cRSqrt);lP3=tD3
xA1
tmp.Rehash(tI1
iN2--iL2;}
if(lP3)tK
cInv);xA1
tree.iN2}
}
tU1
eC2
RootPowerTable{static
const
l03
RootPowers[(1+4)*(1+3)];}
n01
const
l03
t8(1+4)*(1+3)]={l03(1)lS
i83
i83
2*i83
2*2*2)lS
3)lS
3*2)lS
3*iS1
2*iS1
2*2*iS1
3
n43
2
n43
iS1
3*2*iS1
3*2*2*iS1
3*3
n43
3*2
n43
3*iS1
3*3*2*iS1
3*3*2*2*2*2)}
;eC2
PowiResolver{static
const
unsigned
MaxSep=4;static
xG3
i93=5;typedef
int
e73;typedef
long
xM3;typedef
long
tW;eC2
c42{c42():n_int_sqrt(0),n_int_cbrt(0),sep_list(),n11(0){}
int
n_int_sqrt;int
n_int_cbrt;int
sep_list[MaxSep];tW
n11;}
n01
static
c42
CreatePowiResult
tQ2
xT2{c42
x63;e73
t9=FindIntegerFactor(xT2;if(t9==0){
#ifdef DEBUG_POWI
iO2"no factor found for %Lg\n"
,n51
xT2;
#endif
return
x63;}
e83=y81
exponent,t9);xM3
t82=EvaluateFactorCost(t9,0,0,0)+c8
e83);int
iA3=0;int
iB3=0;int
xF3=0;
#ifdef DEBUG_POWI
iO2"orig = %Lg\n"
,n51
xT2;iO2"plain factor = "
l14" %ld\n"
,(int)t9,(long)t82);
#endif
for
tL1
n_s=0;n_s<MaxSep;++n_s){int
x8=0;xM3
yN1=t82;e73
yX1=t9;for(int
s=1;s<i93*4;++s){
#ifdef CBRT_IS_SLOW
if(s>=i93)break;
#endif
int
n_sqrt=s%i93;int
n_cbrt=s/i93;if(n_sqrt+n_cbrt>4)n33
l03
lJ1=exponent;lJ1-=t8
s];iK1=FindIntegerFactor(lJ1
cQ3
y42!=0){tW
xO=y81
lJ1,y42);xM3
cost=EvaluateFactorCost(y42,iA3+n_sqrt,iB3+n_cbrt,xF3+1)+c8
xO);
#ifdef DEBUG_POWI
iO2"Candidate sep %u (%d*sqrt %d*cbrt)factor = "
l14" %ld (for %Lg to %ld)\n"
,s,n_sqrt,n_cbrt,y42,(long)cost,n51
lJ1,(long)xO);
#endif
if(cost<yN1){x8=s;yX1=y42;yN1=cost;}
}
}
if(!x8)break;
#ifdef DEBUG_POWI
iO2"CHOSEN sep %u (%d*sqrt %d*cbrt)factor = "
l14" %ld, exponent %Lg->%Lg\n"
,x8,x8%i93,x8/i93,yX1,yN1,n51(xT2,n51
iY3-t8
x8]));
#endif
iP2
sep_list[n_s]=x8;exponent-=t8
x8];iA3+=x8%i93;iB3+=x8/i93;t82=yN1;t9=yX1;xF3+=1;}
e83=y81
exponent,t9);
#ifdef DEBUG_POWI
iO2"resulting exponent is %ld (from exponent=%Lg, best_factor=%Lg)\n"
,e83,n51
exponent,n51
t9);
#endif
while(t9%2==0){++iP2
n_int_sqrt;t9/=2;}
while(t9%3==0){++iP2
n_int_cbrt;t9/=3;}
return
x63;}
private:static
xM3
c8
tW
xO){static
std::map
e22
iE;if(xO<0){xM3
cost=22
iH
cost+c8-xO);}
std::map
e22::yA3
i=iE.y52
xO
cQ3
i!=iE
c71
xO)return
i
eH2;std::pair
e22
x63(xO,0.0);xM3&cost=iP2
second;while(xO>1){int
y42=0;if(xO<256){y42=xI1::powi_table[xO];if(y42&128)y42&=127;else
y42=0;if(y42&64)y42=-(y42&63)-1;}
if(y42){cost+=c8
y42);xO/=y42;tF2
if(!(xO&1)){xO/=2;cost+=6;}
else{cost+=7;xO-=1;}
}
iE.y63,x63)iH
cost;}
cB1
tW
y81
yO1,iK1)yK
makeLongInteger(value*l03(y42));}
cB1
bool
yP1
yO1,iK1
eZ2
v=value*l03(y42)iH
isLongInteger(v);}
cB1
e73
FindIntegerFactor(yO1){iK1=(2*2*2*2);
#ifdef CBRT_IS_SLOW
#else
y42*=(3*3*3);
#endif
e73
x63=0;if(yP1
value,y42)){x63=y42;while((y42%2)==0&&yP1
value,y42/2))x63=y42/=2;while((y42%3)==0&&yP1
value,y42/3))x63=y42/=3;}
#ifdef CBRT_IS_SLOW
if
iI2==0){if(yP1
value,3
yL3
3;}
#endif
return
x63;}
static
int
EvaluateFactorCost(int
y42,int
s,int
c,int
nmuls){xG3
xH3=6;
#ifdef CBRT_IS_SLOW
xG3
t92=25;
#else
xG3
t92=8;
#endif
int
x63=s*xH3+c*t92;while(y42%2==0){y42/=2;x63+=xH3;}
while(y42%3==0){y42/=3;x63+=t92;}
x63+=nmuls
iH
x63;}
}
;}
t6{xO1
CodeTree
xJ::RecreateInversionsAndNegations(bool
prefer_base2){bool
changed=false
nA3
0
eP3<l71++a)if(lY1.RecreateInversionsAndNegations(prefer_base2))yQ1
if(changed){exit_changed:Mark_Incompletely_Hashed(iG
switch(l42{case
cMul:{eN
nA2
iJ
nQ2,cU1;if(true){bool
nP1=false;l03
xF2=0
nA3
l71
a
nS
c52
0)e42
tL
e52)){nP1=true;xF2=tL
1)eY1
c63}
}
if(nP1
eZ2
immeds=1.0
nA3
l71
a
nS
e62)){immeds*=powgroup
eY1;lR3
a);}
}
for
nD2
a=l71
a-->0;){n41
powgroup=lY1;if(powgroup
c52
0)e42
tL
e52
e5&log2=tL
0);log2.lB1
log2
xC
e93);log2.yB
tF3
immeds,l03(1)/xF2)));log2
xJ2
c73}
}
}
for
nD2
a=l71
a
nS
c52
e52)){xN2
exp_param=tL
1);l03
eJ2
exp_param
eY1;if(e11,l03(-1))){lB1
nA2.push_back(lY1
cG3
lR3
a);}
iQ1
exponent<y71&&tU2
exponent
e5
iR;iR
xC
cPow);iR
cG
tL
0));iR.yB-xT2);iR
xJ2
nA2.push_back(iR)l31}
iQ1
powgroup
e42!nQ2
y12)){nQ2=tL
0)l31
iQ1
powgroup
nC==e93&&!cU1
y12)){cU1=powgroup
l31}
if(!nA2
i13){changed=true
iJ
iC1;iC1
eV3
iC1.SetParamsMove(nA2);iC1
i23
y91
cMul
eL2
SetParamsMove(t1
if
e32
e62)&&fp_equal
e32
eY1
nH2
nB2
cInv)t2
iC1);}
else{if
e32
nY2>=iC1
nY2){nB2
cDiv
iC3
t2
iC1);}
else{nB2
cRDiv)t2
iC1
iC3;}
}
}
if(nQ2
y12
e5
y91
l42;tB2
SetParamsMove(t1
while
e32.RecreateInversionsAndNegations(prefer_base2))tB2
FixIncompleteHashes();nB2
e93)t2
nQ2
iC3;yQ1}
if(cU1
y12
e5
y91
cMul
eL2
yD1
cU1
l8
1)eL2
AddParamsMove(t1
while
e32.RecreateInversionsAndNegations(prefer_base2))tB2
FixIncompleteHashes();DelParams();nB2
e93)t2
cU1
l8
0)iC3;yQ1
cS1
cAdd:{eN
iQ2
nA3
cL1
eA3
cMul){nG2
yA1:iJ&xN3=e92
for
nD2
b=tB2
l71
b-->0;){if
e32
l8
b).lV1
y42=xN3
l8
b)eY1;if
lE3
y42
nL2
yA1;}
tB2
lB1
tB2
lR3
b
eK2
iQ1
nK3
y42,l03(-2)))xB
yA1;}
tB2
lB1
tB2
lR3
b
eL2
yB
l03(2))eK2}
}
if(t3){tB2
tA
xN3
tS1
iQ1
eA3
cDiv&&!IsIntType
xJ::x63){nG2
yB1:iJ&iC1=e92
if(iC1
l8
0).xG2
lE3
iC1
l8
0)eY1
nL2
yB1;}
iC1.lB1
iC1.lR3
0);iC1
xC
cInv
eK2}
if(t3)xB
yB1;}
iC1.tA
iC1
tS1
iQ1
eA3
cRDiv&&!IsIntType
xJ::x63){nG2
xB1:iJ&iC1=e92
if(iC1
l8
1).xG2
lE3
iC1
l8
1)eY1
nL2
xB1;}
iC1.lB1
iC1.lR3
1);iC1
xC
cInv
eK2}
if(t3)xB
xB1;}
iC1.tA
iC1
tS1
if(!iQ2
i13){
#ifdef DEBUG_SUBSTITUTIONS
iO2"Will make a Sub conversion in:\n"
);fflush(stdout);iW
#endif
CodeTree
xJ
c62;c62
xC
cAdd);c62.SetParamsMove(iQ2);c62
i23
cV1;cV1
xC
cAdd);cV1.SetParamsMove(lD2));cV1
xJ2
if(cV1
e62)&&nK3
cV1
eY1,xG1){nB2
cNeg);eG);}
e43
cV1
nY2==1){nB2
cRSub);eG)eB3}
iQ1
c62
nC==cAdd){nB2
cSub)eB3
eG
l8
0))nA3
1
eP3<c62.l71++a)i01
tC2;tC2
xC
cSub);tC2.SetParamsMove(lD2));tC2.yQ2)t2
tC2);eG
l8
a));}
}
else{nB2
cSub)eB3
eG);}
}
#ifdef DEBUG_SUBSTITUTIONS
iO2"After Sub conversion:\n"
);fflush(stdout);iW
#endif
cS1
cPow:{xN2
p0=GetParam(0
l21
p1=GetParam(1
cQ3
p1.xG2
xE1!=y71&&!isInteger
xE1)){eJ
c42
r=eJ
CreatePowiResult(fp_abs
xE1)cQ3
r.n11!=0){bool
lB2=false;if
xE1<y71&&r.sep_list[0]==0&&r.n_int_sqrt>0){lB2=true;}
#ifdef DEBUG_POWI
iO2"Will resolve powi %Lg as powi(chain(%d,%d),%ld)"
,n51
fp_abs
xE1),r.n_int_sqrt,r.n_int_cbrt,r.n11)tK1
n=0;n<eJ
MaxSep;++n){if(r
yE1==0)break;int
n_sqrt=r
yE1%eJ
i93;int
n_cbrt=r
yE1/eJ
i93;iO2"*chain(%d,%d)"
,n_sqrt,n_cbrt);}
iO2"\n"
);
#endif
CodeTree
xJ
eA2=GetParam(0)iJ
c72=eA2;c72.lB1
ChangeIntoRootChain(c72,lB2,r.n_int_sqrt,r.n_int_cbrt);c72
i23
pow;if(r.n11!=1){pow
xC
cPow);pow.yD1
c72);pow.yB
l03(r.n11)));}
else
pow.swap(c72)iJ
mul;mul
eV3
mul.yD1
pow)tK1
n=0;n<eJ
MaxSep;++n){if(r
yE1==0)break;int
n_sqrt=r
yE1%eJ
i93;int
n_cbrt=r
yE1/eJ
i93
iJ
tD2=eA2;tD2.lB1
ChangeIntoRootChain(tD2,false,n_sqrt,n_cbrt);tD2
xJ2
mul.yD1
tD2);}
if
xE1<y71&&!lB2){mul
xJ2
nB2
cInv);nE1
0,mul);lR3
1);}
else{nB2
cMul);SetParamsMove(mul.lD2));}
#ifdef DEBUG_POWI
iW
#endif
yQ1
c73}
}
if(GetOpcode()==cPow&&(!p1.cG2!isLongInteger
xE1)||!IsOptimizableUsingPowi
xJ(makeLongInteger
xE1)))){if(p0
e62)&&p0
eY1>l03(0.0)){if(prefer_base2
eZ2
c82=fp_log2(p0
eY1);if
lE3
c82
nH2
lR3
0);}
else{n0
e31
c82)cA3
cG
p1
n91
cC}
nB2
cExp2);yQ1}
else{l03
c82=fp_log(p0
eY1);if
lE3
c82
nH2
lR3
0);}
else{n0
e31
c82)cA3
cG
p1
n91
cC}
nB2
cExp);yQ1}
}
iQ1
GetPositivityInfo(p0)==l53){if(prefer_base2)i01
log;log
xC
cLog2);log
cG
p0);log
xJ2
n0
p1
cA3.yD1
log
n91
nB2
cExp2);cC
yQ1}
else
i01
log;log
xC
cLog);log
cG
p0);log
xJ2
n0
p1
cA3.yD1
log
n91
nB2
cExp);cC
yQ1}
}
cS1
cDiv:{if(GetParam(0)e62)&&nK3
GetParam(0)eY1
nH2
nB2
cInv);lR3
0);}
c73
yI3
c73
if(changed)goto
exit_changed
iH
changed;}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
lS3
FUNCTIONPARSERTYPES;lS3{using
t6;class
i82{size_t
nQ1
cU3
eK
cU3
eL
cU3
lK1
cU3
t4
cU3
t5
cU3
nA1;eT3
i82():nQ1(0),eK(0),eL(0),lK1(0),t4(0),t5(0),nA1(0){}
void
iH3
OPCODE
op){nQ1+=1;i72
cCos)++eK;i72
cSin)++eL;i72
cSec)++eK;i72
cCsc)++eL;i72
cTan)++lK1;i72
cCot)++lK1;i72
cSinh)++t5;i72
cCosh)++t4;i72
cTanh)++nA1;}
size_t
GetCSEscore()const{size_t
x63=nQ1
iH
x63;}
int
NeedsSinCos()const{bool
yC1=(nQ1==(eK+eL+lK1)cQ3(lK1&&(eL||eK))||(eL&&eK)){if(yC1)return
1
iH
2;}
return
0;}
int
NeedsSinhCosh()const{bool
yC1=(nQ1==(t4+t5+nA1)cQ3(nA1&&(t5||t4))||(t5&&t4)){if(yC1)return
1
iH
2;}
return
0;}
size_t
MinimumDepth()const{size_t
n_sincos=std::min(eK,eL)cU3
n_sinhcosh=std::min(t4,t5
cQ3
n_sincos==0&&n_sinhcosh==0)return
2
iH
1;}
}
n01
class
TreeCountType:public
std::multimap<fphash_t,std::pair<i82,CodeTree
xJ> >{}
xU3
FindTreeCounts(tT1&eT1,xN2
tree,OPCODE
xH2,bool
skip_root=false){cZ
i=eT1.lower_bound
iI.GetHash()cQ3!skip_root){bool
found=false;for(;i!=eT1
c71
tree.GetHash();++i){if
iI
xE
i
eH2
t93)){i
eH2.first.iH3
xH2);found=true
c63}
}
if(!found){i82
count;count.iH3
xH2);eT1.y63,std::make_pair
iI.GetHash(),std::make_pair
eG3,tree)));}
}
lN1
FindTreeCounts(eT1,nB3,tY2);}
eC2
yM{bool
BalanceGood;bool
FoundChild;}
n01
yM
lL1
xN2
root,xN2
child){if(root
xE
child)){yM
x63={true,true}
iH
x63;}
yM
x63={true,false}
;if(root
nC==cIf||root
nC==tC3{yM
cond=lL1
root
l8
0),child);yM
xW=lL1
root
l8
1),child);yM
y5=lL1
root
l8
2),child
cQ3
cond
yN||xW
yN||y5
yN){x63
yN=true;}
x63
eH=((xW
yN==y5
yN)||c91&&(cond
eH||(xW
yN&&y5
yN))&&(xW
eH||c91&&(y5
eH||c91;}
else{bool
iL1=false;bool
nR1=false;for
nD2
b=root
yS1,eQ3<b;++a){yM
tmp=lL1
root
l8
a),child
cQ3
tmp
yN)x63
yN=true;if(tmp
eH==false)iL1=true;iQ1
tmp
yN)nR1=true;}
if(iL1&&!nR1)x63
eH=tD3
return
x63;}
xO1
n13
n03
iD3
xN2
tree,const
xI1
nD3&synth,const
tT1&eT1){for
nD2
b=tree
yS1,eQ3<b;++a){xN2
leaf=nB3;cZ
synth_it;y82
tT1::const_iterator
i=eT1.y93
i!=eT1
lM4;++i){if(i->first!=leaf.GetHash())n33
const
i82&occ
x52
first
cU3
score=occ.GetCSEscore(l21
candidate
x52
second;if(nZ2(candidate))n33
if(leaf
nY2<occ.MinimumDepth())n33
if(score<2)n33
if(lL1
iD3
leaf)eH==false)continue
lV3
if(n13(iD3
leaf,synth,eT1
t41}
return
tD3
xO1
nR2
n03
y73,xN2
expr){for
nC2
y73
n93
expr
t41
for
nC2
nR2(y73
l8
a),expr
yL3
true
iH
tD3
xO1
GoodMomentForCSE
n03
y73,xN2
expr){if(y73
nC==cIf)return
true;for
nC2
y73
n93
expr
t41
size_t
iR2=0;for
nC2
nR2(y73
l8
a),expr))++iR2
iH
iR2!=1;}
}
t6{tU1
size_t
CodeTree
xJ::SynthCommonSubExpressions(xI1::yF1
const{if(GetParamCount()==0)return
0
cU3
stacktop_before
tZ3
GetStackTop();tT1
eT1;FindTreeCounts(eT1,*this,GetOpcode(),true);for(;;){size_t
c92=0;
#ifdef DEBUG_SUBSTITUTIONS_CSE
std::cout<<"Finding a CSE candidate, root is:"
<<std::endl;DumpHashes(*this);
#endif
cZ
cs_it(eT1
lM4);for(cZ
j=eT1.y93
j!=eT1
lM4;){cZ
i(j++);const
i82&occ
x52
first
cU3
score=occ.GetCSEscore(l21
tree
x52
second;
#ifdef DEBUG_SUBSTITUTIONS_CSE
std::cout<<"Score "
<<score<<":\n"
<<std::flush;DumpTreeWithIndent
iI);
#endif
if(nZ2
iI))xS
if
iI
nY2<occ.MinimumDepth())xS
if(score<2)xS
if(lL1*this,tree)eH==false)xS
if(n13(*this,tree,synth,eT1)){tF2
if(!GoodMomentForCSE(*this,tree))xS
score*=tree
nY2;if(score>c92){c92=score;cs_it=i;}
}
if(c92<=0){
#ifdef DEBUG_SUBSTITUTIONS_CSE
std::cout<<"No more CSE candidates.\n"
<<std::flush;
#endif
c73
xN2
tree=cs_it
eH2
t93;
#ifdef DEBUG_SUBSTITUTIONS_CSE
std::cout<<l44"Common Subexpression:"
cL3
xJ
iI)xZ1
std::endl;
#endif
#if 0
int
n21=occ.NeedsSinCos();int
i3=occ.NeedsSinhCosh()iJ
iS2,iT2,cA2,cB2;if(n21){iS2
tE2
iS2
xC
cSin);iS2
xJ2
iT2
tE2
iT2
xC
cCos);iT2
xJ2
if(nZ2(iS2)||nZ2(iT2))iE3==2){eU1
tF2
n21=0;}
}
if(i3){cA2
tE2
cA2
xC
cSinh);cA2
xJ2
cB2
tE2
cB2
xC
cCosh);cB2
xJ2
if(nZ2(cA2)||nZ2(cB2)){if(i3==2){eU1
tF2
i3=0;}
}
#endif
tree.SynthesizeByteCode(synth,false);eU1
#ifdef DEBUG_SUBSTITUTIONS_CSE
synth.template
Dump<0>()xZ1"Done with Common Subexpression:"
cL3
xJ
iI)xZ1
std::endl;
#endif
#if 0
if(n21)iE3==2||i3){tH1}
lY3
cSinCos,1,2
iD1
iS2,1
iD1
iT2,0);}
if(i3)iE3)tH1
if(i3==2){tH1}
lY3
cSinhCosh,1,2
iD1
cA2,1
iD1
cB2,0);}
#endif
}
return
synth.xL
stacktop_before;}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
tU1
lU1
xJ::iU2{using
t6;CopyOnWrite()iJ
tree;tree.GenerateFrom(*mData);FPoptimizer_Optimize::ApplyGrammars
iI);std::vector<unsigned>eC3;std::vector
xJ
immed
cU3
stacktop_max=0;tree.SynthesizeByteCode(eC3,immed,stacktop_max
cQ3
mData->mStackSize!=stacktop_max){mData->mStackSize=unsigned(stacktop_max);
#if !defined(FP_USE_THREAD_SAFE_EVAL) && \
    !defined(FP_USE_THREAD_SAFE_EVAL_WITH_ALLOCA)
mData->mStack.eU3
stacktop_max);
#endif
}
mData->mByteCode.swap(eC3);mData->mImmed.swap(immed);}
#define FUNCTIONPARSER_INSTANTIATE_EMPTY_OPTIMIZE(type) tW1>lU1<type>::iU2{}
#ifdef FP_SUPPORT_MPFR_FLOAT_TYPE
iL3(MpfrFloat)
#endif
#ifdef FP_SUPPORT_GMP_INT_TYPE
iL3(GmpInt)
#endif
#ifdef FP_SUPPORT_COMPLEX_DOUBLE_TYPE
iL3(std::complex<double>)
#endif
#ifdef FP_SUPPORT_COMPLEX_FLOAT_TYPE
iL3(std::complex<float>)
#endif
#ifdef FP_SUPPORT_COMPLEX_LONG_DOUBLE_TYPE
iL3(std::complex<long
double>)
#endif
#define FUNCTIONPARSER_INSTANTIATE_OPTIMIZE(type) template lU1<type>::iU2;
#ifndef FP_DISABLE_DOUBLE_TYPE
iM3(double)
#endif
#ifdef FP_SUPPORT_FLOAT_TYPE
iM3(float)
#endif
#ifdef FP_SUPPORT_LONG_DOUBLE_TYPE
iM3(long
double)
#endif
#ifdef FP_SUPPORT_LONG_INT_TYPE
iM3(long)
#endif
#endif // FP_SUPPORT_OPTIMIZER

#endif
