(* ::Package:: *)

BeginPackage["pesudoRiemannian"]
(*https://mathematica.stackexchange.com/questions/8895/how-to-calculate-scalar-curvature-ricci-tensor-and-christoffel-symbols-in-mathem*)

ChristoffelSymbol1st::usage = ""

ChristoffelSymbol2nd::usage = ""

RicciTensor::usage = ""

RicciScalar::usage = ""

EinsteinTensor::usage = ""

Begin["Private`"]
(*Riemannian formalism*)

(*Intrinsic functions*)
ChristoffelSymbol1st[gLow_,xxx_]:=1/2 (#2-#1[#2,{2,3,1}]+#1[#2,{3,1,2}]&)[Transpose,D[gLow,{xxx}]];
ChristoffelSymbol2nd[gLow_,xxx_]:=TensorContract[Inverse[gLow]\[TensorProduct]ChristoffelSymbol1st[gLow,xxx],{{2,3}}];

RiemannTensor[gLow_,xxx_]:=((Transpose[#,{1,2,4,3}]-#&)[D[#1, {#2}]]+(#1[#2,{1,3,2,4}]-#1[#2,{1,4,2,3}]&)[Transpose,TensorContract[#1\[TensorProduct]#1,{{2,4}}]]&)[ChristoffelSymbol2nd[gLow,xxx],xxx];

RicciTensor[gLow_,xxx_]:=TensorContract[RiemannTensor[gLow,xxx],{{1,3}}];

RicciScalar[gLow_,xxx_]:=TensorContract[Inverse[gLow]\[TensorProduct]RicciTensor[gLow,xxx],{{1,3},{2,4}}];

EinsteinTensor[gLow_,xxx_]:=RicciTensor[gLow,xxx]-2 (RicciScalar[gLow,xxx]gLow)/Length@gLow;
End[] (*"Private`"*)

EndPackage[] (*"pesudoRiemannian"*)
