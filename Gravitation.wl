(* ::Package:: *)

(********1*********2*********3*********4*********5*********6*********7*****
 * GREAT.m = General Relativity, Einstein & All That 4 Mathematica
 * Package written by Yi-Fan Wang: yfwang@thp.uni-koeln.de
 *                                 Universit\[ADoubleDot]t zu K\[ODoubleDot]ln
 *                                 thp.uni-koeln.de/~yfwang
 * Based on
 * https://mathematica.stackexchange.com/questions/8895/how-to-calculate-scalar-curvature-ricci-tensor-and-christoffel-symbols-in-mathem
 * and the package
 * GREAT.m = General Relativity, Einstein & All That 4 Mathematica
 * by Tristan Hubsch: thubsch@howard.edu
 *                    Howard University Physics
 *                    string.howard.edu/~tristan/
 * and the package
 * "EinsteinTensor.m"
 * by Pekka Janhunen: pjanhune@finsun.csc.fi
 *    Finnish Meteorological Institute Geophysics Dept.
 **************************************************************************)

BeginPackage["Gravitation`"]

(*Riemannian formalism*)

(*Functions for intrinsic quantities*)

InverseMetricSeries

ChristoffelSymbol1st::usage="ChristoffelSymbol1st[g,x], with g components of the (0,2) metric tensor (D\[Times]D Matrix) and x the coordinates (D List), gives the Christoffel symbol of the first kind (three lower indices)."

ChristoffelSymbol1stSeries

ChristoffelSymbol2nd::usage="ChristoffelSymbol2nd[g,x], with g components of the (0,2) metric tensor (D\[Times]D Matrix) and x the coordinates (D List), gives the Christoffel symbol of the second kind (1st upper, two lower indices)."

ChristoffelSymbol2ndSeries

LeviCivitaDVector::usage="[gLow_List,xxx_List,vec_List]:=D[vec,xxx]+ChristoffelSymbol2nd[gLow,xxx].vec;"

(*LeviCivitaDVectorSeries*)

LeviCivitaDCovector::usage="[gLow_List,xxx_List,covec_List]:=D[covec,xxx]-covec.ChristoffelSymbol2nd[gLow,xxx];"

(*LeviCivitaDCovectorSeries*)

RiemannTensor::useage="RiemannTensor[g,x], with g components of the (0,2) metric tensor (D\[Times]D Matrix) and x the coordinates (D List), gives the (1,3) Riemann curvature tensor (D\[Times]D\[Times]D\[Times]D List)."

(*RiemannTensorSeries*)

RicciTensor::usage="RicciTensor[g,x], with g components of the (0,2) metric tensor (D\[Times]D Matrix) and x the coordinates (D List), gives the (0,2) Ricci tensor (D\[Times]D Matrix)."

(*RicciTensorSeries*)

RicciScalar::usage="RicciScalar[g,x], with g components of the (0,2) metric tensor (D\[Times]D Matrix) and x the coordinates (D List), gives the Ricci scalar."

(*RicciScalarSeries*)

EinsteinTensor::usage="EinsteinTensor[g,x] with g components of the (0,2) metric tensor (D\[Times]D Matrix) and x the coordinates (D List) gives the (0,2) Einstein tensor (D\[Times]D Matrix)."

(*EinsteinTensorSeries*)

(*Functions for extrinsic quantities*)

(*Lapse etc known*)

MetricTensor4::usage="MetricTensorDecomposed[lapse,shift,hLow] with lapse a number, shift a vector (d List), hLow components of the (0,2) metric tensor (d\[Times]d Matrix), gives a D\[Times]D (0,2) metric (D\[Times]D Matrix, D=d+1)."(* It has only been test for for (d+1) decomposition."*)

MetricTensor3Vector::usage="MetricTensorSurface[g,nVec,s] with g components of the (0,2) metric tensor (D\[Times]D Matrix), n unit normal vector (D List), sC=1 for space-like hypersurface and -1 for time-like hypersurface, gives a d\[Times]d (0,2) metric (D\[Times]D Matrix, D=d+1)."

MetricTensor3Covector::usage="MetricTensorSurface[g,nCvt,sV] with g components of the (0,2) metric tensor (D\[Times]D Matrix), n unit normal covector (D List), sC=1 for space-like hypersurface and -1 for time-like hypersurface, gives a d\[Times]d (0,2) metric (D\[Times]D Matrix, D=d+1)."

MetricTensor3Index::usage="MetricTensor3Index[g,i,sC] with g components of the (0,2) metric tensor (D\[Times]D Matrix), i index (Integer), sC=1 for space-like hypersurface and -1 for time-like hypersurface, gives a d\[Times]d (0,2) metric (D\[Times]D Matrix, D=d+1)."

(*Suppose the metric is d+1 decomposed; nUpp=n^\[Mu]. I believe these are to be improved.*)

NormalCovector::usage="NormalCovector[g,i,sC] with g components of the (0,2) metric tensor (D\[Times]D Matrix), i coordinate index (Integer), sC=1 for space-like hypersurface and -1 for time-like hypersurface, gives the unit covector orthogonal to the hypersurface \!\(\*SuperscriptBox[\(x\), \(i\)]\)==C."

NormalVector::usage="NormalVector[g,i,sC] with g components of the (0,2) metric tensor (D\[Times]D Matrix), i coordinate index (Integer), sC=1 for space-like hypersurface and -1 for time-like hypersurface, gives the unit vector orthogonal to the hypersurface \!\(\*SuperscriptBox[\(x\), \(i\)]\)==C."

KruemmungsTensor4Index::usage="KruemmungsTensor4Index[g,x,sB,i,sC] with g components of the bulk (0,2) metric tensor (D\[Times]D Matrix), x the coordinates (D List), sB the signature of the bulk metric, i coordinate index and sC=1 for space-like hypersurface and -1 for time-like hypersurface, gives the (0,2) second fundamental form of the hypersurface \!\(\*SuperscriptBox[\(x\), \(i\)]\)==C (D\[Times]D Matrix)."

KruemmungsScalar4Index::usage="KruemmungsScalar4Index[g,x,sB,i,sC] with g components of the bulk (0,2) metric tensor (D\[Times]D Matrix), x the coordinates (D List), sB the signature of the bulk metric, i coordinate index and sC=1 for space-like hypersurface and -1 for time-like hypersurface, gives the trace of the second fundamental form of the hypersurface \!\(\*SuperscriptBox[\(x\), \(i\)]\)==C."

(*Gravitation and General Relativity*)

LagrangianScalarEinstein::usage="[gLow_List,xxx_List,coeff_]:=,coeff Tr[Inverse[gLow].(TensorContract[#,{{1,5},{2,4}}]-TensorContract[#,{{1,6},{4,5}}]&[#\[TensorProduct]#&[ChristoffelSymbol2nd[gLow,xxx]]])]"

LagrangianScalarEinsteinHilbert::usage="[gLow_List,xxx_List,coeff_]:=coeff RicciScalar[gLow,xxx];"

LagrangianScalarArnowittDeserMisnerIndex::usage="[gLow_List,xxx_List,sgnBulk_,ind_Integer,sgnCell_,coeff_]:=coeff(Tr[#1.#2.#1.#2]&[Inverse[gLow],KruemmungsTensor4Index[gLow,xxx,sgnBulk,ind,sgnCell]]-KruemmungsScalar4Index[gLow,xxx,sgnBulk,ind,sgnCell\!\(\*SuperscriptBox[\(]\), \(2\)]\)+RicciScalar[Drop[MetricTensor3Index[gLow,ind,sgnCell],{ind},{ind}]]);"

LagrangianScalarArnowittDeserMisnerIndexSeries

Begin["Private`"]

(*Riemannian formalism*)

(*Functions for intrinsic quantities*)

InverseMetricSeries[gLow_List,eps_List]:=Series[Inverse@gLow,eps];

ChristoffelSymbol1st[gLow_List,xxx_List]:=1/2 (#2-#1[#2,{2,3,1}]+#1[#2,{3,1,2}]&)[Transpose,D[gLow,{xxx}]];

ChristoffelSymbol1stSeries[gLow_List,xxx_List,eps_List]:=Series[ChristoffelSymbol1st[gLow,xxx],eps];

ChristoffelSymbol2nd[gLow_List,xxx_List]:=Inverse[gLow].ChristoffelSymbol1st[gLow,xxx];

ChristoffelSymbol2ndSeries[gLow_List,xxx_List,eps_List]:=InverseMetricSeries[gLow,eps].ChristoffelSymbol1stSeries[gLow,xxx,eps];

LeviCivitaDVector[gLow_List,xxx_List,vec_List]:=D[vec,{xxx}]+ChristoffelSymbol2nd[gLow,xxx].vec;

LeviCivitaDVectorSeries[gLow_List,xxx_List,vec_List,eps_List]:=D[vec,{xxx}]+ChristoffelSymbol2ndSeries[gLow,xxx,eps].vec;

LeviCivitaDCovector[gLow_List,xxx_List,cvt_List]:=D[cvt,{xxx}]-cvt.ChristoffelSymbol2nd[gLow,xxx];

LeviCivitaDCovectorSeries[gLow_List,xxx_List,cvt_List,eps_List]:=D[cvt,{xxx}]-cvt.ChristoffelSymbol2ndSeries[gLow,xxx,eps];

RiemannTensor[gLow_List,xxx_List]:=((Transpose[#,{1,2,4,3}]-#&)[D[#1, {#2}]]+(#1[#2,{1,3,2,4}]-#1[#2,{1,4,2,3}]&)[Transpose,TensorContract[#1\[TensorProduct]#1,{{2,4}}]]&)[ChristoffelSymbol2nd[gLow,xxx],xxx];

RiemannTensorSeries[gLow_List,xxx_List,eps_List]:=((Transpose[#,{1,2,4,3}]-#&)[D[#1, {#2}]]+(#1[#2,{1,3,2,4}]-#1[#2,{1,4,2,3}]&)[Transpose,TensorContract[#1\[TensorProduct]#1,{{2,4}}]]&)[ChristoffelSymbol2ndSeries[gLow,xxx,eps],xxx];

RicciTensor[gLow_List,xxx_List]:=TensorContract[RiemannTensor[gLow,xxx],{{1,3}}];

RicciTensorSeries[gLow_List,xxx_List,eps_List]:=TensorContract[RiemannTensorSeries[gLow,xxx],{{1,3}}];

RicciScalar[gLow_List,xxx_List]:=Tr[Inverse[gLow].RicciTensor[gLow,xxx]];

RicciScalarSeries[gLow_List,xxx_List,eps_List]:=Tr[InverseMetricSeries[gLow,eps].RicciTensorSeries[gLow,xxx,eps]];

EinsteinTensor[gLow_List,xxx_List]:=RicciTensor[gLow,xxx]-2RicciScalar[gLow,xxx]gLow/Length@gLow;

EinsteinTensorSeries[gLow_List,xxx_List,eps_List]:=RicciTensorSeries[gLow,xxx,eps]-2RicciScalarSeries[gLow,xxx,eps]gLow/Length@gLow;

(*Functions for extrinsic quantities*)

(*Lapse etc known*)

MetricTensor4[lapse_,shift_,hLow_]:=ArrayFlatten[{{shift.#-lapse^2,{#}},{{#}\[Transpose],hLow}}]&[hLow.shift];

MetricTensor3Covector[gLow_List,norCvt_List,sgnCell_]:=gLow+sgnCell*(norCvt\[TensorProduct]norCvt);

MetricTensor3Vector[gLow_List,norVec_List,sgnCell_]:=MetricTensor3Covector[gLow,gLow.norVec,sgnCell];

MetricTensor3Index[gLow_List,ind_Integer,sgnCell_]:=MetricTensor3Covector[gLow,NormalCovector[gLow,ind,sgnCell],sgnCell];

(*Suppose the metric is d+1 decomposed; nUpp=n^\[Mu]. I believe these are to be improved.*)

NormalCovector[gLow_List,ind_Integer,sgnCell_]:=Array[If[#==ind,-sgnCell/Sqrt[-sgnCell Inverse[gLow][[ind,ind]]],0]&,Length@gLow];

NormalVector[gLow_List,ind_Integer,sgnCell_]:=-sgnCell/Sqrt[-sgnCell*#[[ind,ind]]] #[[ind]]&[Inverse[gLow]];

KruemmungsTensor4Index[gLow_List,xxx_List,sgnBulk_,ind_Integer,sgnCell_]:=MetricTensor3Covector[gLow,#,sgnCell].Inverse[gLow].LeviCivitaDCovector[gLow,xxx,#]&[NormalCovector[gLow,ind,sgnCell]];

(*KruemmungsTensor4IndexSeries[gLow_List,xxx_List,sgnBulk_,ind_Integer,sgnCell_,eps_]:=
	MetricTensor3Covector[gLow,#,sgnCell].InverseMetricSeries[gLow,eps].LeviCivitaDCovectorSeries[gLow,xxx,#,eps]&[NormalCovector[gLow,ind,sgnCell]];*)

KruemmungsScalar4Index[gLow_List,xxx_List,sgnBulk_,ind_Integer,sgnCell_]:=1/#1 Div[#1*#2,xxx]&[Sqrt[sgnBulk*Det@gLow],NormalVector[gLow,ind,sgnCell]];

(*KruemmungsScalar4IndexSeries[gLow_List,xxx_List,sgnBulk_,ind_Integer,sgnCell_]:=1/#1 Div[#1*#2,xxx]&[Sqrt[sgnBulk*Det@gLow],NormalVectorSeries[gLow,ind,sgnCell]];*)

(*General Relativity*)

LagrangianScalarEinstein[gLow_List,xxx_List,coeff_]:=coeff Tr[Inverse[gLow].(TensorContract[#,{{1,5},{2,4}}]-TensorContract[#,{{1,6},{4,5}}]&[#\[TensorProduct]#&[ChristoffelSymbol2nd[gLow,xxx]]])];

LagrangianScalarEinsteinHilbert[gLow_List,xxx_List,coeff_]:=coeff RicciScalar[gLow,xxx];

LagrangianScalarArnowittDeserMisnerIndex[gLow_List,xxx_List,sgnBulk_,ind_Integer,sgnCell_,coeff_]:=coeff
	(Tr[#1.#2.#1.#2]&[Inverse[gLow],KruemmungsTensor4Index[gLow,xxx,sgnBulk,ind,sgnCell]]
		-KruemmungsScalar4Index[gLow,xxx,sgnBulk,ind,sgnCell]^2
		+RicciScalar[Drop[MetricTensor3Index[gLow,ind,sgnCell],{ind},{ind}],Delete[xxx,ind]]);

(*LagrangianScalarArnowittDeserMisnerIndexSeries[gLow_List,xxx_List,sgnBulk_,ind_Integer,sgnCell_,coeff_,eps_]:=coeff
	(Tr[#1.#2.#1.#2]&[InverseMetricSeries[gLow,eps],KruemmungsTensor4IndexSeries[gLow,xxx,sgnBulk,ind,sgnCell,eps]]
		-KruemmungsScalar4IndexSeries[gLow,xxx,sgnBulk,ind,sgnCell,eps]^2
		+RicciScalarSeries[Drop[MetricTensor3Index[gLow,ind,sgnCell],{ind},{ind}],Delete[xxx,ind],eps]);*)

End[] (*"Private`"*)L

Begin["Experimental`"]

(*https://en.wikipedia.org/wiki/Weyl_tensor*)

WeylTensor[gLow_List,xxx_List]:=RiemannTensorLow[gLow,xxx]+(-(1/(#-2))(+#1[#2,{1,4,2,3}]-#1[#2,{1,3,2,4}]+#1[#2,{2,3,1,4}]-#1[#2,{2,4,1,3}]&)[Transpose,RicciTensor[gLow]\[TensorProduct]gLow]+RicciScalar[gLow,xxx]/((#-1)(#-2)) (#1[#2,{1,3,2,4}]-#2&)[Transpose,gLow\[TensorProduct]gLow]&)[Length@gLow];

(*https://en.wikipedia.org/wiki/Ricci_decomposition*)

RicciTensorTraceless[gLow_List,xxx_List]:=RicciTensor[gLow,xxx]-(RicciScalar[gLow,xxx]gLow)/Length@gLow;

RiemannTensorLow[gLow_List,xxx_List]:=gLow.RiemannTensor[gLow,xxx];

RiemannTensorScalar[gLow_List,xxx_List]:=RicciScalar[gLow,xxx]/(#(#-1)) (#1[#2,{1,3,4,2}]-#1[#2,{1,4,3,2}]&)[Transpose,gLow\[TensorProduct]gLow]&[Length@gLow];

RiemannTensorSemiTraceless[gLow_List,xxx_List]:=1/(Length@gLow-2) (#1[#2,{1,3,2,4}]-#1[#2,{1,4,2,3}]+#1[#2,{2,4,1,3}]-#1[#2,{2,3,1,4}]&)[Transpose,gLow\[TensorProduct]RicciTensorTraceless[gLow,xxx]];

End[] (*"Experimental`"*)

EndPackage[] (*"gravitation`"*)
