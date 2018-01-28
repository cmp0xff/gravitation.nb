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

ChristoffelSymbol1st::usage="ChristoffelSymbol1st[g,x], with g a n\[Times]n-matrix
  (metric with lower indices) and x n-vector of coordinates, gives the Christoffel symbol of the first kind
  kind (three lower indices)."

ChristoffelSymbol2nd::usage="ChristoffelSymbol2nd[g,x], with g a n\[Times]n-matrix
  (metric with lower indices) and x
  n-vector of coordinates, gives the Christoffel symbol of the second kind
  kind (1st upper, two lower indices)."

RiemannTensor::useage="RiemannTensor[g,x], with g a n\[Times]n-matrix
  (metric with lower indices) and x n-vector of
  coordinates, gives the Riemann curvature tensor (1st upper, three lower
  indices)."

RicciTensor::usage="RicciTensor[g,x], with g a n\[Times]n-matrix
  (metric with lower indices) and x n-vector of
  coordinates, gives the Ricci tensor (two lower symmetric indices)."

RicciScalar::usage="RicciScalar[g,x], with g a n\[Times]n-matrix
  (metric with lower indices) and x
  n-vector of coordinates, gives the Ricci scalar."

EinsteinTensor::usage="EinsteinTensor[g,x] with g a n\[Times]n-matrix
  (metric with lower indices) and x n-vector (the coordinates)
  gives the Einstein tensor (a n\[Times]n-matrix) with lower indices."

(*Functions for extrinsic quantities*)

(*Lapse etc known*)

MetricTensorDecomposed::usage="MetricTensorDecomposed[lapse,shift,hLow] with lapse a number, shift a d-vector, hLow a
d\[Times]d matrix gives a (d+1)\[Times](d+1) metric. It has only been test for for (d+1) decomposition."

MetricTensorSurface::usage="MetricTensorSurface[g,x,n,sV] with g a n\[Times]n-matrix
  (metric with lower indices), x n-vector (the coordinates), n n-vector (the unit normal vector), sV the square of
  the normal vector"

(*Suppose the metric is d+1 decomposed; nUpp=n^\[Mu]. I believe these are to be improved.*)

NormalVector::usage="NormalVector[g,x,i,sV] with g a n\[Times]n-matrix
  (metric with lower indices), x n-vector (the coordinates), i integer index (the dimension to be separated)... sV the
square of the normal vector... gives the unit normal vector"

KruemmungScalar::usage="KruemmungScalar[g,x,sM,i,sV] with g a n\[Times]n-matrix
  (metric with lower indices), x n-vector (the coordinates), sM the signature of the metric,
i integer index (the dimension to be separated)... sV the square of the normal vector... gives the trace of the second fundamental form"

Begin["Private`"]

(*Riemannian formalism*)

(*Functions for intrinsic quantities*)

ChristoffelSymbol1st[gLow_List,xxx_List]:=1/2 (#2-#1[#2,{2,3,1}]+#1[#2,{3,1,2}]&)[Transpose,D[gLow,{xxx}]];

ChristoffelSymbol2nd[gLow_List,xxx_List]:=TensorContract[Inverse[gLow]\[TensorProduct]ChristoffelSymbol1st[gLow,xxx],{{2,3}}];

RiemannTensor[gLow_List,xxx_List]:=((Transpose[#,{1,2,4,3}]-#&)[D[#1, {#2}]]+(#1[#2,{1,3,2,4}]-#1[#2,{1,4,2,3}]&)[Transpose,TensorContract[#1\[TensorProduct]#1,{{2,4}}]]&)[ChristoffelSymbol2nd[gLow,xxx],xxx];

RicciTensor[gLow_List,xxx_List]:=TensorContract[RiemannTensor[gLow,xxx],{{1,3}}];

RicciScalar[gLow_List,xxx_List]:=TensorContract[Inverse[gLow]\[TensorProduct]RicciTensor[gLow,xxx],{{1,3},{2,4}}];

EinsteinTensor[gLow_List,xxx_List]:=RicciTensor[gLow,xxx]-2 (RicciScalar[gLow,xxx]gLow)/Length@gLow;

(*Functions for extrinsic quantities*)

(*Lapse etc known*)

MetricTensorDecomposed[lapse_,shift_,hLow_]:=ArrayFlatten[{{TensorContract[hLow\[TensorProduct]shift\[TensorProduct]shift,{{1,3},{2,4}}]-lapse^2,{hLow.shift}},{{hLow.shift}\[Transpose],hLow}}];

MetricTensorSurface[gLow_List,xxx_List,norVec_List,sgnVec_]:=gLow-sgnVec*(#\[TensorProduct]#&)[gLow.norVec];

(*Suppose the metric is d+1 decomposed; nUpp=n^\[Mu]. I believe these are to be improved.*)

NormalVector[gLow_List,xxx_List,ind_Integer,sgnVec_]:=-1/Sqrt[sgnVec*#[[ind,ind]]] #[[ind]]&[Inverse[gLow]];

KruemmungScalar[gLow_List,xxx_List,sgnMet_,ind_Integer,sgnVec_]:=1/#1 Div[#1*#2,xxx]&[Sqrt[sgnMet*Det@gLow],NormalVector[gLow,xxx,ind,sgnVec]];

End[] (*"Private`"*)

Begin["Experimental`"]

(*https://en.wikipedia.org/wiki/Weyl_tensor*)

WeylTensor[gLow_List,xxx_List]:=RiemannTensorLow[gLow,xxx]+(-(1/(#-2))(+#1[#2,{1,4,2,3}]-#1[#2,{1,3,2,4}]+#1[#2,{2,3,1,4}]-#1[#2,{2,4,1,3}]&)[Transpose,RicciTensor[gLow]\[TensorProduct]gLow]+RicciScalar[gLow,xxx]/((#-1)(#-2)) (#1[#2,{1,3,2,4}]-#2&)[Transpose,gLow\[TensorProduct]gLow]&)[Length@gLow];

(*https://en.wikipedia.org/wiki/Ricci_decomposition*)

RicciTensorTraceless[gLow_List,xxx_List]:=RicciTensor[gLow,xxx]-(RicciScalar[gLow,xxx]gLow)/Length@gLow;

RiemannTensorLow[gLow_List,xxx_List]:=TensorContract[gLow\[TensorProduct]RiemannTensor[gLow,xxx],{{2,3}}];

RiemannTensorScalar[gLow_List,xxx_List]:=RicciScalar[gLow,xxx]/(#(#-1)) (#1[#2,{1,3,4,2}]-#1[#2,{1,4,3,2}]&)[Transpose,gLow\[TensorProduct]gLow]&[Length@gLow];

RiemannTensorSemiTraceless[gLow_List,xxx_List]:=1/(Length@gLow-2) (#1[#2,{1,3,2,4}]-#1[#2,{1,4,2,3}]+#1[#2,{2,4,1,3}]-#1[#2,{2,3,1,4}]&)[Transpose,gLow\[TensorProduct]RicciTensorTraceless[gLow,xxx]];

End[] (*"Experimental`"*)

EndPackage[] (*"gravitation`"*)
