(* ::Package:: *)

(********1*********2*********3*********4*********5*********6*********7*****
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

(*Holonomic formalism*)

(*Tensor algebra*)

KulkarniNomizu

(*Functions for intrinsic quantities*)

VolumeElement(*[gLow_List?MatrixQ,sgnBulk_Integer?(#==1||#==-1&)]:=\
	Sqrt[sgnBulk Det@gLow];*)

VolumeElementSeries(*[gLow_List?MatrixQ,sgnBulk_Integer?(#==1||#==-1&),eps_List]:=\
	Series[VolumeElement[gLow,sgnBulk],eps];*)

InverseMetricSeries

ChristoffelSymbol1st::usage=\
	"ChristoffelSymbol1st[g,x], with g components of the (0,2) metric tensor \
(D\[Times]D Matrix) and x the coordinates (D List), gives the Christoffel \
symbol of the first kind (three lower indices)."

ChristoffelSymbol1stSeries

ChristoffelSymbol2nd::usage=\
	"ChristoffelSymbol2nd[g,x], with g components of the (0,2) metric tensor \
(D\[Times]D Matrix) and x the coordinates (D List), gives the Christoffel \
symbol of the second kind (1st upper, two lower indices)."

ChristoffelSymbol2ndSeries

LeviCivitaDVector

LeviCivitaDVectorSeries

LeviCivitaDCovector

LeviCivitaDCovectorSeries

RiemannTensor::useage=\
	"RiemannTensor[g,x], with g components of the (0,2) metric tensor \
(D\[Times]D Matrix) and x the coordinates (D List), gives the (1,3) Riemann \
curvature tensor (D\[Times]D\[Times]D\[Times]D List)."

RiemannTensorSeries

RiemannTensorLow

RicciTensor::usage=
	"RicciTensor[g,x], with g components of the (0,2) metric tensor \
(D\[Times]D Matrix) and x the coordinates (D List), gives the (0,2) Ricci \
tensor (D\[Times]D Matrix)."

RicciTensorSeries

RicciScalar::usage=\
	"RicciScalar[g,x], with g components of the (0,2) metric tensor \
(D\[Times]D Matrix) and x the coordinates (D List), gives the Ricci scalar."

RicciScalarSeries

EinsteinTensor::usage=\
	"EinsteinTensor[g,x] with g components of the (0,2) metric tensor \
(D\[Times]D Matrix) and x the coordinates (D List) gives the (0,2) Einstein \
tensor (D\[Times]D Matrix)."

EinsteinTensorSeries

RicciTensorTraceless

RiemannTensorScalar

RiemannTensorSemiTraceless

WeylTensor

(*Functions for extrinsic quantities*)

(*Lapse etc known*)

MetricTensor4ArnowittDeserMisner::usage=\
	"MetricTensorADM[lapse,shift,hLow] with lapse a number, shift a vector \
(d List), hLow components of the (0,2) metric tensor (d\[Times]d Matrix), \
gives a D\[Times]D (0,2) metric (D\[Times]D Matrix, D=d+1)."
(* It has only been test for for (d+1) decomposition."*)

MetricTensor4Direct

MetricTensor34Vector::usage=\
	"MetricTensor34Vector[g,nVec,s] with g components of the (0,2) metric \
tensor (D\[Times]D Matrix), n unit normal vector (D List), sC=1 for space-like \
hypersurface and -1 for time-like hypersurface, gives a d\[Times]d (0,2) \
metric (D\[Times]D Matrix, D=d+1)."

MetricTensor34Covector::usage=\
	"MetricTensorSurface[g,nCvt,sV] with g components of the (0,2) metric \
tensor (D\[Times]D Matrix), n unit normal covector (D List), sC=1 for \
space-like hypersurface and -1 for time-like hypersurface, gives a d\[Times]d \
(0,2) metric (D\[Times]D Matrix, D=d+1)."

MetricTensor34Index::usage=
	"MetricTensor34Index[g,sB,i,sC] with g components of the (0,2) metric tensor \
(D\[Times]D Matrix), i index (Integer), sC=1 for space-like hypersurface and \
-1 for time-like hypersurface, gives a d\[Times]d (0,2) metric \
(D\[Times]D Matrix, D=d+1)."

MetricTensor3Index

(* Suppose the metric is d+1 decomposed; nUpp=n^\[Mu]. I believe these are to be 
 * improved.
*)

NormalCovector::usage=\
	"NormalCovector[g,i,sC] with g components of the (0,2) metric tensor \
(D\[Times]D Matrix), i coordinate index (Integer), sC=1 for space-like \
hypersurface and -1 for time-like hypersurface, gives the unit covector \
orthogonal to the hypersurface \!\(\*SuperscriptBox[\(x\), \(i\)]\)==C."

NormalVector::usage="NormalVector[g,i,sC] with g components of the (0,2) \
metric tensor (D\[Times]D Matrix), i coordinate index (Integer), sC=1 for \
space-like hypersurface and -1 for time-like hypersurface, gives the unit \
vector orthogonal to the hypersurface \!\(\*SuperscriptBox[\(x\), \(i\)]\)==C."

KruemmungsTensor4Covector

KruemmungsTensor4CovectorSeries

KruemmungsTensor4Index::usage=\
	"KruemmungsTensor4Index[g,x,sB,i,sC] with g components of the bulk (0,2) \
metric tensor (D\[Times]D Matrix), x the coordinates (D List), sB the \
signature of the bulk metric, i coordinate index and sC=1 for space-like \
hypersurface and -1 for time-like hypersurface, gives the (0,2) second \
fundamental form of the hypersurface \!\(\*SuperscriptBox[\(x\), \(i\)]\)==C \
(D\[Times]D Matrix)."

KruemmungsTensor4IndexSeries

KruemmungsScalar4Index::usage=\
	"KruemmungsScalar4Index[g,x,sB,i,sC] with g components of the bulk (0,2) \
metric tensor (D\[Times]D Matrix), x the coordinates (D List), sB the \
signature of the bulk metric, i coordinate index and sC=1 for space-like \
hypersurface and -1 for time-like hypersurface, gives the trace of the second \
fundamental form of the hypersurface \!\(\*SuperscriptBox[\(x\), \(i\)]\)==C."

KruemmungsScalar4IndexSeries

KruemmungsTensor3Index

KruemmungsTensor3IndexSeries

KruemmungsScalar4Vector

KruemmungsScalar4Covector

KruemmungsScalar3Index

KruemmungsScalar3IndexSeries

(*Gravitation and General Relativity*)

LagrangianScalarEinstein

LagrangianScalarEinsteinHilbert

LagrangianScalarEinsteinHilbertSeries

LagrangianScalarArnowittDeserMisnerIndex

LagrangianScalarArnowittDeserMisnerIndexSeries

Begin["Private`"]

(*Holonomic formalism*)

(*Tensor algebra*)

(*https://en.wikipedia.org/wiki/Kulkarni%E2%80%93Nomizu_product*)

KulkarniNomizu[h_List?MatrixQ,k_List?MatrixQ]:=\
	#1[#2,{1,3,2,4}]+#1[#2,{2,4,1,3}]-#1[#2,{1,4,2,3}]-#1[#2,{2,3,1,4}]&\
		[Transpose,h\[TensorProduct]k];

(*Functions for intrinsic quantities*)

VolumeElement[gLow_List?MatrixQ,sgnBulk_Integer?(#==1||#==-1&),ass_:True]:=\
	Simplify[Sqrt[sgnBulk Det@gLow],Assumptions->ass];

VolumeElementSeries\
	[gLow_List?MatrixQ,sgnBulk_Integer?(#==1||#==-1&),eps_List,ass_:True]:=\
	Simplify[Series[VolumeElement[gLow,sgnBulk],eps],Assumptions->ass];

InverseMetricSeries[gLow_List?MatrixQ,eps_List,ass_:True]:=\
	Series[Inverse@gLow,eps];

ChristoffelSymbol1st[gLow_List?MatrixQ,xxx_List,ass_:True]:=\
	1/2 (#2-#1[#2,{2,3,1}]+#1[#2,{3,1,2}]&)[Transpose,D[gLow,{xxx}]];

ChristoffelSymbol1stSeries[gLow_List?MatrixQ,xxx_List,eps_List,ass_:True]:=\
	Series[ChristoffelSymbol1st[gLow,xxx],eps];

ChristoffelSymbol2nd[gLow_List?MatrixQ,xxx_List,ass_:True]:=\
	Inverse[gLow].ChristoffelSymbol1st[gLow,xxx,ass];

ChristoffelSymbol2ndSeries[gLow_List?MatrixQ,xxx_List,eps_List,ass_:True]:=\
	InverseMetricSeries[gLow,eps].ChristoffelSymbol1stSeries[gLow,xxx,eps];

LeviCivitaDVector[gLow_List?MatrixQ,xxx_List,vec_List]:=\
	D[vec,{xxx}]+ChristoffelSymbol2nd[gLow,xxx].vec;

LeviCivitaDVectorSeries[gLow_List?MatrixQ,xxx_List,vec_List,eps_List]:=\
	D[Normal[vec],{xxx}]+ChristoffelSymbol2ndSeries[gLow,xxx,eps].vec;

LeviCivitaDCovector[gLow_List?MatrixQ,xxx_List,cvt_List]:=\
	D[cvt,{xxx}]-cvt.ChristoffelSymbol2nd[gLow,xxx];

LeviCivitaDCovectorSeries[gLow_List?MatrixQ,xxx_List,cvt_List,eps_List]:=\
	D[Normal[cvt],{xxx}]-cvt.ChristoffelSymbol2ndSeries[gLow,xxx,eps];

RiemannTensor[gLow_List?MatrixQ,xxx_List]:=\
	(Transpose[#,{1,2,4,3}]-#&)\
		[D[#1, {#2}]]\
	+(#1[#2,{1,3,2,4}]-#1[#2,{1,4,2,3}]&)\
		[Transpose,TensorContract[#1\[TensorProduct]#1,{{2,4}}]]&\
			[ChristoffelSymbol2nd[gLow,xxx],xxx];

RiemannTensorSeries[gLow_List?MatrixQ,xxx_List,eps_List,ass_:True]:=\
	(Transpose[#,{1,2,4,3}]-#&)\
		[D[Normal[#1], {#2}]]\
	+(#1[#2,{1,3,2,4}]-#1[#2,{1,4,2,3}]&)\
		[Transpose,(*Transpose[#1,{1,3,2}].Transpose[#1,{2,1,3}]*)\
				TensorContract[#1\[TensorProduct]#1,{{2,4}}]]&\
			[ChristoffelSymbol2ndSeries[gLow,xxx,eps,ass],xxx];
			
RiemannTensorLow[gLow_List?MatrixQ,xxx_List]:=gLow.RiemannTensor[gLow,xxx];

RicciTensor[gLow_List?MatrixQ,xxx_List]:=\
	TensorContract[RiemannTensor[gLow,xxx],{{1,3}}];

RicciTensorSeries[gLow_List?MatrixQ,xxx_List,eps_List,ass_:True]:=\
	TensorContract[RiemannTensorSeries[gLow,xxx,eps,ass],{{1,3}}];

RicciScalar[gLow_List?MatrixQ,xxx_List]:=\
	Tr[Inverse[gLow].RicciTensor[gLow,xxx]];

RicciScalarSeries[gLow_List?MatrixQ,xxx_List,eps_List,ass_:True]:=\
	Tr[InverseMetricSeries[gLow,eps].RicciTensorSeries[gLow,xxx,eps]];

EinsteinTensor[gLow_List?MatrixQ,xxx_List]:=\
	RicciTensor[gLow,xxx]-2RicciScalar[gLow,xxx]gLow/Length@gLow;

EinsteinTensorSeries[gLow_List?MatrixQ,xxx_List,eps_List,ass_:True]:=\
	RicciTensorSeries[gLow,xxx,eps,ass]\
		-2RicciScalarSeries[gLow,xxx,eps,ass]gLow/Length@gLow;

(*https://en.wikipedia.org/wiki/Ricci_decomposition*)

RicciTensorTraceless[gLow_List?MatrixQ,xxx_List]:=\
	RicciTensor[#1,#2]-(RicciScalar[#1,#2]#1)/Length@#1&\
		[gLow,xxx];

RiemannTensorScalar[gLow_List?MatrixQ,xxx_List]:=\
	KulkarniNomizu[#1,#1]RicciScalar[#1,#2]/(2#3(#3-1))&[gLow,xxx,Length@gLow];
	(*RicciScalar[gLow,xxx]/(#(#-1)) (#1[#2,{1,3,4,2}]-#1[#2,{1,4,3,2}]&)[Transpose,gLow\[TensorProduct]gLow]&[Length@gLow];*)

RiemannTensorSemiTraceless[gLow_List?MatrixQ,xxx_List]:=
	KulkarniNomizu[#1,RicciTensorTraceless[#1,#2]]/(#3-2)&[gLow,xxx,Length@gLow];
	(*1/(Length@gLow-2) (#1[#2,{1,3,2,4}]-#1[#2,{1,4,2,3}]+#1[#2,{2,4,1,3}]-#1[#2,{2,3,1,4}]&)[Transpose,gLow\[TensorProduct]RicciTensorTraceless[gLow,xxx]];*)

(*https://en.wikipedia.org/wiki/Weyl_tensor*)

WeylTensor[gLow_List?MatrixQ,xxx_List]:=\
	RiemannTensorLow[#1,#2]-RiemannTensorSemiTraceless[#1,#2]-RiemannTensorScalar[#1,#2]&\
		[gLow,xxx];
	(*-KulkarniNomizu[(RicciTensor[#1,#2]/(#3-2)+RicciScalar[#1,#2](1/(2#3(#3-1))-1/((#3-2)#3)))\
		#1,#1]&[gLow,xxx,Length@gLow];*)
	(*RiemannTensorLow[gLow,xxx]+\
	(-(1/(#-2))(+#1[#2,{1,4,2,3}]-#1[#2,{1,3,2,4}]+#1[#2,{2,3,1,4}]-#1[#2,{2,4,1,3}]&)\
			[Transpose,RicciTensor[gLow,xxx]\[TensorProduct]gLow]\
		+RicciScalar[gLow,xxx]/((#-1)(#-2)) (#1[#2,{1,3,2,4}]-#2&)\
			[Transpose,gLow\[TensorProduct]gLow]&)[Length@gLow];*)

(*Functions for extrinsic quantities*)

(*Lapse etc known*)

MetricTensor4ArnowittDeserMisner\
	[lapse_,shift_List,hLow_List?MatrixQ]:=\
	ArrayFlatten[{{shift.#-lapse^2,{#}},{{#}\[Transpose],hLow}}]&[hLow.shift];

MetricTensor4Direct[tt_,ts_List,ss_List]:=\
	ArrayFlatten[{{tt,{ts}},{{ts}\[Transpose],ss}}];

MetricTensor34Covector\
	[gLow_List?MatrixQ,norCvt_List,sgnCell_Integer?(#==1||#==-1&)]:=\
	gLow+sgnCell*(norCvt\[TensorProduct]norCvt);

MetricTensor34Vector\
	[gLow_List?MatrixQ,norVec_List,sgnCell_Integer?(#==1||#==-1&)]:=\
	MetricTensor34Covector[gLow,gLow.norVec,sgnCell];

MetricTensor34Index\
	[gLow_List?MatrixQ,sgnBulk_Integer?(#==1||#==-1&),\
		ind_Integer,sgnCell_Integer?(#==1||#==-1&),ass_:True]:=\
	MetricTensor34Covector[gLow,NormalCovector[gLow,sgnBulk,ind,sgnCell,ass],sgnCell];

MetricTensor3Index\
	[gLow_List?MatrixQ,sgnBulk_Integer?(#==1||#==-1&),\
		ind_Integer,sgnCell_Integer?(#==1||#==-1&),ass_:True]:=\
	Drop[MetricTensor34Index[gLow,sgnBulk,ind,sgnCell,ass],{ind},{ind}];

(* Suppose the metric is d+1 decomposed; nUpp=n^\[Mu].
 * I believe these are to be improved. *)

NormalCovector\
	[gLow_List?MatrixQ,sgnBulk_Integer?(#==1||#==-1&),\
		ind_Integer,sgnCell_Integer?(#==1||#==-1&),ass_:True]:=\
	Array[If[#==ind,Simplify[sgnBulk sgnCell*\
			Sqrt[1/Together[sgnBulk sgnCell Inverse[gLow][[ind,ind]]]],\
			Assumptions->ass],0]&,\
		Length@gLow];

NormalVector\
	[gLow_List?MatrixQ,sgnBulk_Integer?(#==1||#==-1&),\
		ind_Integer,sgnCell_Integer?(#==1||#==-1&),ass_:True]:=\
	Simplify[sgnBulk sgnCell Sqrt[1/Together[sgnBulk sgnCell*#[[ind,ind]]]],\
		Assumptions->ass] #[[ind]]&\
		[Inverse[gLow]];

KruemmungsTensor4Covector\
	[gLow_List?MatrixQ,xxx_List,sgnBulk_Integer?(#==1||#==-1&),\
		norCvt_List,sgnCell_Integer?(#==1||#==-1&),ass_:True]:=\
	MetricTensor34Covector[gLow,norCvt,sgnCell].Inverse[gLow].\
		LeviCivitaDCovector[gLow,xxx,norCvt];

KruemmungsTensor4CovectorSeries\
	[gLow_List?MatrixQ,xxx_List,sgnBulk_Integer?(#==1||#==-1&),\
		norCvt_List,sgnCell_Integer?(#==1||#==-1&),eps_List,ass_:True]:=\
	MetricTensor34Covector[gLow,norCvt,sgnCell].InverseMetricSeries[gLow,eps].\
		LeviCivitaDCovector[gLow,xxx,Normal@Series[norCvt,eps]];

KruemmungsTensor4Index\
	[gLow_List?MatrixQ,xxx_List,sgnBulk_Integer?(#==1||#==-1&),\
		ind_Integer,sgnCell_Integer?(#==1||#==-1&),ass_:True]:=\
	KruemmungsTensor4Covector[gLow,xxx,sgnBulk,\
		NormalCovector[gLow,sgnBulk,ind,sgnCell,ass],sgnCell,ass];
	(*MetricTensor34Covector[gLow,#,sgnCell].Inverse[gLow].\
		LeviCivitaDCovector[gLow,xxx,#]&\
		[NormalCovector[gLow,sgnBulk,ind,sgnCell,ass]];*)

KruemmungsTensor4IndexSeries\
	[gLow_List?MatrixQ,xxx_List,sgnBulk_Integer?(#==1||#==-1&),\
		ind_Integer,sgnCell_Integer?(#==1||#==-1&),eps_List,ass_:True]:=\
	KruemmungsTensor4CovectorSeries[gLow,xxx,sgnBulk,\
		NormalCovector[gLow,sgnBulk,ind,sgnCell,ass],sgnCell,eps,ass];
	(*MetricTensor34Covector[gLow,#,sgnCell].InverseMetricSeries[gLow,eps].\
			LeviCivitaDCovectorSeries[gLow,xxx,#,eps]&\
		[Normal@Series[NormalCovector[gLow,sgnBulk,ind,sgnCell,ass],eps]];*)

KruemmungsTensor3Index\
	[gLow_List?MatrixQ,xxx_List,sgnBulk_Integer?(#==1||#==-1&),\
		ind_Integer,sgnCell_Integer?(#==1||#==-1&),ass_:True]:=\
	Drop[KruemmungsTensor4Index[gLow,xxx,sgnBulk,ind,sgnCell,ass],{ind},{ind}];

KruemmungsTensor3IndexSeries\
	[gLow_List?MatrixQ,xxx_List,sgnBulk_Integer?(#==1||#==-1&),\
		ind_Integer,sgnCell_Integer?(#==1||#==-1&),eps_List,ass_:True]:=\
	Drop[MetricTensor34Covector[gLow,#,sgnCell].InverseMetricSeries[gLow,eps].\
			LeviCivitaDCovectorSeries[gLow,xxx,#,eps]&\
		[Normal@Series[NormalCovector[gLow,sgnBulk,ind,sgnCell,ass],eps]],\
	{ind},{ind}];

KruemmungsScalar4Vector\
	[gLow_List?MatrixQ,xxx_List,sgnBulk_Integer?(#==1||#==-1&),\
		norVec_List,ass_:True]:=\
	1/# Div[#*norVec,xxx]&[VolumeElement[gLow,sgnBulk,ass]];

KruemmungsScalar4Covector\
	[gLow_List?MatrixQ,xxx_List,sgnBulk_Integer?(#==1||#==-1&),\
		norCvt_List,ass_:True]:=\
	KruemmungsScalar4Vector[gLow,xxx,sgnBulk,Inverse[gLow].norCvt,ass];

KruemmungsScalar4Index\
	[gLow_List?MatrixQ,xxx_List,sgnBulk_Integer?(#==1||#==-1&),\
		ind_Integer,sgnCell_Integer?(#==1||#==-1&),ass_:True]:=\
	1/#1 Div[#1*#2,xxx]&\
		[VolumeElement[gLow,sgnBulk,ass],NormalVector[gLow,sgnBulk,ind,sgnCell,ass]];

KruemmungsScalar4IndexSeries\
	[gLow_List?MatrixQ,xxx_List,sgnBulk_Integer?(#==1||#==-1&),\
		ind_Integer,sgnCell_Integer?(#==1||#==-1&),eps_List,ass_:True]:=\
	1/#1 Div[Normal[#1*#2],xxx]&\
		[VolumeElementSeries[gLow,sgnBulk,eps,ass],\
			NormalVector[gLow,sgnBulk,ind,sgnCell,ass]];

KruemmungsScalar3Index\
	[gLow_List?MatrixQ,xxx_List,sgnBulk_Integer?(#==1||#==-1&),\
		ind_Integer,sgnCell_Integer?(#==1||#==-1&),ass_:True]:=\
	Tr[Inverse@MetricTensor3Index[gLow,sgnBulk,ind,sgnCell].\
		KruemmungsTensor3Index[gLow,xxx,sgnBulk,ind,sgnCell,ass]];

KruemmungsScalar3IndexSeries\
	[gLow_List?MatrixQ,xxx_List,sgnBulk_Integer?(#==1||#==-1&),\
		ind_Integer,sgnCell_Integer?(#==1||#==-1&),eps_List,ass_:True]:=\
	Tr[Inverse@MetricTensor3Index[gLow,sgnBulk,ind,sgnCell].\
		KruemmungsTensor3IndexSeries[gLow,xxx,sgnBulk,ind,sgnCell,eps,ass]];

(*KruemmungsScalar4IndexSeries\
	[gLow_List?MatrixQ,xxx_List,sgnBulk_Integer?(#==1||#==-1&),ind_Integer,sgnCell_Integer?(#==1||#==-1&)]:=\
		1/#1 Div[#1*#2,xxx]&\
			[Sqrt[sgnBulk*Det@gLow],NormalVectorSeries[gLow,ind,sgnCell]];*)

(*General Relativity*)

LagrangianScalarEinstein[gLow_List?MatrixQ,xxx_List,coeff_]:=\
	coeff Tr[Inverse[gLow].(TensorContract[#,{{1,5},{2,4}}]-\
		TensorContract[#,{{1,6},{4,5}}]&\
			[#\[TensorProduct]#&\
		[ChristoffelSymbol2nd[gLow,xxx]]])];

LagrangianScalarEinsteinHilbert[gLow_List?MatrixQ,xxx_List,coeff_]:=\
	coeff RicciScalar[gLow,xxx];

LagrangianScalarEinsteinHilbertSeries\
	[gLow_List?MatrixQ,xxx_List,coeff_,eps_List,ass_:True]:=\
	coeff RicciScalarSeries[gLow,xxx,eps,ass];

LagrangianScalarArnowittDeserMisnerIndex\
	[gLow_List?MatrixQ,xxx_List,sgnBulk_Integer?(#==1||#==-1&),\
		ind_Integer,sgnCell_Integer?(#==1||#==-1&),coeff_]:=\
	coeff(Tr[#1.#2.#1.#2]&\
		[Inverse@#,
		KruemmungsTensor3Index[gLow,xxx,sgnBulk,ind,sgnCell]]\
	-KruemmungsScalar3Index[gLow,xxx,sgnBulk,ind,sgnCell]^2\
	+RicciScalar[#,Delete[xxx,ind]])&\
		[MetricTensor3Index[gLow,sgnBulk,ind,sgnCell]];

LagrangianScalarArnowittDeserMisnerIndexSeries\
	[gLow_List?MatrixQ,xxx_List,sgnBulk_Integer?(#==1||#==-1&),\
		ind_Integer,sgnCell_Integer?(#==1||#==-1&),coeff_,eps_List,ass_:True]:=\
	coeff(Tr[MatrixPower[#1.#2,2]]&\
		[InverseMetricSeries[#,eps],
		KruemmungsTensor3IndexSeries[gLow,xxx,sgnBulk,ind,sgnCell,eps,ass]]\
	-KruemmungsScalar3IndexSeries[gLow,xxx,sgnBulk,ind,sgnCell,eps,ass]^2\
	+RicciScalarSeries[#,Delete[xxx,ind],eps,ass])&\
		[MetricTensor3Index[gLow,sgnBulk,ind,sgnCell]];

End[] (*"Private`"*)

Begin["Experimental`"]

End[] (*"Experimental`"*)

EndPackage[] (*"gravitation`"*)
