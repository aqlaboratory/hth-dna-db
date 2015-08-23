(* ::Package:: *)

(* ::Text:: *)
(*Returns \[Alpha]-helix closest to DNA.*)
(**)
(*There are three criteria to match. First, the helix of interest has to have the mean of its closest CloseResidues1 residues be less than HelixDNAProximityThreshold1. Distance is measured between the closest atoms in a residue and the closest atom in the DNA molecule. This eliminates HTH-helices that are not the recognition helix. Second, the helix has to have at least CloseResidues2 residues be within HelixDNAProximityThreshold2 of the DNA molecule, with distances measured as above. This in effect insures the right orientation of the helix. Third, the helix has to be within HelixHelixProximityThreshold of two helices, to insure a trihelical bundle.*)


PositionDNABindingHelix::"ShortDNA"="The DNA molecule is too short, rendering results potentially inaccurate.";

Options[PositionDNABindingHelix]={HelixDNAProximityThreshold1->500,HelixDNAProximityThreshold2->650,HelixHelixProximityThreshold->1800,CloseResidues1->5,CloseResidues2->5};
PositionDNABindingHelix[chainProtein_Protein,chainsDNA:{_DNA..},opts___?OptionQ]:=Union@@(PositionDNABindingHelix[chainProtein,#,opts]&/@chainsDNA)
PositionDNABindingHelix[chainProtein_Protein,chainDNA_DNA,opts___?OptionQ]:=Module[{
fNearestDNA,distsAllHelixAtomsToDNAAllAtoms,distsAllHelixAtomsToDNAAllAtomsMeans,distsC\[Alpha]ToDNALine,distsChoice,lineDNA,distsHelixHelix,picksMeanDist,picksNumClose,picks,verticesHelices,verticesPickedHelicesCenters,verticesHelixPairs,
verticesDNA=Flatten[ChainVertices[chainDNA,All],1],
helices=Helix/.Options[chainProtein],
optHelixDNAProximityThreshold1=HelixDNAProximityThreshold1/.{opts}/.Options[PositionDNABindingHelix],
optHelixDNAProximityThreshold2=HelixDNAProximityThreshold2/.{opts}/.Options[PositionDNABindingHelix],
optHelixHelixProximityThreshold=HelixHelixProximityThreshold/.{opts}/.Options[PositionDNABindingHelix],,
optCloseResidues1=CloseResidues1/.{opts}/.Options[PositionDNABindingHelix],
optCloseResidues2=CloseResidues2/.{opts}/.Options[PositionDNABindingHelix]},
fNearestDNA=Nearest[verticesDNA];

If[Length[chainDNA[[1]]]<11,
(* DNA molecule is too short, use less accurate method *)
Message[PositionDNABindingHelix::"ShortDNA"];
distsChoice=distsAllHelixAtomsToDNAAllAtoms,

(* DNA molecule is long enough, use preferred method *)
distsChoice=distsC\[Alpha]ToDNALine;
lineDNA=ChainVertices[chainDNA,{1,11},{{"C","3'"}}];
distsC\[Alpha]ToDNALine=Table[
PointLineDistance[#,lineDNA]&/@ChainVertices[chainProtein,Span[helices[[i]]],{{"C","\[Alpha]"}}],
{i,Length[helices]}]];

distsAllHelixAtomsToDNAAllAtoms=Table[Min[(EuclideanDistance[#,fNearestDNA[#,1][[1]]]&)/@#]&/@ChainVertices[chainProtein,Span[helices[[i]]]],{i,Length[helices]}];
distsAllHelixAtomsToDNAAllAtomsMeans=Mean[PadRight[Sort[#],optCloseResidues1,\[Infinity]]]&/@distsAllHelixAtomsToDNAAllAtoms;
picksMeanDist=(#<=optHelixDNAProximityThreshold1)&/@distsAllHelixAtomsToDNAAllAtomsMeans;
picksNumClose=(Count[#,dist_/;dist<optHelixDNAProximityThreshold2]>=optCloseResidues2)&/@distsAllHelixAtomsToDNAAllAtoms;
picks=MapThread[And,{picksMeanDist,picksNumClose}];
picks=Pick[{#[[1]],#[[2]]+#[[1,1]]-1}&/@Partition[Riffle[helices,MinPosition/@distsChoice],2],picks];

verticesPickedHelicesCenters=ChainVertices[chainProtein,picks[[All,2]],{{"C","\[Alpha]"}}];
verticesHelices=ChainVertices[chainProtein,Span[#],{{"C","\[Alpha]"}}]&/@helices;
verticesHelixPairs=Partition[Append[Riffle[verticesHelices,{#}],#],2]&/@verticesPickedHelicesCenters;
distsHelixHelix=Map[EuclideanDistance[#[[2]],Nearest[#[[1]],#[[2]],1][[1]]]&,verticesHelixPairs,{2}];
Extract[picks[[All,1]],Position[Count[#,dist_/;dist<optHelixHelixProximityThreshold]&/@distsHelixHelix,num_/;num>=3]]]


(* ::Text:: *)
(*Returns a region of windowSize * 2 + 1 of DNA that is closest to an \[Alpha]-helix. The format is *)
(*{{helix1_start, helix1_end}->{{ssDNA_start, ssDNA_end}, ...} or *)
(*{{helix1_start, helix1_end}->{{{dsDNAstrand1_start, dsDNAstrand1_end},{dsDNAstrand2_start, dsDNAstrand2_end}}, ...}*)


PositionBoundDNARegion::OutOfBound="Sought DNA region is out of bound of DNA molecule.";
PositionBoundDNARegion::NoRegion="No protein-bound DNA region has been found.";
PositionBoundDNARegion[chainProtein_Protein,chainDNA_DNA,windowSize_]:=PositionBoundDNARegion[chainProtein,{chainDNA},windowSize]
PositionBoundDNARegion[chainProtein_Protein,chainsDNA:{_DNA..},windowSize_]:=With[
{posHelix=PositionDNABindingHelix[chainProtein,chainsDNA]},
If[posHelix==={},
(* There's no potential helix *)
Message[PositionBoundDNARegion::NoRegion];
{}

,
(* There's a potential helix *)
Module[{
verticesC1DNA=ChainVertices[#,All,{{"C","1'"}}]&/@chainsDNA,
spansHelix=(Span@@#)&/@posHelix,
lengthDNA=Length[chainsDNA[[1,1]]],
fNearestC1DNA},

fNearestC1DNA=Nearest[#->Automatic]&/@verticesC1DNA;
(List@@#->positionBoundDNARegion[chainProtein,#,fNearestC1DNA,lengthDNA,windowSize])&/@spansHelix]]]

positionBoundDNARegion[chainProtein_,spanHelix_,fNearestC1DNA_,lengthDNA_,windowSize_]:=Module[
{verticesC\[Alpha]HTH,posClosestBases,posClosestBase,posDNARegion},

verticesC\[Alpha]HTH=Cases[chainProtein[[1,spanHelix]],{"C","\[Alpha]",coord_}->coord,\[Infinity]];
posClosestBases=Table[SortBy[Tally[fNearestC1DNA[[i]][#,1][[1]]&/@verticesC\[Alpha]HTH],#[[2]]&][[-1]],{i,Length[fNearestC1DNA]}];
posClosestBase=(OrderedQ[posClosestBases[[All,2]]]\[And]Length[fNearestC1DNA]!=1)/.{False->{1,posClosestBases[[1,1]]},True:>{2,posClosestBases[[2,1]]}};
posDNARegion={posClosestBase[[2]]+{-1,1}windowSize,(lengthDNA-(posClosestBase[[2]]-1))+{-1,1}windowSize};
If[Length[fNearestC1DNA]==1,posDNARegion=First[posDNARegion]];
If[Min[posDNARegion]<1\[Or]Max[posDNARegion]>lengthDNA,
Message[PositionBoundDNARegion::OutOfBound];
posDNARegion=Map[Min[Max[#,1],lengthDNA]&,posDNARegion,{-1}]];
If[posClosestBase[[1]]===1,posDNARegion,Reverse[posDNARegion]]]


(* ::Text:: *)
(*Helper function. Generates a sliding window of some fixed length along a given boundary. This takes care of negative boundaries as well, making sure to skip over the 0.*)


SlidingWindow[boundary:{__Integer},n_]:=DeleteCases[Map[
If[Times@@#<=0,#+{0,1},#]&,
Table[{{boundary[[1]]+i,boundary[[1]]+i+n-1}},{i,0,boundary[[2]]-n-boundary[[1]]+1}],{2}],{{0,_}}]
SlidingWindow[boundary:{{__Integer},{__Integer}},n_]:=Partition[Riffle[Flatten[SlidingWindow[boundary[[1]],n],1],Reverse[Flatten[SlidingWindow[boundary[[2]],n],1]]],2]


(* ::Text:: *)
(*Helper function. Gives a window of a certain length n centered at c. Takes into consideration the allowed boundary, and if it spills beyond the boundary, then shifts the center until it can fit the entire window. The Shift option determines whether c is in fact allowed to shift, or if it must return a range centered at c or {} otherwise.*)


Options[CenteredWindow]={Shift->True};
CenteredWindow[boundary_,n_,c_,OptionsPattern[]]:=Module[{ranges=SlidingWindow[boundary,n],devs},
If[ranges=!={},
devs=Abs[Mean/@ranges[[All,1]]-c];
If[OptionValue[Shift]\[Or]MemberQ[devs,0],
ranges[[Ordering[devs][[1]]]],
{}],
{}]]


(* ::Text:: *)
(*Returns whether the pair of DNA chains given form a duplex or not.*)


DNADuplexQ[chainsDNA:{_DNA..},opts:OptionsPattern[ExtractDNADuplex]]:=ExtractDNADuplex[chainsDNA,opts]=!={}
DNADuplexQ[___]=False;


(* ::Text:: *)
(*Returns the region of DNA within the given two strands that forms a duplex. The DNA is cropped so that both strands are of the same length, and are properly lined up. Additionally, it is returned such that both strands run from 5'->3' , even if this were not the form in which they are supplied. The order in which the first strand is given is maintained, and the second may or may not be depending on the directionality.*)
(**)
(*The MaximumOverhang option refers to how much overhang is allowed to be assumed to have existed in the dsDNA passed to it. The MinOverlap option insures that at least that many bases are being aligned, and does so by counting the number thrown out as outliers.*)


ExtractDNADuplex::"ssDNA"="Only a single strand of DNA was passed.";
ExtractDNADuplex::"InvalidOptions"="The minimum overlap requested must exceed the maximum number of outliers.";

Options[ExtractDNADuplex]={MinDistance->600,MaxDistance->1300,MaximumOverhang->4,MaxOutliers->2,MinOverlap->5};
ExtractDNADuplex[chainsDNA:{_DNA..},OptionsPattern[]]:=Module[{vertices,chains,conv,picks,picked,pos,
flagNotReversed=OrderedQ[Length/@chainsDNA[[All,1]]]},

(* Insure that the minimum overlap allowed is longer than the number of outliers to be thrown out. *)
If[OptionValue[MinOverlap]<=OptionValue[MaxOutliers],Message[ExtractDNADuplex::"InvalidOptions"];Return[{}]];

chains=If[flagNotReversed,chainsDNA,Reverse[chainsDNA]];
vertices=ChainVertices[#,All,{{"C","1'"}}]&/@chains;

conv=Join[
ListConvolve[Sequence@@vertices,{-1,1},{},EuclideanDistance,List,1],
Table[Last[
ListConvolve[Sequence@@({Drop[#[[1]],i],Drop[#[[2]],i]}&[vertices]),{-1,1},{},EuclideanDistance,List,1]],
{i,OptionValue[MaximumOverhang]}],
Table[First[
ListConvolve[Sequence@@({Drop[#[[1]],-i],Drop[#[[2]],-i]}&[vertices]),{-1,1},{},EuclideanDistance,List,1]],
{i,OptionValue[MaximumOverhang]}],
ListCorrelate[Sequence@@vertices,{1,-1},{},EuclideanDistance,List,1],
Table[First[
ListCorrelate[Sequence@@({Drop[#[[1]],i],Drop[#[[2]],-i]}&[vertices]),{1,-1},{},EuclideanDistance,List,1]],
{i,OptionValue[MaximumOverhang]}],
Table[Last[
ListCorrelate[Sequence@@({Drop[#[[1]],-i],Drop[#[[2]],i]}&[vertices]),{1,-1},{},EuclideanDistance,List,1]],
{i,OptionValue[MaximumOverhang]}]];
conv=Replace[conv,l_/;Length[l]<OptionValue[MinOverlap]->ConstantArray[\[Infinity],OptionValue[MinOverlap]],{1}];
conv=Delete[#,List/@Ordering[#][[-OptionValue[MaxOutliers];;]]]&/@conv;
picks=Count[#,d_/;(d>OptionValue[MaxDistance]\[Or]d<OptionValue[MinDistance])]===0&/@conv;
picked=Pick[conv,picks];

If[\[Not]MatchQ[picked,{{}...}],
(* Duplex found, return chainsDNA in canonical duplex form *)
pos=MinPosition[Variance/@picked];
pos=Position[picks,True][[pos,1]];
chains=chains[[All,1]];
conv=Join[
ListCorrelate[Sequence@@chains,{1,-1},{},List,List,1],
Table[Last[
ListCorrelate[Sequence@@({Drop[#[[1]],i],Drop[#[[2]],i]}&[chains]),{1,-1},{},List,List,1]],
{i,OptionValue[MaximumOverhang]}],
Table[First[
ListCorrelate[Sequence@@({Drop[#[[1]],-i],Drop[#[[2]],-i]}&[chains]),{1,-1},{},List,List,1]],
{i,OptionValue[MaximumOverhang]}],
ListConvolve[Sequence@@chains,{-1,1},{},List,List,1],
Table[First[
ListConvolve[Sequence@@({Drop[#[[1]],i],Drop[#[[2]],-i]}&[chains]),{-1,1},{},List,List,1]],
{i,OptionValue[MaximumOverhang]}],
Table[Last[
ListConvolve[Sequence@@({Drop[#[[1]],-i],Drop[#[[2]],i]}&[chains]),{-1,1},{},List,List,1]],
{i,OptionValue[MaximumOverhang]}]];
picked=Transpose[conv[[pos]]];
If[\[Not]flagNotReversed,
picked=Reverse[picked],
If[pos>Length[conv]/2,picked=Reverse/@picked]];
DNA/@picked,

(* No duplex found *)
{}]]

ExtractDNADuplex[chainDNA_DNA,___]:=(Message[ExtractDNADuplex::"ssDNA"];chainDNA)


DNAHelixAxes::TooShort="DNA chain is too short to enable reliable extraction of helical axes.";
DNAHelixAxes[protein_Protein,dna:{_DNA..},posFocusBase_Integer,rangeRecHelix_|PatternSequence[]]:=Module[
{vecRecHelix=RecognitionHelixAxis[protein,dna,rangeRecHelix],ranges,vecMajorAxis,vecMinorAxis},
ranges=Table[CenteredWindow[{1,Length[#[[1]]]}&/@dna,11,center],{center,posFocusBase-1,posFocusBase+1}];
If[Flatten[ranges]=!={},
ranges=Flatten[MapAt[Reverse,#,2]&/@ranges,1];
vecMajorAxis=Normalize[Mean[MapThread[
Normalize[Differences[ChainVertices[#1,#2,{{"C","3'"}}]][[1]]]&,
{PadRight[dna,6,dna],ranges}]]];
vecMinorAxis=Normalize[vecRecHelix-(vecRecHelix.vecMajorAxis)vecMajorAxis];
{vecMajorAxis,vecMinorAxis},

Message[DNAHelixAxes::TooShort];
{}]]

DNAHelixAxes[rp_List,posFocusBase_Integer]:=DNAHelixAxes[RPProtein[rp],RPDNA[rp],posFocusBase,RangeProteinMatch[rp]]

DNAHelixAxes[rp_List]:=DNAHelixAxes[RPProtein[rp],RPDNA[rp],\[LeftFloor]Mean[RangeDNAMatch[rp][[1]]]\[RightFloor],RangeProteinMatch[rp]]


RecognitionHelixAxis[protein_Protein,dna:{_DNA..},rangeRecHelix_]:=With[
{vertices=ChainVertices[protein,Span@@rangeRecHelix,{{"C","\[Alpha]"}}]},
Normalize[First[Eigenvectors[Covariance[vertices],1]]]]

RecognitionHelixAxis[protein_Protein,dna:{_DNA..}]:=RecognitionHelixAxis[protein,dna,PositionDNABindingHelix[protein,dna][[1]]]

RecognitionHelixAxis[rp_]:=RecognitionHelixAxis[RPProtein[rp],RPDNA[rp],RangeProteinMatch[rp]]


(* ::Text:: *)
(*Gives the transformation to put the DNA and the HTH into the canonical form, which corresponds to:*)
(**)
(*y is along the major axis of DNA (up/down the molecule)*)
(*x is along the minor axis of DNA (perpendicular to y, parallel to the recognition helix) (right/left along a rung in the molecule)*)
(*z is moving from the DNA to the HTH*)
(**)
(*The center is based on the midpoint of the Subscript[C, 3] atoms of the requested nucleotide (posFocusBase) or the middle nucleotide of the DNA match region if not specified.*)


CanonicalTransformation[protein_Protein,dna:{_DNA..},posFocusBase_Integer,rangeRecHelix_|PatternSequence[]]:=Module[
{axes=DNAHelixAxes[protein,dna,posFocusBase,rangeRecHelix],ptOrigin,basis},
ptOrigin=Mean[Flatten[ChainVertices[#,posFocusBase,{{"C","3'"}}]&/@dna,1]];
basis=CreateBasis[{ptOrigin,ptOrigin+axes[[2]],ptOrigin-axes[[1]]}];
ChangeBasisTransform[{0,0,0},ptOrigin,IdentityMatrix[3],basis]]

CanonicalTransformation[rp_List,posFocusBase_Integer]:=CanonicalTransformation[RPProtein[rp],RPDNA[rp],posFocusBase,RangeProteinMatch[rp]]

CanonicalTransformation[rp_List]:=CanonicalTransformation[RPProtein[rp],RPDNA[rp],\[LeftFloor]Mean[RangeDNAMatch[rp][[1]]]\[RightFloor],RangeProteinMatch[rp]]
