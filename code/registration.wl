(* ::Package:: *)

(* ::Text:: *)
(*Helper function that takes the number of residues and atoms in two alpha helices, and returns all possible correspondences between the two (for in phase and out of phase orientations), where the shorter helix is embedded in the longer one. The lenMinOverlap argument determines how many residues must overlap when some others are hanging off from one helix wrt the other. It is not the minimum overlap in general--i.e. one helix that is very short may only overlap 2 residues so long as it is fully contained, but the lenMinOverlap may be much higher. The default is to allow no overhangs.*)


AlphaHelixCorrespondences[{baseResiduesLength_,basePointsLength_},{modelResiduesLength_,modelPointsLength_},phases_:Both,lenMinOverlap_:\[Infinity],shift_:0]:=Module[
{minResLen,maxResLen,minPtLen,maxPtLen,selPts,corrs,corrsLOverhangs,corrsROverhangs,fMost,fRest,l},
{minResLen,maxResLen}=Sort[{baseResiduesLength,modelResiduesLength}];
{minPtLen,maxPtLen}=Sort[{basePointsLength,modelPointsLength}];

(* Core correspondences *)
selPts=Flatten[Transpose[#]]&/@
Map[
Table[i*maxResLen+#,{i,0,maxPtLen/maxResLen-1}]&,
Flatten[
Table[
{If[phases===Out,Unevaluated[Sequence[]],Range[i,i+minResLen-1]],
If[phases===In,Unevaluated[Sequence[]],Range[i+minResLen-1,i,-1]]},
{i,1,maxResLen-minResLen+1}],
1],
{2}];
corrs=(Range[minPtLen]->#)&/@selPts;

(* Hanging correspondences *)
fMost[l_]:=Delete[l,Table[{i*Length[l]*minResLen/minPtLen},{i,1,minPtLen/minResLen}]];
fRest[l_]:=Delete[l,Table[{1+i*Length[l]*minResLen/minPtLen},{i,0,minPtLen/minResLen-1}]];
corrsLOverhangs=DeleteCases[Flatten[Rest[
NestList[
{If[phases===Out,Null,fRest[#[[1,1]]]->fMost[#[[1,2]]]],
If[phases===In,Null,fMost[#[[2,1]]]->fRest[#[[2,2]]]]}&,
Switch[phases,
Both,corrs[[;;2]],
In,{corrs[[1]],{}},
Out,{{},corrs[[1]]}],
Max[minResLen-lenMinOverlap,0]]],
1],Null];
corrsROverhangs=DeleteCases[Flatten[Rest[
NestList[
{If[phases===Out,Null,fMost[#[[1,1]]]->fRest[#[[1,2]]]],
If[phases===In,Null,fRest[#[[2,1]]]->fMost[#[[2,2]]]]}&,
Switch[phases,
Both,corrs[[-2;;]],
In,{corrs[[-1]],{}},
Out,{{},corrs[[-1]]}],
Max[minResLen-lenMinOverlap,0]]],
1],Null];
corrs=Join[corrsLOverhangs,corrs,corrsROverhangs];

(* Shift and proper placement of base vs. model *)
corrs=Map[shift+#&,corrs,{3}];
If[baseResiduesLength<modelResiduesLength,
Reverse/@corrs,
corrs]]


(* ::Text:: *)
(*Helper function that takes two DNA molecules, either of which can be ssDNA or dsDNA, of equal length, and returns all possible correspondences.*)


DNACorrespondences[
{baseResiduesLength_,basePointsLength_,baseDNAType_},
{modelResiduesLength_,modelPointsLength_,modelDNAType_}]:=Module[{
baseUnpar=Range[basePointsLength],
modelUnpar=Range[modelPointsLength],
basePar=Reverse[Partition[#,\[LeftFloor]baseResiduesLength/2\[RightFloor]]&/@Partition[Range[basePointsLength],baseResiduesLength],2],
modelPar=Reverse[Partition[#,\[LeftFloor]modelResiduesLength/2\[RightFloor]]&/@Partition[Range[modelPointsLength],modelResiduesLength],2],
baseSSDNAQ=(baseDNAType===Single),
modelSSDNAQ=(modelDNAType===Single)},

If[baseSSDNAQ\[And]modelSSDNAQ,
(* Both ssDNA *)
{modelUnpar->baseUnpar},

(* One or both are dsDNA *)
If[\[Not]baseSSDNAQ\[And]\[Not]modelSSDNAQ,
(* Both dsDNA *)
{modelUnpar->baseUnpar,Flatten[modelPar]->baseUnpar},

If[baseSSDNAQ,
(* Base is ssDNA, Model is dsDNA *)
{Flatten[modelPar[[All,1]]]->baseUnpar,Flatten[modelPar[[All,2]]]->baseUnpar},

(* Base is dsDNA, Model is ssDNA *)
{modelUnpar->Flatten[basePar[[All,1]]],modelUnpar->Flatten[basePar[[All,2]]]}]
]]]


(* ::Text:: *)
(*Wrapper function that generates a database of models out of a DNA chain by determining its protein-bound region and creating a series of models out of a given window size around the central/"canonical" base. *)
(**)
(*The return format is:*)
(*{{{dna,...}, protein, database:{{{name, proteinChain # -> dnaChain #, proteinRegion -> dnaRegion, subwindow},{...},{...}}, ... }}, ... }*)
(**)
(*Where the dnaChain # can either be for ssDNA or dsDNA, and the dnaRegion and subwindow accordingly...*)


RegistrationPackage::DNANotBound="Protein is not binding DNA";

RegistrationPackage[name_String,chainsPDB_List,OptionsPattern[RegisterByDNA]]:=Module[
{numPairs,dnaPairs,picks,numDups,numSings,groupsDNA,windows,
dnas=Cases[chainsPDB,_DNA],
proteins=Cases[chainsPDB,_Protein]},

numPairs=Subsets[Flatten[Position[dnas,dna_/;(Length[dna]>=1\[And]Length[dna[[1]]]>=OptionValue[SlideWindow]),{1}]],{2}];
dnaPairs=Part[dnas,#]&/@numPairs;
picks=DNADuplexQ/@dnaPairs;
numDups=Union[Pick[numPairs,picks],SameTest->(#1\[Intersection]#2=!={}&)];
numSings=Complement[Range[Length[dnas]],Flatten[numDups]];
groupsDNA=Join[numDups,numSings];

numPairs=Tuples[{Range[Length[proteins]],groupsDNA}];
picks=With[
{num=Rule@@#,protein=proteins[[#[[1]]]],dna=If[NumberQ[#[[2]]],List,ExtractDNADuplex][dnas[[#[[2]]]]]},
{num,protein,dna,PositionBoundDNARegion[protein,dna,OptionValue[MatchBuffer]]}]&/@numPairs;
picks=DeleteCases[picks,{_,_,_,{}}];
picks=Replace[picks,{e1_,e2_,e3_,e4_}:>Sequence@@({e1,e2,e3,#}&/@e4),{1}];

windows=Map[(Span@@#&),SlidingWindow[#,OptionValue[SlideWindow]]&/@picks[[All,4,2]],{-2}];
MapThread[
Function[
{pick,windowSet},
{pick[[3]],pick[[2]],{name,pick[[1]],pick[[4]],#}&/@windowSet}],
{picks,windows}]]


(* ::Text:: *)
(*This function does the end-to-end registration of two HTHs by aligning their corresponding bound DNA molecules. *)


RegisterByDNA::"NonmatchingAtoms"="The number of atoms of each type are not matching between the two structures.";
RegisterByDNA::"CloseMSDSeparation"="The second best match is within `1` of the best match MSD.";
RegisterByDNA::"NoMatchFound"="No matches found within the requested DNA MSD.";

Options[RegisterByDNA]={DNAVertexTypes->{{"C","1'"|"2'"|"3'"|"4'"|"5'"}},ProteinVertexTypes->{{"C","\[Alpha]"|""}},MatchBuffer->4,SlideWindow->5,MinOverlap->\[Infinity],MSDSeparationWarning->(0.1#&),ResultRank->1,ProteinPhases->In,DNAMSDThreshold->40000,ReturnDNAMSD->False};

RegisterByDNA[{chainsBaseDNA_,chainBaseProtein_,nameBase_String},{chainsModelDNA_,chainModelProtein_,nameModel_String},opts:OptionsPattern[]]:=RegisterByDNA[
{chainsBaseDNA,chainBaseProtein,RegistrationPackage[nameBase,Append[chainsBaseDNA,chainBaseProtein],opts]},
{chainsModelDNA,chainModelProtein,RegistrationPackage[nameModel,Append[chainsModelDNA,chainModelProtein],opts]},
opts]

RegisterByDNA[{chainsBaseDNA_,chainBaseProtein_,dbBaseDNA_List},{chainsModelDNA_,chainModelProtein_,dbModelDNA_List},OptionsPattern[]]:=Module[
{windowsBaseDNA,windowsModelDNA,verticesBaseDNA,verticesModelDNA,labelsBaseDNA,labelsModelDNA,verticesBaseAlphaHelix,verticesModelAlphaHelix,labelsBaseAlphaHelix,labelsModelAlphaHelix,verticesBase,verticesModel,sizesAndLabelsBaseDNA1,sizesAndLabelsModelDNA1,sizesAndLabelsAdjustedBaseDNA1,sizesAndLabelsAdjustedModelDNA1,sizesAndLabelsBaseAlphaHelix,sizesAndLabelsModelAlphaHelix,verticesGroupedByLabelBaseDNA,verticesGroupedByLabelModelDNA,verticesGroupedByLabelBaseAlphaHelix,verticesGroupedByLabelModelAlphaHelix,verticesReorderedByLabelBaseDNA,verticesReorderedByLabelModelDNA,verticesReorderedByLabelBaseAlphaHelix,verticesReorderedByLabelModelAlphaHelix,typeBaseDNA,typeModelDNA,dnaAtomTypesAndNumberSameQ,allAtomTypesInEveryResidueQ,icp,corrsVerticesDNA,corrsResiduesDNA,corrsFinalDNA,corrsChains,corrsVerticesAlphaHelix,corrsResiduesAlphaHelix,corrsFinalAlphaHelix,shift,MSDs,pos,transform,results},

{windowsBaseDNA,windowsModelDNA}=#[[All,4]]&/@{dbBaseDNA,dbModelDNA};

{{verticesBaseDNA,verticesModelDNA},{labelsBaseDNA,labelsModelDNA}}=Table[Flatten[
MapThread[fInfo[#1,#2,OptionValue[DNAVertexTypes]]&,{chainsAndWindows[[1]],#}],1]&/@chainsAndWindows[[2]],
{fInfo,{ChainVertices,ChainAtomLabels}},
{chainsAndWindows,{{chainsBaseDNA,windowsBaseDNA},{chainsModelDNA,windowsModelDNA}}}];

{{verticesBaseAlphaHelix,verticesModelAlphaHelix},{labelsBaseAlphaHelix,labelsModelAlphaHelix}}=Table[
fInfo[chainAndDB[[1]],Span@@chainAndDB[[2]][[1,3,1]],OptionValue[ProteinVertexTypes]],
{fInfo,{ChainVertices,ChainAtomLabels}},
{chainAndDB,{{chainBaseProtein,dbBaseDNA},{chainModelProtein,dbModelDNA}}}];

verticesGroupedByLabelBaseDNA=MapThread[Regroup,{verticesBaseDNA,labelsBaseDNA}];
verticesGroupedByLabelModelDNA=MapThread[Regroup,{verticesModelDNA,labelsModelDNA}];
verticesGroupedByLabelBaseAlphaHelix=Regroup[verticesBaseAlphaHelix,labelsBaseAlphaHelix];
verticesGroupedByLabelModelAlphaHelix=Regroup[verticesModelAlphaHelix,labelsModelAlphaHelix];

verticesReorderedByLabelBaseDNA=Flatten[#[[All,2]],1]&/@verticesGroupedByLabelBaseDNA;
verticesReorderedByLabelModelDNA=Flatten[#[[All,2]],1]&/@verticesGroupedByLabelModelDNA;
verticesReorderedByLabelBaseAlphaHelix=Flatten[verticesGroupedByLabelBaseAlphaHelix[[All,2]],1];
verticesReorderedByLabelModelAlphaHelix=Flatten[verticesGroupedByLabelModelAlphaHelix[[All,2]],1];

typeBaseDNA=Length[chainsBaseDNA]/.{1->Single,2->Double};
typeModelDNA=Length[chainsModelDNA]/.{1->Single,2->Double};

sizesAndLabelsBaseDNA1={#[[1]],Length[#[[2]]]}&/@verticesGroupedByLabelBaseDNA[[1]];
sizesAndLabelsModelDNA1={#[[1]],Length[#[[2]]]}&/@verticesGroupedByLabelModelDNA[[1]];
sizesAndLabelsAdjustedBaseDNA1={#[[1]],#[[2]]*(typeBaseDNA/.{Double->1/2,_->1})}&/@sizesAndLabelsBaseDNA1;
sizesAndLabelsAdjustedModelDNA1={#[[1]],#[[2]]*(typeModelDNA/.{Double->1/2,_->1})}&/@sizesAndLabelsModelDNA1;
sizesAndLabelsBaseAlphaHelix={#[[1]],Length[#[[2]]]}&/@verticesGroupedByLabelBaseAlphaHelix;
sizesAndLabelsModelAlphaHelix={#[[1]],Length[#[[2]]]}&/@verticesGroupedByLabelModelAlphaHelix;

dnaAtomTypesAndNumberSameQ=sizesAndLabelsAdjustedBaseDNA1===sizesAndLabelsAdjustedModelDNA1;
allAtomTypesInEveryResidueQ=And@@(Length[Union[#[[All,2]]]]===1&/@{
sizesAndLabelsBaseDNA1,sizesAndLabelsModelDNA1,sizesAndLabelsBaseAlphaHelix,sizesAndLabelsModelAlphaHelix});

If[dnaAtomTypesAndNumberSameQ\[And]allAtomTypesInEveryResidueQ,
(* Types match, proceed *)
(* ICP *)
corrsVerticesDNA=DNACorrespondences[
{sizesAndLabelsBaseDNA1[[1,2]],Length[verticesReorderedByLabelBaseDNA[[1]]],typeBaseDNA},
{sizesAndLabelsModelDNA1[[1,2]],Length[verticesReorderedByLabelModelDNA[[1]]],typeModelDNA}];
corrsResiduesDNA=DNACorrespondences[
{sizesAndLabelsBaseDNA1[[1,2]],sizesAndLabelsBaseDNA1[[1,2]],typeBaseDNA},
{sizesAndLabelsModelDNA1[[1,2]],sizesAndLabelsModelDNA1[[1,2]],typeModelDNA}];
icp=Table[
ICP[
Join[verticesReorderedByLabelModelDNA[[j]],verticesReorderedByLabelModelAlphaHelix],
verticesReorderedByLabelBaseDNA[[i]],
Correspondences->{corrsVerticesDNA[[k]]}],

{i,Length[windowsBaseDNA]},{j,Length[windowsModelDNA]},{k,Length[corrsVerticesDNA]}];

(* Calculate MSDs of the alpha helices *)
shift=Length[verticesReorderedByLabelModelDNA[[1]]];
corrsVerticesAlphaHelix=AlphaHelixCorrespondences[
{sizesAndLabelsBaseAlphaHelix[[1,2]],Length[verticesReorderedByLabelBaseAlphaHelix]},
{sizesAndLabelsModelAlphaHelix[[1,2]],Length[verticesReorderedByLabelModelAlphaHelix]},
OptionValue[ProteinPhases],
OptionValue[MinOverlap],
shift];
corrsResiduesAlphaHelix=AlphaHelixCorrespondences[
{sizesAndLabelsBaseAlphaHelix[[1,2]],sizesAndLabelsBaseAlphaHelix[[1,2]]},
{sizesAndLabelsModelAlphaHelix[[1,2]],sizesAndLabelsModelAlphaHelix[[1,2]]},
OptionValue[ProteinPhases],
OptionValue[MinOverlap]];
MSDs=Table[
cMSD[
icp[[i,j,l,1]][[corrsVerticesAlphaHelix[[k,1]]]],
verticesReorderedByLabelBaseAlphaHelix[[corrsVerticesAlphaHelix[[k,2]]-shift]]],
{i,Length[windowsBaseDNA]},{j,Length[windowsModelDNA]},{k,Length[corrsVerticesAlphaHelix]},{l,Length[corrsVerticesDNA]}];

(* For all possible results *)
results=Table[
(* Extract, format, and return result *)
pos={i,j,k,l};
transform=Composition[Sequence@@Reverse[TransformationFunction/@icp[[pos[[1]],pos[[2]],pos[[4]],3]]]];

corrsChains={dbModelDNA[[pos[[2]],2]],dbBaseDNA[[pos[[1]],2]]};
corrsFinalDNA=Map[List@@#&,dbModelDNA[[pos[[2]],4]]->dbBaseDNA[[pos[[1]],4]],{2}];
Table[
If[corrsResiduesDNA[[pos[[4]],i,1]]!=1,
corrsFinalDNA=MapAt[Reverse,corrsFinalDNA,i];
corrsChains=MapAt[If[Head[#]===List,Reverse[#],#]&,corrsChains,{i,2}]],
{i,Length[corrsResiduesDNA]}];
If[typeBaseDNA===Single\[Or]typeModelDNA===Single,
corrsFinalDNA=corrsFinalDNA[[All,1]];
corrsChains=MapAt[If[Head[#]===List,First[#],#]&,corrsChains,{{1,2},{2,2}}]];

corrsFinalAlphaHelix=MapThread[Plus,
{List@@corrsResiduesAlphaHelix[[pos[[3]],All,{1,-1}]],
{dbModelDNA[[pos[[2]],3,1,1]],dbBaseDNA[[pos[[1]],3,1,1]]}-1}];

{Extract[MSDs,pos],
{dbModelDNA[[pos[[2]]]],
dbBaseDNA[[pos[[1]]]],
{dbModelDNA[[pos[[2]],1]],corrsChains[[1]],corrsFinalAlphaHelix[[1]]->corrsFinalDNA[[1]]},
{dbBaseDNA[[pos[[1]],1]],corrsChains[[2]],corrsFinalAlphaHelix[[2]]->corrsFinalDNA[[2]]}},
transform,
icp[[i,j,l,2]]},

{i,Length[windowsBaseDNA]},{j,Length[windowsModelDNA]},{k,Length[corrsVerticesAlphaHelix]},{l,Length[corrsVerticesDNA]}];

(* Eliminate results that don't pass DNA MSD threshold *)
results=DeleteCases[Flatten[results,3],{_,_,_,msd_}/;msd>OptionValue[DNAMSDThreshold]];

If[Length[results]>=1,
(* Some matches remain *)
(* Eliminate equivalent results, keeping only the best in each equivalence class, and return flattened and sorted. *)
results=Union[results,
SameTest->(({#1[[2,3,;;2]],#1[[2,3,3,1]],#1[[2,4,;;2]],#1[[2,4,3,1]]}==={#2[[2,3,;;2]],#2[[2,3,3,1]],#2[[2,4,;;2]],#2[[2,4,3,1]]}\[And]((#1[[2,3,3,2]]-#1[[2,4,3,2]])===(#2[[2,3,3,2]]-#2[[2,4,3,2]])))&)];

(* Issue MSD separation warning if appropriate *)
If[Abs[Subtract@@results[[;;2,1]]]<=OptionValue[MSDSeparationWarning][results[[1,1]]],
Message[RegisterByDNA::"CloseMSDSeparation",Abs[Subtract@@results[[;;2,1]]]]];

(* Return sought result *)
If[OptionValue[ReturnDNAMSD],results[[OptionValue[ResultRank]]],results[[OptionValue[ResultRank],1;;3]]],


(* No matches remain *)
Message[RegisterByDNA::"NoMatchFound"];
If[OptionValue[ReturnDNAMSD],{\[Infinity],{},{},\[Infinity]},{\[Infinity],{},{}}]],


(* Atom types and numbers don't match (of DNA molecules), bail. *)
Message[RegisterByDNA::"NonmatchingAtoms"];
If[OptionValue[ReturnDNAMSD],{\[Infinity],{},{},\[Infinity]},{\[Infinity],{},{}}]]]


(* ::Text:: *)
(*Below are helper functions for managing registration packages*)


GetRP[idx_,filesIdxRPs_,OptionsPattern[DefaultFolders]]:=Get[ToFileName[OptionValue[FolderRegistrationPackages],filesIdxRPs[[idx]]]]


IndexDNAChain[rp_]:={Sequence@@rp[[3,1,2,2]]}
RangeDNAMatch[rp_]:=With[{range=rp[[3,1,3,2]]},If[Depth[range]===3,range,{range}]]


IndexProteinChain[rp_]:=rp[[3,1,2,1]]
RangeProteinMatch[rp_]:=rp[[3,1,3,1]]


DNAdsQ[rpOrReg__]:=Length[IndexDNAChain[rpOrReg]]===2


RPProtein[rp_]:=rp[[2]]
RPDNA[rp_]:=rp[[1]]


GetReg[idxes_List,filesIdxRPs_,OptionsPattern[DefaultFolders]]:=Get[ToFileName[OptionValue[FolderRegistrations],StringJoin[Insert[filesIdxRPs[[idxes]]," ",2]]]]


GetClusterReg[cluster_List,filesIdxRPs_,opts:OptionsPattern[DefaultFolders]]:=GetReg[{First[cluster],#},filesIdxRPs,opts]&/@Rest[cluster]


IndexDNAChain[reg_,Base]:={Sequence@@reg[[2,4,2,2]]}
RangeInnerDNAMatch[reg_,Base]:=With[{range=reg[[2,4,3,2]]},If[Depth[range]===3,range,{range}]]
RangeOuterDNAMatch[reg_,Base]:=With[{range=reg[[2,2,3,2]]},If[Depth[range]===3,range,{range}]]
IndexDNAChain[reg_,Model]:={Sequence@@reg[[2,3,2,2]]}
RangeInnerDNAMatch[reg_,Model]:=With[{range=reg[[2,3,3,2]]},If[Depth[range]===3,range,{range}]]
RangeOuterDNAMatch[reg_,Model]:=With[{range=reg[[2,1,3,2]]},If[Depth[range]===3,range,{range}]]


IndexProteinChain[reg_,Base]:=reg[[2,4,2,1]]
RangeInnerProteinMatch[reg_,Base]:=reg[[2,4,3,1]]
RangeOuterProteinMatch[reg_,Base]:=reg[[2,2,3,1]]
IndexProteinChain[reg_,Model]:=reg[[2,3,2,1]]
RangeInnerProteinMatch[reg_,Model]:=reg[[2,3,3,1]]
RangeOuterProteinMatch[reg_,Model]:=reg[[2,1,3,1]]


RegMSD[reg_]:=reg[[1]]
RegTrans[reg_]:=reg[[3]]


RegRPName[reg_,Base]:=StringReplace[StringJoin[Insert[ToString/@{
reg[[2,2,1]],
reg[[2,2,2;;3]]}," ",2]],"->"->"="]
RegRPName[reg_,Model]:=StringReplace[StringJoin[Insert[ToString/@{
reg[[2,1,1]],
reg[[2,1,2;;3]]}," ",2]],"->"->"="]


RegName[reg_]:=RegRPName[reg,Base]~~" "~~RegRPName[reg,Model]


RegProtein[reg_,pick_,OptionsPattern[DefaultFolders]]:=PDBChainsTransform[
RPProtein[Get[ToFileName[OptionValue[FolderRegistrationPackages],RegRPName[reg,pick]]]],
If[pick===Model,RegTrans[reg],TranslationTransform[{0,0,0}]]]


RegDNA[reg_,pick_,OptionsPattern[DefaultFolders]]:=PDBChainsTransform[
RPDNA[Get[ToFileName[OptionValue[FolderRegistrationPackages],RegRPName[reg,pick]]]],
If[pick===Model,RegTrans[reg],TranslationTransform[{0,0,0}]]]


RegFlippedDNAQ[reg_]:=Xor@@(OrderedQ/@{IndexDNAChain[reg,Base],IndexDNAChain[reg,Model]})
