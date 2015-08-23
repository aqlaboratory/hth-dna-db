(* ::Package:: *)

(* ::Text:: *)
(*Accessory functions:*)


PointLineDistance[point_,line_]:=Norm[(line[[2]]-line[[1]])\[Cross](line[[1]]-point)]/Norm[line[[2]]-line[[1]]]


ProjectionPointOntoLine[point_,line_]:=line[[1]]+Norm[point-line[[1]]](Normalize[point-line[[1]]].Normalize[line[[2]]-line[[1]]])Normalize[line[[2]]-line[[1]]]


MinPosition[list_]:=First[MinPosition[list,1]]
MaxPosition[list_]:=First[MaxPosition[list,1]]
MinPosition[list_,n_]:=Ordering[list,n]
MaxPosition[list_,n_]:=Reverse[Ordering[list,-n]]


DNAQ[chain_]:=Or@@((StringLength[#]===2\[And]StringTake[#,1]==="D")&/@chain[[All,1]])
RNAQ[chain_]:=Or@@((StringLength[#]===2\[And]StringTake[#,1]==="R")&/@chain[[All,1]])
ProteinQ[chain_]:=({3}===Union[StringLength/@chain[[All,1]]])\[And]\[Not]MemberQ[chain[[All,1]],"6mp"|"Ade"|"Cpt"|"Gol"|"Gun"|"Hpa"|"Hsl"|"Nh4"|"Ooa"]
BiomoleculeQ[chain_]:=MatchQ[Head[chain],DNA|RNA|Protein]
ChainType[chain_]:=Switch[{DNAQ[chain],RNAQ[chain],ProteinQ[chain]},
{True,False,False},DNA,
{False,True,False},RNA,
{False,False,True},Protein,
{False,False,False},Unidentified]


ChainAtomLabels[chain_,pos_]:=chain[[1,pos,2,All,1;;2]]
ChainAtomLabels[chain_,pos_,vertexTypes_]:=With[
{pattern=Alternatives@@(cap:{Sequence@@#,_}&/@vertexTypes):>Most[cap]},
Cases[chain[[1,pos]],pattern,\[Infinity]]]


ChainVertices[chain_,pos_]:=chain[[1,pos,2,All,3]]
ChainVertices[chain_,pos_,vertexTypes_]:=With[
{pattern=Alternatives@@({Sequence@@#,coord_}&/@vertexTypes)->coord},
Cases[chain[[1,pos]],pattern,\[Infinity]]]


(* ::Text:: *)
(*Returns the file name out of a directory.*)


ExtractFileName[file_String]:=StringDrop[file,StringLength[DirectoryName[file]]]


(* ::Text:: *)
(*Returns an expression of the format {ChainType[{{residue, {atom, atomRole, atomCoords}}..},options]..}, where ChainType can be DNA, RNA, or Protein.*)


PDBChainsImport::"nofile"="File `1` not found.";
PDBChainsImport[file_String]:=If[FileNames[ExtractFileName[file],DirectoryName[file]]==={},
Message[PDBChainsImport::"nofile",file],
PDBChainsImport[Import[file,{"PDB",{"Residues","ResidueAtoms","ResidueRoles","ResidueCoordinates","SecondaryStructure"}}]]]
PDBChainsImport[elemPDB:{_,_,_,_,_}]:=Module[{chains,opts,raw=elemPDB},
opts=Options[raw];
raw=Drop[raw,-1];
chains=Table[{raw[[1,i,j]],Table[raw[[k,i,j,l]],{l,Length[raw[[2,i,j]]]},{k,2,Length[raw]}]},
{i,Dimensions[raw][[2]]},{j,Length[raw[[1,i]]]}];
chains=DeleteCases[#,{"ZN"|"CO",_}]&/@chains;
chains=(ChainType[#][#])&/@chains;
If[("Helix"/.opts)==={},opts={"Helix"->Table[{},{Length[chains]}]}];
chains=Table[Append[chains[[i]],Helix->("Helix"/.opts)[[i]]],{i,Length[chains]}]]


PDBChainsExport[file_String,chains_List]:=Module[{strExport,listExport},
strExport=ExportString[
{chains[[All,1,All,1]],
chains[[All,1,All,2,All,1]],
ReplacePart[chains[[All,1,All,2,All,2]],Position[chains[[All,1,All,2,All,1]],"H"]->""],
SetAccuracy[chains[[All,1,All,2,All,3]],2],
Take[CharacterRange["A","Z"],Length[chains]]},
{"PDB",{"Residues","ResidueAtoms","ResidueRoles","ResidueCoordinates","ResidueChainLabels"}}];
listExport=ImportString[strExport,"List"];
Export[
file,
StringReplace[listExport,RegularExpression["\\-(\\d{3,3})\\.(\\d{3,3})"]:>StringDrop[" -$1.$2",-1]],
"List"]]
PDBChainsExport[file_String,chain_]:=PDBChainsExport[file,{chain}]
