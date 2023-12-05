(* ::Package:: *)

BeginPackage["CVTUtilites`"]


(*Regions: Global variable*)
Clear[regs]; 
regs=<|"Square"->BoundaryDiscretizeRegion@Rectangle[{-1,-1},{1,1}],
"Disk"->BoundaryDiscretizeRegion@Disk[]|>;

(*Distance between arrays of points*)
Clear[arrDistC,arrDistCTot];
arrDistC = Compile[{{arr1,_Real,2},{arr2,_Real,2}},
	(nf=Nearest[arr1];Max[SquaredEuclideanDistance[First@nf@#,#]&/@arr2])];
arrDistCTot = Compile[{{arr1,_Real,2},{arr2,_Real,2}},
	(nf=Nearest[arr1];Total[SquaredEuclideanDistance[First@nf@#,#]&/@arr2])];

(*Graphics*)
Clear[tinyPlots,tinyPlot,tinyRegionPlot,tinyRegionPlots,regionPlots];
tinyPlot[lst_,regID_,size_]:= Graphics[{FaceForm@White,EdgeForm@LightBlue,regs[regID],PointSize->0.07,GrayLevel[.2],Point@lst},
	ImageSize->size]/;MemberQ[Keys[regs],regID];

tinyPlots[lst_,regID_,size_,num_]:= GraphicsGrid[Partition[
Graphics[{FaceForm@White,EdgeForm@LightBlue,regs[regID],PointSize->0.07,GrayLevel[.2],Point@#},
	ImageSize->size]&/@lst,UpTo[num]]]/;MemberQ[Keys[regs],regID];
	
tinyRegionPlot[arr_,regID_,size_]:= With[{mesh=vMesh[arr,regID]},
 Graphics[{EdgeForm@Gray,FaceForm@LightBlue,mesh,
	  PointSize@.04,Point@arr},ImageSize->size]]/;MemberQ[Keys[regs],regID];
	  
tinyRegionPlots[lst_,regID_,size_,num_]:= With[{meshes=vMesh[#,regID]&/@lst},
 GraphicsGrid[Partition[
Table[
                    Graphics[{EdgeForm@Gray,FaceForm@LightBlue,meshes[[k]],
	                           PointSize@.04,Point@lst[[k]]},
ImageSize->size],
{k,Length@lst}],
	UpTo[num]]]]/;MemberQ[Keys[regs],regID];
	
regionPlots[lst_,regID_]:=With[{meshes=vMesh[#,regID]&/@lst,colors=Rescale[Range[Length@lst[[1]]]]},
 GraphicsGrid[Partition[Table[
    Graphics[{EdgeForm@Black,MapThread[{FaceForm[#1],#2}&, {ColorData["DarkBands"]/@colors,meshes[[k]]}],
	PointSize@.04,Point@lst[[k]]},ImageSize->80],{k,Length@lst}],
	UpTo[5]]]]/;MemberQ[Keys[regs],regID];

(*Collect possible centroids from K-Means clustering of given region*)
Clear[dataClust];
dataClust[numSeeds_,\[Delta]_,numTries_,regID_]:=
With[{pts=RandomPoint[regs[regID],\[LeftFloor]8Sqrt[2]/\[Delta]^2\[RightFloor]]},
	SeedRandom[];
	Reap[Monitor[
		Do[Sow[Mean/@FindClusters[pts,numSeeds,
			Method->{"KMeans",
			"InitialCentroids"-> RandomSample[pts,numSeeds]},
			PerformanceGoal->"Speed"]],
		{k,numTries}],
ProgressIndicator[k,{1,numTries}]]][[2,1]]]/;MemberQ[Keys[regs],regID];

(*Dihedral group elements*)
Clear[dihedralTrans,transTypes];
dihedralTrans = Join[
	RotationTransform[#] & /@ {2 \[Pi], \[Pi]/2, \[Pi], -\[Pi]/2}, 
	ReflectionTransform[#] & /@ {{-1, 1}, {1, 0}, {1, 1}, {0, 1}}
	];
transTypes=MapThread[Rule[#1,#2]&,{dihedralTrans, {"\[ScriptE]","\[ScriptA]","\!\(\*SuperscriptBox[\(\[ScriptA]\), \(2\)]\)","\!\(\*SuperscriptBox[\(\[ScriptA]\), \(3\)]\)","\[ScriptX]","\[ScriptA] \[ScriptX]","\!\(\*SuperscriptBox[\(\[ScriptA]\), \(2\)]\)\[ScriptX]","\!\(\*SuperscriptBox[\(\[ScriptA]\), \(3\)]\)\[ScriptX]"}}];
	
(*Ultimate distance with respect to possible dihedral group symmetry*)
Clear[patternDistSquare,patternDistSquareTot];
patternDistSquare[arr1_,arr2_]:=With[{trans = Through[dihedralTrans[arr2]]}, 
	Min[arrDistC[arr1,#]&/@trans]]/;MatrixQ[arr1,NumericQ]&&MatrixQ[arr2,NumericQ];
patternDistSquareTot[arr1_,arr2_]:=With[{trans = Through[dihedralTrans[arr2]]}, 
	Min[arrDistCTot[arr1,#]&/@trans]]/;MatrixQ[arr1,NumericQ]&&MatrixQ[arr2,NumericQ];
	
(*Polar distances*)
Clear[myFPC,polarRotateC,polarDistC,polarRotDist];
myFPC[{r_,\[Theta]_}]:={r Cos[\[Theta]], r Sin[\[Theta]]};
(*Distances between patterns in polar coords*)
polarRotateC = Compile[{{arr,_Real,2},{\[Theta],_Real}},{#[[1]],Mod[#[[2]]+\[Theta],2\[Pi]]}&/@arr];
polarDistC = Compile[{{arr1,_Real,2},{arr2,_Real,2}},
       Max[Min/@DistanceMatrix[arr1,arr2, DistanceFunction->ChessboardDistance]]];
polarDistCTot = Compile[{{arr1,_Real,2},{arr2,_Real,2}},
       Total[Min/@DistanceMatrix[arr1,arr2, DistanceFunction->ChessboardDistance]]];
(*Distance between patterns rotated by angle*)
polarRotDist[arr1_,arr2_,\[Theta]_]:=polarDistC[arr1, polarRotateC[arr2,\[Theta]]]/;
	MatrixQ[arr1,NumericQ]&&MatrixQ[arr2,NumericQ]&&NumericQ[\[Theta]];

(*Ultimate distance with respect to possible rotation symmetry*)
Clear[patternDistPolar,patternDistPolarTot];
patternDistPolar[arr1_,arr2_]:=
  Min[polarDistC[arr1,polarRotateC[arr2,#]]&/@Range[0,2\[Pi],2 \[Delta]]]/;
	MatrixQ[arr1,NumericQ]&&MatrixQ[arr2,NumericQ];
	
patternDistPolarTot[arr1_,arr2_]:=
  Min[polarDistCTot[arr1,polarRotateC[arr2,#]]&/@Range[0,2\[Pi],2 \[Delta]]]/;
	MatrixQ[arr1,NumericQ]&&MatrixQ[arr2,NumericQ];


(*Voronoi mesh*)
Clear[vMesh,cellsOrdered];
vMesh[pts_,regID_]:=With[{mesh=BoundaryDiscretizeRegion/@MeshPrimitives[VoronoiMesh[pts,{{-1,1},{-1,1}}],2],
reg=regs[regID]
},RegionIntersection[#,reg]&/@mesh
];

(*Ordering cells with respect of seeds*)
cellsOrdered[cells_,pts_]:=cells[[SparseArray[Outer[#2@#1&,pts,RegionMember/@cells,1],
Automatic,False]["NonzeroPositions"][[All,2]]]];

(*Lloyd process*)
Clear[lloydStepOrdered,lloydStepFast,lloydProc];
lloydStepOrdered[pts_,regID_]:=With[{cells=cellsOrdered[vMesh[pts,regID],pts]},
RegionCentroid/@cells];
lloydStepFast[pts_,regID_]:=RegionCentroid/@vMesh[pts,regID];
lloydProc[pts_,tol_,maxIt_,regID_,step_]:=NestWhile[step[#,regID]&,
pts,arrDistCTot[#1,#2]>tol&,2,maxIt];

(*Align patterns*)
Clear[alignPattsSquare];
alignPattsSquare[arr_]:=First@MinimalBy[
  Through[dihedralTrans[#]],
  arrDistCTot[#,First@arr]&]&/@arr;

Clear[alignPattsDisk];
alignPattsDisk[arr_]:=Table[
  First@MinimalBy[
    polarRotateC[e,#]&/@ Range[0,2\[Pi],\[Delta]],
    polarDistCTot[#,First@arr]&],
 {e,arr}];



EndPackage[ ]
