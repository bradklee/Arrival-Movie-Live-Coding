fileNames = 
  Import["~/VersionControlled/Arrival/Arrival-Movie-Live-Coding/\
ScriptLogoJpegs/FILELIST"][[All, 1]];

ProcessFilename[filename_] := 
 StringReplace[filename, 
  Join[{".jpg" -> "", "_" -> "", 
    "1" -> ""}, # -> StringJoin[" " <> #] & /@ 
    ToUpperCase[Alphabet[]]]]

ImportIM[n_] := 
  Import["~/VersionControlled/Arrival/Arrival-Movie-Live-Coding/\
ScriptLogoJpegs/" <> fileNames[[n]] ];

ImportPixelData[n_] := ImageData[
    ColorNegate@
     Binarize[
      Import["~/VersionControlled/Arrival/Arrival-Movie-Live-Coding/\
ScriptLogoJpegs/" <> fileNames[[n]] ], .9]
    ][[Range[3300/10] 10, Range[3300/10] 10]];

xoffset = With[{PD = ImportPixelData[#]},
     Min[Cases[#, {x_, Max[#[[All, 2]] ]} :> x ]] &@( 
       Function[{a}, {a, 
          Total[(Dot[#[[1]], Reverse[#[[2]]]] &@
               Partition[RotateRight[#, a][[15 ;; -15]] , 300/2]) & /@
             PD]}] /@ Range[-50, 50])
     ] & /@ Range[38];
yoffset = With[{PD = Transpose@ImportPixelData[#]},
     Min[Cases[#, {x_, Max[#[[All, 2]] ]} :> x ]] &@( 
       Function[{a}, {a, 
          Total[(Dot[#[[1]], Reverse[#[[2]]]] &@
               Partition[RotateRight[#, a][[15 ;; -15]] , 300/2]) & /@
             PD]}] /@ Range[-50, 50])
     ] & /@ Range[38];

imDat[a_] := 
 RotateRight[#, yoffset[[a]] ] & /@ 
  Transpose@(RotateRight[#, xoffset[[a]] ] & /@ ImportPixelData[a])

AllIms = imDat /@ Range[38];
NormIms = Normalize[N@Flatten[#]] & /@ AllIms;
Overlap = Outer[Dot, NormIms, NormIms, 1];

WFOverlaps = ArrayPlot@Overlap
Match = Flatten[
   Position[Overlap[[#]], 
      Max[Complement[Overlap[[#]], Overlap[[#, {#}]]]]] & /@ 
    Range[38]];
MatchPlot = Grid[Partition[Join[MapIndexed[Labeled[
       ArrayPlot[Plus[AllIms[[#1]], AllIms[[#2[[1]] ]]]],
       Column[{ProcessFilename@fileNames[[#2[[1]]]] , 
         ProcessFilename@fileNames[[ #1]] }]
       ] &, Match], {Graphics[], Graphics[]}], 5], Frame -> All]

FourierMask[n_, F_] := 
 With[{dist = 
    Table[N@F[n ArcTan[i - 330/2 + .01, j - 330/2]], {i, 1, 330}, {j, 
      1, 330}]},
  Normalize[Flatten[dist]]
  ]

CosMasks = Flatten[N@FourierMask[#, Cos]] & /@ Range[24];
SinMasks = Flatten[N@FourierMask[#, Sin]] & /@ Range[24];

FourierT[imN_] := {
  (CosMasks[[#]].NormIms[[imN]])^2 & /@ Range[24],
  (SinMasks[[#]].NormIms[[imN]])^2 & /@ Range[24]}
AbsoluteTiming[
 FourierData = Total[FourierT[#]] & /@ Range[38];
 ]

MaxAmp = Flatten[Position[#, Max[#[[2 ;; -1]]]] & /@ FourierData];
ampPartition = (Flatten[Position[MaxAmp, #]] & /@ Union[MaxAmp]);

FourierAnalysis =  Column[{Show[
    MapThread[
     ListLinePlot[FourierData[[#]], PlotRange -> All, PlotStyle -> #2,
        ImageSize -> 250] &, {ampPartition, {Red, Orange, Green, 
       Blend[{Green, Blue}, 1/2], Blue, Purple}}],
    PlotRange -> {{1, 10}, {0, 0.12}}, AxesOrigin -> {1, 0}, 
    ImageSize -> 550],
   Grid[Partition[MapIndexed[
      Labeled[ArrayPlot[AllIms[[#1]]] , 
        StringJoin["N=", ToString[#2[[1]] + 1], 
         ProcessFilename@fileNames[[#1]]]] &, 
      MapThread[#1[[
         Position[#, Max[#]][[1, 1]]  &@
          FourierData[[#1, #2]]]] &, {ampPartition, Union[MaxAmp]}]], 
     3],
    Frame -> All]}]
