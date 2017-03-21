MultiFactorial[n_, nDim_] := 
 Times[n, If[n - nDim > 1, MultiFactorial[n - nDim, nDim], 1]]
GeneralT[n_, m_] := 
 Table[(-m)^(-j) MultiFactorial[i + m (j - 1) + 1, m]/
    MultiFactorial[i + 1, m], {i, 1, n}, {j, 1, i}]
a[n_] := With[{gt = GeneralT[2 n, 2]}, 
  gt[[2 #, Range[#]]] & /@ Range[n]]

c[n_ /; OddQ[n]] := c[n] = 0;
c[n_ /; EvenQ[n]] := c[n] = 2 (n!) (-2)^(n/2)/(n + 2)!;

B2[0, 0] = 1;
B2[n_ /; n > 0, 0] := 0;
B2[0, k_ /; k > 0] := 0;
B2[n_ /; n > 0, k_ /; k > 0] := 
  B2[n, k] = 
   Total[Binomial[n - 1, # - 1] c[#] B2[n - #, k - 1] & /@ 
     Range[1, n - k + 1]];

BasisT[n_] := 
 Table[B2[i, j]/(i!) Q^(i + 2 j), {i, 2, 2 n, 2}, {j, 1, i/2}]
PhaseSpaceExpansion[n_] := 
  Times[Sqrt[2 \[Alpha]], 
   1 + Dot[MapThread[Dot, {BasisT[n], a[n]}], (2 \[Alpha])^
      Range[n]]];
AbsoluteTiming[CES50 = PhaseSpaceExpansion[50];] (* about 2(s)*)

Fast50 = Compile[{{\[Alpha], _Real}, {Q, _Real}}, Evaluate@CES50];


fileNames = 
  Import["~/VersionControlled/Arrival/Arrival-Movie-Live-Coding/\
ScriptLogoJpegs/FILELIST"][[All, 1]];

ImportIM[n_] := 
  Import["~/VersionControlled/Arrival/Arrival-Movie-Live-Coding/\
ScriptLogoJpegs/" <> fileNames[[n]] ];

ImportPixelData[n_] := ImageData[ColorNegate@Binarize[
      Import[
       "~/VersionControlled/Arrival/Arrival-Movie-Live-Coding/\
ScriptLogoJpegs/" <> fileNames[[n]] ], .9]
    ][[10 Range[3300/10], 10 Range[3300/10]]] ;

centers = 
  Import["~/VersionControlled/Arrival/Arrival-Movie-Live-Coding/\
ProcessedData/centers.csv"];

BG = PolarPlot[
   Evaluate[
    CES50 /. {Q -> Cos[\[Phi]], \[Alpha] -> #/10} & /@ 
     Range[9]], {\[Phi], 0, 2 Pi}, Axes -> False, PlotStyle -> Gray];

EK50 = Normal@
   Series[D[
     Expand[CES50^2/2] /. 
      Q^n_ :> (1/2)^n Binomial[n, n/2], \[Alpha]], {\[Alpha], 0, 
     50}];
SameQ[Normal@
  Series[(2/Pi) EllipticK[\[Alpha]], {\[Alpha], 0, 50}], EK50]
Plot[{(2/Pi) EllipticK[\[Alpha]], EK50}, {\[Alpha], .9, 1}, 
 ImageSize -> 500]

tDIY = Mean[
   AbsoluteTiming[Fast50[.9, RandomReal[{0, 1}]] ][[1]] & /@ 
    Range[10000]];
tMma = SetPrecision[
   Mean[AbsoluteTiming[JacobiSN[.9, RandomReal[{0, 1}]] ][[1]] & /@ 
     Range[10000]], 4];
tMma/tDIY

BasisT2[n_] := 
  Table[BellY[i, j, c /@ Range[2 n]]/(i!) Q^(i + 2 j), {i, 2, 2 n, 
    2}, {j, 1, i/2}];
SameQ[BasisT2[20], BasisT[20]]
t1 = AbsoluteTiming[BasisT[#];][[1]] & /@ Range[100]
t2 = AbsoluteTiming[BasisT2[#];][[1]] & /@ Range[25]
ListLinePlot[{t1, t2}, ImageSize -> 500]

ListLinePlot[{t1, t2}, ImageSize -> 500]

H[n_, rep_] := 
 ReplaceRepeated[(1/2) q^2 + (1/2) p^2 + 
   Total[(1/2) c[2 #]/((2 #)!) (q)^(2 + 2 #) & /@ Range[n]], rep]

\[Phi]Dot[n_, rep_] := 
 ReplaceRepeated[
  Divide[-Expand[D[H[n, {}], p] p + D[H[n, {}], q] q], p^2 + q^2], rep]

CES50Squared = Expand[Normal@Series[CES50^2, {\[Alpha], 0, 50}]];
Clear@Pow\[CapitalPsi];
Pow\[CapitalPsi][0] = 1;
Pow\[CapitalPsi][n_] := 
 Pow\[CapitalPsi][n] = 
  Expand[Normal@
    Series[Pow\[CapitalPsi][n - 2]*CES50Squared, {\[Alpha], 0, 50}]]

AbsoluteTiming[Pow\[CapitalPsi][100];]

Expand[H[50, {p^2 -> \[CapitalPsi]^2 - q^2, 
    q -> \[CapitalPsi] Q}] /. {Power[\[CapitalPsi], n_] :> 
    Pow\[CapitalPsi][n]}]

AbsoluteTiming[
 \[Phi]Dot50 = 
   Expand[Expand[\[Phi]Dot[
       50, {p^2 -> \[CapitalPsi]^2 - q^2, 
        q -> \[CapitalPsi] Q}]] /. {Power[\[CapitalPsi], n_] :> 
       Pow\[CapitalPsi][n]}];
 ]

FastCES50 = Compile[{{\[Alpha], _Real}, {Q, _Real}}, Evaluate@CES50];
Fast\[Phi]Dot50 = 
  Compile[{{\[Alpha], _Real}, {Q, _Real}}, Evaluate@\[Phi]Dot50];

 AbsoluteTiming[Fast\[Phi]Dot50[.1, .5]]
 AbsoluteTiming[FastCES50[.9, 0]]

Show[
 ListLinePlot@
  NestList[# + Fast\[Phi]Dot50[.9, Cos[#]] (.1) &, 0, 
   Floor[2 Pi 2/Pi EllipticK[.9] 10]],
 Plot[N[-Pi], {x, 0, 2000}]
 ]


AlienWavefunction[R_, normRad_, Qs_, angles_] := Module[{
   deformedRadii = MapThread[Fast50, {R normRad, Qs}]
   }, Transpose[{deformedRadii, angles}]]

BG = PolarPlot[
   Evaluate[
    CES50 /. {Q -> Cos[\[Phi]], \[Alpha] -> #/10} & /@ 
     Range[9]], {\[Phi], 0, 2 Pi}, Axes -> False, PlotStyle -> Gray];

DepictAlienWF[Coords_] := Show[
  BG, Graphics[Point[#] & /@ Coords],
  ImageSize -> 500, PlotRange -> {{-1.9, 1.9}, {-1.4, 1.4}}]

ExportTimeEvo[tmax_, tnow_, angles_] := 
 Module[{newAngles = angles, tnext = tnow}, Do[
   Export[
    "~/VersionControlled/Arrival/Arrival-Movie-Live-Coding/Brad/\
TimeData/coords" <> StringPadLeft[ToString[tnext], 4, "0"] <> ".csv", 
    AlienWavefunction[.9, normRadii, Cos /@ newAngles, newAngles ]];
   newAngles = MapThread[
     Function[{i\[Alpha], i\[Phi]}, 
      i\[Phi] + Fast\[Phi]Dot50[.9 i\[Alpha], Cos[i\[Phi]]] (.02)],
     {normRadii, newAngles}]; tnext = tnext + 1, tmax]]

Logogram01 = ImportPixelData[32];
cent = centers[[32]];
ArrayPlot@Logogram01;

Positions1 = Position[Logogram01, 1];
onePosCentered = 
  N[With[{cent = {3300/10/2, 3300/10/2} }, # - cent & /@ 
     Positions1]];
radii = Norm /@ onePosCentered;
maxR = Max@radii;
normRadii = radii/maxR;
angles = ArcTan[#[[2]], #[[1]]] & /@ onePosCentered;
Qs = Cos /@ angles;

Length@radii
DepictAlienWF[
 AlienWavefunction[.9, normRadii, Qs, 
   angles] /. {r_, z_} :> {r Cos[z], r Sin[z]}]

(*AbsoluteTiming[
ExportTimeEvo[10,1,angles]
]*)

(*Export["~/VersionControlled/Arrival/Arrival-Movie-Live-Coding/Brad/\
AnimationFrames/im"<> StringPadLeft[ToString[#],4,"0"]<>".gif",
 DepictAlienWF[
Import["~/VersionControlled/Arrival/Arrival-Movie-Live-Coding/Brad/\
TimeData/coords"<> StringPadLeft[ToString[#],4,"0"]<>".csv"]/.{
{r_,z_}\[RuleDelayed]{r Cos[z],r Sin[z]}
}]]&/@Range[10];*)

AbsoluteTiming[WFData =
   Normalize[N[Length /@ BinLists[#, {-1.9, 1.9, .05} ]]] & /@ 
      Transpose[
       Import[
         "~/VersionControlled/Arrival/Arrival-Movie-Live-Coding/Brad/\
TimeData/coords" <> StringPadLeft[ToString[#], 4, "0"] <> 
          ".csv"] /. {r_, \[Phi]_} :> {r Cos[\[Phi]], r Sin[\[Phi]]}
       ] & /@ Range[10];]

(*Export["~/VersionControlled/Arrival/Arrival-Movie-Live-Coding/Brad/\
qWave.csv",WFData[[All,1]] ]
Export["~/VersionControlled/Arrival/Arrival-Movie-Live-Coding/Brad/\
pWave.csv",WFData[[All,2]] ]*)

WFQData = 
  Import["~/VersionControlled/Arrival/Arrival-Movie-Live-Coding/Brad/\
qWave.csv" ];
WFPData = 
  Import["~/VersionControlled/Arrival/Arrival-Movie-Live-Coding/Brad/\
pWave.csv"];

(*With[{bg=ListLinePlot[WFQData[[1]],PlotStyle\[Rule]Gray]},Export[
"~/VersionControlled/Arrival/Arrival-Movie-Live-Coding/Brad/qWave.gif",
Show[bg, ListLinePlot[#],PlotRange\[Rule]{0,0.5},PlotStyle\[Rule]\
Black]&/@WFQData[[Range[1,100]5]]
]]*)

(*With[{bg=ListLinePlot[WFQData[[2000]],PlotStyle\[Rule]Gray]},Export[\

"~/VersionControlled/Arrival/Arrival-Movie-Live-Coding/Brad/qWave2.\
gif",
Show[bg, ListLinePlot[#],PlotRange\[Rule]{0,0.5},PlotStyle\[Rule]\
Black]&/@WFQData[[Range[2000/5,2500/5]5]]
]]*)

QAC = WFQData[[1]].# & /@ WFQData ;
PAC = WFPData [[1]].# & /@ WFPData ;

Row[{
  Show[
   ListLinePlot[PAC[[1 ;; -1]], PlotRange -> {0, 1}, 
    PlotStyle -> Blue],
   ListLinePlot[QAC[[1 ;; -1]], PlotRange -> {0, 1}, 
    PlotStyle -> Green, ImageSize -> 300], ImageSize -> 300
   ],
  Show[
   ListLinePlot[
    Transpose[{Range[2000, 5000] , PAC[[Range[2000, 5000] ]]}], 
    PlotRange -> {.7, 1}, PlotStyle -> Blue],
   ListLinePlot[
    Transpose[{Range[2000, 5000] , QAC[[Range[2000, 5000] ]]}], 
    PlotRange -> {.7, 1}, PlotStyle -> Green],
   ImageSize -> 300
   ]}]

(*Export[
"~/VersionControlled/Arrival/Arrival-Movie-Live-Coding/Brad/\
Autocorrelation.png",
]*)

FourierMetric[data_, n_] := 
 With[{dat = Normalize /@ Partition[data - Mean[data], n]},
  Mean@Flatten[
    Table[Dot[dat[[i]], dat[[j]]], {i, 1, Length[dat]}, {j, 1, i}]]]

Spectrum = Show[
   ListLinePlot[{}, PlotRange -> {{0, 5}, {0, 1}}, 
    PlotStyle -> Blue],
   
   Graphics[
    Line[{{N[2 Pi/(375*0.02)] * #, 0}, {N[2 Pi/(375*0.02)]* #, 
         500}}] & /@ {1, 3, 5}],
   Graphics[{Dashed, 
     Line[{{N[2 Pi/(375*0.02)] * #, 0}, {N[2 Pi/(375*0.02)]* #, 
          500}}] & /@ {2, 4, 6}}],
   
   ListLinePlot[{2 Pi/(#*0.02), 
       FourierMetric[PAC[[2000 ;; -1]], #]} & /@ Range[50, 1000], 
    PlotRange -> All, PlotStyle -> Blue],
   ListLinePlot[{2 Pi/(#*0.02), 
       FourierMetric[QAC[[2000 ;; -1]], #]} & /@ Range[50, 1000], 
    PlotRange -> All, PlotStyle -> Green],
   
   PlotRange -> {{0, 5}, {0, 1}},
   AspectRatio -> 1/2,
   ImageSize -> 550
   ];

SpectrumPartitions = 
  Grid[Partition[
    Labeled[ListLinePlot[Partition[PAC[[2000 ;; -1]], Floor[#]], 
        Axes -> False, ImageSize -> 100],
       SetAccuracy[2 Pi/(#*0.02), 3]] & /@ Flatten[{
       373 + Range[-2, 2],
       (1/2) 373 + Range[-2, 2],
       (1/3) 373 + Range[-2, 2],
       (2/3) 373 + Range[-2, 2]
       }], 5], Frame -> All];

Column[{Spectrum, SpectrumPartitions}]

(*Export[
"~/VersionControlled/Arrival/Arrival-Movie-Live-Coding/Brad/\
TimeFourierAnalysis.png",
Column[{Spectrum,SpectrumPartitions}]]*)
