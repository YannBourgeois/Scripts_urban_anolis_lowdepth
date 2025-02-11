(* ::Package:: *)

(* Copy this package into .Mathematica/Applications/Coatli/ *)
(* Load package with Needs["Coatli`Nt3mom`"] *)
(* Tested with Mathematica 10.0 *)
BeginPackage["Coatli`Nt3mom`"];
E1::usage="E1[n,\[Theta]] returns the first moments E[\!\(\*SubscriptBox[\(\[Xi]\), \(i\)]\)].";
E2::usage="E2[n,\[Theta]] returns the second moments E[\!\(\*SubscriptBox[\(\[Xi]\), \(i\)]\)\!\(\*SubscriptBox[\(\[Xi]\), \(j\)]\)].";
E3::usage="E3[n,\[Theta]] returns the third moments E[\!\(\*SubscriptBox[\(\[Xi]\), \(h\)]\)\!\(\*SubscriptBox[\(\[Xi]\), \(i\)]\)\!\(\*SubscriptBox[\(\[Xi]\), \(j\)]\)].";
\[Tau]1::usage="\[Tau]1[n,\[Theta]] returns terms linear (in \[Theta]) of the third moments.";
\[Tau]2::usage="\[Tau]2[n,\[Theta]] returns terms quadratic (in \[Theta]) of the third moments.";
\[Tau]3::usage="\[Tau]3[n,\[Theta]] returns terms cubic (in \[Theta]) of the third moments.";
\[Mu]3::usage="\[Mu]3[n,\[Theta]] returns the third central moments \!\(\*SubscriptBox[\(\[Mu]\), \(3\)]\)[\!\(\*SubscriptBox[\(\[Xi]\), \(h\)]\)\!\(\*SubscriptBox[\(\[Xi]\), \(i\)]\)\!\(\*SubscriptBox[\(\[Xi]\), \(j\)]\)]."; 
E3S::usage="E3S[n,\[Theta]] returns the third moment E[\!\(\*SuperscriptBox[\(S\), \(3\)]\)] of the number S of segregating sites.";
\[Mu]3S::usage="\[Mu]3S[n,\[Theta]] returns the third central moment \!\(\*SubscriptBox[\(\[Mu]\), \(3\)]\)[S] of the number S of segregating sites.";
ta::usage="ta[n,i,j]";
tb::usage="tb[n,i,j]";
taa::usage="taa[n,h,i,j]";
tab::usage="tab[n,h,i,j]";
tba::usage="tba[n,h,i,j]";
tbb::usage="tbb[n,h,i,j]";
write3dArray::usage="write3dArray[file,m,p] writes m to file with p decimal places.";
testEquality::usage="testEquality[n] checks if simplified txx yields the same as unsimplified txx.";
Begin["`Private`"];
(* functions for the second moments *)
a[n_]:=HarmonicNumber[n-1];
\[Alpha][n_,i_]:=(a[n]-a[i])/(n-i);
\[Beta][n_,i_]:=If[i<n,(2n)/((n-i+1)(n-i)) (a[n+1]-a[i])-2/(n-i),0];
(* additional functions for the third moments *) 
\[Alpha]2[n_,i_,j_]:=\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(k = 2\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(t = 1\), \(k - 1\)]\((
\*FractionBox[\(Binomial[i - 1, t - 1] Binomial[n - i - j, k - t - 1]\), \(Binomial[n - 1, k - 1]\)] 
\*FractionBox[\(1\), \(k \((k - 1)\)\)] \[Alpha][k, t])\)\)\);
\[Beta]2[n_,i_,j_]:=\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(k = 2\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(t = 1\), \(k - 1\)]\((
\*FractionBox[\(Binomial[i - 1, t - 1] Binomial[n - i - j, k - t - 1]\), \(Binomial[n - 1, k - 1]\)] 
\*FractionBox[\(1\), \(k \((k - 1)\)\)] 
\*FractionBox[\(\[Beta][k, t]\), \(2\)])\)\)\);
\[Alpha]3[n_,h_,i_,j_]:=(h+1)\[Alpha]2[n,i,j]-2h \[Alpha]2[n,i,j+1]+(h-1)\[Alpha]2[n,i,j+2];
\[Beta]3[n_,h_,i_,j_]:=(h+1)\[Beta]2[n,i,j]-2h \[Beta]2[n,i,j+1]+(h-1)\[Beta]2[n,i,j+2];
\[Alpha]4[n_,h_,i_,j_]:=(h+1)\[Alpha]2[n,i+1,j]-2h \[Alpha]2[n,i,j+1]+(h-1)\[Alpha]2[n,i-1,j+2];
\[Beta]4[n_,h_,i_,j_]:=(h+1)\[Beta]2[n,i+1,j]-2h \[Beta]2[n,i,j+1]+(h-1)\[Beta]2[n,i-1,j+2];

(* Tensor with terms linear in \[Theta] *)
\[Tau]1[n_,\[Theta]_ : 1]:=\[Theta] Table[Table[Boole[h==i&&i==j]1/i,{j,1,n-1},{i,1,n-1}],{h,1,n-1}];

(* Vector of the first moments *)
E1[n_,\[Theta]_ : 1]:=Table[1/i \[Theta],{i,1,n-1}];

(* Second moments *)
ta[n_,i_,j_]:=\[Piecewise]{
 {((\[Beta][n,j]-\[Beta][n,j+1])/2), j<i},
 {(\[Beta][n,j]/2), i==j}
}
tb[n_,i_,j_]:=\[Piecewise]{
 {1/(i j)-1/(i(i+j))-(\[Beta][n,j]-\[Beta][n,j+1])/2, i+j<n},
 {\[Alpha][n,j]-\[Beta][n,j]/2, i+j==n}
}
\[Tau]2[n_,i_,j_]:=ta[n,i,j]+ta[n,j,i]+tb[n,i,j]+tb[n,j,i];

(* Tensor with terms quadratic in \[Theta] *)
\[Tau]2[n_,\[Theta]_ : 1]:=\[Theta]^2 Table[Table[Table[Boole[i==j]\[Tau]2[n,h,i]+Boole[h==i]\[Tau]2[n,j,i]+Boole[j==h]\[Tau]2[n,i,h],{j,1,n-1}],{i,1,n-1}],{h,1,n-1}];

(* Matrix of the second moments *)
E2[n_, \[Theta]_ : 1]:=Table[\[Tau]2[n,i,j]+Boole[i==j]1/i \[Theta],{i,1,n-1},{j,1,n-1}];

(* Formulas of Klassmann and Ferretti 2016 for the third moments *)
(* The unsimplified formulas *)
myBinomial[n_,k_]:=If[n<0||k<0,If[n==-1 && k==-1,1,0],Binomial[n,k]];
taaSlow[n_,h_,i_,j_]:=\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(kpp = 2\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(kp = 2\), \(kpp\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(k = 2\), \(kp\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(t = 1\), \(kp - 1\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(t1 = 0\), \(kpp - 2\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(t2 = 1\), \(kpp - t1\)]\((k\ kp\ kpp 
\*FractionBox[\(myBinomial[kp - t - 1, k - 2]\), \(Binomial[kp - 1, k - 1]\)] 
\*FractionBox[\(t\), \(kp\)] 
\*FractionBox[\(myBinomial[t1 - 1, t - 2] myBinomial[kpp - t1 - t2 - 1, kp - t - 1]\), \(Binomial[kpp - 1, kp - 1]\)] 
\*FractionBox[\(t2\), \(kpp\)] 
\*FractionBox[\(myBinomial[h - i - 1, t1 - 1] myBinomial[i - j - 1, t2 - 2] myBinomial[n - h - 1, kpp - t1 - t2 - 1]\), \(Binomial[n - 1, kpp - 1]\)] 
\*FractionBox[\(1\), \(k \((k - 1)\) kp \((kp - 1)\) kpp \((kpp - 1)\)\)])\)\)\)\)\)\)\);
tabSlow[n_,h_,i_,j_]:=\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(kpp = 2\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(kp = 2\), \(kpp\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(k = 2\), \(kp\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(t = 2\), \(kp - 1\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(t1 = 1\), \(kpp - 2\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(t2 = 1\), \(kpp - t1\)]\((k\ kp\ kpp 
\*FractionBox[\(myBinomial[kp - t - 1, k - 2]\), \(Binomial[kp - 1, k - 1]\)] 
\*FractionBox[\(t\), \(kp\)] 
\*FractionBox[\(myBinomial[t1 - 1, t - 2] myBinomial[kpp - t1 - t2 - 1, kp - t - 1]\), \(Binomial[kpp - 1, kp - 1]\)] 
\*FractionBox[\(t1\), \(kpp\)] 
\*FractionBox[\(myBinomial[h - i - j - 1, t1 - 2] myBinomial[i - 1, t2 - 1] myBinomial[n - h - 1, kpp - t1 - t2 - 1]\), \(Binomial[n - 1, kpp - 1]\)] 
\*FractionBox[\(1\), \(k \((k - 1)\) kp \((kp - 1)\) kpp \((kpp - 1)\)\)])\)\)\)\)\)\)\);
tba3Slow[n_,h_,i_,j_]:=\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(kpp = 2\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(kp = 2\), \(kpp\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(k = 2\), \(kp\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(t = 1\), \(kp - 1\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(t1 = 0\), \(kpp - 2\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(t2 = 1\), \(kpp - t1\)]\((k\ kp\ kpp 
\*FractionBox[\(myBinomial[kp - t - 1, k - 2]\), \(Binomial[kp - 1, k - 1]\)] 
\*FractionBox[\(t\), \(kp\)] 
\*FractionBox[\(myBinomial[t1 - 1, t - 2] myBinomial[kpp - t1 - t2 - 1, kp - t - 1]\), \(Binomial[kpp - 1, kp - 1]\)] 
\*FractionBox[\(kpp - t1 - t2\), \(kpp\)] 
\*FractionBox[\(myBinomial[h - i - 1, t1 - 1] myBinomial[i - 1, t2 - 1] myBinomial[n - h - j - 1, kpp - t1 - t2 - 2]\), \(Binomial[n - 1, kpp - 1]\)] 
\*FractionBox[\(1\), \(k \((k - 1)\) kp \((kp - 1)\) kpp \((kpp - 1)\)\)])\)\)\)\)\)\)\);
tba2Slow[n_,h_,i_,j_]:=\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(kpp = 2\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(kp = 2\), \(kpp\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(k = 2\), \(kp\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(t = 1\), \(kp - 1\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(t1 = 1\), \(kpp - 1\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(t2 = 1\), \(kpp - t1\)]\((k\ kp\ kpp 
\*FractionBox[\(myBinomial[kp - t - 1, k - 2]\), \(Binomial[kp - 1, k - 1]\)] 
\*FractionBox[\(kp - t\), \(kp\)] 
\*FractionBox[\(myBinomial[t1 - 1, t - 1] myBinomial[kpp - t1 - t2 - 1, kp - t - 2]\), \(Binomial[kpp - 1, kp - 1]\)] 
\*FractionBox[\(t1\), \(kpp\)] 
\*FractionBox[\(myBinomial[h - j - 1, t1 - 2] myBinomial[i - 1, t2 - 1] myBinomial[n - h - i - 1, kpp - t1 - t2 - 1]\), \(Binomial[n - 1, kpp - 1]\)] 
\*FractionBox[\(1\), \(k \((k - 1)\) kp \((kp - 1)\) kpp \((kpp - 1)\)\)])\)\)\)\)\)\)\);
tba1Slow[n_,h_,i_,j_]:=\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(kpp = 2\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(kp = 2\), \(kpp\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(k = 2\), \(kp\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(t = 1\), \(kp - 1\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(t1 = 1\), \(kpp - 1\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(t2 = 1\), \(kpp - t1\)]\((k\ kp\ kpp 
\*FractionBox[\(myBinomial[kp - t - 1, k - 2]\), \(Binomial[kp - 1, k - 1]\)] 
\*FractionBox[\(kp - t\), \(kp\)] 
\*FractionBox[\(myBinomial[t1 - 1, t - 1] myBinomial[kpp - t1 - t2 - 1, kp - t - 2]\), \(Binomial[kpp - 1, kp - 1]\)] 
\*FractionBox[\(t2\), \(kpp\)] 
\*FractionBox[\(myBinomial[h - 1, t1 - 1] myBinomial[i - j - 1, t2 - 2] myBinomial[n - h - i - 1, kpp - t1 - t2 - 1]\), \(Binomial[n - 1, kpp - 1]\)] 
\*FractionBox[\(1\), \(k \((k - 1)\) kp \((kp - 1)\) kpp \((kpp - 1)\)\)])\)\)\)\)\)\)\);
tbaSlow[n_,h_,i_,j_]:=tba3Slow[n,i,j,h]+tba2Slow[n,i,h,j]+tba1Slow[n,h,i,j];
tbbSlow[n_,h_,i_,j_]:=\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(kpp = 2\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(kp = 2\), \(kpp\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(k = 2\), \(kp\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(t = 1\), \(kp - 1\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(t1 = 1\), \(kpp - 2\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(t2 = 1\), \(kpp - t1\)]\((k\ kp\ kpp 
\*FractionBox[\(myBinomial[kp - t - 1, k - 2]\), \(Binomial[kp - 1, k - 1]\)] 
\*FractionBox[\(kp - t\), \(kp\)] 
\*FractionBox[\(myBinomial[t1 - 1, t - 1] myBinomial[kpp - t1 - t2 - 1, kp - t - 2]\), \(Binomial[kpp - 1, kp - 1]\)] 
\*FractionBox[\(kpp - t1 - t2\), \(kpp\)] 
\*FractionBox[\(myBinomial[h - 1, t1 - 1] myBinomial[i - 1, t2 - 1] myBinomial[n - h - i - j - 1, kpp - t1 - t2 - 2]\), \(Binomial[n - 1, kpp - 1]\)] 
\*FractionBox[\(1\), \(k \((k - 1)\) kp \((kp - 1)\) kpp \((kpp - 1)\)\)])\)\)\)\)\)\)\);

(* the simplified formulas *)
taa[n_,h_,i_,j_]:=\[Piecewise]{
 {\[Beta]4[n,i-j,i-j,j]-\[Beta]4[n,i-j,i-j+1,j], j<i && i<h},
 {\[Beta]4[n,i-j,i-j,j], j<i &&i==h},
 {\[Beta]2[n,1,j]-\[Beta]2[n,2,j], j==i && i<h},
 {\[Beta]2[n,1,j], j==i && i==h}
}
tab[n_,h_,i_,j_]:=\[Piecewise]{
 {\[Beta]3[n,h-i-j,i,j]-\[Beta]3[n,h-i-j,i+1,j], i+j<h},
 {\[Beta]2[n,i,j]-\[Beta]2[n,i+1,j], i+j==h}
}
tba[n_,h_,i_,j_]:=\[Piecewise]{
 {1/(2h) (\[Beta][n,j]-\[Beta][n,j+1]-\[Beta][n,h+j]+\[Beta][n,h+j+1])-\[Beta]4[n,i-j,i-j,j]+\[Beta]4[n,i-j,i-j+1,j]-\[Beta]3[n,i-j,h,j]+\[Beta]3[n,i-j,h+1,j]+\[Beta]3[n,n-h-i,j,h]-\[Beta]3[n,n-h-i,j+1,h], j<i && h+i<n},
 {\[Alpha]4[n,n-h-j,n-h-j,j]-\[Beta]4[n,n-h-j,n-h-j,j]+\[Alpha]3[n,n-h-j,h,j]-\[Beta]3[n,n-h-j,h,j]+\[Beta]2[n,j,h]-\[Beta]2[n,j+1,h], j<i && h+i==n},
 {1/(2h) (\[Beta][n,j]-\[Beta][n,h+j])+\[Beta]3[n,n-h-j,j,h]-\[Beta]2[n,h,j]+\[Beta]2[n,h+1,j]-\[Beta]2[n,1,j]+\[Beta]2[n,2,j], j==i && h+i<n},
 {1 /2 (\[Alpha]2[n,n-j,j]+\[Alpha]2[n,j,n-j])+\[Alpha]2[n,1,j]-\[Beta]2[n,1,j], j==i && h+i==n}
}
tbb[n_,h_,i_,j_]:=\[Piecewise]{
 {1/((h+i+j)(h+i)h)+1 /(i j(h+i))-1 /(i h(i+j))-1/(2i) (\[Beta][n,j]-\[Beta][n,j+1]-\[Beta][n,i+j]+\[Beta][n,i+j+1])-\[Beta]3[n,n-h-i-j,i,j]+\[Beta]3[n,n-h-i-j,i+1,j], h+i+j<n},
 {1 /i (\[Alpha][n,j]-\[Alpha][n,i+j])-1/(2i) (\[Beta][n,j]-\[Beta][n,i+j])-\[Beta]2[n,i,j]+\[Beta]2[n,i+1,j], h+i+j==n}
}

t[n_,h_,i_,j_]:=taa[n,h,i,j]+tab[n,h,i,j]+tba[n,h,i,j]+tbb[n,h,i,j];

(* Tensor with terms cubic in \[Theta] *)
\[Tau]3[n_,\[Theta]_:1]:=\[Theta]^3 Table[Table[t[n,h,i,j]+t[n,i,j,h]+t[n,j,h,i]+t[n,h,j,i]+t[n,j,i,h]+t[n,i,h,j],{i,1,n-1},{j,1,n-1}],{h,1,n-1}];

(* Third moments *)
E3[n_,\[Theta]_ : 1]:=\[Tau]1[n,\[Theta]]+\[Tau]2[n,\[Theta]]+\[Tau]3[n,\[Theta]];

(* Third central moments *)
\[Mu]3[n_,\[Theta]_ : 1]:=Module[{e3=E3[n,\[Theta]],e2=E2[n,\[Theta]],e1=E1[n,\[Theta]]},Table[Table[e3[[h,i,j]]-e1[[h]]e2[[i,j]]-e1[[i]]e2[[h,j]]-e1[[j]]e2[[h,i]]+2e1[[h]]e1[[i]]e1[[j]],{h,1,n-1},{i,1,n-1}],{j,1,n-1}]];

(* Third moment of segregating sites S *)
E3S[n_,\[Theta]_ : 1]:=HarmonicNumber[n-1,1]\[Theta]+3(HarmonicNumber[n-1,1]^2+HarmonicNumber[n-1,2])\[Theta]^2+(HarmonicNumber[n-1,1]^3+3HarmonicNumber[n-1,1]HarmonicNumber[n-1,2]+2HarmonicNumber[n-1,3])\[Theta]^3;

(* Third central moment of segregating sites S *)
\[Mu]3S[n_, \[Theta]_ : 1]:=HarmonicNumber[n-1,1]\[Theta]+3HarmonicNumber[n-1,2]\[Theta]^2+2HarmonicNumber[n-1,3]\[Theta]^3;

testEquality[n_]:=Module[{},{Table[taa[n,h,i,j],{h,1,n-1},{i,1,n-1},{j,1,n-1}]===Table[taaSlow[n,h,i,j],{h,1,n-1},{i,1,n-1},{j,1,n-1}],Table[tab[n,h,i,j],{h,1,n-1},{i,1,n-1},{j,1,n-1}]===Table[tabSlow[n,h,i,j],{h,1,n-1},{i,1,n-1},{j,1,n-1}],Table[tba[n,h,i,j],{h,1,n-1},{i,1,n-1},{j,1,n-1}]===Table[tbaSlow[n,h,i,j],{h,1,n-1},{i,1,n-1},{j,1,n-1}],Table[tbb[n,h,i,j],{h,1,n-1},{i,1,n-1},{j,1,n-1}]===Table[tbbSlow[n,h,i,j],{h,1,n-1},{i,1,n-1},{j,1,n-1}]}];

(* Print cubic array *)
write3dArray[file_,m__,p_:4]:=Module [{n=Length[m]},
f=OpenWrite[file];
Do[Do[Do[WriteString[f,NumberForm[m[[h,i,j]]//N,{p+1,p}]];If[i<n,WriteString[f," "],WriteString[f,"\n"]],{i,1,n}],{j,1,n}];If[h<n,WriteString[f,"\n"]],{h,1,n}];
Close[f]];

End[];
EndPackage[];
