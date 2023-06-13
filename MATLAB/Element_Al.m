clc
clear all
disp(' ');
disp('------------------------------------');
disp('SFE CALCULATION OF HESA FeNiCrAlCu');
disp('UPDATED ON FEB 2023');
disp('------------------------------------');
disp(' ');

%Input
    atAl = 5:1:40 ; %Rentang Komposisi
    T = 300 ; %Suhu (Kelvin)
    Tcelcius = T-273; %Suhu (Celcius)

%DATA FIX
%Atomic Percent
    ArFe=55.85; ArNi=58.69; ArCr=52; ArAl=26.98; ArCu=63.546; 
    atFe = (100-(atAl))./4;
    atNi = (100-(atAl))./4;
    atCr = (100-(atAl))./4;
    atCu = (100-(atAl))./4;

    XFe=atFe/100 ; XNi=atNi/100 ; XCr=atCr/100; XAl=atAl/100; XCu=atCu/100;
    T=Tcelcius+273;

%Pure Element Contribution
GFe = -2243.38 + 4.309.*T;
GNi = 1046 + 1.255.*T;
GCr = -2846-0.163.*T;
GAl = 2800 + 5.*T;
GCu = 600 + 0.2.*T;
deltaG1 = XFe.*GFe + XNi.*GNi + XCr.*GCr + XAl.*GAl + XCu.*GCu;

%Excess Free Energy Contribution
GFeCr = 2095;
GFeNi = 2095;
GFeAl = 3326.28;
GCrNi = 27+12.1.*T+(-19895+16.38.*T).*(XCr-XNi);
deltaG2 = XFe.*XCr.*GFeCr + XFe.*XNi.*GFeNi + XFe.*XAl.*GFeAl + XCr.*XNi.*GCrNi;

%Magnetic Contribution

%Fe
tcFeFCC=-201.*XFe; %Kenapa di kali -3 Pak? ini tulisannya Tn
BFeFCC=-2.1.*XFe;


%Ni
tcNiFCC=633.*XNi; %ini ga dikali, tulisannya Tc dan Bo
tcNiHCP=633.*XNi;
BNiFCC=0.52.*XNi;
BNiHCP=0.52.*XNi;

%Cr
tcCrFCC=-1109.*XCr; %ini tulisannya Tn, dikali -3
tcCrHCP=-1109.*XCr;
BCrFCC=-2.46.*XCr;
BCrHCP=-2.46.*XCr;

%Al

%Cu

%NiCr
tcNiCrFCC = -3605.*XCr.*XNi;
BNiCrFCC = -1.91.*XCr.*XNi;

%FeCr
tcFeCrFCC = -1009.*XCr-201.*XFe;
tcFeCrHCP = -1009.*XCr-201.*XFe;
BFeCrFCC = -2.46.*XCr -2.1.*XFe;
BFeCrHCP = -2.46.*XCr -2.1.*XFe;

%Calculate tc and B, also check if it's below zero, if it's below zero,
%divide by 3, so it won't cause error in the next function

tcFCC = tcFeFCC + tcNiFCC + tcCrFCC + tcNiCrFCC + tcFeCrFCC ;
tcFCC(tcFCC<0) = tcFCC(tcFCC<0)./-3;

tcHCP = tcNiHCP + tcCrHCP + tcFeCrHCP;
tcHCP(tcHCP<0) = tcHCP(tcHCP<0)./-3;

BFCC = BFeFCC + BNiFCC + BCrFCC + BNiCrFCC + BFeCrFCC;
BFCC(BFCC<0) = BFCC(BFCC<0)./-3;

BHCP = BNiHCP + BCrHCP + BFeCrHCP;
BHCP(BHCP<0) = BHCP(BHCP<0)./-3;

tauHCP = T./ tcHCP;
tauFCC = T./ tcFCC;

D=2.34;
p=0.28;
R=8.314;

tauHCP(tauHCP<=1) = 1-( 79.*(tauHCP(tauHCP<=1)).^-1./(140.*p)+ 474./497.*(1./p-1).*( (tauHCP(tauHCP<=1)).^3./6 + (tauHCP(tauHCP<=1)).^9./135 + (tauHCP(tauHCP<=1)).^15./600 ) )./D;
tauHCP(tauHCP>1) = -( (tauHCP(tauHCP>1)).^-5./10 + (tauHCP(tauHCP>1)).^-15./315 + (tauHCP(tauHCP>1)).^-25./1500 )./D;
tauFCC(tauFCC<=1) = 1-( 79.*(tauFCC(tauFCC<=1)).^-1./(140.*p)+ 474./497.*(1./p-1).*( (tauFCC(tauFCC<=1)).^3./6 + (tauFCC(tauFCC<=1)).^9./135 + (tauFCC(tauFCC<=1)).^15./600 ) )./D;
tauFCC(tauFCC>1) = -( (tauFCC(tauFCC>1)).^-5./10 + (tauFCC(tauFCC>1)).^-15./315 + (tauFCC(tauFCC>1)).^-25./1500 )./D;

GmgHCP = R.*T.*log(BHCP+1).*tauHCP;
GmgFCC = R.*T.*log(BFCC+1).*tauFCC;
deltaGmg = GmgHCP - GmgFCC ; 

%SFE
a = 0.36e-9;
N = 6.022e23;
rho = 4./(sqrt(3)).*1./(a.^2.*N);
sigma = 8; %tanya pak tria kenapa sigma 8, emang 8
Gchem = deltaG1+deltaG2;
deltaGtotal = Gchem+deltaGmg;
SFE = 2.*rho.*deltaGtotal*1000+2.*sigma; %kenapa kali 1000? biar jadi mili joule

%SFE tiap Gchemical
SFE_Fe = GFe.*XFe.*2.*rho*1000;
SFE_Ni = GNi.*XNi.*2.*rho*1000;
SFE_Cr = GCr.*XCr.*2.*rho*1000;
SFE_Al = GAl.*XAl.*2.*rho*1000;
SFE_Cu = GCu.*XCu.*2.*rho*1000;

SFE_FeCr = XFe.*XCr.*GFeCr.*2.*rho*1000;
SFE_FeNi = XFe.*XNi.*GFeNi.*2.*rho*1000;
SFE_FeAl = XFe.*XAl.*GFeAl.*2.*rho*1000;
SFE_CrNi = XCr.*XNi.*GCrNi.*2.*rho*1000;

SFE_mg = deltaGmg.*2.*rho*1000;

%dSMix
dSmix = -(XFe.*log(XFe)+XNi.*log(XNi)+XCr.*log(XCr)+XAl.*log(XAl)+XCu.*log(XCu));

%Biar nilai konstant bisa masuk di table
Tcelcius = Tcelcius .*XAl./XAl;
T = T .*XAl./XAl;

atFe = atFe.*XAl./XAl;
atNi = atNi.*XAl./XAl;
atCr = atCr.*XAl./XAl;
atAl = atAl.*XAl./XAl;
atCu = atCu.*XAl./XAl;

XFe = XFe.*XAl./XAl;
XNi = XNi.*XAl./XAl;
XCr = XCr.*XAl./XAl;
XAl = XAl.*XAl./XAl;
XCu = XCu.*XAl./XAl;

GFe = GFe.*XAl./XAl;
GNi = GNi.*XAl./XAl;
GCr = GCr.*XAl./XAl;
GAl = GAl.*XAl./XAl;
GCu = GCu.*XAl./XAl;

GFeCr = GFeCr.*XAl./XAl;
GFeNi = GFeNi.*XAl./XAl;
GFeAl = GFeAl.*XAl./XAl;
GCrNi = GCrNi.*XAl./XAl;

SFE_Fe = SFE_Fe.*XAl./XAl;
SFE_Ni = SFE_Ni.*XAl./XAl;
SFE_Cr = SFE_Cr.*XAl./XAl;
SFE_Al = SFE_Al.*XAl./XAl;
SFE_Cu = SFE_Cu.*XAl./XAl;

SFE_FeCr = SFE_FeCr.*XAl./XAl;
SFE_FeNi = SFE_FeNi.*XAl./XAl;
SFE_FeAl = SFE_FeAl.*XAl./XAl;
SFE_CrNi = SFE_CrNi.*XAl./XAl;
SFE_mg = SFE_mg.*XAl./XAl;

%Output
disp(' ');
disp('---RESULT---');
format bank; %ini apa? format penulisan yang menggunakan 2 desimal dibelakang koma
    %Output 1-1
    %Menulis Hasil ke dalam Excel
    col_header = {'T(C)','T(K)','atFe','atNi','atCr','atAl','atCu','SFE_Fe','SFE_Ni','SFE_Cr','SFE_Al','SFE_Cu','SFE_FeNi','SFE_FeCr','SFE_FeAl','SFE_CrNi','GChemical','dgMagnetic','SFE_mg','dGTotal','SFE','SMix'};
    result = [Tcelcius;T;atFe;atNi;atCr;atAl;atCu;SFE_Fe;SFE_Ni;SFE_Cr;SFE_Al;SFE_Cu;SFE_FeNi;SFE_FeCr;SFE_FeAl;SFE_CrNi;Gchem;deltaGmg;SFE_mg;deltaGtotal;SFE;dSmix]';
    %col_header={'T(C)','XFe', 'GFe', 'SFE_Fe', 'XNi', 'GNi', 'SFE_Ni', 'XCr', 'GCr', 'SFE_Cr', 'XAl', 'GAl', 'SFE_Al', 'XCu', 'GCu', 'SFE_Cu', 'GFeNi','SFE_FeNi', 'GFeCr','SFE_FeCr', 'GFeAl','SFE_FeAl', 'GCrNi','SFE_CrNi', 'Gchemical', 'dGmg','SFE_mg', 'dGTotal', 'SFE', 'SMix'};
    %result = [Tcelcius; XFe; GFe; SFE_Fe; XNi; GNi; SFE_Ni; XCr; GCr; SFE_Cr; XAl; GAl; SFE_Al; XCu; GCu; SFE_Cu; GFeNi; SFE_FeNi; GFeCr; SFE_FeCr; GFeAl; SFE_FeAl; GCrNi; SFE_CrNi; Gchem; deltaGmg; SFE_mg; deltaGtotal; SFE; dSmix]' ;
    xlswrite('masterdata3.xlsx', col_header, 'Sheet4', 'A1');
    xlswrite('masterdata3.xlsx', result, 'Sheet4', 'A1082');


