clc
clear all
disp(' ');
disp('------------------------------------');
disp('SFE CALCULATION OF HESA FeNiCrAlCu');
disp('UPDATED ON MAR 2023');
disp('------------------------------------');
disp(' ');

%Input
    T = [1, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000, 1025, 1050, 1075, 1100, 1125, 1150, 1175, 1200]; %Kelvin
    Tcelcius = T-273; %Celcius

%DATA FIX
%Atomic Percent
    ArFe=55.85; ArNi=58.69; ArCr=52; ArAl=26.98; ArCu=63.546; 
    atFe = 20;
    atNi = 20;
    atCr = 20;
    atAl = 20;
    atCu = 20;

    XFe=atFe/100 ; XNi=atNi/100 ; XCr=atCr/100; XAl=atAl/100; XCu=atCu/100;

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
GFeAl = 3328;
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
Tcelcius = Tcelcius .*T./T;
dSmix = dSmix .*T./T;

atFe = atFe.*T./T;
atNi = atNi.*T./T;
atCr = atCr.*T./T;
atAl = atAl.*T./T;
atCu = atCu.*T./T;

XFe = XFe.*T./T;
XNi = XNi.*T./T;
XCr = XCr.*T./T; 
XAl = XAl.*T./T; 
XCu = XCu.*T./T; 

GFe = GFe.*T./T;
GNi = GNi.*T./T;
GCr = GCr.*T./T;
GAl = GAl.*T./T;
GCu = GCu.*T./T;

GFeCr = GFeCr.*T./T;
GFeNi = GFeNi.*T./T;
GFeAl = GFeAl.*T./T;
GCrNi = GCrNi.*T./T;

SFE_Fe = SFE_Fe.*T./T;
SFE_Ni = SFE_Ni.*T./T;
SFE_Cr = SFE_Cr.*T./T;
SFE_Al = SFE_Al.*T./T;
SFE_Cu = SFE_Cu.*T./T;

SFE_FeCr = SFE_FeCr.*T./T;
SFE_FeNi = SFE_FeNi.*T./T;
SFE_FeAl = SFE_FeAl.*T./T;
SFE_CrNi = SFE_CrNi.*T./T;
SFE_mg = SFE_mg.*T./T;

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
    xlswrite('masterdataT.xlsx', col_header, 'Sheet1', 'A1');
    xlswrite('masterdataT.xlsx', result, 'Sheet1', 'A2');  


