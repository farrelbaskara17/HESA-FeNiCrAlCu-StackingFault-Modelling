clc
clear all
disp(' ');
disp('------------------------------------');
disp('SFE CALCULATION OF HESA FeNiCrAlCu');
disp('UPDATED ON FEB 2023');
disp('------------------------------------');
disp(' ');

%Input
    atFe = 5:5:40 ; 
    atNi = 5:5:40 ;
    T = 300 ; %Suhu (Kelvin)
    Tcelcius = T-273; %Suhu (Celcius)

%DATA FIX
%Atomic Percent
    ArFe=55.85; ArNi=58.69; ArCr=52; ArAl=26.98; ArCu=63.546; 

    atAl = (100-(atFe)-(atNi))./3;
    atCr = (100-(atFe)-(atNi))./3;
    atCu = (100-(atFe)-(atNi))./3;
    
    [atFe,atNi]=meshgrid(atFe,atNi);

    XFe=atFe/100 ; XNi=atNi/100 ; XCr=atCr/100; XAl=atAl/100; XCu=atCu/100;

%Pure Element Contribution
GFe = -2243.38 + 4.309.*T;
GNi = 1046 + 1.255.*T;
GCr = -2846-0.163.*T;
GAl = 2800+5.*T;
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
Tcelcius = Tcelcius .*XFe./XFe;
T = T .*XFe./XFe;

atFe = atFe.*XFe./XFe;
atNi = atNi.*XFe./XFe;
atCr = atCr.*XFe./XFe;
atAl = atAl.*XFe./XFe;
atCu = atCu.*XFe./XFe;

XFe = XFe.*XFe./XFe;
XNi = XNi.*XFe./XFe;
XCr = XCr.*XFe./XFe; 
XAl = XAl.*XFe./XFe; 
XCu = XCu.*XFe./XFe; 

GFe = GFe.*XFe./XFe;
GNi = GNi.*XFe./XFe;
GCr = GCr.*XFe./XFe;
GAl = GAl.*XFe./XFe;
GCu = GCu.*XFe./XFe;

GFeCr = GFeCr.*XFe./XFe;
GFeNi = GFeNi.*XFe./XFe;
GFeAl = GFeAl.*XFe./XFe;
GCrNi = GCrNi.*XFe./XFe;

SFE_Fe = SFE_Fe.*XFe./XFe;
SFE_Ni = SFE_Ni.*XFe./XFe;
SFE_Cr = SFE_Cr.*XFe./XFe;
SFE_Al = SFE_Al.*XFe./XFe;
SFE_Cu = SFE_Cu.*XFe./XFe;

SFE_FeCr = SFE_FeCr.*XFe./XFe;
SFE_FeNi = SFE_FeNi.*XFe./XFe;
SFE_FeAl = SFE_FeAl.*XFe./XFe;
SFE_CrNi = SFE_CrNi.*XFe./XFe;
SFE_mg = SFE_mg.*XFe./XFe;

%Output
disp(' ');
disp('---RESULT---');
format bank; %ini apa? format penulisan yang menggunakan 2 desimal dibelakang koma
    %Output 1-1
    %Result
    result = [atFe(:)'; atNi(:)'; Gchem(:)'; deltaGmg(:)'; deltaGtotal(:)'; SFE(:)']' ;
    title = '       atFe        atNi      Gchemical(J/mol)     Gmagnetic(J/mol)      deltaG(J/mol)           SFE(mJ/m2)';  
    disp(title), disp(result);

    figure(1);
    contourf(atFe,atNi,Gchem);
    xlabel('atFe'); ylabel('atNi'); zlabel('Gchemical (J/mol)');
    set(gca,'YDir','reverse');
        colorbar;
    hold on;
    figure(2);
    contourf(atFe,atNi,deltaGmg);
    xlabel('atFe'); ylabel('atNi'); zlabel('Gmagnetic (J/mol)');
    set(gca,'YDir','reverse');
        colorbar;
    hold on;
    figure(3);
    contourf(atFe,atNi,deltaGtotal)
    xlabel('atFe'); ylabel('atNi'); zlabel('deltaG(J/mol)');
    set(gca,'YDir','reverse');
        colorbar;
    hold on;
    figure(4);
    contourf(atFe,atNi,SFE)
    xlabel('at%-Fe'); ylabel('at%-Ni'); zlabel('SFE(mJ/m2)');
    set(gca,'YDir');
        colorbar;
    hold on;

