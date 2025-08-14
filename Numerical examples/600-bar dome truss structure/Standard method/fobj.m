%-----------------------------------------------------------------
% The MATLAB source code of the 600-bar dome truss structure
% Written by Kiarash Biabani Hamedani - Postdoctoral researcher at Iran University of Science and Technology
% Supervisor: Distinguished professor Ali Kaveh
%-----------------------------------------------------------------

function [X,W,Frequencies,PFIT,e1,e2,e3] = fobj(nodes,X,B,f,VarMin,VarMax,se1,se2)

X = max(min(X,VarMax),VarMin);
E = 20e10;
Density = 7850;
Mass = 100;
minfrq = [5 7];
FF = zeros(216*3,216*3);
MM = zeros(216*3,216*3);
W = 0;
AA = zeros(600,1);

for j = 1:24
        AA(25*(j-1)+1,1) = X(1,1);
    AA(25*(j-1)+2,1) = X(1,2);
    AA(25*(j-1)+3,1) = X(1,3);
    AA(25*(j-1)+4,1) = X(1,4);
    AA(25*(j-1)+5,1) = X(1,5);
    AA(25*(j-1)+6,1) = X(1,6);
    AA(25*(j-1)+7,1) = X(1,7);
    AA(25*(j-1)+8,1) = X(1,8);
    AA(25*(j-1)+9,1) = X(1,9);
    AA(25*(j-1)+10,1) = X(1,10);
    AA(25*(j-1)+11,1) = X(1,11);
    AA(25*(j-1)+12,1) = X(1,12);
    AA(25*(j-1)+13,1) = X(1,13);
    AA(25*(j-1)+14,1) = X(1,14);
    AA(25*(j-1)+15,1) = X(1,15);
    AA(25*(j-1)+16,1) = X(1,16);
    AA(25*(j-1)+17,1) = X(1,17);
    AA(25*(j-1)+18,1) = X(1,18);
    AA(25*(j-1)+19,1) = X(1,19);
    AA(25*(j-1)+20,1) = X(1,20);
    AA(25*(j-1)+21,1) = X(1,21);
    AA(25*(j-1)+22,1) = X(1,22);
    AA(25*(j-1)+23,1) = X(1,23);
    AA(25*(j-1)+24,1) = X(1,24);
    AA(25*(j-1)+25,1) = X(1,25);
end

for i = 1:size(AA,1)
    XTOOL = nodes(B(i,2),1)-nodes(B(i,1),1);
    YTOOL = nodes(B(i,2),2)-nodes(B(i,1),2);
    ZTOOL = nodes(B(i,2),3)-nodes(B(i,1),3);
    ELEMTOOL = sqrt(XTOOL^2+YTOOL^2+ZTOOL^2);
    XT = XTOOL/ELEMTOOL;
    YT = YTOOL/ELEMTOOL;
    ZT = ZTOOL/ELEMTOOL;
    se = [XT^2,XT*YT,XT*ZT;XT*YT,YT^2,YT*ZT;XT*ZT,YT*ZT,ZT^2];
    I2 = (B(i,1)-1)*3;
    H2 = (B(i,2)-1)*3;
    EAL = AA(i,1)*E/ELEMTOOL;
    APL6 = AA(i,1)*Density*ELEMTOOL/6;
    Weight = AA(i,1)*ELEMTOOL*Density;
    W = W+Weight;
    SE = EAL*[se -se;-se se];
    elementdof = [I2+1 I2+2 I2+3 H2+1 H2+2 H2+3];
    FF(elementdof,elementdof) = FF(elementdof,elementdof)+SE;
    ME = APL6*[2 0 0 1 0 0;0 2 0 0 1 0;0 0 2 0 0 1;1 0 0 2 0 0;0 1 0 0 2 0;0 0 1 0 0 2];
    MM(elementdof,elementdof) = MM(elementdof,elementdof)+ME;
end

for j = 1:216
    MM(3*j-2,3*j-2) = MM(3*j-2,3*j-2)+Mass;
    MM(3*j-1,3*j-1) = MM(3*j-1,3*j-1)+Mass;
    MM(3*j,3*j) = MM(3*j,3*j)+Mass;
end

FF = FF(f,f);
MM = MM(f,f);
sw = eig(FF,MM);
Fr = sqrt(sw)/(2*pi);
Frequencies = Fr(1:5,1);
penalty = 0;
penalty = penalty+max(0,(1-(Frequencies(1,1)/minfrq(1,1))))+max(0,(1-(Frequencies(2,1)/minfrq(1,2))));
e1 = 1;
e2 = 2;
e3 = 2.5;

ppenalty = (1+e1*penalty)^(e2+e3*((se1-1)/(se2-1)));
PFIT = W*ppenalty;

end