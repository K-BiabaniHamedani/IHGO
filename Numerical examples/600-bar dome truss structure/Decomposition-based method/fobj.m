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
FF = zeros(81,81);
MM = zeros(81,81);
W = 0;
AA = [X';X'];

for i = 1:50
    XTOOL = nodes(B(i,2),1)-nodes(B(i,1),1);
    YTOOL = nodes(B(i,2),2)-nodes(B(i,1),2);
    ZTOOL = nodes(B(i,2),3)-nodes(B(i,1),3);
    ELEMTOOL = sqrt(XTOOL^2+YTOOL^2+ZTOOL^2);
    roi = sqrt(nodes(B(i,1),1)^2+nodes(B(i,1),2)^2);
    roj = sqrt(nodes(B(i,2),1)^2+nodes(B(i,2),2)^2);
    loi = nodes(B(i,1),1)/roi;
    moi = nodes(B(i,1),2)/roi;
    loj = nodes(B(i,2),1)/roj;
    moj = nodes(B(i,2),2)/roj;
    Roi = [loi -moi 0; moi loi 0; 0 0 1];
    Roj = [loj -moj 0; moj loj 0; 0 0 1];
    R = [Roi zeros(3,3);zeros(3,3) Roj];
    XT = XTOOL/ELEMTOOL;
    YT = YTOOL/ELEMTOOL;
    ZT = ZTOOL/ELEMTOOL;
    se = [XT^2,XT*YT,XT*ZT;XT*YT,YT^2,YT*ZT;XT*ZT,YT*ZT,ZT^2];
    I2 = (B(i,1)-1)*3;
    H2 = (B(i,2)-1)*3;
    EAL = AA(i,1)*E/ELEMTOOL;
    APL6 = AA(i,1)*Density*ELEMTOOL/6;
    Weight = AA(i,1)*ELEMTOOL*Density;
    W = W+12*Weight;
    SE = EAL*[se -se;-se se];
    K_cylindrical_coordinate = R'*SE*R;
    elementdof = [I2+1 I2+2 I2+3 H2+1 H2+2 H2+3];
    FF(elementdof,elementdof) = FF(elementdof,elementdof)+K_cylindrical_coordinate;
    ME = APL6*[2 0 0 1 0 0;0 2 0 0 1 0;0 0 2 0 0 1;1 0 0 2 0 0;0 1 0 0 2 0;0 0 1 0 0 2];
    M_cylindrical_coordinate = R'*ME*R;
    MM(elementdof,elementdof) = MM(elementdof,elementdof)+M_cylindrical_coordinate;
end

for j = 1:27
    MM(3*j-2,3*j-2) = MM(3*j-2,3*j-2)+Mass;
    MM(3*j-1,3*j-1) = MM(3*j-1,3*j-1)+Mass;
    MM(3*j,3*j) = MM(3*j,3*j)+Mass;
end

MM = MM(f,f);
FF = FF(f,f);
AK = FF(25:48,25:48);
AM = MM(25:48,25:48);
BK = FF(1:24,25:48);
BM = MM(1:24,25:48);
sw = [];

for J = 1:11
    LANDA1 = complex(cos(2*J*pi/24),sin(2*J*pi/24));
    LANDA2 = complex(cos(2*J*pi/24),-sin(2*J*pi/24));
    kj = AK+LANDA1*BK+LANDA2*BK';
    mj = AM+LANDA1*BM+LANDA2*BM';
    sw = [sw;sqrt(real(eig(kj,mj)))/(2*pi)'];
end

sw = [sw;sw];
LANDA1 = -1;
kj = AK+LANDA1*BK+LANDA1*BK';
mj = AM+LANDA1*BM+LANDA1*BM';
sw = [sw;sqrt(real(eig(kj,mj)))/(2*pi)'];
LANDA1 = 1;
kj = AK+LANDA1*BK+LANDA1*BK';
mj = AM+LANDA1*BM+LANDA1*BM';
sw = [sw;sqrt(real(eig(kj,mj)))/(2*pi)'];
Fr = sort(sw);
Frequencies = Fr(1:5,1);
penalty = 0;
penalty = penalty+max(0,(1-(Frequencies(1,1)/minfrq(1,1))))+max(0,(1-(Frequencies(3,1)/minfrq(1,2))));
e1 = 1;
e2 = 2;
e3 = 2.5;

ppenalty = (1+e1*penalty)^(e2+e3*((se1-1)/(se2-1)));
PFIT = W*ppenalty;

end