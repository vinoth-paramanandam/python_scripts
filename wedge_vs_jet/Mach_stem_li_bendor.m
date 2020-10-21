% Mach stem height in Weak Mach domain iusing Li and Ben-Dor Method 


clc 
clear all
close all
global M0 M1 Theta1 P0 T0 R Rho0 P1 g Theta3 M2 P2 a  Nu_M1 Nu_Md  M_ A Mc Md  Phi1 Phi2 Phic Mc_d W H L 


g=1.4;% gama
%across incident shock
P0=101325;
T0=300;
R=288;
Rho0=P0/(R*T0);

M0=  input('flow Mach number = ');
W= input('length of wedge, in meter = ');
H= 1; %input('height between symmetric line and leading edge of the wedge,in meter= ');
a=(g+1)/(g-1);

% X0=[1.6 45*(pi/180) 20*(pi/180) 200000 1.4 45*(pi/180) 20*(pi/180) 300000];
% X=fsolve(@vonNew_,X0);
% Theta_von=X(7)*(180/pi);

Theta1_= input ('wedge angle = '); 
Theta1=Theta1_*(pi/180);
L= W*cos(Theta1)  ;%input('')

%% PRESSURE POLAR FOR THE INITIAL CONDITION


mu=asin(1/M0);
Phi2=mu:.001:pi/2;


for i=1:length(Phi2)
    M1n(i)=M0*sin(Phi2(i));
    Theta_right(i)= atan((2.*cot(Phi2(i)).*(M1n(i).^2-1))/(M0^2*(g+cos(2.*Phi2(i)))+2));
    P10_r(i)= 1+((2*g*(M1n(i).^2-1))/(g+1));
    i=i+1;
end

Phi2=pi/2:-0.001:mu;

for i=1:length(Phi2)
    M1n(i)=M0*sin(Phi2(i));
    Theta_left(i)= -atan((2.*cot(Phi2(i)).*(((M1n(i).^2)-1))/(M0^2*(g+cos(2.*Phi2(i)))+2)));
    P10_l(i)= 1+((2*g*(M1n(i).^2-1))/(g+1));
    i=i+1;
end
Theta_1=[Theta_right*(180/pi), Theta_left*(180/pi)];
P10=[P10_r, P10_l];

while P10==1
    Theta_1=0;
end

%reflection shock polar

X0=40*(pi/180);
Phi=fsolve(@T_B, X0);
P10_ref=1+(((M0*sin(Phi))^2-1)*(2*g/(g+1)));
M1=(1/(sin(Phi-Theta1)))*sqrt((1+(((g-1)/2)*(M0*sin(Phi))^2))/((g*(M0*sin(Phi))^2)-((g-1)/2)));

mu=asin(1/M1);
Phi2=mu:.001:pi/2;


for i=1:length(Phi2)
    M1n(i)=M1*sin(Phi2(i));
    Theta2_right(i)= atan((2.*cot(Phi2(i)).*(M1n(i).^2-1))/(M1^2*(g+cos(2.*Phi2(i)))+2));
    P21_r(i)= 1+((2*g*(M1n(i).^2-1))/(g+1));
    P20_r(i)=P21_r(i).*P10_ref;
    i=i+1;
end

Phi2=pi/2:-0.001:mu;

for i=1:length(Phi2)
    M1n(i)=M1*sin(Phi2(i));
    Theta2_left(i)= -atan((2.*cot(Phi2(i)).*(M1n(i).^2-1))/(M1^2*(g+cos(2.*Phi2(i)))+2));
    P21_l(i)= 1+((2*g*(M1n(i).^2-1))/(g+1));
    P20_l(i)=P21_l(i).*P10_ref;
    i=i+1;
end

Theta_2r=[(Theta2_right+Theta1)*(180/pi),(Theta2_left+Theta1)*(180/pi)];
Theta_2l=[-(Theta2_right+Theta1)*(180/pi),-(Theta2_left+Theta1)*(180/pi)];
P21=[P21_r, P21_l];
P20=[P20_r, P20_l];
while P21==1
    Theta_2r=0;
    Theta_2l=0;
end

%ploting shock polar
plot(Theta_1,P10,'b','Linewidth',1.5)
hold on
plot(Theta_2r,P20,'r','Linewidth',1.5)
hold on 
% plot(Theta_2l,P20,'g','Linewidth',1.5)
% hold on
legend('i polar','r polar +ve');
legend('location','southeast')
xlabel('Theta')
ylabel('P/P0')
hold on
grid on


%%
% finding initial conditions for  solving three shock theory

M0_=M0;
Phii = fsolve(@T_B,pi/4);
M1i = sqrt(((1+((g-1)*M0^2*sin(Phii)*sin(Phii)))+((((g+1)^2/4)-g*sin(Phii)*sin(Phii))*M0^4*(sin(Phii))^2))/(((g*M0^2*(sin(Phii))^2)-((g-1)/2))*(((M0^2*(g-1)/2)*(sin(Phii))^2)+1)));
P1i = P0*((2/(g+1))*((g*M0^2*(sin(Phii))^2)-(g-1)/2));
Rho1i= Rho0*(((g+1)*M0^2*(sin(Phii))^2)/(((g-1)*M0^2*(sin(Phii))^2)+2));
T1i=(P1i/(R*Rho1i));
M0=M1i;
Phii_ = fsolve(@T_B,pi/4);
M2i = sqrt(((1+((g-1)*M0^2*sin(Phii_)*sin(Phii_)))+((((g+1)^2/4)-g*sin(Phii_)*sin(Phii_))*M0^4*(sin(Phii_))^2))/(((g*M0^2*(sin(Phii_))^2)-((g-1)/2))*(((M0^2*(g-1)/2)*(sin(Phii_))^2)+1)));
P2i = P0*((2/(g+1))*((g*M0^2*(sin(Phii_))^2)-(g-1)/2));
Rho2i= Rho0*(((g+1)*M0^2*(sin(Phii_))^2)/(((g-1)*M0^2*(sin(Phii_))^2)+2));
T2i=(P2i/(R*Rho2i));

M0=M0_;
M1n = sqrt(((1+((g-1)*M0^2))+((((g+1)^2/4)-g)*M0^4))/(((g*M0^2)-((g-1)/2))*(((M0^2*(g-1)/2))+1)));
P1n = P0*((2/(g+1))*((g*M0^2)-(g-1)/2));
Rho1n= Rho0*(((g+1)*M0^2)/(((g-1)*M0^2)+2));
T1n=(P1n/(R*Rho1n));

% X0=[M1i Phii 15*(pi/180) P1i M2i Phii_ 15*(pi/180) P2i];
% % X0=[1.6 45*(pi/180) 20*(pi/180) 200000 1.4 45*(pi/180) 20*(pi/180) 300000];
% X=fsolve(@vonNew_,X0);
% Theta_von=X(3)*(180/pi);





%% solving of three shock theory

X0=[M1i Phii P1i Rho1i T1i M2i Phii_ Theta1 (P2i) (Rho2i) (T2i) M1n 90*(pi/180) 5*(pi/180) P1n Rho1n T1n];
% opts=optimoptions(@fsolve,'Display','off','Algorithm','trust-region-dogleg')
X=fsolve(@MRMach,X0);

% defining parameters obtained from three shock theory

M1=X(1);
Phi1=X(2);
P1=X(3);
Rho1=X(4);
T1=X(5);
M2=X(6);
Phi2=X(7);
Theta2=X(8);
P2=X(9);
Rho2=X(10);
T2=X(11);
M3=X(12);
Phi3=X(13);
Theta3=X(14);
P3=X(15);
Rho3=X(16);
T3=X(17);

%%    %% for curved Mach stem
% for foot of the Mach stem is point G


Mg=sqrt(((1+((g-1)*M0^2*sin(pi/2)*sin(pi/2)))+((((g+1)^2/4)-g*sin(pi/2)*sin(pi/2))*M0^4*(sin(pi/2))^2))/(((g*M0^2*(sin(pi/2))^2)-((g-1)/2))*(((M0^2*(g-1)/2)*(sin(pi/2))^2)+1)));
Pg=P0*((2/(g+1))*((g*M0^2*(sin(pi/2))^2)-(g-1)/2));
Rhog= Rho0*(((g+1)*M0^2*(sin(pi/2))^2)/(((g-1)*M0^2*(sin(pi/2))^2)+2));
Tg=(Pg/(R*Rhog));

a3=sqrt(g*R*T3);
ag=sqrt(g*R*Tg);
Ug=Mg*ag;
U3=M3*a3;

% Average values

Rho_= (Rho3+Rhog)/2;
a_=(a3+ag)/2;
U_=(1/(2*Rho_))*((Rho3*U3*cos(Theta3))+(Rhog*Ug));
M_=U_/a_; 

%%              %Expansion fan interactions

Nu_M1=(sqrt(a))*(atan(sqrt((M1^2-1)/a)))-(atan(sqrt(M1^2-1)));
Nu_M2=(sqrt(a))*(atan(sqrt((M2^2-1)/a)))-(atan(sqrt(M2^2-1)));
Nu_Md=Nu_M2+Theta3;


% Y0=[45*(pi/180) 2.0 5*(pi/180) 200000 2.2 45*(pi/180) 5*(pi/180) 100000 2.0 100000];
Y0=[Nu_M1 M1 Theta1 P1 M1 Phi2 Theta1 P1 M2 P2];
Y=fsolve(@Exp,Y0);

Nu_Mc = Y(1);
   Mc = Y(2);
    A = Y(3);
   Pc = Y(4);
 Mc_d = Y(5);
 Phic = Y(6);
Thetac_d = Y(7);
 Pc_d = Y(8);
   Md = Y(9);
   Pd = Y(10);
   
   
   
   
   
   %% Geometrical relations

R0=[H L H L H H L H L H L H L H];
R=fsolve(@Geom_relation,R0);
   
yB =R(1);
xB =R(2);
yC =R(3);
xC =R(4);
yF =R(5);
Hs =R(6);
xF =R(7);
yE =R(8);
xE =R(9);
yD =R(10);
xD =R(11);
yT =R(12);
xT =R(13);
Hm =R(14);
% disp(Theta_von) 
disp('Hm/H');
disp(Hm/H);


mub = asin(1/M1);
Htmax_rr= (W*sin(mub+Theta1)*sin(Phi1-Theta1))/(sin(mub+Theta1-Phi1));
Htmin_rr= (W*sin(Phi2-Theta1)*sin(Phi1-Theta1))/(sin(Phi2-Theta1+Phi1));
Htmax_mr= Hm + (W*sin(mub+Theta1)*sin(Phi1-Theta1))/(sin(mub+Theta1-Phi1));
Htmin_mr = Hm + (W*sin(Phi2-Theta1)*sin(Phi1-Theta1))/(sin(Phi2-Theta1+Phi1));
disp('Htmin_rr, Htmax_rr')
disp(Htmin_rr)
disp(Htmax_rr)
disp('Htmin_mr, Htmax_mr')
disp(Htmin_mr)
disp(Htmax_mr)

%%  FUNCTIONS for the code


function F=vonNew_(X)
global M0 g P0 

M1=X(1);
Phi1=X(2);
Theta1=X(3);
P1=X(4);

M2=X(5);
Phi2=X(6);
Theta2=X(7);
P2=X(8);


F(1)=M1-sqrt(((1+((g-1)*M0^2*sin(Phi1)*sin(Phi1)))+((((g+1)^2/4)-g*sin(Phi1)*sin(Phi1))*M0^4*(sin(Phi1))^2))/(((g*M0^2*(sin(Phi1))^2)-((g-1)/2))*(((M0^2*(g-1)/2)*(sin(Phi1))^2)+1)));
F(2)=Theta1-atan((2*cot(Phi1)*(M0^2*(sin(Phi1))^2-1))/(M0^2*(g+cos(2*Phi1))+2));
F(3)=P1-P0*((2/(g+1))*((g*M0^2*(sin(Phi1))^2)-(g-1)/2));

P10=((2/(g+1))*((g*M0^2)-(g-1)/2));
%across reflected shock

F(4)=M2-sqrt(((1+((g-1)*M1^2*sin(Phi2)*sin(Phi2)))+((((g+1)^2/4)-g*sin(Phi2)*sin(Phi2))*M1^4*(sin(Phi2))^2))/(((g*M1^2*(sin(Phi2))^2)-((g-1)/2))*(((M1^2*(g-1)/2)*(sin(Phi2))^2)+1)));
F(5)=Theta1-atan((2*cot(Phi2)*(M1^2*(sin(Phi2))^2-1))/(M1^2*(g+cos(2*Phi2))+2));
F(6)=P2-P1*((2/(g+1))*((g*M1^2*(sin(Phi2))^2)-(g-1)/2));
F(7)=Theta1-Theta2;
F(8)=(P2/P0)-P10;
end


function L=TBM_(Phi)
global M0 Theta1 g

L=Theta1-atan(2*cot(Phi)*(((M0*sin(Phi))^2-1)/(M0^2*(g+cos(2*Phi))+2)));
end

%% function for Three shock theory


function F=MRMach(X)
global M0 Theta1 g P0 R Rho0

M1=X(1);
Phi1=X(2);
P1=X(3);
Rho1=X(4);
T1=X(5);


M2=X(6);
Phi2=X(7);
Theta2=X(8);
P2=X(9);
Rho2=X(10);
T2=X(11);

M3=X(12);
Phi3=X(13);
Theta3=X(14);
P3=X(15);
Rho3=X(16);
T3=X(17);

%across incident shock wave

F(1)=M1-sqrt(((1+((g-1)*M0^2*sin(Phi1)*sin(Phi1)))+((((g+1)^2/4)-g*sin(Phi1)*sin(Phi1))*M0^4*(sin(Phi1))^2))/(((g*M0^2*(sin(Phi1))^2)-((g-1)/2))*(((M0^2*(g-1)/2)*(sin(Phi1))^2)+1)));
F(2)=Theta1-atan((2*cot(Phi1)*(M0^2*(sin(Phi1))^2-1))/(M0^2*(g+cos(2*Phi1))+2));
F(3)=P1-P0*((2/(g+1))*((g*M0^2*(sin(Phi1))^2)-(g-1)/2));
F(4)=Rho1- Rho0*(((g+1)*M0^2*(sin(Phi1))^2)/(((g-1)*M0^2*(sin(Phi1))^2)+2));
F(5)=T1-(P1/(R*Rho1));

%across reflected shock

F(6)=M2-sqrt(((1+((g-1)*M1^2*sin(Phi2)*sin(Phi2)))+((((g+1)^2/4)-g*sin(Phi2)*sin(Phi2))*M1^4*(sin(Phi2))^2))/(((g*M1^2*(sin(Phi2))^2)-((g-1)/2))*(((M1^2*(g-1)/2)*(sin(Phi2))^2)+1)));
F(7)=Theta2-atan((2*cot(Phi2)*(M1^2*(sin(Phi2))^2-1))/(M1^2*(g+cos(2*Phi2))+2));
F(8)=P2-P1*((2/(g+1))*((g*M1^2*(sin(Phi2))^2)-(g-1)/2));
F(9)=Rho2- Rho1*(((g+1)*M1^2*(sin(Phi2))^2)/(((g-1)*M1^2*(sin(Phi2))^2)+2));
F(10)=T2-(P2/(R*Rho2));

%across Mach stem

F(11)=M3-sqrt(((1+((g-1)*M0^2*sin(Phi3)*sin(Phi3)))+((((g+1)^2/4)-g*sin(Phi3)*sin(Phi3))*M0^4*(sin(Phi3))^2))/(((g*M0^2*(sin(Phi3))^2)-((g-1)/2))*(((M0^2*(g-1)/2)*(sin(Phi3))^2)+1)));
F(12)=Theta3-atan((2*cot(Phi3)*(M0^2*(sin(Phi3))^2-1))/(M0^2*(g+cos(2*Phi3))+2));
F(13)=P3-P0*((2/(g+1))*((g*M0^2*(sin(Phi3))^2)-(g-1)/2));
F(14)=Rho3- Rho0*(((g+1)*M0^2*(sin(Phi3))^2)/(((g-1)*M0^2*(sin(Phi3))^2)+2));
F(15)=T3-(P3/(R*Rho3));
F(16)=Theta1-Theta2-Theta3;
F(17)=P3-P2;
end



%%  function for Expansion fan interaction

function Q=Exp(Y)

global M1 a P1 g M2 P2 Theta1 Nu_M1 Nu_Md 

Nu_Mc = Y(1);
   Mc = Y(2);
    A = Y(3);
   Pc = Y(4);
 Mc_d = Y(5);
 Phic = Y(6);
Thetac_d = Y(7);
 Pc_d = Y(8);
   Md = Y(9);
   Pd = Y(10);

Q(1)= Nu_Mc-((sqrt(a)*atan(sqrt((Mc^2-1)/a)))-atan(sqrt(Mc^2-1)));
Q(2)= Nu_Mc-Nu_M1-Theta1+A;
Q(3)=Pc-P1*(((2+((g-1)*(M1^2)))/(2+((g-1)*(Mc^2))))^(g/(g-1)));
Q(4)= Mc_d-sqrt(((1+((g-1)*Mc^2*sin(Phic)*sin(Phic)))+((((g+1)^2/4)-g*sin(Phic)*sin(Phic))*Mc^4*(sin(Phic))^2))/(((g*Mc^2*(sin(Phic))^2)-((g-1)/2))*(((Mc^2*(g-1)/2)*(sin(Phic))^2)+1)));
Q(5)= Thetac_d-atan((2*cot(Phic)*(Mc^2*(sin(Phic))^2-1))/(Mc^2*(g+cos(2*Phic))+2));
Q(6)= Pc_d-Pc*((2/(g+1))*((g*Mc^2*(sin(Phic))^2)-(g-1)/2));
Q(7)= A-Thetac_d;
Q(8)= Nu_Md-(sqrt(a)*atan(sqrt((Md^2-1)/a)))+atan(sqrt(Md^2-1));
Q(9)= Pd-P2*(((2+(g-1)*M2^2)/(2+(g-1)*Md^2))^(g/(g-1)));
Q(10)= Pc_d-Pd;
end


 %%  Function for geometrical relations
   
function Z=Geom_relation(R)
global A M1 Mc M2 Md Theta3 Theta1 W Phi1 Phi2 Phic Mc_d M_ H L g 

mub = asin(1/M1);
muc = asin(1/Mc);
mu2 = asin(1/M2);
mud = asin(1/Md);
muc_=asin(1/Mc_d);
xR = L ; %W*cos(Theta1);
yR = H-W*sin(Theta1);

yB =R(1);
xB =R(2);
yC =R(3);
xC =R(4);
yF =R(5);
Hs =R(6);
xF =R(7);
yE =R(8);
xE =R(9);
yD =R(10);
xD =R(11);
yT =R(12);
xT =R(13);
Hm =R(14);
% mu_c_d=asin(1/Mc_d);
Z(1)= yB-yR+(xB-xR)*tan(mub+Theta1);
Z(2)= yC-yR+(tan(muc+A))*(xC-xR);
Z(3)= yF-yB+(tan(mu2+Theta3))*(xF-xB);
Z(4)= yE-Hs;
Z(5)= yE-yD+tan(mud)*(xE-xD);
Z(6)= yB-yT-tan(Phi2-Theta1)*(xB-xT);
Z(7)= yT-Hm;
Z(8)= xT-(H-Hm)*cot(Phi1);
Z(9)= yF-yT+tan(Theta3)*(xF-xT);
Z(10)= yB-yC-(xB-xC)*tan(delta((Phi2-Theta1),(Phic-A)));
% Z(11)= yC-yD-(xC-xD)*tan(delta((-),(-asin(1/Md))));
Z(11)= yC-yD-(xC-xD)*tan(delta((-muc_),(-asin(1/Md))));
Z(12)= yB-yD-(xB-xD)*tan(delta((-Theta3),(0)));
Z(13)= yF-yE-(xF-xE)*tan(delta((-Theta3),(0)));
Z(14)= Hm-(Hs/M_)*((2+((g-1)*M_^2))/(g+1))^((g+1)/(2*(g-1)));
end

function d = delta(d1,d2)
d = atan((2*tan(d1) + tan(d2-d1))/(2 - tan(d1)*tan(d2-d1)));
end
function F=T_B(Phi)
global M0 Theta1 g

F=Theta1-atan(2.*cot(Phi)*(((M0.*sin(Phi)).^2-1)/(M0^2*(g+cos(2.*Phi))+2)));
end