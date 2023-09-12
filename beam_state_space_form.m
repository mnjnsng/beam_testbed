clear all
close all
clc
% beam_state_space_form.m, Vivek Natarajan

% This file sets up the physical parameters of the mold oscillation test 
% bed, for use with Simulink models Beam_simulation.mdl and 
% Testbed_simulation.mdl.  This should be run before either of those
% models.

Pgain=-0.6;
Tgain=1; % make it zero if you want to stop tracking controller
d=1.27e-6;
Area=0.0046;
ps=3000*6894.75;
pt=30*6894.75;
c=1.5e-4*2; 
V_A=(2.375+0.5)*0.0254^3;
V_B=(3.8+0.5)*0.0254^3;
piston_mass=2;
beta=1500e6;
L=0.015;
fric=1000;
L_1=0.88;

Pa0=50*pt;
Pb0=50*pt;
x0=0;
xdot0=0;
t0=0;

dx=L_1/35;

%material properties and dimensions
G=77e9;
rho=7870;
E=2e11;
kd=0.83;
CSArea=0.0088/3;
I=2.2085e-5/1.25;
mass=rho*CSArea;
k=sqrt(I/CSArea);
gamma1=10;
gamma2=10;
gamma3=1;

n=L_1/dx+1;                 %number of rod nodes


C1=kd*G*CSArea/mass/dx^2;
C2=kd*G*CSArea/mass/dx/2;
C3=E*I/k^2/mass/dx^2;
C4=kd*G*CSArea/k^2/mass/dx/2;
C5=kd*G*CSArea/k^2/mass;
D1=gamma1/mass;
D2=gamma2/k^2/mass;
mold=5000*0.45;


A1=[0 0 0 0 0 1 0 0 0 0 0 0; C1 0 C2 0 -2*C1 -D1 0 0 C1 0 -C2 0; 0 0 0 0 0 0 0 1 0 0 0 0; -C4 0 C3 0 0 0 -2*C3-C5 -D2 C4 0 C3 0];
A1_prev=A1(1:4,1:4);
A1_curr=A1(1:4,5:8);
A1_next=A1(1:4,9:12);

a=ones(n-2,1);
b=ones(n-3,1);

A_left_prev=A1_prev;
A_left_curr=A1_curr;
A_left_next=A1_next;

for i=1:n-4
    A_left_prev=blkdiag(A_left_prev,A1_prev);
    A_left_curr=blkdiag(A_left_curr,A1_curr);
    A_left_next=blkdiag(A_left_next,A1_next);
end
A_left_curr=blkdiag(A_left_curr,A1_curr);
A_left=horzcat(vertcat(zeros(4,4*(n-3)),A_left_prev),zeros(4*(n-2),4))+A_left_curr+vertcat(horzcat(zeros(4*(n-3),4),A_left_next),zeros(4,4*(n-2)));
 
A_right=A_left;
%left BC for left beam
A_left_BCleft=blkdiag([0 0 0 0; 0 0  C2 0; 0 0 0 0;0 0 C3 0;],zeros(4*(n-3)));


A_left=A_left+A_left_BCleft;

%right BC for right beam

A_right_BCright=blkdiag(zeros(4*(n-3)),[0 0 0 0; 0 0 -C2 0; 0 0 0 0;0 0 C3 0;]);

MoldDyn=[zeros(1,4*(n-2)) 0 1; zeros(1,4*(n-3)) kd*G*CSArea/dx/mold 0 kd*G*CSArea/mold 0 -kd*G*CSArea/dx/mold -gamma3];

A_right=A_right+A_right_BCright;

A_right=vertcat(horzcat(A_right,[zeros(4*(n-3),2); 0 0; C1 0; 0 0; C4 0]),MoldDyn);

%middle BC
A=blkdiag(A_left,A_right);

A_left_middleBC=horzcat(zeros(4,4*(n-4)),[0 0 0 0; 0 0 -C2*(-1/6) 0; 0 0 0 0; 0 0 C3*(-1/6) 0],[0 0 0 0; 0 0 -C2*(2/3) 0; 0 0 0 0; 0 0 C3*(2/3) 0],[0 0 0 0; 0 0 -C2*(2/3) 0; 0 0 0 0; 0 0 C3*(2/3) 0],[0 0 0 0; 0 0 -C2*(-1/6) 0; 0 0 0 0; 0 0 C3*(-1/6) 0],zeros(4,4*(n-4)+2));
A_right_middleBC=horzcat(zeros(4,4*(n-4)),[0 0 0 0; 0 0 C2*(-1/6) 0; 0 0 0 0; 0 0 C3*(-1/6) 0],[0 0 0 0; 0 0 C2*(2/3) 0; 0 0 0 0; 0 0 C3*(2/3) 0],[0 0 0 0; 0 0 C2*(2/3) 0; 0 0 0 0; 0 0 C3*(2/3) 0],[0 0 0 0; 0 0 C2*(-1/6) 0; 0 0 0 0; 0 0 C3*(-1/6) 0],zeros(4,4*(n-4)+2));
A=A+vertcat(zeros(4*(n-3),4*(n-2+n-2)+2), A_left_middleBC, A_right_middleBC,zeros(4*(n-3)+2,4*(n-2+n-2)+2));

B=[0 C1 0 -C4 zeros(1,4*(n-3)+4*(n-2)+2); [repmat([0 1 0 0],1,(n-2+n-2)), 0, 1]*(-9.8)]';


C=[ 1 0 0 0 0 0 zeros(1,4*(n-2+n-2)-4); 0 0 1 0 zeros(1,4*(n-2+n-2)-2); zeros(1,4*(n-2+n-2)),1,0;zeros(1,4*(n-2+n-2)),0,1;];

D=[0, 0; 0 0; 0 0;0 0; ];
   
a1=fric/piston_mass;
a2=1/piston_mass;
a3=Area/piston_mass;
a4=4*beta*Area/(V_A+V_B);
a5=4*beta*c/(V_A+V_B);