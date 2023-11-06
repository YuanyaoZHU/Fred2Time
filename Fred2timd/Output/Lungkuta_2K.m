clear;
clc;
%% 一阶
dt=0.01;
T = 100;
N=T/dt;
Xn=zeros(1,N);
Xn(1)=5;
t=zeros(1,N);
for i=1:N
   t(i)=(i-1)*dt; 
end
for i=1:N-1
   X0=Xn(i);
   y= lgkt4k(t(i),dt,X0);
   Xn(i+1)=y;
end
figure(1);
plot(t,Xn);

%% 2阶常微分
dt=0.01;
T = 100;
N=T/dt;
Vn=zeros(1,N);
Sn=zeros(1,N);
Vn(1)=0;
Sn(1)=0;
t=zeros(1,N);
for i=1:N
   t(i)=(i-1)*dt; 
   X(i)=-cos(t(i))+1;
end
for i=1:N-1
   S0=Sn(i);
   V0=Vn(i);
   [Sn(i+1),Vn(i+1)]= lgkt4k2k(t(i),dt,S0,V0);
end


figure(2);
subplot(2,1,1);
plot(t,Sn);
subplot(2,1,2);
plot(t,X);


figure(3);
plot(t,Vn);



%% 2阶常微分 矩阵
clear;
clc;
dt=0.01;
T = 100;
N=T/dt;
Vn=zeros(6,N);
Sn=zeros(6,N);
Vn(1,1)=0;
Sn(3,1)=100;
t=zeros(1,N);
for i=1:N
   t(i)=(i-1)*dt; 
end
for i=1:N-1
   S0(1:6,1)=Sn(:,i);
   V0(1:6,1)=Vn(:,i);
   [Sn(:,i+1),Vn(:,i+1)]= lgkt4k2kM(t(i),dt,S0,V0);
end
figure(4);
subplot(2,3,1);
plot(t,Sn(1,:));
title('Surge');

subplot(2,3,2);
plot(t,Sn(2,:));
title('Sway');

subplot(2,3,3);
plot(t,Sn(3,:));
title('Heave');

subplot(2,3,4);
plot(t,Sn(4,:));
title('Roll');

subplot(2,3,5);
plot(t,Sn(5,:));
title('Pitch');

subplot(2,3,6);
plot(t,Sn(6,:));
title('Yaw');


%% 1阶常微分
function X1 = lgkt4k(t,dt,X0)
    k1=f(t,X0);
    k2=f(t+dt/2,X0+dt/2*k1);
    k3=f(t+dt/2,X0+dt/2*k2);
    k4=f(t+dt,X0+dt*k3);
    X1=X0+dt/6*(k1+2*k2+2*k3+k4);
end

function y=f(t,x)
    m=60;
    g=9.81;
    f0=100;
    f1=m;

    %omega=2*pi*0.1;
    F=-m*g;  %sin(omega*t);
    y=(F-f0*x)/f1;
end

%% 2阶常微分
function [a1,b1]=lgkt4k2k(t,dt,a0,b0)
    k1a=ff(t,a0,b0);
    k1b=gg(t,a0,b0);
    k2a=ff(t+dt/2,a0+dt/2*k1a,b0+dt/2*k1b);
    k2b=gg(t+dt/2,a0+dt/2*k1a,b0+dt/2*k1b);
    k3a=ff(t+dt/2,a0+dt/2*k2a,b0+dt/2*k2b);
    k3b=gg(t+dt/2,a0+dt/2*k2a,b0+dt/2*k2b);
    k4a=ff(t+dt,a0+dt*k3a,b0+dt*k3b);
    k4b=gg(t+dt,a0+dt*k3a,b0+dt*k3b);
    a1=a0+dt/6*(k1a+2*k2a+2*k3a+k4a);
    b1=b0+dt/6*(k1b+2*k2b+2*k3b+k4b);
end

function y1=ff(t,a,b)
y1=b;
end

function y2=gg(t,a,b)
m=60;
g=9.81;
f0=0;
f1=0;
f2=1;
F=cos(t); 
   
y2=(F-(f1*b+f0*a))/f2;
end
%% 2阶常微分方程 矩阵形式
function [A1,B1]=lgkt4k2kM(t,dt,A0,B0)
    K1A = FF(t,A0,B0);
    K1B = GG(t,A0,B0);
    K2A = FF(t+dt/2,A0+dt/2*K1A,B0+dt/2*K1B);
    K2B = GG(t+dt/2,A0+dt/2*K1A,B0+dt/2*K1B);
    K3A = FF(t+dt/2,A0+dt/2*K2A,B0+dt/2*K2B);
    K3B = GG(t+dt/2,A0+dt/2*K2A,B0+dt/2*K2B);
    K4A = FF(t+dt,A0+dt*K3A,B0+dt*K3B);
    K4B = GG(t+dt,A0+dt*K3A,B0+dt*K3B);
    A1=A0+dt/6*(K1A+2*K2A+2*K3A+K4A);
    B1=B0+dt/6*(K1B+2*K2B+2*K3B+K4B);    
end

function Y1 = FF(t,A,B)
    Y1 = B;
end

function Y2 = GG(t,A,B)
mass = 1000;
g = 9.81;
a=1;
b=1;
c=1;
I11 = 1/12*mass*(b^2+c^2);
I22 = 1/12*mass*(a^2+c^2);
I33 = 1/12*mass*(a^2+b^2);
M11=mass;
M22=mass;
M33=mass;
M44=I11;
M55=I22;
M66=I33;

M_ = [M11 0 0 0 0 0;
     0 M22 0 0 0 0;
     0 0 M33 0 0 0;
     0 0 0 M44 0 0;
     0 0 0 0 M55 0;
     0 0 0 0 0 M66;];

B11 = 10;
B22 = 10;
B33 = 5;
B44 = 50;
B55 = 50;
B66 = 50;

B_ = [B11 0 0 0 0 0;
     0 B22 0 0 0 0;
     0 0 B33 0 0 0;
     0 0 0 B44 0 0;
     0 0 0 0 B55 0;
     0 0 0 0 0 B66;];
 
 C11 = 100;
 C22 = 100;
 C33 = 100;
 C44 = 100;
 C55 = 100;
 C66 = 100;
 
 C_ = [C11 0 0 0 0 0;
     0 C22 0 0 0 0;
     0 0 C33 0 0 0;
     0 0 0 C44 0 0;
     0 0 0 0 C55 0;
     0 0 0 0 0 C66;];
 Fn = [0;0;0;0;0;0];%-mass*g
 
 F2 = M_;
 F1 = B_;
 F0 = C_;
 F2inv = inv(F2);
 
 Y2 = F2inv*(Fn-(F1*B+F0*A)); %F2\ 与inv(F2)* 是等效的
 
end



