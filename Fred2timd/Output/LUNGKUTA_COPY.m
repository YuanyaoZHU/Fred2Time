h=0.001;

t0=0:h:400;

w=5;

ep=0.02;

Fm=0.1;

Fah=0.05;

u(1)=0;u(2)=0;

for i=1:length(t0) % 进行多次迭代

  tao=t0(i);

  u=RK(u,tao,h,ep,w,Fm,Fah);

  Result(i,:)=u; % 将每次迭代的位移和速度保存

end

figure(1)

subplot(2,1,1)

plot(t0,Result(:,1)) % 绘制位移图

xlabel('Time')

ylabel('displacement')

subplot(2,1,2)

plot(t0,Result(:,2))% 绘制速度图

xlabel('Time')

ylabel('velocity')


function u=RK(u,tao,h,ep,w,Fm,Fah)

KK1=u(2);%u1是位移；u2是速度

KK2=u(2)+h/2*KK1;

KK3=u(2)+h/2*KK2;

KK4=u(2)+h*KK3;

u(1)=u(1)+h/6*(KK1+2*KK2+2*KK3+KK4);

K1=-2*ep*u(2)-u(1)+Fm+Fah*cos(w*tao);

K2=-2*ep*(u(2)+h/2*K1)-u(1)-h/2+Fm+Fah*cos(w*tao);

K3=-2*ep*(u(2)+h/2*K2)-u(1)-h/2+Fm+Fah*cos(w*tao);

K4=-2*ep*(u(2)+h*K3)-u(1)-h+Fm+Fah*cos(w*tao);

u(2)=u(2)+h/6*(K1+2*K2+2*K3+K4);

end