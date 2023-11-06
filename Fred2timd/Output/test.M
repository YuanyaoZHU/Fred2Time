%% use for case a
clc;
clear;
timestep=0.;
RtA=180/pi;%弧度变角度
AtR=pi/180;%角度变弧度
DISP = load('ODISP.txt');
VLCT = load('OVLCT.txt');
MOORINGFORCE = load('ZHUMOORFORCE.txt');
MOORINGFORCE0 = load('ZHUMOORFORCE0.txt');


MOO1 = load('ZHUWMOOR1.txt');
MOO2 = load('ZHUWMOOR2.txt');
MOO3 = load('ZHUWMOOR3.txt');
%Ldamp = load('ZHULDAMPING.txt');
%COOR_Fair1 = load('ZHUFAIR1.txt');
%COOR_Fair2 = load('ZHUFAIR2.txt');
%COOR_Fair3 = load('ZHUFAIR3.txt');
%LOCAL_Fair1 = load('ZHUFAIRLOCAL1.txt');
%LOCAL_Fair2 = load('ZHUFAIRLOCAL2.txt');
%LOCAL_Fair3 = load('ZHUFAIRLOCAL3.txt');
%DAMP=load('ZHUDAMPINGFORCE.txt');
%HSRF=load('ZHUHSFORCE.txt');
%RETA=load('ZHUREATFORCE.txt');
%WAVE=load('ZHUWAVEFORCE.txt');
%WIND=load('ZHUWINDFORCE.txt');

%%
%DISP = load('ODISP.txt');
figure(1);% 显示位移曲线
for i= 1:3
    subplot(2,3,i)
    plot(DISP(:,1),DISP(:,i+1));% MAXDISP(i+1)
    title('6-degree');
end
for i=4:6
   subplot(2,3,i)
   plot(DISP(:,1),DISP(:,i+1)*RtA);
end
%%
figure(2); %显示系泊力合力
MoringTitle = ["x" "y" "z" "rx" "ry" "rz"];
N=length(MOORINGFORCE(:,1));


for i= 1:6
    subplot(2,3,i)
    plot(1:N,MOORINGFORCE(1:end,i+1));% MAXDISP(i+1)
    title(MoringTitle(i));
end
%%
figure(3); %显示系泊力合力
MoringTitle = ["x" "y" "z" "rx" "ry" "rz"];
for i= 1:6
    subplot(2,3,i)
    plot(DISP(:,1),MOO1(1:end,i));% MAXDISP(i+1)
    title(MoringTitle(i));
end
%%
figure(2);
for(i= 1:6)
    subplot(2,3,i)
    plot(DISP(:,1),Ldamp(:,i));% MAXDISP(i+1)
    title('LDAMP')
end

figure(3);
for(i= 1:6)
    subplot(2,3,i)
    plot(VLCT(:,1),VLCT(:,i+1));% MAXDISP(i+1)
    title('VLCT')
end

%%
MAXDISP = max(abs(DISP));
MAXMOOR = max(abs(MOORINGFORCE));
MAXDAMP = max(abs(DAMP));
MAXHSRF = max(abs(HSRF));
MAXRETA = max(abs(RETA));
MAXWAVE = max(abs(WAVE));
MAXWIND = max(abs(WIND));

figure(1);
for(i= 1:6)
    subplot(2,3,i)
    plot(DISP(:,1),DISP(:,i+1));% MAXDISP(i+1)
    title('6-degree')
end
for(i=1:6) 
    hold on;
    plot(DISP(:,1),MOORINGFORCE(:,i)/MAXMOOR(i));
    plot(DISP(:,1),DAMP(:,i)/MAXDAMP(i));
    plot(DISP(:,1),HSRF(:,i)/MAXHSRF(i));
    plot(DISP(:,1),RETA(:,i)/MAXRETA(i));
    plot(DISP(:,1),WAVE(:,i)/MAXWAVE(i));
    plot(DISP(:,1),WIND(:,i)/MAXWIND(i));
    hold off;
    legend('DISP','MOOR','DAMP','HSRF','RETA','WAVE','WIND','location','best');    
    set(gca,'LooseInset', get(gca,'TightInset'))
end


for(i= 1:6)
    figure(7);
    subplot(2,3,i)
    plot(DISP(:,1),MOORINGFORCE(:,i));
    title('MOOR')
end
for(i= 1:6)
    set(gca,'LooseInset', get(gca,'TightInset'))
    figure(3)
    subplot(2,3,i)
    plot(DISP(:,1),DAMP(:,i));
    title('damp')
    set(gca,'LooseInset', get(gca,'TightInset'))
    figure(4)
    subplot(2,3,i)
    plot(DISP(:,1),HSRF(:,i));
    title('HSRF') 
    set(gca,'LooseInset', get(gca,'TightInset'))
    figure(5)
    subplot(2,3,i)
    plot(DISP(:,1),RETA(:,i));
    title('RETA')  
    set(gca,'LooseInset', get(gca,'TightInset'))
    figure(6)
    subplot(2,3,i)
    plot(DISP(:,1),WAVE(:,i));
    title('WAVE')  
    set(gca,'LooseInset', get(gca,'TightInset'))
end





for(i=1:length(DAMP(:,1)))
    if(MOORINGFORCE(i)==0)
        BI(i) = 0;
    else
        BI(i) = DAMP(i)/MOORINGFORCE(i);
    end
end

figure(7);
plot(DISP(:,1),BI);

figure(8);
plot(DISP(:,1),DAMP(:,2));
hold on;
plot(DISP(:,1),MOORINGFORCE(:,2));

hold off;
legend('damp','moor');

figure(9);
hold on;
for(i=1:3)
    plot(DISP(:,1),DISP(:,i+1));
end
hold off;
%%
parfor(i=1:length(MOO1(:,1))-1)
    MOOR1(i) = sqrt(MOO1(i,1)^2+MOO1(i,2)^2)
    MOOR2(i) = sqrt(MOO2(i,1)^2+MOO2(i,2)^2)
    MOOR3(i) = sqrt(MOO3(i,1)^2+MOO3(i,2)^2)
end
figure(10);
subplot(3,1,1)
plot((1:length(MOOR1)),MOOR1/10e6);
subplot(3,1,2)
plot((1:length(MOOR2)),MOOR2/10e6);
subplot(3,1,3)
plot((1:length(MOOR3)),MOOR3/10e6);
%%
figure(11);
subplot(3,2,1)
plot((1:length(MOO1(:,1)))/4*0.5,MOO1(:,1)/10e6);
subplot(3,2,2)
plot((1:length(MOO1(:,1)))/4*0.5,MOO1(:,2)/10e6);
subplot(3,2,3)
plot((1:length(MOO2(:,1)))/4*0.5,MOO2(:,1)/10e6);
subplot(3,2,4)
plot((1:length(MOO2(:,1)))/4*0.5,MOO2(:,2)/10e6);
subplot(3,2,5)
plot((1:length(MOO3(:,1)))/4*0.5,MOO3(:,1)/10e6);
subplot(3,2,6)
plot((1:length(MOO3(:,1)))/4*0.5,MOO3(:,2)/10e6);

%%
figure(12);
subplot(3,2,1)
plot((1:length(LOCAL_Fair1(:,1)))/4*0.5,LOCAL_Fair1(:,1));
subplot(3,2,2)
plot((1:length(LOCAL_Fair1(:,1)))/4*0.5,LOCAL_Fair1(:,2));
subplot(3,2,3)
plot((1:length(LOCAL_Fair2(:,1)))/4*0.5,LOCAL_Fair2(:,1));
subplot(3,2,4)
plot((1:length(LOCAL_Fair2(:,1)))/4*0.5,LOCAL_Fair2(:,2));
subplot(3,2,5)
plot((1:length(LOCAL_Fair3(:,1)))/4*0.5,LOCAL_Fair3(:,1));
subplot(3,2,6)
plot((1:length(LOCAL_Fair3(:,1)))/4*0.5,LOCAL_Fair3(:,2));


