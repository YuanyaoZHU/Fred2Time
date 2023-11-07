clc;
clear;
filename = 'Spar.1';
num_freedom=10;
fileID = fopen(filename);
dot = textscan(fileID,'%f %f %f %f %f','headerlines',num_freedom*2); 
dot3 = load('Spar.3');
dot1 = cell2mat(dot);
hst = load('Spar.hst');

g=9.80665;%重力加速度
rho = 1025;%海水密度
rstMatrix = zeros(6,6);
for i=1:36
  rstMatrix(hst(i,1),hst(i,2)) = hst(i,3)*rho*g;
end

line = ones(6,6);
for i=1:length(dot1(:,1))
    dot1(i,1) = 1/dot1(i,1)*2*pi;
end
NumOmega = length(dot1(:,1))/num_freedom;
HamsAddMass1 = zeros(NumOmega,8);
HamsDamping1 = zeros(NumOmega,8);
HamsAddMass2 = zeros(NumOmega,8);
HamsDamping2 = zeros(NumOmega,8);
HamsAddMass3 = zeros(NumOmega,8);
HamsDamping3 = zeros(NumOmega,8);
HamsAddMass4 = zeros(NumOmega,8);
HamsDamping4 = zeros(NumOmega,8);
HamsAddMass5 = zeros(NumOmega,8);
HamsDamping5 = zeros(NumOmega,8);
HamsAddMass6 = zeros(NumOmega,8);
HamsDamping6 = zeros(NumOmega,8);

HamsAddMass1(:,2) = dot1(1:num_freedom:end,1);


for i = 1:NumOmega
   HamsAddMass1(i,1)=HamsAddMass1(i,2)^2/g; 
end

HamsAddMass2(:,1:2) = HamsAddMass1(:,1:2);
HamsAddMass3(:,1:2) = HamsAddMass1(:,1:2);
HamsAddMass4(:,1:2) = HamsAddMass1(:,1:2);
HamsAddMass5(:,1:2) = HamsAddMass1(:,1:2);
HamsAddMass6(:,1:2) = HamsAddMass1(:,1:2);
HamsDamping1(:,1:2) = HamsAddMass1(:,1:2);
HamsDamping2(:,1:2) = HamsAddMass1(:,1:2);
HamsDamping3(:,1:2) = HamsAddMass1(:,1:2);
HamsDamping4(:,1:2) = HamsAddMass1(:,1:2);
HamsDamping5(:,1:2) = HamsAddMass1(:,1:2);
HamsDamping6(:,1:2) = HamsAddMass1(:,1:2);


for i=1:length(dot1(:,1))
    if dot1(i,2) == 1
        switch dot1(i,3)
            case 1
                HamsAddMass1(line(1,1),3) = dot1(i,4)*rho;
                HamsDamping1(line(1,1),3) = dot1(i,5)*rho*HamsDamping1(line(1,1),2);
                line(1,1)=line(1,1)+1;
            case 2
                HamsAddMass1(line(1,2),4) = dot1(i,4)*rho;
                HamsDamping1(line(1,2),4) = dot1(i,5)*rho*HamsDamping1(line(1,2),2);
                line(1,2)=line(1,2)+1;
            case 3
                HamsAddMass1(line(1,3),5) = dot1(i,4)*rho;
                HamsDamping1(line(1,3),5) = dot1(i,5)*rho*HamsDamping1(line(1,3),2);
                line(1,3)=line(1,3)+1;
            case 4
                HamsAddMass1(line(1,4),6) = dot1(i,4)*rho;
                HamsDamping1(line(1,4),6) = dot1(i,5)*rho*HamsDamping1(line(1,4),2);
                line(1,4)=line(1,4)+1;
            case 5
                HamsAddMass1(line(1,5),7) = dot1(i,4)*rho;
                HamsDamping1(line(1,5),7) = dot1(i,5)*rho*HamsDamping1(line(1,5),2);
                line(1,5)=line(1,5)+1;
            case 6
                HamsAddMass1(line(1,6),8) = dot1(i,4)*rho;
                HamsDamping1(line(1,6),8) = dot1(i,5)*rho*HamsDamping1(line(1,6),2);
                line(1,6)=line(1,6)+1;
        end
    elseif dot1(i,2) == 2
        switch dot1(i,3)
            case 1
                HamsAddMass2(line(2,1),3) = dot1(i,4)*rho;
                HamsDamping2(line(2,1),3) = dot1(i,5)*rho*HamsDamping2(line(2,1),2);
                line(2,1)=line(2,1)+1;
            case 2
                HamsAddMass2(line(2,2),4) = dot1(i,4)*rho;
                HamsDamping2(line(2,2),4) = dot1(i,5)*rho*HamsDamping2(line(2,2),2);
                line(2,2)=line(2,2)+1;
            case 3
                HamsAddMass2(line(2,3),5) = dot1(i,4)*rho;
                HamsDamping2(line(2,3),5) = dot1(i,5)*rho*HamsDamping2(line(2,3),2);
                line(2,3)=line(2,3)+1;
            case 4
                HamsAddMass2(line(2,4),6) = dot1(i,4)*rho;
                HamsDamping2(line(2,4),6) = dot1(i,5)*rho*HamsDamping2(line(2,4),2);
                line(2,4)=line(2,4)+1;
            case 5
                HamsAddMass2(line(2,5),7) = dot1(i,4)*rho;
                HamsDamping2(line(2,5),7) = dot1(i,5)*rho*HamsDamping2(line(2,5),2);
                line(2,5)=line(2,5)+1;
            case 6
                HamsAddMass2(line(2,6),8) = dot1(i,4)*rho;
                HamsDamping2(line(2,6),8) = dot1(i,5)*rho*HamsDamping2(line(2,6),2);
                line(2,6)=line(2,6)+1;
        end
    elseif dot1(i,2) == 3
        switch dot1(i,3)
            case 1
                HamsAddMass3(line(3,1),3) = dot1(i,4)*rho;
                HamsDamping3(line(3,1),3) = dot1(i,5)*rho*HamsDamping3(line(3,1),2);
                line(3,1)=line(3,1)+1;
            case 2
                HamsAddMass3(line(3,2),4) = dot1(i,4)*rho;
                HamsDamping3(line(3,2),4) = dot1(i,5)*rho*HamsDamping3(line(3,2),2);
                line(3,2)=line(3,2)+1;
            case 3
                HamsAddMass3(line(3,3),5) = dot1(i,4)*rho;
                HamsDamping3(line(3,3),5) = dot1(i,5)*rho*HamsDamping3(line(3,3),2);
                line(3,3)=line(3,3)+1;
            case 4
                HamsAddMass3(line(3,4),6) = dot1(i,4)*rho;
                HamsDamping3(line(3,4),6) = dot1(i,5)*rho*HamsDamping3(line(3,4),2);
                line(3,4)=line(3,4)+1;
            case 5
                HamsAddMass3(line(3,5),7) = dot1(i,4)*rho;
                HamsDamping3(line(3,5),7) = dot1(i,5)*rho*HamsDamping3(line(3,5),2);
                line(3,5)=line(3,5)+1;
            case 6
                HamsAddMass3(line(3,6),8) = dot1(i,4)*rho;
                HamsDamping3(line(3,6),8) = dot1(i,5)*rho*HamsDamping3(line(3,6),2);
                line(3,6)=line(3,6)+1;
        end
    elseif dot1(i,2) == 4
        switch dot1(i,3)
            case 1
                HamsAddMass4(line(4,1),3) = dot1(i,4)*rho;
                HamsDamping4(line(4,1),3) = dot1(i,5)*rho*HamsDamping4(line(4,1),2);
                line(4,1)=line(4,1)+1;
            case 2
                HamsAddMass4(line(4,2),4) = dot1(i,4)*rho;
                HamsDamping4(line(4,2),4) = dot1(i,5)*rho*HamsDamping4(line(4,2),2);
                line(4,2)=line(4,2)+1;
            case 3
                HamsAddMass4(line(4,3),5) = dot1(i,4)*rho;
                HamsDamping4(line(4,3),5) = dot1(i,5)*rho*HamsDamping4(line(4,3),2);
                line(4,3)=line(4,3)+1;
            case 4
                HamsAddMass4(line(4,4),6) = dot1(i,4)*rho;
                HamsDamping4(line(4,4),6) = dot1(i,5)*rho*HamsDamping4(line(4,4),2);
                line(4,4)=line(4,4)+1;
            case 5
                HamsAddMass4(line(4,5),7) = dot1(i,4)*rho;
                HamsDamping4(line(4,5),7) = dot1(i,5)*rho*HamsDamping4(line(4,5),2);
                line(4,5)=line(4,5)+1;
            case 6
                HamsAddMass4(line(4,6),8) = dot1(i,4)*rho;
                HamsDamping4(line(4,6),8) = dot1(i,5)*rho*HamsDamping4(line(4,6),2);
                line(4,6)=line(4,6)+1;
        end
    elseif dot1(i,2) == 5
        switch dot1(i,3)
            case 1
                HamsAddMass5(line(5,1),3) = dot1(i,4)*rho;
                HamsDamping5(line(5,1),3) = dot1(i,5)*rho*HamsDamping5(line(5,1),2);
                line(5,1)=line(5,1)+1;
            case 2
                HamsAddMass5(line(5,2),4) = dot1(i,4)*rho;
                HamsDamping5(line(5,2),4) = dot1(i,5)*rho*HamsDamping5(line(5,2),2);
                line(5,2)=line(5,2)+1;
            case 3
                HamsAddMass5(line(5,3),5) = dot1(i,4)*rho;
                HamsDamping5(line(5,3),5) = dot1(i,5)*rho*HamsDamping5(line(5,3),2);
                line(5,3)=line(5,3)+1;
            case 4
                HamsAddMass5(line(5,4),6) = dot1(i,4)*rho;
                HamsDamping5(line(5,4),6) = dot1(i,5)*rho*HamsDamping5(line(5,4),2);
                line(5,4)=line(5,4)+1;
            case 5
                HamsAddMass5(line(5,5),7) = dot1(i,4)*rho;
                HamsDamping5(line(5,5),7) = dot1(i,5)*rho*HamsDamping5(line(5,5),2);
                line(5,5)=line(5,5)+1;
            case 6
                HamsAddMass5(line(5,6),8) = dot1(i,4)*rho;
                HamsDamping5(line(5,6),8) = dot1(i,5)*rho*HamsDamping5(line(5,6),2);
                line(5,6)=line(5,6)+1;
        end
    elseif dot1(i,2) == 6
        switch dot1(i,3)
            case 1
                HamsAddMass6(line(6,1),3) = dot1(i,4)*rho;
                HamsDamping6(line(6,1),3) = dot1(i,5)*rho*HamsDamping6(line(6,1),2);
                line(6,1)=line(6,1)+1;
            case 2
                HamsAddMass6(line(6,2),4) = dot1(i,4)*rho;
                HamsDamping6(line(6,2),4) = dot1(i,5)*rho*HamsDamping6(line(6,2),2);
                line(6,2)=line(6,2)+1;
            case 3
                HamsAddMass6(line(6,3),5) = dot1(i,4)*rho;
                HamsDamping6(line(6,3),5) = dot1(i,5)*rho*HamsDamping6(line(6,3),2);
                line(6,3)=line(6,3)+1;
            case 4
                HamsAddMass6(line(6,4),6) = dot1(i,4)*rho;
                HamsDamping6(line(6,4),6) = dot1(i,5)*rho*HamsDamping6(line(6,4),2);
                line(6,4)=line(6,4)+1;
            case 5
                HamsAddMass6(line(6,5),7) = dot1(i,4)*rho;
                HamsDamping6(line(6,5),7) = dot1(i,5)*rho*HamsDamping6(line(6,5),2);
                line(6,5)=line(6,5)+1;
            case 6
                HamsAddMass6(line(6,6),8) = dot1(i,4)*rho;
                HamsDamping6(line(6,6),8) = dot1(i,5)*rho*HamsDamping6(line(6,6),2);
                line(6,6)=line(6,6)+1;
        end        
    end
    
    
end
%% 生成波浪激励力
line1 = ones(6);
HamsExcit1 = zeros(NumOmega,4);
HamsExcit2 = zeros(NumOmega,4);
HamsExcit3 = zeros(NumOmega,4);
HamsExcit4 = zeros(NumOmega,4);
HamsExcit5 = zeros(NumOmega,4);
HamsExcit6 = zeros(NumOmega,4);

HamsExcit1(:,1:2) = HamsAddMass1(:,1:2);
HamsExcit2(:,1:2) = HamsAddMass1(:,1:2);
HamsExcit3(:,1:2) = HamsAddMass1(:,1:2);
HamsExcit4(:,1:2) = HamsAddMass1(:,1:2);
HamsExcit5(:,1:2) = HamsAddMass1(:,1:2);
HamsExcit6(:,1:2) = HamsAddMass1(:,1:2);
for i=1:length(dot3(:,1))
   if dot3(i,2) == 0
       switch dot3(i,3)
           case 1
               HamsExcit1(line1(1),3) = dot3(i,6)*rho*g;
               HamsExcit1(line1(1),4) = dot3(i,7)*rho*g;
               line1(1)=line1(1)+1;
           case 2 
               HamsExcit2(line1(2),3) = dot3(i,6)*rho*g;
               HamsExcit2(line1(2),4) = dot3(i,7)*rho*g;
               line1(2)=line1(2)+1;
           case 3 
               HamsExcit3(line1(3),3) = dot3(i,6)*rho*g;
               HamsExcit3(line1(3),4) = dot3(i,7)*rho*g;
               line1(3)=line1(3)+1;
           case 4 
               HamsExcit4(line1(4),3) = dot3(i,6)*rho*g;
               HamsExcit4(line1(4),4) = dot3(i,7)*rho*g;
               line1(4)=line1(4)+1;
           case 5 
               HamsExcit5(line1(5),3) = dot3(i,6)*rho*g;
               HamsExcit5(line1(5),4) = dot3(i,7)*rho*g;
               line1(5)=line1(5)+1;
           case 6 
               HamsExcit6(line1(6),3) = dot3(i,6)*rho*g;
               HamsExcit6(line1(6),4) = dot3(i,7)*rho*g;
               line1(6)=line1(6)+1;
       end
   end
end
%%
HamsAddMass1_E = zeros(5*NumOmega,8);
HamsAddMass2_E = zeros(5*NumOmega,8);
HamsAddMass3_E = zeros(5*NumOmega,8);
HamsAddMass4_E = zeros(5*NumOmega,8);
HamsAddMass5_E = zeros(5*NumOmega,8);
HamsAddMass6_E = zeros(5*NumOmega,8);

HamsDamping1_E = zeros(5*NumOmega,8);
HamsDamping2_E = zeros(5*NumOmega,8);
HamsDamping3_E = zeros(5*NumOmega,8);
HamsDamping4_E = zeros(5*NumOmega,8);
HamsDamping5_E = zeros(5*NumOmega,8);
HamsDamping6_E = zeros(5*NumOmega,8);

HamsExcit1_E = zeros(5*NumOmega,4);
HamsExcit2_E = zeros(5*NumOmega,8);
HamsExcit3_E = zeros(5*NumOmega,8);
HamsExcit4_E = zeros(5*NumOmega,8);
HamsExcit5_E = zeros(5*NumOmega,8);
HamsExcit6_E = zeros(5*NumOmega,8);

for i=1:NumOmega
   for j = 1:5
       for k = 1:8
         if i==1
             if k==1 || k==2
               HamsAddMass1_E(5*(i-1)+j,k) = HamsAddMass1(i,k)*j*0.2;
               HamsAddMass2_E(5*(i-1)+j,k) = HamsAddMass2(i,k)*j*0.2;
               HamsAddMass3_E(5*(i-1)+j,k) = HamsAddMass3(i,k)*j*0.2;
               HamsAddMass4_E(5*(i-1)+j,k) = HamsAddMass4(i,k)*j*0.2;
               HamsAddMass5_E(5*(i-1)+j,k) = HamsAddMass5(i,k)*j*0.2;
               HamsAddMass6_E(5*(i-1)+j,k) = HamsAddMass6(i,k)*j*0.2;
               
               HamsDamping1_E(5*(i-1)+j,k) = HamsDamping1(i,k)*j*0.2;
               HamsDamping2_E(5*(i-1)+j,k) = HamsDamping2(i,k)*j*0.2;
               HamsDamping3_E(5*(i-1)+j,k) = HamsDamping3(i,k)*j*0.2;
               HamsDamping4_E(5*(i-1)+j,k) = HamsDamping4(i,k)*j*0.2;
               HamsDamping5_E(5*(i-1)+j,k) = HamsDamping5(i,k)*j*0.2;
               HamsDamping6_E(5*(i-1)+j,k) = HamsDamping6(i,k)*j*0.2;
             else
                 HamsAddMass1_E(5*(i-1)+j,k) = HamsAddMass1(i,k);
                 HamsAddMass2_E(5*(i-1)+j,k) = HamsAddMass2(i,k);
                 HamsAddMass3_E(5*(i-1)+j,k) = HamsAddMass3(i,k);
                 HamsAddMass4_E(5*(i-1)+j,k) = HamsAddMass4(i,k);
                 HamsAddMass5_E(5*(i-1)+j,k) = HamsAddMass5(i,k);
                 HamsAddMass6_E(5*(i-1)+j,k) = HamsAddMass6(i,k);
                 
                 HamsDamping1_E(5*(i-1)+j,k) = HamsDamping1(i,k);
                 HamsDamping2_E(5*(i-1)+j,k) = HamsDamping2(i,k);
                 HamsDamping3_E(5*(i-1)+j,k) = HamsDamping3(i,k);
                 HamsDamping4_E(5*(i-1)+j,k) = HamsDamping4(i,k);
                 HamsDamping5_E(5*(i-1)+j,k) = HamsDamping5(i,k);
                 HamsDamping6_E(5*(i-1)+j,k) = HamsDamping6(i,k);
             end
         else
           HamsAddMass1_E(5*(i-1)+j,k) = (HamsAddMass1(i,k)-HamsAddMass1(i-1,k))*j*0.2+HamsAddMass1(i-1,k);
           HamsAddMass2_E(5*(i-1)+j,k) = (HamsAddMass2(i,k)-HamsAddMass2(i-1,k))*j*0.2+HamsAddMass2(i-1,k);
           HamsAddMass3_E(5*(i-1)+j,k) = (HamsAddMass3(i,k)-HamsAddMass3(i-1,k))*j*0.2+HamsAddMass3(i-1,k);
           HamsAddMass4_E(5*(i-1)+j,k) = (HamsAddMass4(i,k)-HamsAddMass4(i-1,k))*j*0.2+HamsAddMass4(i-1,k);
           HamsAddMass5_E(5*(i-1)+j,k) = (HamsAddMass5(i,k)-HamsAddMass5(i-1,k))*j*0.2+HamsAddMass5(i-1,k);
           HamsAddMass6_E(5*(i-1)+j,k) = (HamsAddMass6(i,k)-HamsAddMass6(i-1,k))*j*0.2+HamsAddMass6(i-1,k);
           
           HamsDamping1_E(5*(i-1)+j,k) = (HamsDamping1(i,k)-HamsDamping1(i-1,k))*j*0.2+HamsDamping1(i-1,k);
           HamsDamping2_E(5*(i-1)+j,k) = (HamsDamping2(i,k)-HamsDamping2(i-1,k))*j*0.2+HamsDamping2(i-1,k);
           HamsDamping3_E(5*(i-1)+j,k) = (HamsDamping3(i,k)-HamsDamping3(i-1,k))*j*0.2+HamsDamping3(i-1,k);
           HamsDamping4_E(5*(i-1)+j,k) = (HamsDamping4(i,k)-HamsDamping4(i-1,k))*j*0.2+HamsDamping4(i-1,k);
           HamsDamping5_E(5*(i-1)+j,k) = (HamsDamping5(i,k)-HamsDamping5(i-1,k))*j*0.2+HamsDamping5(i-1,k);
           HamsDamping6_E(5*(i-1)+j,k) = (HamsDamping6(i,k)-HamsDamping6(i-1,k))*j*0.2+HamsDamping6(i-1,k);
         end
       end
       for k = 1:4
         if i==1
             if k==1 || k==2
               HamsExcit1_E(5*(i-1)+j,k) = HamsExcit1(i,k)*j*0.2;
               HamsExcit2_E(5*(i-1)+j,k) = HamsExcit2(i,k)*j*0.2;
               HamsExcit3_E(5*(i-1)+j,k) = HamsExcit3(i,k)*j*0.2;
               HamsExcit4_E(5*(i-1)+j,k) = HamsExcit4(i,k)*j*0.2;
               HamsExcit5_E(5*(i-1)+j,k) = HamsExcit5(i,k)*j*0.2;
               HamsExcit6_E(5*(i-1)+j,k) = HamsExcit6(i,k)*j*0.2;
               
             else
                 HamsExcit1_E(5*(i-1)+j,k) = HamsExcit1(i,k); 
                 HamsExcit2_E(5*(i-1)+j,k) = HamsExcit2(i,k); 
                 HamsExcit3_E(5*(i-1)+j,k) = HamsExcit3(i,k); 
                 HamsExcit4_E(5*(i-1)+j,k) = HamsExcit4(i,k); 
                 HamsExcit5_E(5*(i-1)+j,k) = HamsExcit5(i,k); 
                 HamsExcit6_E(5*(i-1)+j,k) = HamsExcit6(i,k); 
             end
         else
           HamsExcit1_E(5*(i-1)+j,k) = (HamsExcit1(i,k)-HamsExcit1(i-1,k))*j*0.2+HamsExcit1(i-1,k);
           HamsExcit2_E(5*(i-1)+j,k) = (HamsExcit2(i,k)-HamsExcit2(i-1,k))*j*0.2+HamsExcit2(i-1,k);
           HamsExcit3_E(5*(i-1)+j,k) = (HamsExcit3(i,k)-HamsExcit3(i-1,k))*j*0.2+HamsExcit3(i-1,k);
           HamsExcit4_E(5*(i-1)+j,k) = (HamsExcit4(i,k)-HamsExcit4(i-1,k))*j*0.2+HamsExcit4(i-1,k);
           HamsExcit5_E(5*(i-1)+j,k) = (HamsExcit5(i,k)-HamsExcit5(i-1,k))*j*0.2+HamsExcit5(i-1,k);
           HamsExcit6_E(5*(i-1)+j,k) = (HamsExcit6(i,k)-HamsExcit6(i-1,k))*j*0.2+HamsExcit6(i-1,k);
           
         end
       end
       
   end
end

HamsAddMass1_E(end,:)=[];
HamsAddMass2_E(end,:)=[];
HamsAddMass3_E(end,:)=[];
HamsAddMass4_E(end,:)=[];
HamsAddMass5_E(end,:)=[];
HamsAddMass6_E(end,:)=[];

HamsDamping1_E(end,:)=[];
HamsDamping2_E(end,:)=[];
HamsDamping3_E(end,:)=[];
HamsDamping4_E(end,:)=[];
HamsDamping5_E(end,:)=[];
HamsDamping6_E(end,:)=[];

HamsExcit1_E(end,:)=[];
HamsExcit2_E(end,:)=[];
HamsExcit3_E(end,:)=[];
HamsExcit4_E(end,:)=[];
HamsExcit5_E(end,:)=[];
HamsExcit6_E(end,:)=[];


%%
fid = fopen('OAMASS1.txt','wt');
fprintf(fid,'%7.3f %7.3f %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e\n',HamsAddMass1_E.');
fclose(fid);
fid = fopen('OAMASS2.txt','wt');
fprintf(fid,'%7.3f %7.3f %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e\n',HamsAddMass2_E.');
fclose(fid);
fid = fopen('OAMASS3.txt','wt');
fprintf(fid,'%7.3f %7.3f %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e\n',HamsAddMass3_E.');
fclose(fid);
fid = fopen('OAMASS4.txt','wt');
fprintf(fid,'%7.3f %7.3f %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e\n',HamsAddMass4_E.');
fclose(fid);
fid = fopen('OAMASS5.txt','wt');
fprintf(fid,'%7.3f %7.3f %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e\n',HamsAddMass5_E.');
fclose(fid);
fid = fopen('OAMASS6.txt','wt');
fprintf(fid,'%7.3f %7.3f %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e\n',HamsAddMass6_E.');
fclose(fid);

fid = fopen('ODAMPING1.txt','wt');
fprintf(fid,'%7.3f %7.3f %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e\n',HamsDamping1_E.');
fclose(fid);
fid = fopen('ODAMPING2.txt','wt');
fprintf(fid,'%7.3f %7.3f %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e\n',HamsDamping2_E.');
fclose(fid);
fid = fopen('ODAMPING3.txt','wt');
fprintf(fid,'%7.3f %7.3f %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e\n',HamsDamping3_E.');
fclose(fid);
fid = fopen('ODAMPING4.txt','wt');
fprintf(fid,'%7.3f %7.3f %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e\n',HamsDamping4_E.');
fclose(fid);
fid = fopen('ODAMPING5.txt','wt');
fprintf(fid,'%7.3f %7.3f %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e\n',HamsDamping5_E.');
fclose(fid);
fid = fopen('ODAMPING6.txt','wt');
fprintf(fid,'%7.3f %7.3f %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e\n',HamsDamping6_E.');
fclose(fid);

fid = fopen('OEXFOR1.txt','wt');
fprintf(fid,'%7.3f %7.3f %14.5e %14.5e\n',HamsExcit1_E.');
fclose(fid);
fid = fopen('OEXFOR2.txt','wt');
fprintf(fid,'%7.3f %7.3f %14.5e %14.5e\n',HamsExcit2_E.');
fclose(fid);
fid = fopen('OEXFOR3.txt','wt');
fprintf(fid,'%7.3f %7.3f %14.5e %14.5e\n',HamsExcit3_E.');
fclose(fid);
fid = fopen('OEXFOR4.txt','wt');
fprintf(fid,'%7.3f %7.3f %14.5e %14.5e\n',HamsExcit4_E.');
fclose(fid);
fid = fopen('OEXFOR5.txt','wt');
fprintf(fid,'%7.3f %7.3f %14.5e %14.5e\n',HamsExcit5_E.');
fclose(fid);
fid = fopen('OEXFOR6.txt','wt');
fprintf(fid,'%7.3f %7.3f %14.5e %14.5e\n',HamsExcit6_E.');
fclose(fid);

%save('OAMASS1.txt','HamsAddMass1','-ascii');