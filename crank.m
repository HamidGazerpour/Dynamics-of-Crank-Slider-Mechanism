% Crank movement 
clc; close all; clear all;
% required force for cutting:
r_wire = 0.006/2; A_wire = pi*r_wire^2; sigma_wire = 415 *10^6; 
F = sigma_wire * A_wire; % [N]
% Data
omega = 90*2*pi/60;  % rad/s

a = 0.02; % length of crank [m]
b = 0.045;
Mass = 0.4 + 1.8;
Jtot = 0.155;

link=line;
set(link,'XData',[],'Color','b' ,'LineWidth',1.5,'YData',[],'marker','o');
grid;
axis('equal',[-a*1.5 (a+b)*1.1 -a*1.5 a*1.5]);

i=1;
alpha0=rad2deg(pi/2+pi/6);
beta0   = asind(-a*sind((alpha0))/b);
beta(1) = asind(-a*sind((alpha0))/b);
alpha_crank(1)= rad2deg(pi/2+pi/6); beta_dot(1)=0;

for j = 0:1:360

    res1 = MC_03(j);
    t(i) = j;
    time(i) = deg2rad(j)/omega;
    p3(i) = res1.pos;
    v3(i) = res1.vel*omega;
    a3(i) = res1.acc*omega^2;

    x(i) = b*cosd(beta0)+a*cosd(alpha0)  +  p3(i);
    x_dot(i) = v3(i); x_2dot(i) = a3(i);
    
     alpha_crank(i+1) = acosd((x(i)^2+a^2-b^2)/(2*x(i)*a));
    
     beta(i+1) = asind((-a*sind(alpha_crank(i+1))/b));
%     beta(i+1)= beta(i) + beta_dot(i)*0.1*60/90;
    alpha_crank_dot(i+1) = -x_dot(i)*cosd(beta(i+1))/(a*sind(alpha_crank(i+1)-beta(i+1)));
    beta_dot(i+1) = -deg2rad(alpha_crank_dot(i+1))*a*cosd(alpha_crank(i+1))/(b*cosd(beta(i+1)));

BB=-deg2rad(alpha_crank_dot(i+1))^2*a*(cosd(alpha_crank(i+1))*cosd(beta(i+1)) + sind(alpha_crank(i+1))*sind(beta(i+1)))  -deg2rad(beta_dot(i+1))^2*b;
AA=a* (sind(alpha_crank(i+1))*cosd(beta(i+1))-cosd(alpha_crank(i+1))*sind(beta(i+1)));
alpha_crank_2dot(i+1) = (-cosd(beta(i+1))*x_2dot(i) +BB)/AA;


    Ax(i)=a*cosd(alpha_crank(i));
    Ay(i)=a*sind(alpha_crank(i));
    Bx(i)=Ax(i)+b*sqrt(1-(a/b*sind(alpha_crank(i)))^2);


    set(link,'XData',[0 Ax(i) Bx(i)],'YData',[0 Ay(i) 0]);
    pause(0.001)
    if x_dot(i) > 0
        Cr(i) = (F*x_dot(i) + Mass*x_2dot(i)*x_dot(i)) / alpha_crank(i) ;
    else
        Cr(i) = (Mass*x_2dot(i)*x_dot(i)) / alpha_crank(i);
    end
    CJ(i) = - Jtot*alpha_crank_2dot(i+1); % minus because of opposite directions of x and alpha
    Crs(i) = Cr(i) + CJ(i);
    i=i+1;
    
end

figure;
subplot(3,1,1); plot([0 time],alpha_crank,'r','LineWidth',1);grid;
title('Angular position of crank'); xlabel('time [s]'); ylabel('${\alpha}$ [deg]','interpreter','latex');

subplot(3,1,2); plot([0 time],alpha_crank_dot,'b','LineWidth',1);grid;
title('Angular Velocity of crank'); xlabel('time [s]'); ylabel('$\dot{\alpha}$ [deg/s]','interpreter','latex');

subplot(3,1,3); plot([0 time],alpha_crank_2dot,'g','LineWidth',1);grid;
title('Angular Acceleration of crank'); xlabel('time [s]'); ylabel('$\ddot{\alpha}$ [deg/s]','interpreter','latex');

%%
figure;
subplot(3,1,1); plot(time,Cr,'r','LineWidth',1);grid; 
title('Required Torque because of external force'); xlabel('time [s]'); ylabel('C_{r} [N.m]');

subplot(3,1,2); plot(time,CJ,'g','LineWidth',1);grid; 
title('Required Torque because of moment of inertia'); xlabel('time [s]'); ylabel('C_{rJ} [N.m]');

subplot(3,1,3); plot(time,Crs,'b','LineWidth',1);grid; 
title('Total required Torque '); xlabel('time [s]'); ylabel('C_{rs} [N.m]');

figure;
plot(alpha_crank_dot(2:end)/(2*pi),Crs,'k','linewidth',1);
title('Motor torque');xlabel('rpm'); ylabel('Motor torque');


% the variable to use in simscape model:
alpha_crank(2)=0;
time_alpha = [10*time' -(alpha_crank(2:end))']; %the var which is used in simscape model
save 'time_alpha' time_alpha 
