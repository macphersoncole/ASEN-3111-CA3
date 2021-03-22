%% ASEN 3111 - Computational Assignment 3 - Main
% Script to model a NACA airfoil using the votex panel method to find the 
% sectional coeffecient of lift and the coeffecient of pressure across the 
% airfoil
%
% Author: Cole MacPherson
% Collaborators: R. Block, Z. Lesan, S. Mansfield, A. Uprety
% Date: 26st Feb 2021

%% Housekeeping

clc;
clear;
close all;
fprintf('Code takes about 12 seconds to run on my laptop\n\n');
tic

%% Problem 2 Error
% define the NACA airfoil desired to model
NACA_num = '0012';
% calculate important airfoil values that can be found form the NACA 
% airfoil number
m = str2double(NACA_num(1))/100; % maximum camber
p = str2double(NACA_num(2))/10; % location of maximum camber
t = str2double(NACA_num(3:4))/100; % thickness

% define other characteristics of the airfoil and the flow around it
c = 1; % chord length [m]
V_inf = 50; % free-stream velocity [m/s]
alpha = 0; % angle of attack [deg]

%{

 Test case to ensure the Vortex_Panel function is working correctly by
 comparing the results of the Vortex_Panel function to the reults found in
 Kuethe and Chow textbook section 5.10

 alpha = 8;
 x = [1, 0.933, 0.75, 0.5, 0.25, 0.067, 0, 0.067, 0.25, 0.5, 0.75, 0.933, 1]'; % temp
 y = [0, -0.005, -0.017, -0.033, -0.042, -0.033, 0, 0.045, 0.076, 0.072, 0.044, 0.013, 0]'; % temp
 [c_l, C_p] = Vortex_Panel(x,y,V_inf,alpha);

%}

%{ 

 Find the number of panels required to get the coefficient of pressure
 within 1% difference of the previous number of panels. This will be
 done by running a while loop to find the sectional coefficient of lift and
 the coefficient of pressure at a specified number of panels and then
 finding sectional coefficient of lift and the coefficient of pressure
 agian for a greater number of panels. Once this is done, the percent
 difference between the two will be found and compared to the desired
 error. If the error found is greater than the error desired, the while loop
 will run again with an increasing number of panels. If the error found is
 less than the desired error, that number of panels will be taken as
 accurate and be used for the rest of the analysis.

%}

% initialization
err = 1; % initialize error for while loop
N_0 = 10; % intialize number of employed panels to model the airfoil
N = N_0;
i = 1; % initialize the vector index

% find the x and y values of the surface of the airfoil for the first
% number of panels
[x_prev,y_prev] = NACA_Airfoils(m,p,t,c,N_0); % (x,y) coordinates of the airfoil surface

% find sectional coefficient of lift and coefficient of pressure for
% the first number of panels
[c_l_prev, C_p_prev] = Vortex_Panel(x_prev,y_prev,V_inf,alpha); % c_l and C_p of the airfoil

fprintf('Preforming panel analysis for vortex panel method...\n\n');

while err > 0.01 % 1015 panels of 0.5% ----- 515 panels of 1% 
    % NOTE: occillation occurs for error less than 1%
    
    % increase the number of panels
    N = N + N_0/2;
    N_vec(i) = N;
    
    % find the x and y values of the surface of the airfoil for the second
    % number of panels
    [x_next,y_next] = NACA_Airfoils(m,p,t,c,N); % (x,y) coordinates of the airfoil surface

    % find sectional coefficient of lift and coefficient of pressure for
    % the second number of panels
    [c_l_next, C_p_next] = Vortex_Panel(x_next,y_next,V_inf,alpha); % c_l and C_p of the airfoil

    % calculate error
    
    err = 1 - abs(c_l_prev-c_l_next/c_l_next);
    err_vec(i) = err;
    
    % define the previous sectional coefficient of lift and coefficient of 
    % pressure as the next sectional coefficient of lift and coefficient of
    % pressure if the error was not small enough. Else define the sectional
    % coefficient of lift, coefficient of pressure, x-locations, and 
    % y-locations from the final values calculated
    if err > 0.01
        C_p_prev = C_p_next;
        c_l_prev = c_l_next;
    else
        C_p = C_p_next;
        c_l = c_l_next;
        x = x_next;
        y = y_next;
    end
    
    i = i + 1; % increment the vector index
    
end

% print the number of panels needed to the comand window
fprintf('For a %0.2f%% relative error %i panels are needed to model the airfoil\n\n',err*100,(N-1)*2);

% plot airfoil for ideal number of panels calculated
figure
plot(x,y,'linewidth',2);
axis equal
axis([-c/4 c+c/4 min(y)+min(y)/2 max(y)+max(y)/2]);
xlabel('x');
ylabel('y');
title(['NACA ' NACA_num ' with ' num2str((N-1)*2) ' panels']);
grid on

% plot error vs number of panels
figure
plot(2*(N_vec-1),err_vec*100,'linewidth',2);
xlabel('Number of Panels');
ylabel('Relative Error [%]');
title(['Relative Error vs Number of Panels (NACA ' NACA_num ') [#2]']);
grid on

%% Problem 2 Angle of Attack Variation

fprintf('Preforming angle of attack variation analysis on the coefficinet of pressure and coefficient of lift...\n\n');

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of -5 degrees
[c_l_neg5, C_p_neg5] = Vortex_Panel(x,y,V_inf,-5); % c_l and C_p of the airfoil

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of 5 degrees
[c_l_5, C_p_5] = Vortex_Panel(x,y,V_inf,5); % c_l and C_p of the airfoil

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of 10 degrees
[c_l_10, C_p_10] = Vortex_Panel(x,y,V_inf,10); % c_l and C_p of the airfoil

% plot coefficient of pressure vs angle of attack
figure
subplot(2,2,1)
plot(x(1:end-1)./c,C_p_neg5,'linewidth',2); % -5 deg AoA
set(gca, 'YDir','reverse')
xlabel('x/c');
ylabel('C_{P}');
title('AoA = -5');
xlim([0 1]);
hold on
grid on
subplot(2,2,2)
plot(x(1:end-1)./c,C_p,'linewidth',2); % 0 deg AoA
set(gca, 'YDir','reverse')
xlabel('x/c');
ylabel('C_{P}');
title('AoA = 0');
xlim([0 1]);
hold on
grid on
subplot(2,2,3)
plot(x(1:end-1)./c,C_p_5,'linewidth',2); % 5 deg AoA
set(gca, 'YDir','reverse')
xlabel('x/c');
ylabel('C_{P}');
title('AoA = 5');
xlim([0 1]);
hold on
grid on
subplot(2,2,4)
plot(x(1:end-1)./c,C_p_10,'linewidth',2); % 10 deg AoA
set(gca, 'YDir','reverse')
xlabel('x/c');
ylabel('C_{P}');
title('AoA = 10');
xlim([0 1]);
hold on
grid on
sgtitle(['C_{p} vs Percent Chord at varying AoA (NACA ' NACA_num ')  [#2]']); % title the whole subplot
hold off

figure
plot(x(1:end-1)./c,C_p_neg5,'linewidth',2,'DisplayName','AoA = -5^{o}'); % -5 deg AoA
set(gca, 'YDir','reverse')
xlabel('x/c');
ylabel('C_{P}');
title(['C_{p} vs Percent Chord at varying AoA (NACA ' NACA_num ')  [#2]']);
xlim([0 1]);
hold on
plot(x(1:end-1)./c,C_p,'linewidth',2,'DisplayName','AoA = 0^{o}');
hold on
plot(x(1:end-1)./c,C_p_5,'linewidth',2,'DisplayName','AoA = 5^{o}');
hold on
plot(x(1:end-1)./c,C_p_10,'linewidth',2,'DisplayName','AoA = 10^{o}');
legend
grid on

% plot sectional coefficient of lift vs angle of attack
figure
plot([-5, 0, 5, 10],[c_l_neg5, c_l, c_l_5, c_l_10],'linewidth',2);
xlabel('AoA [deg]');
ylabel('C_{l}');
title(['Sectional Coefficient of lift vs Angle of Attack (NACA ' NACA_num ') [#2]']);
grid on

%% Problem 3

fprintf('Preforming lift slope determination calculations...\n\n');

% Caclulate lift slope
m_L = (c_l_10-c_l_neg5)/(deg2rad(10)-deg2rad(-5)); % lift slope

% Calculate zero-lift angle of attack
AoA_L0 = -(c_l/m_L);

% print lift slope and zero-lift angle of attack values
fprintf('For the NACA %s airfoil\n\tLift slope: %0.3f [1/rad]\n\tZero-lift AoA: %0.3f [deg]\n\n',NACA_num,m_L,rad2deg(AoA_L0));

% plot NACA 0012 sectional coefficient of lift vs angle of attack
figure
plot([-5, 0, 5, 10],[c_l_neg5, c_l, c_l_5, c_l_10],'linewidth',2,'DisplayName',['NACA ' NACA_num]);
hold on

% redefine NACA airfoil
NACA_num = '2412';
% calculate important airfoil values that can be found form the NACA 
% airfoil number
m = str2double(NACA_num(1))/100; % maximum camber
p = str2double(NACA_num(2))/10; % location of maximum camber
t = str2double(NACA_num(3:4))/100; % thickness

% find the x and y values of the surface of the airfoil
[x,y] = NACA_Airfoils(m,p,t,c,N); % (x,y) coordinates of the airfoil surface

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of -5 degrees
[c_l_neg5, ~] = Vortex_Panel(x,y,V_inf,-5); % c_l and C_p of the airfoil

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of -5 degrees
[c_l, ~] = Vortex_Panel(x,y,V_inf,0); % c_l and C_p of the airfoil

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of 5 degrees
[c_l_5, ~] = Vortex_Panel(x,y,V_inf,5); % c_l and C_p of the airfoil

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of 10 degrees
[c_l_10, ~] = Vortex_Panel(x,y,V_inf,10); % c_l and C_p of the airfoil

% Caclulate lift slope
m_L = (c_l_10-c_l_neg5)/(deg2rad(10)-deg2rad(-5)); % lift slope

% Calculate zero-lift angle of attack
AoA_L0 = -(c_l/m_L);

% print lift slope and zero-lift angle of attack values
fprintf('For the NACA %s airfoil\n\tLift slope: %0.3f [1/rad]\n\tZero-lift AoA: %0.3f [deg]\n\n',NACA_num,m_L,rad2deg(AoA_L0));

% plot NACA 2412 sectional coefficient of lift vs angle of attack
plot([-5, 0, 5, 10],[c_l_neg5, c_l, c_l_5, c_l_10],'linewidth',2,'DisplayName',['NACA ' NACA_num]);
hold on

% redefine NACA airfoil
NACA_num = '4412';
% calculate important airfoil values that can be found form the NACA 
% airfoil number
m = str2double(NACA_num(1))/100; % maximum camber
p = str2double(NACA_num(2))/10; % location of maximum camber
t = str2double(NACA_num(3:4))/100; % thickness

% find the x and y values of the surface of the airfoil
[x,y] = NACA_Airfoils(m,p,t,c,N); % (x,y) coordinates of the airfoil surface

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of -5 degrees
[c_l_neg5, ~] = Vortex_Panel(x,y,V_inf,-5); % c_l and C_p of the airfoil

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of -5 degrees
[c_l, ~] = Vortex_Panel(x,y,V_inf,0); % c_l and C_p of the airfoil

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of 5 degrees
[c_l_5, ~] = Vortex_Panel(x,y,V_inf,5); % c_l and C_p of the airfoil

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of 10 degrees
[c_l_10, ~] = Vortex_Panel(x,y,V_inf,10); % c_l and C_p of the airfoil

% Caclulate lift slope
m_L = (c_l_10-c_l_neg5)/(deg2rad(10)-deg2rad(-5)); % lift slope

% Calculate zero-lift angle of attack
AoA_L0 = -(c_l/m_L);

% print lift slope and zero-lift angle of attack values
fprintf('For the NACA %s airfoil\n\tLift slope: %0.3f [1/rad]\n\tZero-lift AoA: %0.3f [deg]\n\n',NACA_num,m_L,rad2deg(AoA_L0));

% plot NACA 4412 sectional coefficient of lift vs angle of attack
plot([-5, 0, 5, 10],[c_l_neg5, c_l, c_l_5, c_l_10],'linewidth',2,'DisplayName',['NACA ' NACA_num]);
hold on

% redefine NACA airfoil
NACA_num = '2424';
% calculate important airfoil values that can be found form the NACA 
% airfoil number
m = str2double(NACA_num(1))/100; % maximum camber
p = str2double(NACA_num(2))/10; % location of maximum camber
t = str2double(NACA_num(3:4))/100; % thickness

% find the x and y values of the surface of the airfoil
[x,y] = NACA_Airfoils(m,p,t,c,N); % (x,y) coordinates of the airfoil surface

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of -5 degrees
[c_l_neg5, ~] = Vortex_Panel(x,y,V_inf,-5); % c_l and C_p of the airfoil

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of -5 degrees
[c_l, ~] = Vortex_Panel(x,y,V_inf,0); % c_l and C_p of the airfoil

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of 5 degrees
[c_l_5, ~] = Vortex_Panel(x,y,V_inf,5); % c_l and C_p of the airfoil

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of 10 degrees
[c_l_10, ~] = Vortex_Panel(x,y,V_inf,10); % c_l and C_p of the airfoil

% Caclulate lift slope
m_L = (c_l_10-c_l_neg5)/(deg2rad(10)-deg2rad(-5)); % lift slope

% Calculate zero-lift angle of attack
AoA_L0 = -(c_l/m_L);

% print lift slope and zero-lift angle of attack values
fprintf('For the NACA %s airfoil\n\tLift slope: %0.3f [1/rad]\n\tZero-lift AoA: %0.3f [deg]\n\n',NACA_num,m_L,rad2deg(AoA_L0));

% plot NACA 2424 sectional coefficient of lift vs angle of attack
plot([-5, 0, 5, 10],[c_l_neg5, c_l, c_l_5, c_l_10],'linewidth',2,'DisplayName',['NACA ' NACA_num]);
xlabel('AoA [deg]');
ylabel('C_{l}');
title('NACA Airfoils Lift Slopes [#3]');
hold on
grid on
legend('location','northwest');
hold off

%% End of Housekeeping
toc
