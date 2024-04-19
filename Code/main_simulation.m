%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pressure calculation from skin friction for simulation cases
%
% 08/10/2022
% By Zemin Cai
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clc;
clear all;
addpath(GetAbsolutePath('..\Data'));
addpath(GetAbsolutePath('..\Results'));

%% Data preparing
% simulation data
CaseNum = 1;
alpha = [0.01 0.08 0.1 0.2 0.3 0.5];
switch CaseNum
    case 1
        %%%%%%%%%%%%%%% Falkner-Skan flow simulation %%%%%%%%%%%%%%%%%%
        %p_gt = load('..\Data\Falkner_Skan_Flow\dp_Falkner_Skan_m0p5.dat');   % m = 1/2
        %p_gt = load('..\Data\Falkner_Skan_Flow_m0\dp_Falkner_Skan_m0.dat');
        %p_gt = load('..\Data\Falkner_Skan_Flow_m1d3\dp_Falkner_Skan_m1d3.dat');   % m = 1/3
        p_gt = load('..\Data\Falkner_Skan_Flow_m2d3\dp_Falkner_Skan_m2d3.dat');   % m = 2/3
        
        %tau_x = load('..\Data\Falkner_Skan_Flow\tor_x_Falkner_Skan_m0p5.dat');
        %tau_x = load('..\Data\Falkner_Skan_Flow_m0\tor_x_Falkner_Skan_m0.dat');
        %tau_x = load('..\Data\Falkner_Skan_Flow_m1d3\tor_x_Falkner_Skan_m1d3.dat');
        tau_x = load('..\Data\Falkner_Skan_Flow_m2d3\tor_x_Falkner_Skan_m2d3.dat');
        
        %tau_y = load('..\Data\Falkner_Skan_Flow\tor_y_Falkner_Skan_m0p5.dat');
        %tau_y = load('..\Data\Falkner_Skan_Flow_m0\tor_y_Falkner_Skan_m0.dat');
        %tau_y = load('..\Data\Falkner_Skan_Flow_m1d3\tor_y_Falkner_Skan_m1d3.dat');
        tau_y = load('..\Data\Falkner_Skan_Flow_m2d3\tor_y_Falkner_Skan_m2d3.dat');
        
        alpha = 0.01;
        %FolderName = 'Falkner_Skan_Flow';
        %FolderName = 'Falkner_Skan_Flow_m0';
        %FolderName = 'Falkner_Skan_Flow_m1d3';
        FolderName = 'Falkner_Skan_Flow_m2d3';
        if ~exist(['..\Results\' FolderName], 'dir')
            mkdir(['..\Results\' FolderName]);
        end
        ResultsFileName = 'PressureFromSkinFriction_Falkner_Skan_Flow.mat';
        N=0.7       
    case 2
        %%%%%% 70 degree delta wing at the AoA of 20 deg. simulation %%%%%
        p_gt = load('..\Data\Delta_wing_simu\dp_Falkner_Skan_m0p5_delta.dat');
        tau_x = load('..\Data\Delta_wing_simu\tor_x_delta.dat');
        tau_y = load('..\Data\Delta_wing_simu\tor_y_delta.dat');
        alpha = 0.01;
        
        FolderName = 'Delta_wing_simu';
        if ~exist(['..\Results\' FolderName], 'dir')
            mkdir(['..\Results\' FolderName]);
        end
        ResultsFileName = 'PressureFromSkinFriction_Delta_wing_simu.mat';
        N=0.1
end

% downsmapling to accelarate the processing     
[height_full, width_full] = size(tau_x);
p_gt_downsampled = imresize(p_gt, N);
[p_x_downsampled, p_y_downsampled] = gradient(p_gt_downsampled);
tau_x_downsampled = imresize(tau_x, N);
tau_y_downsampled = imresize(tau_y, N);

%Phi_m_downsampled = tau_x_downsampled.*p_x_downsampled + tau_y_downsampled.*p_y_downsampled;
Phi_m_downsampled = -0.1*ones(size(tau_x_downsampled));   % you can change the factor. e.g. 0.01

[height_downsampled, width_downsampled] = size(tau_x_downsampled);
tau_div_downsampled = divergence(tau_x_downsampled, tau_y_downsampled);
[Phi_m_x_downsampled, Phi_m_y_downsampled] = gradient(Phi_m_downsampled);

%% Calculating
alpha_x = alpha;
alpha_y = alpha;

A = tau_x_downsampled.^2 + alpha_x;
B = tau_x_downsampled.*tau_y_downsampled;
C = tau_y_downsampled.^2 + alpha_y;
[tau_x_x_downsampled, tau_x_y_downsampled] = gradient(tau_x_downsampled); 
D = tau_x_downsampled.*tau_x_x_downsampled + tau_y_downsampled.*tau_x_y_downsampled + tau_x_downsampled.*tau_div_downsampled;
[tau_y_x_downsampled, tau_y_y_downsampled] = gradient(tau_y_downsampled); 
E = tau_x_downsampled.*tau_y_x_downsampled + tau_y_downsampled.*tau_y_y_downsampled + tau_y_downsampled.*tau_div_downsampled;
F = Phi_m_downsampled.*tau_div_downsampled + tau_x_downsampled.*Phi_m_x_downsampled+ tau_y_downsampled.*Phi_m_y_downsampled;
h = 1;

P_s1 = zeros(height_downsampled, width_downsampled);
P_s1(:, 1) = p_gt_downsampled(:, 1);
P_s1(1, :) = p_gt_downsampled(1, :);
P_s1(height_downsampled, :) = p_gt_downsampled(height_downsampled, :);
P_s1(:, width_downsampled) = p_gt_downsampled(:, width_downsampled);

% Creating the coefficient matrix
disp(sprintf('Now, creating the coeffcient matrix...'));
E_CoefMatrix = CoeffMatrix(height_downsampled, width_downsampled, A, B, C, D, E, h);
size(E_CoefMatrix)

% Creating the right hand side vector
disp(sprintf('Now, creating the right hand side vector...'));
T = RHS_Vec(P_s1, A, B, C, D, E, F, h);

% solving the equations system
disp(sprintf('Now, solving the linear equations system...'));
tic
sol = linsolve(E_CoefMatrix, T);
toc

% merging
for j = 1:1:(height_downsampled-2)
    for i = 1:1:(width_downsampled-2)
        P_s1(j+1, i+1)= sol((i-1)*(height_downsampled-2)+j);
    end
end

% converting to the pressure
p_s1 = P_s1;

% resize the pressure field to the true size
p_s1_fullsize = imresize(p_s1, [height_full, width_full]);

RMSE = (sum(sum((p_gt_downsampled-p_s1).^2))/(height_downsampled*width_downsampled)).^0.5;
disp(sprintf('The Mean Square Error: %2.12f', RMSE));

%% Data saving
Phi_m = -0.1;
path = ['../Results/' FolderName '/' ResultsFileName];
save(path, 'p_gt', 'p_s1_fullsize', 'tau_x', 'tau_y', 'Phi_m', 'alpha', 'RMSE');

%% Displaying
figure(1);
s=[min(min(p_gt))*0.9999 max(max(p_gt))*1.0001];
imagesc(p_gt,s);
colormap(pink);
xlabel('x (pixels)', 'FontSize', 14);
ylabel('y (pixels)', 'FontSize', 14);
axis image;
colorbar;
title('The ground truth of pressure field', 'FontSize', 14);

figure(2);
s=[min(min(p_s1_fullsize))*0.9999 max(max(p_s1_fullsize))*1.0001];
imagesc(p_s1_fullsize,s);
colormap(pink);
xlabel('x (pixels)', 'FontSize', 14);
ylabel('y (pixels)', 'FontSize', 14);
axis image;
colorbar;
title('The calculated pressure field', 'FontSize', 14);

figure(3);
imagesc(A);
colormap(pink);
xlabel('x (pixels)');
ylabel('y (pixels)');
axis image;
colorbar;
title('A matrix');

figure(4);
imagesc(B);
colormap(pink);
xlabel('x (pixels)');
ylabel('y (pixels)');
axis image;
colorbar;
title('B matrix');

figure(5);
imagesc(C);
colormap(pink);
xlabel('x (pixels)');
ylabel('y (pixels)');
axis image;
colorbar;
title('C matrix');

figure(6);
imagesc(D);
colormap(pink);
xlabel('x (pixels)');
ylabel('y (pixels)');
axis image;
colorbar;
title('D matrix');

figure(7);
imagesc(E);
colormap(pink);
xlabel('x (pixels)');
ylabel('y (pixels)');
axis image;
colorbar;
title('E matrix');

figure(8);
imagesc(F);
colormap(pink);
xlabel('x (pixels)');
ylabel('y (pixels)');
axis image;
colorbar;
title('F matrix');

figure(9);
s=[min(min(tau_x_downsampled))*1.01 max(max(tau_x_downsampled))*1.01];
imagesc(tau_x_downsampled,s);
colormap(pink);
xlabel('x (pixels)', 'FontSize',14);
ylabel('y (pixels)', 'FontSize',14);
axis image;
colorbar;
title('$m_{x}$ component of the skin friction vector', 'Interpreter', 'latex',  'FontSize',13,'FontWeight','bold');

p_recovered = p_s1_fullsize;
[rows, cols] = size(p_s1);
x_vec = linspace(1, cols, cols);
y_vec = linspace(1, rows, rows);

[rows_fullsize, cols_fullsize] = size(p_s1_fullsize);
x_vec_fullsize = linspace(1, cols_fullsize, cols_fullsize);
y_vec_fullsize = linspace(1, rows_fullsize, rows_fullsize);

figure(11);
s=[min(min(p_s1_fullsize))*0.9999 max(max(p_s1_fullsize))*1.0001];
imagesc(imresize(p_s1_fullsize,0.1),s);
colormap(pink);
xlabel('x (pixels)', 'FontSize', 14);
ylabel('y (pixels)', 'FontSize', 14);
axis image;
hold on;
figure(11);
d=10;
h=quiver(tau_x(1:d:end,1:d:end), tau_y(1:d:end,1:d:end),3);
set(h,'Color','red')
xlabel('x (pixels)', 'FontSize', 14);
ylabel('y (pixels)', 'FontSize', 14);
axis image;
set(gca,'YDir','reverse');
%title('The Ground Truth of Skin Friction and Pressure of AoA = 10 [deg]', 'FontSize', 14);

figure(12);
s=[min(min(p_s1_fullsize))*0.9999 max(max(p_s1_fullsize))*1.0001];
imagesc(p_s1_fullsize,s);
colormap(pink);
colorbar; 
xlabel('x (pixels)', 'FontSize', 14);
ylabel('y (pixels)', 'FontSize', 14);
axis image;
hold on;
figure(12);
[x,y]=meshgrid(1:width_full,1:height_full);
x30=(x-0);
y30=(y-0);
dn=5;
dm=5;
[sx,sy]=meshgrid(1:dn:width_full,1:dm:height_full);
h=streamslice(x30, y30, tau_x, tau_y, 20);
set(h, 'Color', 'red');
xlabel('x (pixels)', 'FontSize', 14);
ylabel('y (pixels)', 'FontSize', 14);
axis image;
set(gca,'YDir','reverse');
%title('Skin Friction Lines and Extracted Pressure of AoA = 10 [deg]', 'FontSize', 14);