%% ========================================================================
% SF2565 Program construction in C++ for Scientific Computing
% Project 4: Caroline Eriksson
% Date: December 2019
%% ========================================================================
clear all; close all; clc

% Load interior files
load interior_x.txt -ascii
load interior_y.txt -ascii

% Load differentials
load differential_x.txt -ascii
load differential_y.txt -ascii

% Load second order differentials
load differential2_x.txt -ascii
load differential2_y.txt -ascii

% Load Laplacian
load laplace.txt -ascii

% ========================================================================

% Plot of the differential with regards to x.
% Creating vertical grid lines
figure(1)
for i = 1:50 
   plot3(interior_x(i,:), interior_y(i,:), differential_x(i,:), 'm')
   hold on
end

% Creating horizontal grid lines
for i = 1:20 
   plot3(interior_x(:,i), interior_y(:,i), differential_x(:,i), 'm')
   hold on
end

plot3(interior_x(:,:), interior_y(:,:), differential_x(:,:), 'k.')
grid on
title('First order differential - x')
xlabel('x')
ylabel('y')
zlabel('du/dx')
axis([-10 5 0 3 -0.4 0.28]);
hold off

% Plot of the differential with regards to y.
% Creating vertical grid lines
figure(2)
for i = 1:50 
   plot3(interior_x(i,:), interior_y(i,:), differential_y(i,:), 'm')
   hold on
end

% Creating horizontal grid lines
for i = 1:20 
   plot3(interior_x(:,i), interior_y(:,i), differential_y(:,i), 'm')
   hold on
end

plot3(interior_x(:,:), interior_y(:,:), differential_y(:,:), 'k.')
grid on
title('First order differential - y')
xlabel('x')
ylabel('y')
zlabel('du/dy')
axis([-10 5 0 3 0.5 1.5]);
hold off


% Creating a plot of the domain

X = [-10:0.306122:5];
Y = [0:0.1578947368:3];
U = zeros(length(X), length(Y));

for i = 1:50
    for j = 1:20
        U(i,j) = sin((interior_x(i,j)/10)^2)*cos(interior_x(i,j)/10) + interior_y(i,j);
    end
end

figure(3)
for i = 1:50 
   plot3(interior_x(i,:), interior_y(i,:), U(i,:), 'm')
   hold on
end

% Creating horizontal grid lines
for i = 1:20 
   plot3(interior_x(:,i), interior_y(:,i), U(:,i), 'm')
   hold on
end
plot3(interior_x(:,:), interior_y(:,:), U(:,:), 'm')
grid on
xlabel('x')
ylabel('y')
zlabel('u')
title('Domain')
hold off


% Computing the "exact" differentials with Matlab's diff function
syms x y
u = sin((x/10)^2)*cos(x/10) + y;
Diffx = diff(u,x);
Diffy = diff(u,y);

% Plots of Matlab's diff functions over our domain
Diff_x = zeros(length(X), length(Y));

for i = 1:50
    for j = 1:20
        Diff_x(i,j) = (interior_x(i,j)*cos(interior_x(i,j)/10)*cos(interior_x(i,j)^2/100))/50 - (sin(interior_x(i,j)/10)*sin(interior_x(i,j)^2/100))/10;
    end
end

% Creating vertical grid lines
figure(4)
for i = 1:50 
   plot3(interior_x(i,:), interior_y(i,:), Diff_x(i,:), 'm')
   hold on
end

% Creating horizontal grid lines
for i = 1:20 
   plot3(interior_x(:,i), interior_y(:,i), Diff_x(:,i), 'm')
   hold on
end
plot3(interior_x(:,:), interior_y(:,:), Diff_x(:,:),'k')
grid on
xlabel('x')
ylabel('y')
zlabel('du/dx')
title("Matlab's exact differentiation - x")
axis([-10 5 0 3 -0.4 0.28]);
hold off


Diff_y = ones(length(X), length(Y));

% Creating vertical grid lines
figure(5)
for i = 1:50 
   plot3(interior_x(i,:), interior_y(i,:), Diff_y(i,:), 'm')
   hold on
end

% Creating horizontal grid lines
for i = 1:20 
   plot3(interior_x(:,i), interior_y(:,i), Diff_y(:,i), 'm')
   hold on
end
plot3(interior_x(:,:), interior_y(:,:), Diff_y(:,:),'k')
grid on
xlabel('x')
ylabel('y')
zlabel('du/dy')
title("Matlab's exact differentiation - y")
axis([-10 5 0 3 0.5 1.5]);
hold off

% ========================================================================
% Computing the error for du/dx

error_x = zeros(length(X),length(Y));
for i = 1:length(X)
    for j = 1:length(Y) 
        error_x(i,j) = abs(Diff_x(i,j) - differential_x(i,j));
    end
end

% Computing the error for du/dy

error_y = zeros(length(X),length(Y));
for i = 1:length(X)
    for j = 1:length(Y) 
        error_y(i,j) = abs(Diff_y(i,j) - differential_y(i,j));
    end
end

% Error plot for du/dx
% Creating vertical grid lines
figure(6)
for i = 1:50 
   plot3(interior_x(i,:), interior_y(i,:), error_x(i,:), 'k')
   hold on
end
 
% Creating horizontal grid lines
for i = 1:20 
   plot3(interior_x(:,i), interior_y(:,i), error_x(:,i), 'k')
   hold on
end
plot3(interior_x(:,:), interior_y(:,:), error_x(:,:), 'k.')
xlabel('x')
ylabel('y')
zlabel('Error')
grid on
title('Error - x')
hold off


% Error plot for du/dy
% Creating vertical grid lines
figure(7)
for i = 1:50 
   plot3(interior_x(i,:), interior_y(i,:), error_y(i,:), 'k')
   hold on
end

% Creating horizontal grid lines
for i = 1:20 
   plot3(interior_x(:,i), interior_y(:,i), error_y(:,i), 'k')
   hold on
end
plot3(interior_x(:,:), interior_y(:,:), error_y(:,:), 'k.')
grid on
xlabel('x')
ylabel('y')
zlabel('Error')
title('Error - y')
hold off

% ========================================================================
% Laplace

% Creating vertical grid lines
figure(8)
for i = 1:50 
   plot3(interior_x(i,2:19), interior_y(i,2:19), laplace(i,2:19), 'm')
   hold on
end
 
% Creating horizontal grid lines
for i = 2:19 
   plot3(interior_x(:,i), interior_y(:,i), laplace(:,i), 'm')
   hold on
end

plot3(interior_x(:,2:19), interior_y(:,2:19), laplace(:,2:19), 'k.')
grid on
title('The computed Laplacian')
xlabel('x')
ylabel('y')
zlabel('Laplace')
%axis([-10 5 0 3 -0.06 0.025]);
hold off

% Computing the "exact" Laplacian with Matlab's diff functions
syms x y
u = sin((x/10)^2)*cos(x/10) + y;
Diffx = diff(u,x);
Diff2x = diff(Diffx,x);
Diffy = diff(u,y);
Diff2y = diff(Diffy,y);

Laplace_m = zeros(length(X),length(Y));
Diff2_x = zeros(length(X),length(Y));

for i = 1:50
    for j = 1:20
        Diff2_x(i,j) = (cos(interior_x(i,j)/10)*cos(interior_x(i,j)^2/100))/50 - (cos(interior_x(i,j)/10)*sin(interior_x(i,j)^2/100))/100 - (interior_x(i,j)^2*cos(interior_x(i,j)/10)*sin(interior_x(i,j)^2/100))/2500 - (interior_x(i,j)*sin(interior_x(i,j)/10)*cos(interior_x(i,j)^2/100))/250;
        Laplace_m(i,j) = Diff2_x(i,j);
    end
end

% Plot of Matlab's Laplacian over our domain
figure(9)
%fsurf(Laplace_m, [-10 5 0 3])
for i = 1:50 
   plot3(interior_x(i,:), interior_y(i,:), Laplace_m(i,:), 'm')
   hold on
end
 
% Creating horizontal grid lines
for i = 1:20 
   plot3(interior_x(:,i), interior_y(:,i), Laplace_m(:,i), 'm')
   hold on
end
plot3(interior_x(:,:), interior_y(:,:), Laplace_m(:,:), 'k.')
grid on
title("Matlab's exact Laplacian")
xlabel('x')
ylabel('y')
zlabel('Laplace')
%axis([-10 5 0 3 -0.06 0.02])


% Computing the error for the Laplacian

error_laplace = zeros(length(X),length(Y));
for i = 1:length(X)
    for j = 1:length(Y)
        
    error_laplace(i,j) = abs(Laplace_m(i,j) - laplace(i,j));
    end
end

% Error plot for laplace
% Creating vertical grid lines
figure(10)
for i = 1:50 
   plot3(interior_x(i,2:19), interior_y(i,2:19), error_laplace(i,2:19), 'k')
   hold on
end
 
% Creating horizontal grid lines
for i = 2:19
   plot3(interior_x(:,i), interior_y(:,i), error_laplace(:,i), 'k')
   hold on
end
plot3(interior_x(:,2:19), interior_y(:,2:19), error_laplace(:,2:19), 'k.')
grid on
title('Error of the Laplacian')
xlabel('x')
ylabel('y')
zlabel('Error')
hold off
