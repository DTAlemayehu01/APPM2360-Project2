clear;
clc;
close all;
%
% PART 1: Incomplete SVD Decomposition
%
A = readmatrix('mariana_depth.csv');
dataLat = readmatrix('mariana_latitude.csv');
dataLong = readmatrix('mariana_longitude.csv');

A_T = A';

V = zeros(1440, 50);

n = 1440;
u = randn(n,1);
u = u / norm(u);

for i = 1:50
    while 1
        u_new = A_T*A*u;
        summation = 0;
        for j = 1:i
            summation = summation + (u_new'*V(j))*V(j);
        end
        u_new = u_new - summation';
        u_new = u_new/norm(u_new);
        if norm(u_new - u) > 0.002
            break;
        end
        u = u_new;
    end
    V(:,i) = u_new';
end

[sigma] = eig(A*A_T, "matrix");
sigma = sqrt(sigma);
sigma = trimdata(sigma, [50 50]);

U = [];
for i = 1:50
    U(:,i) = A*V(:,i)/sigma(i,i);
end

figure(11)
spy(sigma)
grid
figure(12)
spy(U)
grid
figure(13)
spy(V')
grid

%
% PART 2: SVD Efficacy for 50 columns of U, V
%
SVD_space = numel(sigma) + numel(U) + numel(V);
A_space = numel(A);
space_efficiency = SVD_space/A_space;
space_efficiency

SVD_nzs = nnz(sigma) + nnz(U) + nnz(V);
A_nzs = nnz(A);
nnz_efficiency = SVD_nzs/A_nzs;
nnz_efficiency

%
% PART 3: Image/Contour Map
%
figure(1)
surf(dataLong,dataLat,(U*sigma*V'/1000)','EdgeColor','none')
colorbar
title('Mariana Depth along Latitude and Longitude, (Incomplete SVD n=50)');
xlabel('Longitude (degrees)');
ylabel('Latitude (degrees)');
zlabel('Elevation (meters)');
grid

figure(2)
Levels = [-11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 11];
%contour(dataLong,dataLat,(U*sigma*V'/1000)', Level)
contour(dataLong,dataLat,(U*sigma*V'/1000)')
title('Mariana Depth along Latitude and Longitude, (Incomplete SVD n=50)');
xlabel('Longitude (degrees)');
ylabel('Latitude (degrees)');
cb = contourcbar("eastoutside");
cb.XLabel.String = "Elevation (km)";
grid

%
% Part 4: U V column counts?
%

%
% U, V with 10 columns
%
V10 = zeros(1440, 10);

n = 1440;
u = randn(n,1);
u = u / norm(u);

for i = 1:10
    while 1
        u_new = A_T*A*u;
        summation = 0;
        for j = 1:i
            summation = summation + (u_new'*V10(j))*V10(j);
        end
        u_new = u_new - summation';
        u_new = u_new/norm(u_new);
        if norm(u_new - u) > 0.002
            break;
        end
        u = u_new;
    end
    V10(:,i) = u_new';
end

[sigma10] = eig(A*A_T, "matrix");

sigma10 = sqrt(sigma10);
sigma10 = trimdata(sigma,[10,10]);

U10 = [];
for i = 1:10
    U10(:,i) = A*V(:,i)/sigma10(i,i);
end

figure(3)
surf(dataLong,dataLat,(U10*sigma10*V10'/1000)','EdgeColor','none')
colorbar
title('Mariana Depth along Latitude and Longitude, (Incomplete SVD n=10)');
xlabel('Longitude');
ylabel('Latitude');
zlabel('Elevation (meters)');
grid

figure(4)
Levels = [-11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 11];
%contour(dataLong,dataLat,(U10*sigma10*V10'/1000)', Levels)
contour(dataLong,dataLat,(U10*sigma10*V10'/1000)')
title('Mariana Depth along Latitude and Longitude, (Incomplete SVD n=10)');
xlabel('Longitude (degrees)');
ylabel('Latitude (degrees)');
cb = contourcbar("eastoutside");
cb.XLabel.String = "Elevation (km)";
grid

% 
% U, V with 100 columns
% 
V100 = zeros(1440, 100);

n = 1440;
u = randn(n,1);
u = u / norm(u);

for i = 1:100
    while 1
        u_new = A_T*A*u;
        summation = 0;
        for j = 1:i
            summation = summation + (u_new'*V100(j))*V100(j);
        end
        u_new = u_new - summation';
        u_new = u_new/norm(u_new);
        if norm(u_new - u) > 0.002
            break;
        end
        u = u_new;
    end
    V100(:,i) = u_new';
end

[sigma100] = eig(A*A_T, "matrix");

sigma100 = sqrt(sigma100);
sigma100 = trimdata(sigma100,[100 100]);

U100 = [];
for i = 1:100
    U100(:,i) = A*V100(:,i)/sigma100(i,i);
end

figure(5)
surf(dataLong,dataLat,(U100*sigma100*V100'/1000)','EdgeColor','none')
colorbar
title('Mariana Depth along Latitude and Longitude, (Incomplete SVD (n=100)');
xlabel('Longitude');
ylabel('Latitude');
zlabel('Elevation (meters)');
grid

figure(6)
Levels = [-11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 11];
%contour(dataLong,dataLat,(U100*sigma100*V100'/1000)', Levels)
contour(dataLong,dataLat,(U100*sigma100*V100'/1000)')
title('Mariana Depth along Latitude and Longitude, (Incomplete SVD n=100)');
xlabel('Longitude (degrees)');
ylabel('Latitude (degrees)');
cb = contourcbar("eastoutside");
cb.XLabel.String = "Elevation (km)";
grid