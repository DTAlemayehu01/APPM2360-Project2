clc
clear
close all

%% Part 2.1 - Samuel H. Meyn

%% 2.1.1
%
% Data Plots
%

%Import all data
Z = importdata('mariana_depth.csv');
Y = importdata('mariana_longitude.csv');
X = importdata('mariana_latitude.csv');
Z = transpose(Z);
Zkm = Z .* 1/1000;

figure(1)

Levels = [-11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 11];
contour(Y, X, Zkm, Levels)
title('Mariana Depth along Latitude and Longitude');
xlabel('Longnitude (deg)');
ylabel('Latitude (deg)');
cb = contourcbar("eastoutside");
cb.XLabel.String = "Elevation (km)";

figure(2)

surf(Y,X,Zkm, 'EdgeColor','none')
title('Mariana Depth along Latitude and Longitude');
xlabel('Longnitude (deg)');
ylabel('Latitude(deg)');
zlabel('Elevation (meters)');
colorbar;

%% 2.1.2

% Finding Maximal Depth
Min= min(Z, [], 'all');
for i=1:length(Z)
    for j=1:width(Z)
        if Z(i,j) == Min
            break
        end
    end
end

maxDepthLatitude = j;
maxDepthLongitude = i;
maxDepthLongitude
maxDepthLatitude
Min

%% 2.1.3

% Mean depth sub -6km
meanNum = 0;
N = 0;

for i=1:length(Z)
    for j=1:width(Z)
        if Z(i,j) < -6000
            meanNum = meanNum + Z(i,j);
            N = N + 1;
        end
    end
end

meanTrenchDepth = meanNum/N;
meanTrenchDepth/1000
