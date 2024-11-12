clc
clear
close all

%% Part 2.1 - Samuel H. Meyn

%% 2.1.1

%Import all data
Z = importdata('mariana_depth.csv');
Y = importdata('mariana_longitude.csv');
X = importdata('mariana_latitude.csv');
Z = transpose(Z);
Zkm = Z .* 1/1000;

figure(1)

contour(Y, X, Zkm);
xlabel('Longnitude [deg]');
ylabel('Latitude[deg]');
colorbar;

figure(2)

surf(Y,X,Z, 'EdgeColor','none')
xlabel('Longnitude [deg]');
ylabel('Latitude[deg]');
colorbar

%% 2.1.2

M= min(Z, [], 'all');
for i=1:length(Z)
    for j=1:width(Z)
        if Z(i,j) == M
            break
        end
    end
end

maxDepthLatitude = j;
maxDepthLongitude = i;

%% 2.1.3

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
