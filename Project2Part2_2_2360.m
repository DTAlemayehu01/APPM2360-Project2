clear;
clc;
close all;

A = readmatrix('mariana_depth.csv');
%dataLat = readmatrix('mariana_latitude.csv');
%dataLong = readmatrix('mariana_longitude.csv');


%% Part 1
tic
%Transpose A
A_T = A';
%Create Unit vector
n = 1;
x = randn(n,1);
x = x / norm(x);

%Iterate until the unit vector for n+1 minus the previous unit vector
%doesn't change
for i = 1:10
unit_x = A_T*A*x;
unit_x = unit_x / (norm(unit_x));

if  norm(unit_x - x) <= 0
break;
end

x = unit_x;

end
%Show final eigenvalue and eigenvector
eigen_final = x;
k = eig(A'*A);
x_n_vector = sqrt(sum(x^2));
N = length(x);
hold on;
plot(1:N, x_n_vector, 'k-');
toc
%% Part 2

v_50 = zeros(50,1);

%Initialize a unit vector
n = 1;
u = randn(n,1);
u = u / norm(u);
% 
for i = 1:50
u_new_asterisk = A'*A*u;
u_new = u_new_asterisk - sum((u_new_asterisk'*v_50(i)*v_50(i)));
u_new = u_new / norm(u_new);
if norm(u_new - u) < 0.002
break;
end

u = u_new;
end

v_i = u_new;

figure(2)
semilogy(1:50, trimdata(k, [50 1]), 'k-');
grid