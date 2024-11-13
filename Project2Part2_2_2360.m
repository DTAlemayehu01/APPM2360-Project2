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
x = randn(length(A'*A),1);
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
eigen_first = x;
x_n_vector = sqrt(sum(x.^2));
N = length(x);
hold on;
plot(1:N, x_n_vector, 'k-');
toc
%% Part 2
v_matrix = zeros(length(A'*A),50);
v_vals50 = zeros(50,1);
%Initialize a unit vector
% 

for i = 1:50
u = randn(length(A'*A),1);
u = u / norm(u);
    while true
    u_new_asterisk = A'*A*u;
        for j = 1:(i-1)
          u_new_asterisk = u_new_asterisk - (u_new_asterisk'*v_matrix(:,j))*v_matrix(:,j);
        end
        u_new = u_new_asterisk/ norm(u_new_asterisk);
        if norm(u_new - u) < 0.002
            break;
        end
    u = u_new;

    end
v_matrix(:,i) = u_new;
v_vals50(i) = (u_new')*(A')*(A)*(u_new);
end
v_i = u_new;

figure(2)
semilogy(1:50, trimdata(v_vals50, [50 1]), 'k-');
grid
