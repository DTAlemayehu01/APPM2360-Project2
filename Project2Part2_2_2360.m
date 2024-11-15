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
n = (1440);
x = randn(n,1);
x = x / norm(x);

%{
Iterate until the unit vector for n+1 minus the previous unit vector
doesn't change
%}
for i = 1:10

%Initialize the new unit vector as x_new
x_new = A_T*A*x;
x_new = x_new / (norm(x_new));
    
    %Check x_new - x until the value is small
    if  norm(x_new - x) <= 0
        break;
    end
    
x = x_new;

end
%Show final eigenvalue and eigenvector
eigen_first = x'*A_T*A*x;
%x_n_vector = sqrt(sum(x^2));
N = length(x);
%plot of the final vector against the length of the matrix from 1 to N
hold on;
plot(1:N,x, 'k-');
xlabel('nth Component')
ylabel('nth Component Magnitude')
title('Vector Components as a Function of the nth Vector Component')
eigen_first
toc
%% Part 2
tic
%Initializes vectors as matrices of zeros
v_matrix = zeros(length(A_T*A),50);
v_vals50 = zeros(50,1);

%for loop to find eigenvalues and eigenvectors
for i = 1:50

%Initialize a unit vector
u = randn(length(A_T*A),1);
u = u / norm(u);

%checks if 
    while true
    u_new_asterisk = A_T*A*u;
        %iterated to compute eigenvalues for each column
        for j = 1:(i-1)
          u_new_asterisk = u_new_asterisk - (u_new_asterisk'*v_matrix(:,j))*v_matrix(:,j);
        end

        %Reinitialize unit vector so it can be checked against the previous
        u_new = u_new_asterisk/ norm(u_new_asterisk);

        %Checks eigenvalue for a tolerance less than 0.002
        if norm(u_new - u) < 0.002
            break;
        end
    %Iterates the new unit vector
    u = u_new;

    end
%{
Saves eigen vectors and values to their corresponding 1440 x 50 and 50 x 1
matrices
%}
v_matrix(:,i) = u_new;
v_vals50(i) = (u_new')*(A_T)*(A)*(u_new);

end
%final U_n+1 as v_i
v_i = u_new;

%eigenvectors as columns of 1440 x 50 matrix
V = v_matrix;

%semilog plot of the 50 eigenvalues computed 
figure(2)
semilogy(1:50, trimdata(v_vals50, [50 1]), 'k-');
xlabel('nth Largest Eigenvalue')
ylabel('Eigenvalues (Natural Log Scaled)')
title('Eigenvalues as a function of the nth Largest Eigenvalue')
grid

figure(3)
semilogx(1:50, trimdata(v_vals50, [50 1]), 'k-');
xlabel('nth Largest Eigenvalue (Natural Log Scaled)')
ylabel('Eigenvalues')
title('Eigenvalues as a function of the nth Largest Eigenvalue')
grid

toc
