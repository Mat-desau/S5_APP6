clc
clear all
close all

%% Problème 6
%E17
x = [-2 -1.5 -1 -0.5 0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0];
y = [14 8.75 5.0 2.75 2.0 2.75 5.0 8.75 14 20.75 29 38.75 50];

h = (x(2) - x(1));

F_trap(1) = y(1);
for n = 2:length(y)
    F_trap(n) = (y(1) + y(n) + 2*sum(y(2:(n-1))))*(h/2);
end

Diff1(1) = (y(2) - y(1))/h;
Diff1(2) = (y(end) - y(end-1))/h;

Erreur = ((h^2)/12)*(Diff1(end)-Diff1(1))

%E18
F_simp(1) = 0;
x_simp(1) = x(1);
for n = 3:2:length(y)
    F_simp(((n-1)/2)+1) = (y(1) + y(n) + 4*sum(y(2:2:(n-1))) + 2*sum(y(3:2:(n-1))))*(h/3);
    x_simp(((n-1)/2)+1) = x(n);
end

Diff3(1) = (y(4) - 3*y(3) + 3*y(2) - y(1)) / (h^3);
Diff3(2) = (y(end) - 3*y(end-1) + 3*y(end-2) - y(end-3)) / (h^3);

Erreur = ((h^4)/180)*(Diff3(2) - Diff3(1));

% figure
% scatter(x, F_trap)
% 
% figure
% scatter(x_simp, F_simp)

%% Problème 7
%E20
clc
clear all
close all

t = [0 1 2 3 4 5 6];
a = [0.7 0.9 0.5 0.1 0.3 1.7 4.9];

v_ini = 4;

h = (t(2) - t(1));

F_trap(1) = v_ini;
for n = 2:length(t)
    F_trap(n) = v_ini + (a(1) + a(n) + 2*sum(a(2:(n-1))))*(h/2);
    
end

Diff1(1) = (a(2) - a(1))/h;
Diff1(2) = (a(end) - a(end-1))/h;

Erreur = ((h^2)/12)*(Diff1(end)-Diff1(1))

%E18
F_simp(1) = v_ini;
x_simp(1) = t(1);
for n = 3:2:length(t)
    F_simp(((n-1)/2)+1) = v_ini + (a(1) + a(n) + 4*sum(a(2:2:(n-1))) + 2*sum(a(3:2:(n-1))))*(h/3);
    x_simp(((n-1)/2)+1) = t(n);
end

Diff3(1) = (a(4) - 3*a(3) + 3*a(2) - a(1)) / (h^3);
Diff3(2) = (a(end) - 3*a(end-1) + 3*a(end-2) - a(end-3)) / (h^3);

Erreur = ((h^4)/180)*(Diff3(2) - Diff3(1))

%% Problème 8
clc
close all
clear all

x_ini = 0.85;
%Pour pas que ce soit 0 en partant
x = x_ini;
F = (x.^3) - (6*x.^2) + (7*x) + 2;
D = (3*x.^2) - (12*x) + (7);

tol = 1*10^(-8);
iter = 0;
Max_iter = 100;

while (iter < Max_iter) && (abs(F) > tol)
    x = x - F/D;
    F = (x.^3) - (6*x.^2) + (7*x) + 2;
    D = (3*x.^2) - (12*x) + (7);
    iter = iter + 1;
end

%% Problème 9
clc 
clear all
close all

P = [262.615 195.318 174.949 155.430 150.345 153.575 188.786 221.194 242.943 280.956 332.294];
v = [20.000 20.628 22.553 25.894 30.862 37.769 47.048 59.284 75.244 95.931 122.646];
h = [0 100 200 300 400 500 600 700 800 900 1000];

Y = log((2.*P)./(v.^2));
X = h;

A = [length(v) sum(X);
     sum(X)    sum(X.^2)];
C = [sum(Y);
     sum(Y.*X)];

rep = pinv(A) * C; 

b = rep(1);
m = rep(2);

P0 = exp(b)
hs = -1/m