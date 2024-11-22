clc 
clear all
close all

%% Problème 11 (Erreur local erreur globale)
e = [0.08443 0.02603 0.01048 0.00319 0.00040];
t = [0.05 0.04 0.03 0.02 0.01];

for n = 2:length(e)
    m = n-1;
    P(n-1) = log(e(n)/e(m)) / log(t(n)/t(m));
end

% figure
% plot(t(1:length(P)), P)

%% Problème 12 
load()