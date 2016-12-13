%% This solves a non-linear heat equations of the form ut-triangle u = 1 using Matlab pde toolbox
close all; clear all; clc;

R1 = [3;4;-1;1;1;-1;-1;-1;1;1];
C1 = [1;0;0;0.4];
C1 = [C1;zeros(length(R1) - length(C1),1)];
gd = [R1,C1];
sf = 'R1+C1';
ns = char('R1','C1')';
g = decsg(gd,sf,ns);

numberOfPDE = 1;
pdem = createpde(numberOfPDE);
geometryFromEdges(pdem,g);

figure
pdegplot(pdem,'EdgeLabels','on','EdgeLabels','on')
axis([-1.1 1.1 -1.1 1.1]);
axis equal
title 'Geometry With Edge and Subdomain Labels'

% Solution is zero at all four outer edges of the square
applyBoundaryCondition(pdem,'Edge',(1:4),'u',0);
specifyCoefficients(pdem,'m',0,...
                         'd',1,...
                         'c',1,...
                         'a',0,...
                         'f',1);
                     
setInitialConditions(pdem,1,'Face',2);
msh = generateMesh(pdem);
figure;
pdemesh(pdem);
axis equal