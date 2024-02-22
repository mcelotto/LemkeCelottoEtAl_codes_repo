%
%  This source code is part of:
%  NIT - Neuroscience Information Toolbox
%  Copyright (C) 2020  Roberto Maffulli, Miguel Angel Casal Santiago
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%%% # Mixed-vine Copula tutorials
%%% ## Tutorial 1
%%% Tutorial 1 uses the MVC to reproduce figure 1 of the paper: *Onken, Arno, and Stefano Panzeri. "Mixed vine copulas as joint models of spike counts and local field potentials." Advances in Neural Information Processing Systems. 2016.*

clear all
close all
clc

% Construct 3D mixed copula vine
disp('Constructing mixed copula vine...');
d = 3; % Dimension
vine.type = 'c-vine'; % Canonical vine type
% Set margins
vine.margins = cell(d,1);
% Standard normal margin
vine.margins{1}.dist = 'norm';
vine.margins{1}.theta = [0;1];
vine.margins{1}.iscont = true; % Continuous margin
% Poisson margin
vine.margins{2}.dist = 'poiss';
vine.margins{2}.theta = 5;
vine.margins{2}.iscont = false; % Discrete margin
% Gamma margin
vine.margins{3}.dist = 'gam';
vine.margins{3}.theta = [2;4];
vine.margins{3}.iscont = true; % Continuous margin
% Set copula families
vine.families = cell(d);
vine.theta = cell(d);
% Gaussian copula family
vine.families{1,2} = 'gaussian';
vine.theta{1,2} = 0.5;
% Student copula family
vine.families{1,3} = 'student';
vine.theta{1,3} = [0.5;2];
% Clayton copula family
vine.families{2,3} = 'clayton';
vine.theta{2,3} = 5;
fprintf('\n');
% Construct argument to specify which margins are continuous
iscont = false(d,1);
for i = 1:d
    iscont(i) = vine.margins{i}.iscont;
end

% Test probability density function
disp('Calculating probability density function on a grid...');
% Calculate probability density function on lattice
x1gv = linspace(-3,3,100);
x2gv = 0:10;
x3gv = linspace(0.5,25,100);
[x1,x2,x3] = ndgrid(x1gv,x2gv,x3gv);
dx = [abs(x1gv(1)-x1gv(2))...
    abs(x2gv(1)-x2gv(2))...
    abs(x3gv(1)-x3gv(2))];
p = mixedvinepdf(vine,[x1(:),x2(:),x3(:)]);
p = reshape(p,size(x1));
% Plot 2D margins
figure('Name','2D margins of mixed vine copula PDF','Position',[0,0,1600,1000]);

margin12 = marginalize_pdf(p,[1,2],iscont,x1gv,x2gv,x3gv);
subplot(3,3,1);
imagesc(x1gv,x2gv,margin12');
colorbar
set(gca,'YDir','normal');
xlabel('Margin 1');
ylabel('Margin 2');

margin13 = marginalize_pdf(p,[1,3],iscont,x1gv,x2gv,x3gv);
subplot(3,3,2);
imagesc(x1gv,x3gv,margin13');
colorbar
set(gca,'YDir','normal');
xlabel('Margin 1');
ylabel('Margin 3');

margin23 = marginalize_pdf(p,[2,3],iscont,x1gv,x2gv,x3gv);
subplot(3,3,3);
imagesc(x2gv,x3gv,margin23');
colorbar
set(gca,'YDir','normal');
xlabel('Margin 2');
ylabel('Margin 3');

% Test sampling
disp('Sampling from mixed copula vine...');
% Draw samples
cases = 200;
x = mixedvinernd(vine,cases);
% Plot samples in 2D
subplot(3,3,4);
scatter(x(:,1),x(:,2),20,[0 0 0],'filled');
xlabel('Margin 1');
ylabel('Margin 2');
xlim([-3 3])
ylim([0 10])
subplot(3,3,5);
scatter(x(:,1),x(:,3),20,[0 0 0],'filled');
xlabel('Margin 1');
ylabel('Margin 3');
xlim([-3 3])
ylim([0 25])
subplot(3,3,6);
scatter(x(:,2),x(:,3),20,[0 0 0],'filled');
xlabel('Margin 2');
ylabel('Margin 3');
xlim([0 10])
ylim([0 25])
fprintf('\n');

% Fit copula on draw data
vine_fit = mixedvinefit(x,vine.type,iscont);
p = mixedvinepdf(vine_fit,[x1(:),x2(:),x3(:)]);
p = reshape(p,size(x1));
% Plot 2D margins
margin12 = marginalize_pdf(p,[1,2],iscont,x1gv,x2gv,x3gv);
subplot(3,3,7);
imagesc(x1gv,x2gv,margin12');
colorbar
set(gca,'YDir','normal');
xlabel('Margin 1');
ylabel('Margin 2');
margin13 = marginalize_pdf(p,[1,3],iscont,x1gv,x2gv,x3gv);
subplot(3,3,8);
imagesc(x1gv,x3gv,margin13');
colorbar
set(gca,'YDir','normal');
xlabel('Margin 1');
ylabel('Margin 3');
margin23 = marginalize_pdf(p,[2,3],iscont,x1gv,x2gv,x3gv);
subplot(3,3,9);
imagesc(x2gv,x3gv,margin23');
colorbar
set(gca,'YDir','normal');
xlabel('Margin 2');
ylabel('Margin 3');

