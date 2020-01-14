% Phase-consistent DMD
% Aditya Nair, Benjamin Strom, Bingni Brunton, Steven Brunton

%%%%%%%%%%%%%%%%%%%%%%%%Description%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates and plots concatenated and phase-consistent DMD modes for
% data collected from overlapping windows.
%%%%%%%%%%%%%%%%%%%%%%%%Description%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;close all;

%% Import libraries for Optimized DMD and phase-consistency analysis

addpath('../optdmd-master/src');
addpath('../optdmd-master/examples');
addpath('../Functions');

%% Load Dataset

load('../Data/cylinder_dataset.mat');      % Dataset
% d (Data in each overlapping window),
% mv (Mean flow data in each overlapping window)
% fullCoords (Full domain)

%% Set Parameters

comp           = {'VORT'};         % vorticity dataset
% comp           = {'U','V','W'};  % if data contains u,v,w data
temp_resln     = 229;              % Number of temporal snapshots
NO_FOV         = 2;                % number of FOV's
r              = 10;               % number of modes to keep
dt             = 0.08;             % time step between modes
t              = 0:dt:(temp_resln-1)*dt;  % snapshot timestamps
flag           = 1;                % 1 - Exact DMD, 2 - Optimized DMD

%% Create Data Matrix, X 

X = [];
m = zeros(NO_FOV,1);n = zeros(NO_FOV,1);
for ii = 1:NO_FOV
    if (length(comp) == 3)
        m(ii) = size(d(ii).U,1);
        n(ii) = size(d(ii).U,2);
        X = [X; reshape(d(ii).U,m(ii)*n(ii),temp_resln); ...
            reshape(d(ii).V,m(ii)*n(ii),temp_resln); ...
            reshape(d(ii).W,m(ii)*n(ii),temp_resln)];
    elseif (length(comp) == 1)
        m(ii) = size(d(ii).VORT,1);
        n(ii) = size(d(ii).VORT,2);
        X = [X; reshape(d(ii).VORT,m(ii)*n(ii),temp_resln)];
    end
end
X = double(X);

% Pull out the F.O.V. coordinates
for ii = 1:NO_FOV
    coords(ii).x = d(ii).X;
    coords(ii).y = d(ii).Y;
end

% needed to re-combine overlapping fields-of-view
overlap = FOVoverlap(coords);

%% Perform Multi-domain DMD analysis

if flag == 1
    modes          = 2:3;              % mode to use in phase calc
    [~,omega,Phi,b,~,~] = Exact_DMD(X(:,1:end-1),X(:,2:end),r,dt);    
elseif flag == 2
    modes          = 3:4;              % mode to use in phase calc
    %function to restrict DMD eigen values to imaginary only
    proxfun = @(alpha)(1i.*sign(imag(alpha)).*max(abs(imag(alpha)),0.01));
    opts = varpro_opts('maxiter',400); % increase opt dmd max iterations
    [Phi,omega,b] = optdmd(X,t,r,2,opts,[],[],[],[],proxfun);% run opt dmd
    [~,inds] = sort(abs(b),'descend');% sort modes by magnitude
    omega = omega(inds);
    Phi = Phi(:,inds);
    b = b(inds);
end

%% Get data for one mode set

sind = 1;
eind = m(1)*n(1);
phiFOV = [];
for ii = 1:NO_FOV
    for j = 1:length(comp)
        for k = 1:length(modes)
            phiTemp = reshape(Phi(sind:eind,modes(k)),m(ii),n(ii));
            a = real(phiTemp);
            c = imag(phiTemp);
            A = abs(phiTemp); % Amplitude          
            bkAk = b(modes(k)).*A;
            wk = imag(omega(modes(k)));
            phi0 = 2*atan2(c,a+A);         
            % Field-of-view specific amplitudes, phases, and frequencies
            phiFOV(ii).(comp{j})(k).phi0 = phi0;
            phiFOV(ii).(comp{j})(k).bkAk = bkAk;
            phiFOV(ii).(comp{j})(k).wk = wk;
        end
        sind = eind + 1;
        eind = eind + m(ii)*n(ii);
    end
    % interpolant objects, to be used later
    if (length(comp) == 3)
        F{ii} = griddedInterpolant(d(ii).X',d(ii).Y',phiFOV(ii).U(1).bkAk');
    elseif (length(comp) == 1)
        F{ii} = griddedInterpolant(d(ii).X',d(ii).Y',phiFOV(ii).VORT(1).bkAk');
    end
end

%% Phase-consistency analysis

T2 = 2*pi/abs(imag(omega(modes(1)))); % Time Period
t2 = [1:32]/32*T2;                    % Time vector used for phase alignment
opts = optimset('MaxIter',10*200*2,'MaxFunEvals',10*200*2);
phiAll = fminsearch(@(phiV)phaseconsDMD(comp,phiFOV,F,t2,...
    phiV,overlap,coords,mv),rand(1,NO_FOV).*T2,opts);

%% Reconstruct full field, no phase correction

paUnC = DMD_reconstruct(comp,phiFOV,zeros(1,NO_FOV),[],t2);
for ii = 1:NO_FOV
    paUnC(ii).X = d(ii).X;
    paUnC(ii).Y = d(ii).Y;
end
stitchUnC = combineFOVs(paUnC,fullCoords.X,fullCoords.Y);

%% Reconstruct full field-of-view, phase corrected

paC = DMD_reconstruct(comp,phiFOV,phiAll,[],t2);
for ii = 1:NO_FOV
    paC(ii).X = d(ii).X;
    paC(ii).Y = d(ii).Y;
end
stitchC = combineFOVs(paC,fullCoords.X,fullCoords.Y);

%% Plotting data

load('../Data/colormap_cylinder.mat');
r1 = 0.5;
q = linspace(0,2*pi,361);
sx = r1*cos(q); sy = r1*sin(q);

% At a particular time
% Plot reconstructed concatenated mode

FH = figure;hold on;
h=pcolor(fullCoords.X,fullCoords.Y,stitchUnC.VORT(:,:,20));shading interp;
set(h,'EdgeColor','None');
caxis([-1 1]);
fill(sx,sy,'k');
colormap(cmap);
axis off;axis equal;
xlim([-1 11.1]);ylim([-2.05 2.05]);
set(gcf, 'Units', 'Inches','PaperSize', [1.832, 0.64],'Position',...
    [0,0,1.832*1.4, 0.64*1.4])
print(FH, '-depsc','-loose',['concatenated_mode.eps']);close;

% Plot reconstructed phase-consistent mode

FH = figure;hold on;
h=pcolor(fullCoords.X,fullCoords.Y,stitchC.VORT(:,:,1));shading interp;
set(h,'EdgeColor','None');
caxis([-1 1]);
fill(sx,sy,'k');
colormap(cmap);
axis off;axis equal;
xlim([-1 11.1]);ylim([-2.05 2.05]);
set(gcf, 'Units', 'Inches','PaperSize', [1.832, 0.64],'Position',...
    [0,0,1.832*1.4, 0.64*1.4])
print(FH, '-depsc','-loose',['phase_consistent_mode.eps']);close;
