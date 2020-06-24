clear all
close all
%%
%{
    This script is to test the instantaneous wavenumber in a simple 2D
    example.

    author: Jannis.Bulling@bam.de

    This script is to test the instantaneous wavenumber in a simple 2D
    example.
    Copyright (C) 2020  Jannis Bulling

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
%}
%%
% choose one of tree defect geometries.
% CsStr = 'circle';
% CsStr = 'diamond';
CsStr = 'rectangle';

% half the window size for the median filter
N_FILTER = 2^4;
% the window size for the local wavenumber map
N_WINDOW = 2^5;
% the number for the padded fft inside the window of the local wavenumber
% approach
N_FFT_WIN = 2^8;

% length of the signal
Lx   = 15;
Ly   = 10;
% number of point in both directions
N_X  = 2^8;
N_Y  = 2^8;
%% signal parameters
% 1st wavenumber for char(x) = 0
k0 = 05;
% 2nd wavenumber for char(x) = 1
k1 = 08; 
% 1st angle for char(x) = 0
phi0 = pi/3;
% 1st angle for char(x) = 0
phi1 = pi/2; 

% additional parameters to def. the signal
kv0 = k0*[cos(phi0), sin(phi0)];
kv1 = k1*[cos(phi1), sin(phi1)];

% def. x and y
x    = linspace(-Lx/2, Lx/2, N_X + 1).'; x(end) = [];
dx   = Lx / (N_X-1);
y    = linspace(-Ly/2, Ly/2, N_Y + 1); y(end) = [];
dy   = Ly / (N_Y-1);

% get a grid
[X,Y] = ndgrid(x,y);
R = sqrt( X.^2 + Y.^2 );

% define the characteristic funktion for the different wavenumbers
% interval boundaries
ell = [2, 4];
switch CsStr
    case 'circle'
        char = (R < ell(1)); 
        phi = linspace(0, 2*pi, 101);
        PlotDefect = @() plot3(ell(1)*cos(phi), ell(1)*sin(phi), k1+0*phi, 'r-');
    case 'rectangle'
        char = ( abs(X) < ell(1) ).*( abs(Y) < ell(2) );
        PlotDefect = @() plot3(ell(1)*[1 -1 -1 1 1], ell(2)*[1 1 -1 -1 1], k1*ones(1,5), 'r-');
    case 'diamond'
        char = (abs(X) + abs(Y) < ell(1));
        PlotDefect = @() plot3(ell(1)*[1 0 -1 0 1], ell(1)*[0 1 0 -1 0], k1*ones(1,5), 'r-');
    otherwise
        error('Unknown Case!')
end

% signal for char(x) = 0
signal0 = real(exp(1i*kv0(1)*x)*exp(1i*kv0(2)*y));
% signal for char(x) = 1
signal1 = real(exp(1i*kv1(1)*x)*exp(1i*kv1(2)*y));

% define the final signal and the reference
Img = signal0 .* (1 - char) + signal1 .* char;
ref =      k0 .* (1 - char) +      k1 .* char;

% get insnstantaneous wavenumber
disp('unfiltered insnstantaneous wavenumber')
tic
IW_unfiltered = GetInstantaneousWavenumber(Img, [dx, dy]);
toc
disp('median filtered insnstantaneous wavenumber')
IW_median = GetMedianFilter(IW_unfiltered, N_FILTER);
toc
disp('local wavenumber point-wise moving window')
tic
[local_k_pw, k_window] = GetMaxLocalWavenumber_PointwiseWindow(Img, [dx, dy], N_WINDOW, N_FFT_WIN);
toc
disp('local wavenumber line-wise moving window')
tic
[local_k_lw, ~] = GetMaxLocalWavenumber_LinewiseWindow(Img, [dx, dy], N_WINDOW, N_FFT_WIN);
toc

% compute the errors
error_IW_unfiltered = norm(IW_unfiltered - ref)/norm(ref);
error_IW_median     = norm(IW_median     - ref)/norm(ref);
error_local_k_pw    = norm(local_k_pw    - ref)/norm(ref);
error_local_k_lw    = norm(local_k_lw    - ref)/norm(ref);

% print errors
disp('unfiltered insnstantaneous wavenumber')
fprintf('error: %5.2e \n', error_IW_unfiltered)
disp('median filtered insnstantaneous wavenumber')
fprintf('error: %5.2e \n', error_IW_median)
disp('local wavenumber point-wise moving window')
fprintf('error: %5.2e \n', error_local_k_pw)
disp('local wavenumber line-wise moving window')
fprintf('error: %5.2e \n', error_local_k_lw)
%%
figure
    subplot(2,3,1)
        hold on
        surf(x, y, Img.', 'LineStyle', 'none'); view(2);
        axis equal
        axis tight
        colorbar
        PlotDefect();
        title('wave field')
    subplot(2,3,4)
        hold on
        surf(x, y, ref.', 'LineStyle', 'none'); view(2);
        set(gca, 'CLim', [0 10])
        axis equal
        axis tight
        colorbar
        PlotDefect();
        title('reference')
    subplot(2,3,2)
        hold on
        surf(x, y, IW_unfiltered.', 'LineStyle', 'none'); view(2);
        set(gca, 'CLim', [0 10])
        axis equal
        axis tight
        colorbar
        PlotDefect();
        title({'unfiltered','instantaneous wavenumber'})
    subplot(2,3,5)
        hold on
        surf(x, y, IW_median.', 'LineStyle', 'none'); view(2);
        set(gca, 'CLim', [0 10])
        axis equal
        axis tight
        colorbar
        PlotDefect();
        title({'median filtered','instantaneous wavenumber'})
    subplot(2,3,3)
        hold on
        surf(x, y, local_k_pw.', 'LineStyle', 'none'); view(2);
        set(gca, 'CLim', [0 10])
        axis equal
        axis tight
        colorbar
        PlotDefect();
        title({'local k','point-wise moving window'})
    subplot(2,3,6)
        hold on
        surf(x, y, local_k_lw.', 'LineStyle', 'none'); view(2);
        set(gca, 'CLim', [0 10])
        axis equal
        axis tight
        colorbar
        PlotDefect();
        title({'local k','line-wise moving window'})

