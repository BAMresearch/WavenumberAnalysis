clear all
close all
%%
%{
    This script is to test the instantaneous wavenumber and the local 
    wavenumber approach in a simple 1d example.

    author: Jannis.Bulling@bam.de

    This script is to test the instantaneous wavenumber and the local 
    wavenumber approach in a simple 1d example.
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
% number of image points
N_X       = 2^9;
% half the window size for the median filter
N_FILTER  = 2^4;
% the window size for the local wavenumber map
N_WINDOW  = 2^5;
%%
% length of the image
Lx  = 20.3*pi;

% define the gird and step size
x  = linspace(0, Lx, N_X + 1); x(end) = [];
dx = ( x(end) - x(1) )/(N_X - 1);

% define the characteristic funktion for the different wavenumbers
% interval boundaries
ell  = [8, 12]*pi;
char = (ell(1) < x) .* (x < ell(2));

% 1st wavenumber for char(x) = 0
k0 = 01;
% 2nd wavenumber for char(x) = 1
k1 = 04; 

% signal for char(x) = 0
signal0 = real(exp(1i*k0*x));
% signal for char(x) = 1
signal1 = real(exp(1i*k1*x));

% define the final signal and the reference
Img = signal0 .* (1 - char) + signal1 .* char;
ref =      k0 .* (1 - char) +      k1 .* char;

% get insnstantaneous wavenumber
IW_unfiltered = GetInstantaneousWavenumber(Img, dx);
% mean filter for insnstantaneous wavenumber
IW_median     = GetMedianFilter(IW_unfiltered, N_FILTER);
% get insnstantaneous wavenumber
[local_k, k]  = GetMaxLocalWavenumber_1d(Img, dx, N_WINDOW);

% compute the errors
error_IW_unfiltered = norm(IW_unfiltered - ref)/norm(ref);
error_IW_median     = norm(IW_median     - ref)/norm(ref);
error_local_k       = norm(local_k       - ref)/norm(ref);

% print errors
disp('unfiltered insnstantaneous wavenumber')
fprintf('error: %5.2e \n', error_IW_unfiltered)
disp('median filtered insnstantaneous wavenumber')
fprintf('error: %5.2e \n', error_IW_median)
disp('local wavenumber')
fprintf('error: %5.2e \n', error_local_k)
%% Plot and evaluation
figure(1)
    subplot(2, 2, 1)
        plot(x, Img, '.k-')
        ylabel('')
        xlabel('x [m]')
        title('signal')
    subplot(2, 2, 2)
        hold on
        plot(x, ref, '.g-')
        plot(x, local_k, '.k-')
        legend('reference','approximation')
        ylabel('[1/m]')
        xlabel('x [m]')
        title('local wavenumber')
    subplot(2, 2, 3)
        hold on
        plot(x, ref, '.g-')
        plot(x, IW_unfiltered, '.k-')
        legend('reference','approximation')
        ylabel('[1/m]')
        xlabel('x [m]')
        title({'unfiltered','insnstantaneous wavenumber'})
    subplot(2, 2, 4)
        hold on
        plot(x, ref, '.g-')
        plot(x, IW_median, '.k-')
        legend('reference','approximation')
        title({'median filtered','insnstantaneous wavenumber'})
        ylabel('[1/m]')
        xlabel('x [m]')