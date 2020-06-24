clear all
close all
%%
%{
    This script gives a simple example, for a thickness map.

    author: Jannis.Bulling@bam.de

    This script gives a simple example, for a thickness map.
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
%% Load data
example = load('data.mat');
%% Define parameters
% give frequency
frq  = example.simulation.frequency.value;
% spectrum
Spec = example.simulation.spectrum.value;
% step size
dx   = example.simulation.dx.value;
% number of points
N_POINT = example.simulation.N_POINT.value;
% grid
x = dx(1)*(0:N_POINT(1)-1);
y = dx(2)*(0:N_POINT(2)-1);

% corners of the defect
corner = example.defect.position.value;

PlotDefect = @() plot3(corner(1,[1:4,1]),...
                       corner(2,[1:4,1]),...
                       1e3*ones(1,5),...
                       'w-', 'LineWidth', 1.5);

% relation between wavenumber and thickness for this frequency
rk = example.relation.wavenumber.value;
rt = example.relation.thickness.value;

% tolerance for the auto filter
tol = 5e-2;

% half the window size for the median filter
N_FILTER = 6;
% the window size for the local wavenumber map
N_WINDOW = 2^3;
% the number for the padded fft inside the window of the local wavenumber
% approach
N_FFT_WIN = 2^8;
%%
k             = GetShiftedWavenumberByDx(dx, size(Spec));
Filter        = GetRadialFilter(k, [300 900], 180 + [-60 60]);
Spec_filtered = Spec.*Filter;
%% Plot
close all

k_lim = 1000;

figure
    surf(k{1}, k{2}, abs(Spec).', 'LineStyle', 'none'); 
    view(2);
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('Displacement (mm)');
    xlim(k_lim*[-1,+1])
    ylim(k_lim*[-1,+1])
    title('unfiltered spectrum')
figure
    surf(k{1}, k{2}, Filter.', 'LineStyle', 'none'); 
    view(2);
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('Displacement (mm)');
    xlim(k_lim*[-1,+1])
    ylim(k_lim*[-1,+1])
    title('filter')
figure
    surf(k{1}, k{2}, abs(Spec_filtered).', 'LineStyle', 'none'); 
    view(2);
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('Displacement (mm)');
    xlim(k_lim*[-1,+1])
    ylim(k_lim*[-1,+1])
    title('filtered spectrum')
%% IFFT of one frequency slice with and without filter
% shift back and ifft
Img = ifft2( ifftshift(Spec_filtered) );
% at this point, the signal is still complex, because the time direction is
% not transferred back
Img = real(Img);
% cut to the original size
Img = Img(1:N_POINT(1), 1:N_POINT(2));
%%
tic
IW        = GetInstantaneousWavenumber(Img, dx);
IW_Median = GetMedianFilter(IW, N_FILTER);
toc
tic
[local_k_pw, wink_pw] = GetMaxLocalWavenumber_PointwiseWindow(Img, dx, N_WINDOW, N_FFT_WIN);
toc
PW_Median = GetMedianFilter(local_k_pw, N_FILTER);
tic
[local_k_lw, wink_lw] = GetMaxLocalWavenumber_LinewiseWindow(Img, dx, N_WINDOW, N_FFT_WIN);
toc
LW_Median = GetMedianFilter(local_k_lw, N_FILTER);

[TM_IW, CmTicks_IW] = GetThicknessMap(IW_Median, rk, rt);
[TM_pw, CmTicks_pw] = GetThicknessMap(PW_Median, rk, rt);
[TM_lw, CmTicks_lw] = GetThicknessMap(LW_Median, rk, rt);
%%
[k_vec,II]  = sort(IW_Median(:));
t_vec       = TM_IW(II);
figure
    hold on
    plot(k_vec, t_vec, 'r')
    set(gca, 'Box','on', 'TickDir','out', 'FontName', 'arial', 'FontSize', 16)
    stem(rk, rt, '.k')
    xlabel('Wavenumber (rad/m)')
    ylabel('Thickness (mm)')
    title({'relation between','wavenumber and thickness'})
%%

Img_norm = (2*(Img - min(Img(:))) / (max(Img(:)) - min(Img(:)))) - 1;  

figure
    hold on
    surf(x, y, Img_norm.', 'LineStyle', 'none'); 
    set(gca, 'Box','on', 'TickDir','out', 'FontName', 'arial', 'FontSize', 16)
    view(2);
    clrbar = colorbar;
    clrbar.Label.String = 'Out-of-plane displacement (arb. u.)';
    clrbar.Label.FontSize = 16;
    PlotDefect();
    xlabel('x (m)');
    ylabel('y (m)');
    xlim([0 x(end)])
    ylim([0 y(end)])
    %savefig('Output\wavefield')
%%
figure
    hold on
    surf(x, y, abs(IW_Median).', 'LineStyle', 'none');
    set(gca, 'Box','on', 'TickDir','out', 'FontName', 'arial', 'FontSize', 16, 'CLim', [500 900])
    set(gca, 'Xtick', [0 .1 .2 .3])
    set(gca, 'Ytick', [0 .1 .2 .3])
    clrbar = colorbar;
    clrbar.Label.String = 'Wavenumber (rad/m)';
    clrbar.Label.FontSize = 16;
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('IW (1/m)')
    PlotDefect();
    view(2);
    title({'Instantaneous', 'wavenumber'})
    savefig('Output\WM_iw')
figure
    hold on
    surf(x, y, abs(PW_Median).', 'LineStyle', 'none');
    set(gca, 'Box','on', 'TickDir','out', 'FontName', 'arial', 'FontSize', 16, 'CLim', [500 900])
    set(gca, 'Xtick', [0 .1 .2 .3])
    set(gca, 'Ytick', [0 .1 .2 .3])
    clrbar = colorbar;
    clrbar.Label.String = 'Wavenumber (rad/m)';
    clrbar.Label.FontSize = 16;
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('IW (1/m)')
    PlotDefect();
    view(2);
    title({'Local wavenumber','with a point-wise moving window'})
    savefig('Output\WM_pw')
figure
    hold on
    surf(x, y, abs(LW_Median).', 'LineStyle', 'none');
    set(gca, 'Box','on', 'TickDir','out', 'FontName', 'arial', 'FontSize', 16, 'CLim', [500 900])
    set(gca, 'Xtick', [0 .1 .2 .3])
    set(gca, 'Ytick', [0 .1 .2 .3])
    clrbar = colorbar;
    clrbar.Label.String = 'Wavenumber (rad/m)';
    clrbar.Label.FontSize = 16;
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('IW (1/m)')
    PlotDefect();
    view(2);
    title({'Local wavenumber','with a line-wise moving window'})
    savefig('Output\WM_lw')
%%
figure
    hold on
    surf(x, y, TM_IW.', 'LineStyle', 'none'); view(2);
    set(gca, 'Box','on', 'TickDir','out', 'FontName', 'arial', 'FontSize', 16)
    set(gca, 'Xtick', [0 .1 .2 .3])
    set(gca, 'Ytick', [0 .1 .2 .3])
    clrbar = colorbar;
    clrbar.Label.String = 'Effective Thickness (mm)';
    clrbar.Label.FontSize = 14;
    clrbar.Ticks = [1 1.5 2];
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('Effective Thickness (mm)');
    PlotDefect();
    savefig('Output\TM_iw')
figure
    hold on
    surf(x, y, TM_pw.', 'LineStyle', 'none'); view(2);
    set(gca, 'Box','on', 'TickDir','out', 'FontName', 'arial', 'FontSize', 16)
    set(gca, 'Xtick', [0 .1 .2 .3])
    set(gca, 'Ytick', [0 .1 .2 .3])
    clrbar = colorbar;
    clrbar.Label.String = 'Effective Thickness (mm)';
    clrbar.Label.FontSize = 14;
    clrbar.Ticks = [1 1.5 2];
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('Effective thickness (mm)');
    PlotDefect();
    savefig('Output\TM_pw')
figure
    hold on
    surf(x, y, TM_lw.', 'LineStyle', 'none'); view(2);
    set(gca, 'Box','on', 'TickDir','out', 'FontName', 'arial', 'FontSize', 16)
    set(gca, 'Xtick', [0 .1 .2 .3])
    set(gca, 'Ytick', [0 .1 .2 .3])
    clrbar = colorbar;
    clrbar.Label.String = 'Effective Thickness (mm)';
    clrbar.Label.FontSize = 14;
    clrbar.Ticks = [1 1.5 2];
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('Effective thickness (mm)');
    PlotDefect();
    savefig('Output\TM_lw')