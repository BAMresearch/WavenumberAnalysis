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
sv.fname  = 'Data\Data_Simulation.mat'; 
OutputDir = 'Output\Simulation\';
example   = load(sv.fname);
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

PlotDefect = @(h) plot3(corner(1,[1:4,1]),...
                        corner(2,[1:4,1]),...
                        h*ones(1,5),...
                        'w-', 'LineWidth', 1.5);

% relation between wavenumber and thickness for this frequency
rk = example.relation.wavenumber.value;
rt = example.relation.thickness.value;

% half the window size for the median filter
sv.N_FILTER = 6;
% the window size for the local wavenumber map
sv.N_WINDOW = 10;
% the number for the padded fft inside the window of the local wavenumber
% approach
sv.N_FFT_WIN = 2^9;
% wavenumber filter 
sv.FILTER.K     = [300 900];
sv.FILTER.Angle = 180 + [-60 +60];
%% All options for the surf-Command and axis-Command
TLIM = [min(rt), max(rt)];
KLIM = [500 900];

XTICK = [0 .1 .2 .3]; % Ticks for the axis
YTICK = [0 .1 .2 .3]; % Ticks for the axis

SurfOpt{1} = 'LineStyle';
SurfOpt{2} = 'none';

PlotOpt{1} = 'LineWidth';
PlotOpt{2} = 2;

AxisOpt{1,01} = 'Box';
AxisOpt{2,01} = 'on';
AxisOpt{1,02} = 'TickDir';
AxisOpt{2,02} = 'out';
AxisOpt{1,03} = 'FontName';
AxisOpt{2,03} = 'arial';
AxisOpt{1,04} = 'FontSize';
AxisOpt{2,04} = 16;
AxisOpt{1,05} = 'Xtick';
AxisOpt{2,05} = XTICK;
AxisOpt{1,06} = 'Ytick';
AxisOpt{2,06} = YTICK;
AxisOpt{1,07} = 'xlim';
AxisOpt{2,07} = [XTICK(1), XTICK(end)];
AxisOpt{1,08} = 'ylim';
AxisOpt{2,08} = [YTICK(1), YTICK(end)];

ColIW = [.8 0  0]; StyleIW = '-';
ColPW = [0  0 .8]; StylePW = '--';
ColLW = [0 .6 .6]; StyleLW = '-.';

iy = 151;
%%
k             = GetShiftedWavenumberByDx(dx, size(Spec));
Filter        = GetRadialFilter(k, sv.FILTER.K, sv.FILTER.Angle);
Spec_filtered = Spec.*Filter;
%% Plot
k_lim = 1000;

figure
    surf(k{1}, k{2}, abs(Spec).', SurfOpt{:});  
    view(2);
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('Displacement (mm)');
    xlim(k_lim*[-1,+1])
    ylim(k_lim*[-1,+1])
    title('unfiltered spectrum')
    savefig([OutputDir,'Spec'])
figure
    surf(k{1}, k{2}, Filter.', SurfOpt{:});  
    view(2);
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('Displacement (mm)');
    xlim(k_lim*[-1,+1])
    ylim(k_lim*[-1,+1])
    title('filter')
    savefig([OutputDir,'Filter'])
figure
    surf(k{1}, k{2}, abs(Spec_filtered).', SurfOpt{:});  
    view(2);
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('Displacement (mm)');
    xlim(k_lim*[-1,+1])
    ylim(k_lim*[-1,+1])
    title('filtered spectrum')
    savefig([OutputDir,'FilteredSpec'])
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
IW = GetInstantaneousWavenumber(Img, dx);
sv.cTime_IW = toc;
tic
[local_k_pw, wink_pw] = GetMaxLocalWavenumber_PointwiseWindow(Img, dx, sv.N_WINDOW, sv.N_FFT_WIN);
sv.cTime_PW = toc;
tic
[local_k_lw, wink_lw] = GetMaxLocalWavenumber_LinewiseWindow(Img, dx, sv.N_WINDOW, sv.N_FFT_WIN);
sv.cTime_LW = toc;
%%
IW_Median = GetMedianFilter(IW, sv.N_FILTER);
PW_Median = GetMedianFilter(local_k_pw, sv.N_FILTER);
LW_Median = GetMedianFilter(local_k_lw, sv.N_FILTER);
%%
TM_IW = GetThicknessMap(IW_Median, rk, rt);
TM_pw = GetThicknessMap(PW_Median, rk, rt);
TM_lw = GetThicknessMap(LW_Median, rk, rt);
%% save meta data
save([OutputDir,'MetaData.mat'], '-struct', 'sv')
%%
[k_vec, II] = sort(IW_Median(:));
t_vec       = TM_IW(II);

figure
    hold on
    plot(k_vec, t_vec, 'r')
    set(gca, AxisOpt{:,1:4})
    stem(rk, rt, '.k')
    xlabel('Wavenumber (rad/m)')
    ylabel('Thickness (mm)')
    title({'relation between','wavenumber and thickness'})
    savefig([OutputDir,'Fit'])
%%
Img_norm = (2*(Img - min(Img(:))) / (max(Img(:)) - min(Img(:)))) - 1;  

figure
    hold on
    surf(x, y, Img_norm.', SurfOpt{:}); 
    set(gca, AxisOpt{:,1:8})
    clrbar = colorbar;
    clrbar.Label.String   = 'Out-of-plane displacement (arb. u.)';
    clrbar.Label.FontSize = 16;
    xlabel('x (m)');
    ylabel('y (m)');
    view(2)
    xlim([0 x(end)])
    ylim([0 y(end)])
    title('wavefield')
    savefig([OutputDir,'Wavefield'])
%%
figure
    hold on
    surf(x, y, abs(IW_Median).', SurfOpt{:}); 
    set(gca, AxisOpt{:})
    set(gca, 'Clim', KLIM)
    clrbar = colorbar;
    clrbar.Ticks = [500 600 700 800 900];
    clrbar.Label.String   = 'Wavenumber (rad/m)';
    clrbar.Label.FontSize = 16;
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('IW (1/m)')
    view(2);
%     title({'Instantaneous', 'wavenumber'})
    PlotDefect(1e3);
    savefig([OutputDir,'WM_iw'])
%%
figure
    hold on
    surf(x, y, abs(PW_Median).', SurfOpt{:}); 
    set(gca, AxisOpt{:})
    set(gca, 'Clim', KLIM)
    clrbar = colorbar;
    clrbar.Ticks = [500 600 700 800 900];
    clrbar.Label.String   = 'Wavenumber (rad/m)';
    clrbar.Label.FontSize = 16;
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('IW (1/m)')
    view(2);
%     title({'Local wavenumber','with a point-wise moving window'})
    PlotDefect(1e3);
    savefig([OutputDir,'WM_pw'])
%%
figure
    hold on
    surf(x, y, abs(LW_Median).', SurfOpt{:}); 
    set(gca, AxisOpt{:})
    set(gca, 'Clim', KLIM)
    clrbar = colorbar;
    clrbar.Ticks = [500 600 700 800 900];
    clrbar.Label.String   = 'Wavenumber (rad/m)';
    clrbar.Label.FontSize = 16;
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('IW (1/m)')
    view(2);
%     title({'Local wavenumber','with a line-wise moving window'})
    PlotDefect(1e3);
    savefig([OutputDir,'WM_lw'])
%%
LineIW = TM_IW(:,iy);
figure
    hold on
    surf(x, y, TM_IW.', SurfOpt{:});
    plot3(x, y(iy)*ones(size(x)), LineIW, StyleIW, PlotOpt{:}, 'Color', ColIW)
    set(gca, AxisOpt{:})
    set(gca, 'Clim', TLIM)
    clrbar = colorbar;
    clrbar.Label.String   = 'Effective Thickness (mm)';
    clrbar.Label.FontSize = 14;
    xlabel('x (m)');
    ylabel('y (m)');
    view(2);
    zlabel('Effective Thickness (mm)');
    PlotDefect(3);
    savefig([OutputDir,'TM_iw'])
%%
LinePW = TM_pw(:,iy);
figure
    hold on
    surf(x, y, TM_pw.', SurfOpt{:});
    plot3(x, y(iy)*ones(size(x)), LinePW, StylePW, PlotOpt{:}, 'Color', ColPW)
    set(gca, AxisOpt{:})
    set(gca, 'Clim', TLIM)
    clrbar = colorbar;
    clrbar.Label.String   = 'Effective Thickness (mm)';
    clrbar.Label.FontSize = 14;
    xlabel('x (m)');
    ylabel('y (m)');
    view(2);
    zlabel('Effective thickness (mm)');
    PlotDefect(3);
    savefig([OutputDir,'TM_pw'])
%%
LineLW = TM_lw(:,iy);
figure
    hold on
    surf(x, y, TM_lw.', SurfOpt{:});
    plot3(x, y(iy)*ones(size(x)), LineLW, StyleLW, PlotOpt{:}, 'Color', ColLW)
    set(gca, AxisOpt{:})
    set(gca, 'Clim', TLIM)
    clrbar = colorbar;
    clrbar.Label.String   = 'Effective Thickness (mm)';
    clrbar.Label.FontSize = 14;
    xlabel('x (m)');
    ylabel('y (m)');
    view(2);
    zlabel('Effective thickness (mm)');
    PlotDefect(3);
    savefig([OutputDir,'TM_lw'])
%%
x1 = example.defect.position.value(1,1);
x2 = example.defect.position.value(1,3);

ref = 2 - (x1 < y & y < x2);

figure
    hold on
    set(gca, AxisOpt{:,1:4})
    set(gca, 'XTick', XTICK)
    set(gca, 'YTick', [1,2])
    plot(x, ref, ':', PlotOpt{:}, 'Color', [.5 .5 .5])
    plot(x, LineIW, StyleIW, PlotOpt{:}, 'Color', ColIW)
    plot(x, LinePW, StylePW, PlotOpt{:}, 'Color', ColPW)
    plot(x, LineLW, StyleLW, PlotOpt{:}, 'Color', ColLW)
    legend({'REF', 'IW', 'LW pw','LW lw'}, 'Location', 'southeast')
    xlabel('x (m)');
    ylabel('Effective thickness (mm)');
    ylim([.9 2.1])
    savefig([OutputDir,'Cut'])