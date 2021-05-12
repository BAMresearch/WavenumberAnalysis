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
sv.fname  = 'Data\Data_AlSide_0_8mm_100kHz.mat'; 
OutputDir = 'Output\ExperimentALUside\';
example   = load(sv.fname);
%% Define parameters
% give frequency
frq  = example.experiment.frequency.value;
% spectrum
Spec = example.experiment.spectrum.value;
% step size
dx   = example.experiment.dx.value;
% number of points
N_POINT = example.experiment.N_POINT.value;
% grid
x = dx(1)*(0:N_POINT(1)-1);
y = dx(2)*(0:N_POINT(2)-1);

% relation between wavenumber and thickness for this frequency
rk = example.relation.wavenumber.value([1,end]);
rt = example.relation.thickness.value([1,end]);

% half the window size for the median filter
sv.N_FILTER = 8;
% the window size for the local wavenumber map
sv.N_WINDOW = 21; % 0.8 mm case
% the number for the padded fft inside the window of the local wavenumber
% approach
sv.N_FFT_WIN = 2^8;
% wavenumber filter 
sv.FILTER.K = [300 1100];
%% All options for the surf-Command and axis-Command
TLIM = [min(rt), max(rt)];
KLIM = [300 900];

XTICK = [0 0.1 0.2]; % Ticks for the axis
YTICK = [0 0.1 0.2]; % Ticks for the axis

SurfOpt{1} = 'LineStyle';
SurfOpt{2} = 'none';

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
AxisOpt{2,07} = [x(1), x(end)];
AxisOpt{1,08} = 'ylim';
AxisOpt{2,08} = [y(1), y(end)];
%%
k             = GetShiftedWavenumberByDx(dx, size(Spec));
Filter        = GetRadialFilter(k, sv.FILTER.K);
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
    clrbar.Label.String   = 'Wavenumber (rad/m)';
    clrbar.Label.FontSize = 16;
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('IW (1/m)')
    view(2);
    title({'Instantaneous', 'wavenumber'})
    savefig([OutputDir,'WM_iw'])
%%
figure
    hold on
    surf(x, y, abs(PW_Median).', SurfOpt{:}); 
    set(gca, AxisOpt{:})
    set(gca, 'Clim', KLIM)
    clrbar = colorbar;
    clrbar.Label.String   = 'Wavenumber (rad/m)';
    clrbar.Label.FontSize = 16;
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('IW (1/m)')
    view(2);
    title({'Local wavenumber','with a point-wise moving window'})
    savefig([OutputDir,'WM_pw'])
%%
figure
    hold on
    surf(x, y, abs(LW_Median).', SurfOpt{:}); 
    set(gca, AxisOpt{:})
    set(gca, 'Clim', KLIM)
    clrbar = colorbar;
    clrbar.Label.String   = 'Wavenumber (rad/m)';
    clrbar.Label.FontSize = 16;
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('IW (1/m)')
    view(2);
    title({'Local wavenumber','with a line-wise moving window'})
    savefig([OutputDir,'WM_lw'])
%%
figure
    hold on
    surf(x, y, TM_IW.', SurfOpt{:});
    set(gca, AxisOpt{:})
    set(gca, 'Clim', TLIM)
    clrbar = colorbar;
    clrbar.Label.String   = 'Effective Thickness (mm)';
    clrbar.Label.FontSize = 14;
    clrbar.Ticks = rt;
    xlabel('x (m)');
    ylabel('y (m)');
    view(2);
    zlabel('Effective Thickness (mm)');
    savefig([OutputDir,'TM_iw'])
%%
figure
    hold on
    surf(x, y, TM_pw.', SurfOpt{:}); 
    set(gca, AxisOpt{:})
    set(gca, 'Clim', TLIM)
    clrbar = colorbar;
    clrbar.Label.String   = 'Effective Thickness (mm)';
    clrbar.Label.FontSize = 14;
    clrbar.Ticks = rt;
    xlabel('x (m)');
    ylabel('y (m)');
    view(2);
    zlabel('Effective thickness (mm)');
    savefig([OutputDir,'TM_pw'])
%%
figure
    hold on
    surf(x, y, TM_lw.', SurfOpt{:});
    set(gca, AxisOpt{:})
    set(gca, 'Clim', TLIM)
    clrbar = colorbar;
    clrbar.Label.String   = 'Effective Thickness (mm)';
    clrbar.Label.FontSize = 14;
    clrbar.Ticks = rt;
    xlabel('x (m)');
    ylabel('y (m)');
    view(2);
    zlabel('Effective thickness (mm)');
     savefig([OutputDir,'TM_lw'])