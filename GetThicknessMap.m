function [TM, CmTicks] = GetThicknessMap(WM, k, d)
% GetThicknessMap
%   This function defines a thickness map for a wavenumber field.
%   Input:
%      WM: wavenumber field
%       k: ticks for the wavenumber
%       d: corresponding thickness for a wavenumber in k
%   Output:
%      TM: thickness map
% CmTicks: sorted ticks for the colorbar
%
%   See also GetInstantaneousWavenumber

% author: Jannis.Bulling@bam.de

%     This function defines a thickness map for a wavenumber field.
%     Copyright (C) 2020  Jannis Bulling
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
    %% 
    % reshape for a simpler vector concatenation of arrays
    d = reshape(d, 1, []);
    k = reshape(k, 1, []);
    
    % interpolate all wavenumbers smaller than the minimum thickness
    % to the minimum thickness
    WM_min         = min( WM(:) );
    [k_min, i_min] = min(k);
    if(WM_min<k_min)
        k = [  WM_min, k];
        d = [d(i_min), d];
    end
    % interpolate all wavenumbers larger than the maximum thickness 
    % to the maximum thickness
    WM_max         = max( WM(:) );
    [k_max, i_max] = max(k);
    if(WM_max>k_max)
        k = [k, WM_max];
        d = [d, d(i_max)];
    end
    %% define ticks for the colorbar
    CmTicks = uniquetol(d);
    %% interpolation of the wave numbers to the thickness
    TM = interp1(k, d, WM(:), 'pchip', NaN);
    TM = reshape(TM, size(WM) );
end

