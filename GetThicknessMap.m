function TM = GetThicknessMap(WM, k, d)
% GetThicknessMap
%   This function defines a thickness map for a wavenumber field.
%   Input:
%      WM: wavenumber field
%       k: ticks for the wavenumber
%       d: corresponding thickness for a wavenumber in k
%   Output:
%      TM: thickness map
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
    
    if    (numel(d) > 2)
        % quadratic polynomial fit for the depth
        p = polyfit(k, d, 2);
    elseif(numel(d) == 2)
        % linear polynomial fit for the depth
        p = polyfit(k, d, 1);
    else
        error('d is to small!')
    end
    
    
    % interpolate all wavenumbers smaller than the minimum thickness
    % to the minimum thickness
    k_min = min(k);
    % interpolate all wavenumbers larger than the maximum thickness 
    % to the maximum thickness
    k_max = max(k);
    %% interpolation of the wave numbers to the thickness
    TM = polyval(p,  k_min).*(WM(:) <= k_min)...
       + polyval(p, WM(:) ).*(k_min < WM(:) & WM(:) < k_max)...
       + polyval(p,  k_max).*(k_max <= WM(:));
    TM = reshape(TM, size(WM) );
end

