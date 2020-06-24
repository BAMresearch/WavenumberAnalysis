function Filter = GetRadialFilter(k, r, theta)
% GetRadialFilter
%   This function computes a radial filter for the wavenumbers given in the 
%   cell k. The parameters r and theta define the radius and angle in 
%   degrees of the filter, respectively.
%
%   See also GetWavenumberShift, GetTukeyWindow

%     This function computes a radial filter.
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

    % default value for r
    if(nargin < 2 || isempty(r))
        km = sqrt(k{1}(end)^2 + k{2}(end)^2);
        r = [.25, .5]*km;
    end
    % default value for theta
    if(nargin < 3 || isempty(theta))
        theta = [-1 -1 361 361];
    end
    
    % sorting
    r     = sort(r);
    theta = sort(theta);
    
    % default behaviour
    if(numel(r)==1)
        r(4) = r(1);
        r(3) = .9*r(4);
        r(2) = -1;
        r(1) = -1;
    elseif(numel(r)==2)
        r([1,4]) = r;
        r(2) = r(1) + .1*(r(4) - r(1));
        r(3) = r(1) + .9*(r(4) - r(1));
    end
    if(numel(theta)==2)
        theta([1,4]) = theta;
        theta(2) = theta(1) + 10;
        theta(3) = theta(4) - 10;
    end
    
    % degree to radian
    theta = theta*pi/180;
    
    % wavenumber grid
    [K1, K2] = ndgrid(k{1},k{2});
    % polar grid
    [TH, RR] = cart2pol(K1, K2);
    
    % turn the angle such that theta(1) becomes 0
    TH = TH - theta(1) + pi;
    
    % for periodic behaviour
    TH(TH > + pi) = TH(TH > + pi) - 2 * pi;
    TH(TH < - pi) = TH(TH < - pi) + 2 * pi;

    % turn the angle back
    theta = theta - theta(1) - pi;
    
    % get the filter mask
    Filter = GetTukeyWindow(TH, theta) .* GetTukeyWindow(RR, r);
end

