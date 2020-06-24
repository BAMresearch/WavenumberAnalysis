function k = GetShiftedWavenumberByDx(dx, N_X)
% GetShiftedWavenumberByDx
%   This function defines a cell with the wavenumber-vectors.
%
%   Input:
%       dx:  vector with the step size in all directions
%       N_X: number of points in the grid in all directions
%   Output:
%       k: cell with the wavenumber-vectors
%
%   See also GetRadialFilter

% author: Jannis.Bulling@bam.de

%     This function defines a cell with the wavenumber-vectors.
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
    % number of directions
    N_DIR = numel(dx);
    % initialize cell
    k = cell(N_DIR, 1);
    for i_dir = 1:N_DIR % loop for all directions
        % step size for the wavenumber 
        dk       = ( 2*pi ) / ( dx(i_dir)*N_X(i_dir) );
        % number of wavenumbers
        NH       = floor(N_X(i_dir)/2);
        % handle odd and even case
        oe       = 1 - mod(N_X(i_dir), 2);
        % define the wavenumber-vector
        k{i_dir} = (-NH:NH-oe)*dk;
    end
end

