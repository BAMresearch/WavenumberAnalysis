function [local_k, k_window] = GetMaxLocalWavenumber_1d(Img, dx, N_BIN, N_FFT)
% GetMaxLocalWavenumber_1d
%   This function computes the maximal local wavenumber k(x) moving window.
%   The signal must be given for a equidistant grid with step-size dx. 
%
%   Input:
%         Img: vector with the wavefield in space
%          dx: scalar value
%       N_BIN: half the window width in bins
%       N_FFT: number for the padded fft in the window
%   Output:
%       local_k: vector or matrix with the local wavenumber in space.
%       k_window: cell with the wavenumbers of the window
%
%   See also GetThicknessMap

% author: Jannis.Bulling@bam.de  

%     This function computes the maximal local wavenumber.
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
    if( isvector(Img) )
        % errors
        if(nargin < 4)
            N_FFT = max(size(Img));
        end
        
        N_ROW   = numel(Img);      % get size
        local_k = zeros(N_ROW, 1); % output
        
        % wavenumbers for the window
        k_window    = GetShiftedWavenumberByDx(dx, N_FFT);
        % shift the wavenumbers
        k_window{1} = abs(ifftshift(k_window{1}));
        
        for i_row = 1:N_ROW
            % get window index
            idx_row = i_row-N_BIN:i_row+N_BIN;  
            % correct to the periodic index
            idx_row = mod(idx_row - 1, N_ROW) + 1;
            
            % amplitude values in the window
            Amp_win = abs(fft(Img(idx_row), N_FFT));

            % get index of the max value
            [~, i_max] = max(Amp_win);

            % assign the wave number that corresponds to the maximum value
            local_k(i_row) = k_window{1}(i_max);
        end
        % reshape so that input and output have the same shape
        local_k = reshape(local_k, size(Img));
    else
        error('Img is not a vector!')
    end
end

