function [local_k, k_window] = GetMaxLocalWavenumber_LinewiseWindow(Img, dx, N_BIN, N_FFT)
% GetMaxLocalWavenumber_PointwiseWindow
%   This function computes the maximal local wavenumber k(x) or k(x,y) in a
%   line-wise moving window.
%   The signal must be given for a equidistant grid with step-size dx. 
%
%   Input:
%       1d case:
%             Img: vector with the wavefield in space
%              dx: scalar value
%           N_BIN: half the window width in bins
%           N_FFT: number for the padded fft in the window
%       2d case:
%             Img: matrix the wavefield in space
%              dx: vector of length 2
%           N_BIN: half the window width in bins
%           N_FFT: number for the padded fft in the window
%   Output:
%        local_k: vector or matrix with the local wavenumber in space.
%       k_window: cell with the wavenumbers of the window
%
%   See also GetThicknessMap

%     This function computes the maximal local wavenumber using a 
%     line-wise moving window.
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
    %% default value
    if(nargin < 4)
        N_FFT = max(size(Img));
    end
    if    ( isvector(Img) )
        % call 1d case
        [local_k, k_window] = GetMaxLocalWavenumber_1d(Img, dx, N_BIN, N_FFT);
    elseif( ismatrix(Img) )
        [N_ROW, N_COL] = size(Img); % get size
        local_k1          = zeros(N_ROW, N_COL); % local wavenumber for x
        local_k2          = zeros(N_ROW, N_COL); % local wavenumber for y
        
        k_window    = GetShiftedWavenumberByDx(dx, [N_FFT, N_FFT]);
        % shift the wavenumbers
        k_window{1} = abs(ifftshift(k_window{1}));
        k_window{2} = abs(ifftshift(k_window{2}));
                  
        for i_row = 1:N_ROW
            % get window index
            idx_row = i_row-N_BIN:i_row+N_BIN;
            % correct to the periodic index
            idx_row = mod(idx_row - 1, N_ROW) + 1;

            % amplitude values in the window
            Amp_win = abs( fft(Img(idx_row, :), N_FFT, 1) );

            % get index of the max value
            [~, i_max] = max(Amp_win, [], 1);

            % assign the wave number that corresponds to the maximum value
            local_k1(i_row, :) = k_window{1}(i_max);
        end
%         clear i_row idx_row i_max
        for i_col = 1:N_COL % loop over columns
            % get window index
            idx_col = i_col-N_BIN:i_col+N_BIN;
            % correct to the periodic index
            idx_col = mod(idx_col - 1, N_COL) + 1;

            % amplitude values in the window
            Amp_win = abs( fft(Img(:, idx_col), N_FFT, 2) );

            % get index of the max value
            [~, j_max] = max(Amp_win, [], 2);

            % assign the wave number that corresponds to the maximum value
            local_k2(:,i_col) = k_window{2}(j_max);
        end
        % compute the length of the wavenumber vector
        local_k = sqrt(local_k1.^2 + local_k2.^2);
    else
        warning('IMGin has to be a vector or matrix.')
        local_k  = [];
        k_window = [];
    end
end

