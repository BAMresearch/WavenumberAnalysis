function [local_k, k_window] = GetMaxLocalWavenumber_PointwiseWindow(Img, dx, N_BIN, N_FFT)
% GetMaxLocalWavenumber_PointwiseWindow
%   This function computes the maximal local wavenumber k(x) or k(x,y) 
%   using a point-wise moving window.
%   The signal must be given for a equidistant grid with step-size dx. 
%
%   Input:
%       1d case:
%             Img: vector of the wavefield in space
%              dx: scalar value
%           N_BIN: half the window width in bins
%           N_FFT: number for the padded fft in the window
%       2d case:
%             Img: matrix of the wavefield in space
%              dx: vector of length 2 for both x- and y-directions
%           N_BIN: half the window width in bins
%           N_FFT: number for the padded fft in the window
%   Output:
%        local_k: vector or matrix with the local wavenumber in space.
%       k_window: cell with the wavenumbers of the window
%
%   See also GetThicknessMap

% author: Jannis.Bulling@bam.de  

%     This function computes the maximal local wavenumber using a 
%     point-wise moving window.
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
    %% default values

    if    ( isvector(Img) )
        % call 1d case
        [local_k, k_window] = GetMaxLocalWavenumber_1d(Img, dx, N_BIN, N_FFT);
    elseif( ismatrix(Img) )
        %% 2d case
        
        if(nargin < 4)
            N_FFT = max(size(Img));
        end
    
        [N_ROW, N_COL] = size(Img);           % get size
        local_k        = zeros(N_ROW, N_COL); % output
        
        % wavenumbers for the window
        k_window    = GetShiftedWavenumberByDx(dx, [N_FFT, N_FFT]);
        % shift the wavenumbers
        k_window{1} = abs(ifftshift(k_window{1}));
        k_window{2} = abs(ifftshift(k_window{2}));
                  
        for i_col = 1:N_COL     % loop over columns
            % get window index
            idx_col = i_col-N_BIN:i_col+N_BIN;
            % correct to the periodic index
            idx_col = mod(idx_col - 1, N_COL) + 1;
            for i_row = 1:N_ROW % loop over rows
                % get window index
                idx_row = i_row-N_BIN:i_row+N_BIN;
                % correct to the periodic index
                idx_row = mod(idx_row - 1, N_ROW) + 1;
                
                % amplitude values in the window
                Amp_win = abs( fft2( Img(idx_row,idx_col), N_FFT, N_FFT) );
                
                % get index of the max value
                [~,    ij_max] = max( Amp_win(:) );
                % get subindex
                [i_max, j_max] = ind2sub([N_FFT, N_FFT], ij_max);
                
                % compute the length of the wavenumber vector
                local_k(i_row, i_col) = sqrt(k_window{1}(i_max)^2 + k_window{2}(j_max)^2);
            end
        end
    else
        warning('IMGin has to be a vector or matrix.')
        local_k  = [];
        k_window = [];
    end
end

