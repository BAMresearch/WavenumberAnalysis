function ImgOut = GetMedianFilter(ImgIn, N_BIN)
% GetMedianFilter
%     Input: 
%            ImgIn: one color 1D or 2D image as vector or matrix,
%                   respectively.
%            N_BIN: an integer with half the filter window width 
%                   for example: N_BIN = 0 will give the original image
%                                N_BIN = 1 will use the three points for 1D 
%                                          and nine points for 2D
% 
%     Output: 
%             ImgOut: a median filtered image with the same size as ImgIn   
% 
%   See also median

% author: Jannis.Bulling@bam.de    
    
%     This function filters an image using a median filter. 
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
    
    if    ( isvector(ImgIn) )
        %% 1d case
        % initialized variables
        N_ROW  = numel(ImgIn);     % get size of grayscale image
        ImgOut = zeros(N_ROW, 1);  % output / median filtered image
        
        for i_row = 1 : N_ROW % loop over rows
            % get window index
            idx_row = i_row-N_BIN:i_row+N_BIN;
            % correct to the periodic index
            idx_row = mod(idx_row - 1, N_ROW) + 1;
            % get values for the window
            win = ImgIn(idx_row);
            
            % find median
            ImgOut(i_row, 1) = median( win(:) );   
        end
        % reshape so that input and output have the same shape
        ImgOut = reshape(ImgOut, size(ImgIn));
    elseif( ismatrix(ImgIn) )
        %% 2d case
        % initialized variables
        [N_ROW, N_COL] = size(ImgIn); % get size of grayscale image
        ImgOut = zeros(N_ROW, N_COL); % output / median filtered image
        
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
                
                % get values for the window
                win = ImgIn(idx_row, idx_col);
                
                % find median
                ImgOut(i_row,i_col) = median( win(:) );
            end
        end
    else
        warning('ImgIn has to be a vector or matrix.')
        ImgOut = []; 
    end
end