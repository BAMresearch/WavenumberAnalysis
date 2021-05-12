function H = GetHilbert(Img)
% GetHilbert
%     A simple implementation of a Hilbert transformation by FFT.
%
%   Input:
%       1d case:
%             Img: vector with the wavefield in space
%       2d case:
%             Img: matrix the wavefield in space
%   Output:
%       H: vector or matrix with the hilbert transformation for the rows

% author: Jannis.Bulling@bam.de

%     A simple implementation of a Hilbert transformation by FFT.
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
    % switch to wavenumber domain
    H_hat = fft(Img);
    N_H   = ceil(size(Img,1)/2)+1;
    
    % hilbert transformation
    if( isvector(Img) )
        H_hat(2:N_H)       = 2*H_hat(2:N_H);
        H_hat(N_H + 1:end) = 0;
    else
        H_hat(2:N_H,:)       = 2*H_hat(2:N_H,:);
        H_hat(N_H + 1:end,:) = 0;
    end
    % switch back to space domain
    H = ifft(H_hat);
end

