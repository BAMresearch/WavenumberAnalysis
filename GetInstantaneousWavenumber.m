function IW = GetInstantaneousWavenumber(Img, dx)
% GetInstantaneousWavenumber
%   This function computes the instantaneous wavenumber k(x) or k(x,y). 
%   The instantaneous wavenumber k is given for a equidistant grid with 
%   step-size dx. 
%
%   Input:
%       1d case:
%          Img: vector with the wavefield in space
%           dx: scalar value
%       2d case:
%          Img: matrix the wavefield in space
%           dx: vector of length 2
%   Output:
%       IW: vector or matrix with the instantaneous wavenumber in space.
%
%   Inspired by the following work:
%       Mesnil, Olivier, Cara AC Leckey, and Massimo Ruzzene. 
%       "Instantaneous and local wavenumber estimations for damage 
%       quantification in composites." 
%       Structural Health Monitoring 14.3 (2015): 193-204.
%
%   See also GetHilbert, GetThicknessMap

% author: Jannis.Bulling@bam.de

%     This function computes the instantaneous wavenumber.
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
    if(~ismatrix(Img))
        error('ndims ~= 2')
    end
    
    if( isvector(Img) )
        % hilbert-transformation
        h = GetHilbert(Img);
        % unwrap for a continuous phase
        angX = unwrap( angle(h) );
        
        % derivative of the phase
        IW = abs( gradient(angX, dx(1)) );
    else
        % hilbert-transformation in x-direction
        H = GetHilbert( Img );
        % unwrap for a continuous phase in x-direction
        angX = unwrap( angle(H) );
        
        % hilbert-transformation in y-direction
        H = GetHilbert( Img.' );
        % unwrap for a continuous phase in y-direction
        angY = unwrap( angle(H) ).';
        
        % derivative of the phase
        [kX_x, ~] = gradient(angX.', dx(1), dx(2));
        [~, kY_y] = gradient(angY.', dx(1), dx(2));
        
        % compute the length of the wavenumber vector
        IW = sqrt( kX_x.^2 + kY_y.^2 ).';
    end
end

