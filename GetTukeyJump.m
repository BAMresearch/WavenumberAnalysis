function j = GetTukeyJump(r, p)
% GetTukeyJump
%   This function computes a continuous junp j(r) define by Tukey with
%   two parameters saved in p.
%   p(1) is the starting point and p(2) is the end point of the jump.
%
%      p(2)
%       _____
%   ___/
%    p(1)
%
%   See also GetTukeyWindow

%     This function computes a continuous junp
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

    % default input
    if(nargin < 1)
        r = linspace(0, 1, 100);
    end
    if(nargin < 2)
        p = [0.5 0.5];
    end
    
    if(p(1)<p(2)) 
        % continuous junp
        j = ( ( p(1)  < r(:) ) & ( r(:) < p(2) ) )...
          .* (1 - cos( pi*( ( r(:) - p(1) ) / ( p(2) - p(1) ) ) ) )/2 ...
          + ( p(2) <= r(:) );
    else
        % discontinuous junp
        j = ( p(2) <= r(:) );
    end
  
    j = reshape(j(:), size(r));
end

