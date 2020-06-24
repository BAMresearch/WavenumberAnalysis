function w = GetTukeyWindow(r, p)
% GetTukeyWindow
%   This function computes a continuous window w(r) define by Tukey with
%   four parameters saved in p.
%   p(1) is the starting point and p(2) is the end point of the 1st jump.
%   p(3) is the starting point and p(4) is the end point of the 2nd jump.
%
%     p(2)  p(3)
%       _____
%   ___/     \___
%     p(1)   p(4)
%
%   See also GetTukeyJump

%     This function computes a continuous window.
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
        p = [0.2 0.3 0.6 0.8];
    end
    % define the window by two jumps
    w = GetTukeyJump(r, p(1:2)).*(1 - GetTukeyJump(r, p(3:4)));
    
    w = reshape(w(:), size(r));
end

