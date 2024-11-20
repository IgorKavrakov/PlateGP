function LP = GP_HypsPrior(hyps, flagUniformD)
% Prior model for GP parameters

% By Gledson Rodrigo Tondo
%%%%%%%%% COPYRIGHT NOTICE %%%%%%%%% 
%  This file is part of PlateGP.
%  PlateGP is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  PlateGP is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with PlateGP.  If not, see <https://www.gnu.org/licenses/>.

% Copyright (c) Igor Kavrakov, Gledson Rodrigo Tondo, Guido Morgenthal 2024
    
    % Start log-prior vector
    LP = (-Inf)*ones(size(hyps));    
    
    % Convert hyps back from log form
    thisHyps = exp(hyps);
    
    % Number of fixed hyps, everything after this are noise (D, A, lx, ly)
    numFixedParams = 4;
    
    % 01 - Prior on A: uniform distribution
    upperLim = 1e20;
    lowerLim = 1e-3;
    if thisHyps(1) <= upperLim && thisHyps(1) >= lowerLim
        LP(1) = log(1/(upperLim-lowerLim));
    end
   
    % 02 - Prior on LX: uniform distribution
    upperLim = 50;
    lowerLim = 1e-1;
    if thisHyps(2) <= upperLim && thisHyps(2) >= lowerLim
        LP(2) = log(1/(upperLim-lowerLim));
    end
    
    % 03 - Prior on LY: uniform distribution
    upperLim = 50;
    lowerLim = 1e-1;
    if thisHyps(3) <= upperLim && thisHyps(3) >= lowerLim
        LP(3) = log(1/(upperLim-lowerLim));
    end
    
    % 04 - Prior on D: uniform/normal distribution
    if nargin == 1; flagUniformD = true; end
    if flagUniformD
        upperLim = 30;
        lowerLim = 1;
        if thisHyps(4) <= upperLim && thisHyps(4) >= lowerLim
            LP(4) = log(1/(upperLim-lowerLim));
        end
    else
        mu = 20; stdev = 1.5;
        LP(4) = - log(stdev^2) - 0.5*log(2*pi) - (thisHyps(4)-mu)^2/(2*stdev^2);
    end

    % 05 - Priors on noise
    upperLim = 1e-1;
    lowerLim = 0;
    for i = numFixedParams+1:length(thisHyps)
        if thisHyps(i) <= upperLim && thisHyps(i) >= lowerLim
            LP(i) = log(1/(upperLim-lowerLim));
        end
    end    

    % Final value of log prior
    LP = sum(LP);
end