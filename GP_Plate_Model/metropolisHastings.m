function [hypsOpt,lik,accRatio] = metropolisHastings(hyps, nsamples, logpdf, proprnd, logproppdf, burnin, thin, verbose)
% Generate samples from a distribution using a Markov Chain
% Implemented with the Metropolis-Hastings algorithm

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

    warning off
    if ~exist('verbose', 'var'); verbose = true; end
    if isa(logproppdf, 'function_handle'); isSymmetric = false; else; isSymmetric = true; end
    if ~exist('burnin', 'var'); burnin = 1; end
    if ~exist('thin', 'var'); thin = 1; end
    
    if verbose; fprintf(1,'Running Metropolis-Hastings: %3d%%\n', 0); end
    
    % Make sure params are a column vector
    hyps = hyps(:);
    
    % Allocate outputs
    N = burnin + thin*nsamples;
    hypsOpt = zeros(length(hyps), N);
    hypsOpt(:,1) = hyps;  
    accRatio = 0;
    lik = zeros(N,1);
    
    % Allocate random acceptance numbers
    alpha = log(rand(N,1));
    
    % MH loop
    for i = 2:N
        % Get current params and sample proposals
        thisX = hypsOpt(:,i-1);
        propX= proprnd(thisX); 
        
        % Always assume non-simmetric proposal dist
        if isSymmetric
            q1 = 0; q2 = 0;
        else
            q1 = logproppdf(thisX,propX);
            q2 = logproppdf(propX,thisX);
        end
        
        % MH ratio
        newPost = logpdf(propX);
        oldPost = logpdf(thisX);
        
        % MH ratio in log-form
        rho = (q1+newPost)-(q2+oldPost); 
        
        % Accept or reject the proposal
        acc = alpha(i) <= rho;
        if acc
            hypsOpt(:,i) = propX;
            lik(i) = newPost;
            accRatio = accRatio + 1;
        else
            hypsOpt(:,i) = thisX;
            lik(i) = oldPost;
        end
    
        % Verbose
        if verbose; fprintf(1,'\b\b\b\b\b%3.0f%%\n', 100*(i/N)); end
    end
    
    % Acceptance rate
    accRatio = accRatio/N; 

    % Burnin and thin
    hypsOpt = hypsOpt(:, burnin:end)';
    hypsOpt = hypsOpt(1:thin:end,:);
    
    % Restore warning
    warning on
end