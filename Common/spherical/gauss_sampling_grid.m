function [ alpha, beta, weights ] = gauss_sampling_grid( L )
%GAUSS_SAMPLING_GRID Creates a Gaussian sampling grid with 2*L points 
% along the azimuth and L points along the colatitude; all 2*L^2 data
% points are returned in row vectors; L should be N+1 to avoid aliasing
% requires Chebfun
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This work is supplementary material for the book                        %
%                                                                         %
% Jens Ahrens, Analytic Methods of Sound Field Synthesis, Springer-Verlag %
% Berlin Heidelberg, 2012, http://dx.doi.org/10.1007/978-3-642-25743-8    %
%                                                                         %
% It has been downloaded from http://soundfieldsynthesis.org and is       %
% licensed under a Creative Commons Attribution-NonCommercial-ShareAlike  %
% 3.0 Unported License. Please cite the book appropriately if you use     %
% these materials in your own work.                                       %
%                                                                         %
% (c) 2015 by Jens Ahrens                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Alpha             = linspace( 0, 2*pi, 2*L+1 );
Alpha             = Alpha( 1:end-1 );
[ Beta, Weights ] = legpts( L, [ 0 1 ] );
Beta              = acos( Beta );

[ alpha, beta ] = meshgrid( Alpha, Beta );
[ ~ , weights ] = meshgrid( Alpha, Weights );

alpha   = alpha(:).';
beta    = beta(:).';
weights = weights(:).';

end

