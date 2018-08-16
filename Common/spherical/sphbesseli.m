function [ out ] = sphbesseli( nu, z )
%function [ out ] = sphbesseli( nu, z )
%
% SPHBESSELI modified spherical bessel function of first type of order nu 
% and argument z

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
% (c) 2012 by Jens Ahrens                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out = zeros( size( z ) );

% avoid division by "0"
if ( nu==0 )
    out( find( z==0 ) ) = 1;
elseif ( nu~=0 )
    out( find( z==0 ) ) = 0;
end

% finally, evaluate for z~=0
out( find( z~=0 ) ) = sqrt( pi ./ ( 2 .* z( find( z~=0 ) ) ) ) .* ...
                                    besseli( nu + 0.5, z( find( z~=0 ) ) );
