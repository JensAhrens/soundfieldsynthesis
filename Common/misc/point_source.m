function [ S ] = point_source( x, y, x_s, y_s, z_s, k )
%function [ S ] = point_source( x, y, x_s, y_s, z_s, k )
%
% calculates the sound field at the grid points given by x and y with z = 0
% of a monochromatic point source located at [ x_s y_s z_s ];
% k = 2*pi*f/c
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
% (c) 2012 by Jens Ahrens                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = sqrt( ( x - x_s ).^2 + ( y - y_s ).^2 + ( z_s ).^2 );

S = 1 ./ r .* exp( -1i .* k .* r );