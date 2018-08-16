function h = fig_size( x, y );
% function h = fig_size( x, y );
% Creates a new figure with x and y dimensions in cm;
% Bases on a script by Sascha Spors

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set( gcf, 'Position', [ 20, 45, 44*x, 44*y ] );
x = x + 2;
y = y + 2;
set( gca, 'FontName', 'Times' );
set( gca, 'FontSize', 8 );
set( gcf, 'PaperUnits', 'centimeters' );
set( gcf, 'Papertype', 'A4' )
tmp = get( gcf,'Papersize' );
set( gcf, 'PaperPosition', [ ( tmp(1) - x ) / 2,( tmp(2) - y ) / 2, x, y ] );
set( gcf, 'Color', [ 1, 1, 1 ] );

