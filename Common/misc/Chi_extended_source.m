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

function [ out ] = Chi_extended_source( m, eta )
% evaluates Eq.(E.37)

out = 0;

for l = 0 : 2*eta - 1
    
    alpha_2 = (l+1) * pi / eta;
    alpha_1 =  l    * pi / eta;
        
    if ( m == 0 )        
        integral = alpha_2 - alpha_1;
    else
        integral = 1i / m * ( exp( -1i * m * alpha_2 ) - exp( -1i * m * alpha_1 ) );
    end
    
    out = out + (-1)^l * integral;

end

end