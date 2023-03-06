function [exp_refl] = generateexp_refl(wavelengths, x, y, wave_start, wave_end)
%GENERATEEXP_REFL This function generates the experimental reflectance
%given reference and spectrum data. 
%   This function takes inputs of spectrum, reference and "dark" data to
%   generate plots for the experimental reflectance, and the first and
%   second derivatives of these. 
% 
% 
% ******* INPUTS *********
% 
% spectrum 
% reference 
% dark 
% 
% 
% ****** OUTPUTS *******
% 
% experimental reflectance 
% 1st diff of exp refl 
% 2nd diff of exp refl 


%load data from miniCYRIL 
% wavelengths = Data.Wavelengths; 
% x = Data.Reference; 
% y = Data.Spectrum; 

% create a vector that has the wavelenghts in the first column, the
% reference data in the second column and the spectrum data in the third
% column 
% horzcat is horizontal concatenation 
vector = horzcat(wavelengths', x', y');

waves=650:915; 
interp_data = interp1(wavelengths, vector, waves); 

%then here we need to find the wavelengths which correspond to 680:850
[p] = find(round(interp_data) == wave_start, 1);
[q] = find(round(interp_data) == wave_end, 1); 


% between the wavelengths of p and q, divide the reference by the spectrum
% data. 
r = interp_data(p:q, 2);
s = interp_data(p:q, 3); 
refl = r./(s.*100000); 

% this plots the experimental reflectance
exp_refl = smooth(refl); 

%take the first and second derivatives of the experimental reflectance 
diff_1_exp_refl = diff(smooth(refl),1); 
diff_2_exp_refl = diff(smooth(diff_1_exp_refl));

end 
