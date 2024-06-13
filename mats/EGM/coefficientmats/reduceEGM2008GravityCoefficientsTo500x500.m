% 
% Cut back on filesize for the EGM 2008 gravity model from it's full
% 2190x2190 representation (115.4 MB) to 500x500 (3.1 MB).
% 
% Matt Werner - June 12, 2024
% 

clear

% Cut the tide-free coefficients back from 2190x2190 to 500x500
load('EGM2008GravityCoeff_TideFree.mat')
% THE FOLLOWING LINE...
% %% find(EGM2008GravityCoeff_TideFree(:,1) == 500 & EGM2008GravityCoeff_TideFree(:,2) == 500)
% RESULTS IN THE INDEX: 125748
% The 125748+1 index corresponds to degree 501, order 0
EGM2008GravityCoeff_TideFree(125748+1:end,:) = [];
save('EGM2008GravityCoeff_TideFree_500x500.mat', 'EGM2008GravityCoeff_TideFree');

clear

% Cut the zero-tide coefficients back from 2190x2190 to 500x500
load('EGM2008GravityCoeff_ZeroTide.mat')
% The index to find the degree=order=500 element is the same, since both
% coefficients are exactly the same except for the C20 coefficient.
EGM2008GravityCoeff_ZeroTide(125748+1:end,:) = [];
save('EGM2008GravityCoeff_ZeroTide_500x500.mat', 'EGM2008GravityCoeff_ZeroTide');
