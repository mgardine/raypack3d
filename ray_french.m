function h = ray_french()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function ray_french
%
% This function generates the 'french' colormap with values going from blue
% to red through white (standard tomography colormap)
%
% Usage: 
%   h = ray_french()
%       OR
%   colormap(french) in plotting routine
%
% Output:
%   h:      A m-by-3 matrix of real numbers between 0 and 1. Each row is an
%           RGB vector that defines one color.
%
% Author:
% Matt Gardine
% June 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = size(get(gcf,'colormap'),1); 

n3 = fix(m/2);

r = [ones(n3,1); (sqrt(1-((1:n3)/n3).^.7))'];
g = [flipud( (sqrt(1-((1:n3)/n3).^.7))'); (sqrt(1-((1:n3)/n3).^.7))'];
b = [flipud((sqrt(1-((1:n3)/n3).^.7))'); ones(n3,1)];

h = [r g b];

if size(h,1) < m
    h(ceil(m/2)+1:m,:) = h(ceil(m/2):end,:);
    h(ceil(m/2),:) = 1;
end