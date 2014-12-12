function varargout = getOctantSubs(szx,mode)
%Get the indices of the elements of each octant/quadrant for a dataset of
%given size. Called by JJL2POLY()
% A = getOctantSubs(szx,mode)
%
%   sze - size of the dataset.  Can be 1, 2, or 3D.
%   mode - [('reg') 'minmax'] In REG mode, the function outputs a cell
%       array with entries for each element that is contained in each
%       octant/quadrant.  In MINMAX mode, the ouputs are the entries in the
%       cell array are the minimum and maximum bounding limits for each
%       octant/quadrant
if nargin < 2 || isempty(mode)
    mode = 'reg';
end
if length(szx)==2
    szx(szx<=1) = [];
end

switch length(szx(:))
    case 1
        Oct = cell(2,1);
        if strcmpi(mode,'reg')
            Oct{1} = (1:szx(1)/2)';
            Oct{2} = (szx(1)/2+1:szx(1))';
        elseif strcmpi(mode,'minmax')
            Oct{1} = [1,szx(1)/2];
            Oct{2} = [szx(1)/2+1,szx(1)];
        end
    case 2
        Oct = cell(2,2);
        if strcmpi(mode,'reg')
            Oct{1} = {1:szx(1)/2,       1:szx(2)/2};
            Oct{2} = {szx(1)/2+1:szx(1),1:szx(2)/2};
            Oct{3} = {1:szx(1)/2,       (szx(2)/2+1):szx(2)};
            Oct{4} = {szx(1)/2+1:szx(1),(szx(2)/2+1):szx(2)};
        elseif strcmpi(mode,'minmax')
            Oct{1} = {[1,szx(1)/2],     [1,szx(2)/2]};
            Oct{2} = {[szx(1)/2+1,szx(1)],[1,szx(2)/2]};
            Oct{3} = {[1,szx(1)/2],       [(szx(2)/2+1),szx(2)]};
            Oct{4} = {[szx(1)/2+1,szx(1)],[(szx(2)/2+1),szx(2)]};
        end
    case 3
        Oct = cell(2,2,2);
        if strcmpi(mode,'reg')
            Oct{1} = {1:szx(1)/2,       1:szx(2)/2,         1:szx(3)/2};
            Oct{2} = {szx(1)/2+1:szx(1),1:szx(2)/2,         1:szx(3)/2};
            Oct{3} = {1:szx(1)/2,       (szx(2)/2+1):szx(2),1:szx(3)/2};
            Oct{4} = {szx(1)/2+1:szx(1),(szx(2)/2+1):szx(2),1:szx(3)/2};
            Oct{5} = {1:szx(1)/2,       1:szx(2)/2,         (szx(3)/2+1):szx(3)};
            Oct{6} = {szx(1)/2+1:szx(1),1:szx(2)/2,         (szx(3)/2+1):szx(3)};
            Oct{7} = {1:szx(1)/2,       (szx(2)/2+1):szx(2),(szx(3)/2+1):szx(3)};
            Oct{8} = {szx(1)/2+1:szx(1),(szx(2)/2+1):szx(2),(szx(3)/2+1):szx(3)};
        elseif strcmpi(mode,'minmax')
            Oct{1} = {[1,szx(1)/2],       [1,szx(2)/2],         [1,szx(3)/2]};
            Oct{2} = {[szx(1)/2+1,szx(1)],[1,szx(2)/2],         [1,szx(3)/2]};
            Oct{3} = {[1,szx(1)/2],       [(szx(2)/2+1),szx(2)],[1,szx(3)/2]};
            Oct{4} = {[szx(1)/2+1,szx(1)],[(szx(2)/2+1),szx(2)],[1,szx(3)/2]};
            Oct{5} = {[1,szx(1)/2],       [1,szx(2)/2],         [(szx(3)/2+1),szx(3)]};
            Oct{6} = {[szx(1)/2+1,szx(1)],[1,szx(2)/2],         [(szx(3)/2+1),szx(3)]};
            Oct{7} = {[1,szx(1)/2],       [(szx(2)/2+1),szx(2)],[(szx(3)/2+1),szx(3)]};
            Oct{8} = {[szx(1)/2+1,szx(1)],[(szx(2)/2+1),szx(2)],[(szx(3)/2+1),szx(3)]};
        end
end

if nargout > 1
    for i = 1:nargout
        varargout{i} = Oct{i};
    end
else
    varargout = {Oct};
end
