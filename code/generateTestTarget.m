function target = generateTestTarget(n,R)
%target = generateTestTarget(n,R)
%target = generateTestTarget(s,R)
%target = generateTestTarget()
%
%Generates a 2D or 3D resolution test target.
%   n --> [(2) 3] Number of dimensions to make the target.
%   s --> Seed.  Serves as a template for creating the target. Default is a
%       [16 x 16] or [16 x 16 x 16] predefined template. The 2D target has
%       vertical, horizontal, and checker-board patterns of multiple sizes,
%       beginning with lines of 1 pixel width.  The 3D default template
%       has vertical, horizontal and z-dimension planes, as well as 3D
%       checker-board patterns.
%   R --> Recusion Level.  Number of times to recursively grow the seed
%       target. Each time the target is grown it's size doubles (Default =
%       4)
if nargin == 0
    n = 2;
    s1 =[];
else
    if isscalar(n)
        s1 = [];
    else
        s1 = n;
        n = ndims(s1);
    end
end
if nargin < 2
    R = 5;
end
if isempty(s1)
    q1 = [
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        1 1 1 1 1 1 1 0;
        0 0 0 0 0 0 0 0;
        1 1 1 1 1 1 1 0;
        0 0 0 0 0 0 0 0;
        1 1 1 1 1 1 1 0;
        0 0 0 0 0 0 0 0];
    q2 = [
        0 0 0 0 0 0 0 0 ;
        0 0 1 0 1 0 1 0 ;
        0 0 1 0 1 0 1 0 ;
        0 0 1 0 1 0 1 0 ;
        0 0 1 0 1 0 1 0 ;
        0 0 1 0 1 0 1 0 ;
        0 0 1 0 1 0 1 0 ;
        0 0 1 0 1 0 1 0 ];
    q3 = [
        0 0 0 0 0 0 0 0;
        0 1 0 1 0 1 0 1;
        0 0 1 0 1 0 1 0;
        0 1 0 1 0 1 0 1;
        0 0 1 0 1 0 1 0;
        0 1 0 1 0 1 0 1;
        0 0 1 0 1 0 1 0;
        0 0 0 0 0 0 0 0;
        ];
    q4 = [
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        ];
    if n==2
        s1 = [q2 q1;
            q3 q4 ];
    elseif n == 3
        q1 = repmat(q1,[1,1,16]);
        q1(:,:,1) = 0; q1(:,:,16)=0;
        q2 = repmat(q2,[1,1,16]);
        q2(:,:,1) = 0; q2(:,:,16)=0;
        q3 = repmat(q3,[1,1,16]);
        q3(2:7,2:8,2:2:16) = 1-q3(2:7,2:8,2:2:16);
        q3(:,:,1) = 0; q3(:,:,16)=0;
        q4 = repmat(q4,[1,1,16]);
        q4(2:8,2:7,1:2:16) = 1;

        s1 = [q2 q1;
            q3 q4 ];
    end
end
target = pad(s1,R);
if nargout == 0
    if n==2
        imshow(target,[]);
        axis on;
    else
        i = find(target);
        [x,y,z] = ind2sub(size(target),i);
        scatter3(x,y,z,1,'filled');
        axis on;
    end
end

function s = pad(s1,T)
s = s1;

%Break recursion
if T >= 1
    sz = size(s1);
    for i = 1:ndims(s1)
        %For each dimension we expand by 2
        order = circshift(1:ndims(s1),[1,(i-1)]);
        m = floor(1:.5:(sz(i)+.5));
        s = permute(s,order);
        s = s(m,:,:);
        s = ipermute(s,order);
    end
    if ndims(s1)==2
        s((end/2+1):end,(end/2+1):end,floor(end/2+1):end) = s1;
    elseif ndims(s1)==3
        s((end/2+1):end,(end/2+1):end,(end/2+1):end) = s1;
    end

    s = pad(s,T-1);
end