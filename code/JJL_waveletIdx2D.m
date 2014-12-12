function idx = JJL_waveletIdx2D(SIZE,SCALES,origRC)
%EXPERIMENTAL -- DO NOT USE
if nargin ==0
    SIZE =[ 480 720];
    SCALES = 1:3;
    origRC = [188 431];
end

idx = [];
for i = 1:length(SCALES)
    sc = SIZE./2^SCALES(i);
    for j = 1:size(origRC,1)
        oc = origRC(j,:)./2^SCALES(i);
        idx(end+1,:) = oc+sc;
        idx(end+1,:) = oc+sc.*[0 1];
        idx(end+1,:) = oc+sc.*[1 0];
    end
end
for j = 1:size(origRC,1)
    idx(end+1,:) = origRC(j,:)./2^(SCALES(end)+1);
end

if nargout == 0
    im = zeros(SIZE);
    i = sub2ind(SIZE,round(idx(:,1)),round(idx(:,2)));
    im(i) = 1;
    imshow(im,[]);
end
