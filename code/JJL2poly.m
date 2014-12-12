function out = JJL2poly(w,szx,J,type)
%out = JJL2poly(w,J,type)
%type - ['dwt2D',
%       'dwt3D',
%       'dwt1D',
%       'dtdwt2D',
%       'dtdwt3D',
%       'dtdwt1D',
%       'cmplxdtdwt2D',
%       'cmplxdtdwt3D']

out = [];
switch type
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Non - Expansive dwt's
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'dwt1D',
        for j = 1:J
            Oct = getOctantSubs(szx);
            out{j} = w(Oct{2});
            szx = szx/2;
        end
        out{J+1} = w(Oct{1});
    case 'dwt2D',
        for j = 1:J
            Oct = getOctantSubs(szx);
            for i = 1:3
                out{j}{i} = w(Oct{i+1}{1},Oct{i+1}{2});
            end
            szx = szx/2;
        end
        out{J+1} = w(Oct{1}{1},Oct{1}{2});
    case 'dwt3D',
        for j = 1:J
            Oct = getOctantSubs(szx);
            for i = 1:3
                out{j}{i} = w(Oct{i+1}{1},Oct{i+1}{2},Oct{i+1}{3});
            end
            szx = szx/2;
        end
        out{J+1} = w(Oct{1}{1},Oct{1}{2},Oct{1}{3});

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Expansive DWT's
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'dtdwt1D',
        for j = 1:J
            Oct = getOctantSubs(szx);
            for k = 1:2
                out{j}{k} = w{k}(Oct{2});
            end
            szx = szx/2;
        end
        for k = 1:2
            out{J+1}{k} = w{k}(Oct{1});
        end
    case 'dtdwt2D',
        for j = 1:J
            Oct = getOctantSubs(szx);
            for i = 1:3
                for k = 1:2
                    out{j}{k}{i} = w{k}(Oct{i+1}{1},Oct{i+1}{2});
                end
            end
            szx = szx/2;
        end
        for k = 1:2
            out{J+1}{k} = w{k}(Oct{1}{1},Oct{1}{2});
        end
    case 'dtdwt3D',
        for j = 1:J
            Oct = getOctantSubs(szx);
            for i = 1:7
                for k = 1:4
                    out{j}{k}{i} = w{k}(Oct{i+1}{1},Oct{i+1}{2},Oct{i+1}{3});
                end
            end
            szx = szx/2;
        end
        for k = 1:4
            out{J+1}{k} = w{k}(Oct{1}{1},Oct{1}{2},Oct{1}{3});
        end
    case 'cplxdtdwt2D',
        for j = 1:J
            Oct = getOctantSubs(szx);
            for i = 1:3
                for k = 1:2
                    for m = 1:2
                        out{j}{k}{m}{i} = w{k,m}(Oct{i+1}{1},Oct{i+1}{2});
                    end
                end
            end
            szx = szx/2;
        end
        for k = 1:2
            for m = 1:2
                out{J+1}{k}{m} = w{k,m}(Oct{1}{1},Oct{1}{2});
            end
        end
    case 'cplxdtdwt3D',
        for j = 1:J
            Oct = getOctantSubs(szx);
            for i = 1:7
                for k = 1:4
                    for m = 1:2
                        for n = 1:2
                            out{j}{k}{m}{n}{i} = ...
                                w{k,m,n}(Oct{i+1}{1},Oct{i+1}{2},...
                                Oct{i+1}{3});
                        end
                    end
                end
            end
            szx = szx/2;
        end
        for k = 1:4
            for m =1:2
                for n = 1:2
                    out{J+1}{k}{m}{n} = ...
                        w{k,m,n}(Oct{1}{1},Oct{1}{2},Oct{1}{3});
                end
            end
        end
    otherwise
end