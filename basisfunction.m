function [N] = basisfunction(n, npts, u, t)
% Evaluate basis function at specified u.
%
% Input arguments:
% n:
%    NURBS order (2 for linear, 3 for quadratic, 4 for cubic, etc.)
% npts:
%    number of control points
% u:
%    point to evaluate
% t:
%    knot vector
%
% Output arguments:
% N:
%   vector with size npts containing value of the basis function at u

%Written by Graziano Fuccio, email: g.fuccio359@gmail.com
nplusc = npts + n;
N = zeros(npts);

for i = 1: nplusc - 1
    if u >= t(i) && u <= t(i+1)
        N(i) = 1;
    else
        N(i) = 0;
    end
end

%application of the formula that you can find on the NURBS book
for k = 2 : n
    for i = 1 : nplusc - k
        if N(i) ~= 0
            Length1 = t(i + k - 1) - t(i);
            d = ((u - t(i)) * N(i)) / Length1;
            if Length1 == 0
                d = 0;
            end
        else
            d = 0;
        end
        if(N(i + 1) ~= 0)
            Length2 = t(i + k) - t(i + 1);
            e = ((t(i + k) - u) * N(i + 1)) / Length2;
            if Length2 == 0
                e = 0;
            end
        else
            e = 0;
        end
        N(i) = d + e;
    end
end

N = N(:, 1);
N = N';

end

