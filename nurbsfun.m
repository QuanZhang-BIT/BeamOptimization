function [C, N, R, S, U] = nurbsfun(n, t, w, P, U)
% Evaluate NURBS at specified locations.
%
% Input arguments:
% n:
%    NURBS order (2 for linear, 3 for quadratic, 4 for cubic, etc.)
% t:
%    knot vector
% w:
%    weight vector
% P:
%    control points, typically 2-by-m
% U (optional):
%    values where the NURBS is to be evaluated, or a positive
%    integer to set the number of points to automatically allocate
%
% Output arguments:
% C:
%    points of the NURBS curve
% N:
%    basis function value for each element in U
% R:
%    rational function value for each element in U
% S:
%    column containing the value of the summarize of basis function per
%    weight
% U:
%    points where the NURBS curve is evaluated

%Written by Graziano Fuccio, email: g.fuccio359@gmail.com
    validateattributes(n, {'numeric'}, {'positive','integer','scalar'});
    assert(all( t(2:end)-t(1:end-1) >= 0 ), 'nurbs:InvalidArgumentValue', ...
        'Knot vector values should be nondecreasing.');
    validateattributes(P, {'numeric'}, {'real','2d'});
    nctrl = numel(t)- n;
    assert(size(P,2) == nctrl, 'nurbs:DimensionMismatch', ...
        'Invalid number of control points, %d given, %d required.', size(P,2), nctrl);
    assert(size(P,2) == numel(w), 'nurbs:DimensionMismatch', ...
        'Invalid number of weight points, %d given, %d required.', numel(w), size(P,2));
    if nargin < 5
        U = linspace(t(1), t(end), 15*size(P,2));% allocate points uniformly,
    elseif isscalar(U) && U > 1
        validateattributes(U, {'numeric'}, {'positive','integer','scalar'});
        U = linspace(t(1), t(end), U);  % allocate points uniformly
    end

    totalU = numel(U);%total elements to evaluate
    nc = size(P, 2);
    N = zeros(totalU, nc);
    
    %calcolate basis function
    for i = 1 : totalU
        u = U(i);
        N(i, :) = basisfunction(n, nc, u, t);
    end
    
    %calcolate denominator of rational basis functions
    S = zeros(totalU, 1);
    for i = 1 : totalU
        tmp = N(i, :);
        for j = 1 : size(tmp, 2)
            S(i) = S(i) + (tmp(j) * w(j));
        end
    end
    
    %calcolate rational basis functions
    R = zeros(totalU, nc);
    for i = 1 : totalU
        tmp = N(i, :);
        for j = 1 : size(tmp, 2)
            if(S(i) ~= 0)
                R(i, j) = (tmp(j) * w(j)) / S(i);
            else
                R(i, j) = 0;
            end
        end
    end
    
    %calcolate curve points
    C = zeros(2, totalU);
    for i = 1 : totalU
        tmp = R(i, :);
        sum = [0; 0];
        for j = 1 : size(tmp, 2)
            sum = sum + (P(:, j) * tmp(j));
        end
        C(:, i) = sum;
    end
end

