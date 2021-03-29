function [x0, okay, feas, margin] = findinteriorpoint(A, b, Aeq, beq, tol, maxinside)
% [x0, okay, feas, margin] = findinteriorpoint(A, b, [Aeq], [beq], [tol])
%
% Find a strictly feasible point x0 so that A*x0 <= b - tol. If no such point
% can be found, okay is set to False.
%
% If there is at least a feasible point (but not necessarily on the interior),
% then feas is true, and x0 gives that point. If both okay and feas are false,
% then x0 is meaningless.
%
% margin gives the value e such that A*x0 <= b - e.
%
% The origin is always checked first. If it does not work, an LP is solved
% to find a valid point.
%
% tol is used to decide how much on the interior we need to be. If not
% supplied, the default value is 1e-8*max(1, norm(b)/N). Note that this may
% be too large if b is poorly scaled.
if nargin() < 2
    error('A and b are mandatory.')
elseif nargin() < 5
    tol = 1e-8*max(1, norm(b)/length(b));
end
if nargin() < 6
    maxinside = 1;
else
    maxinside = max(tol, maxinside);
end
[m, n] = size(A);
if nargin() < 4 || isempty(Aeq)
    Aeq = zeros(0, n);
    beq = [];
end
meq = size(Aeq, 1);
okay = false();

% Check whether the origin is on the inside.
if all(abs(beq) < tol) && all(b > tol)
    x0 = zeros(n, 1);
    okay = true();
    feas = true();
    margin = min(b);
end

% Try to use fminsearch if there are no equality constraints. Doesn't work
% well if the number of dimensions is very high, so we cap it at 10.
if ~okay && meq == 0 && m <= 10
    options = optimset('display', 'off');
    [x0, maxr] = fminsearch(@(x) max(max(A*x - b), -1e5*tol), A\b, options);
    okay = (maxr < -tol);
    feas = okay;
    margin = -maxr;
end

% Solve LP otherwise.
if ~okay
    c = [zeros(n, 1); -1];
    AA = [A, ones(m, 1)];
    AAeq = [Aeq, zeros(meq,1)];
    lb = [-inf(n, 1); 0];
    ub = [inf(n, 1); maxinside];
    if isOctave()
        ctype = [repmat('U', m, 1); repmat('S', meq, 1)];
        [xtilde, ~, err, extra] = glpk(c, [AA; AAeq], [b; beq], lb, ub,  ...
            ctype, repmat('C', n + 1, 1), 1, struct('msglev', 0));
        okay = (err == 0 && extra.status == 5);
    else
        options = optimoptions('linprog', 'display', 'off', 'algorithm', 'dual-simplex');
        [xtilde, ~, exitflag] = linprog(c, AA, b, AAeq, beq, lb, ub, [], options);
        okay = (exitflag == 1);
    end
    if isempty(xtilde)
        margin = -inf();
    else
        margin = xtilde(end); % Amount constraints could be tightened.
    end
    okay = okay && (margin >= tol);
    if isempty(xtilde)
        x0 = zeros(n, 1);
        okay = false();
    else
        x0 = xtilde(1:end-1);
    end

    % Check feasibility of x0.
    feas = all(A*x0 - b < tol);
    if feas && ~isempty(Aeq)
        feas = all(abs(Aeq*x0 - beq) < tol);
    end
    okay = okay && feas;
end
end%function