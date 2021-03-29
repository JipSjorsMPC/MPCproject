function [P, K] = riciter(A, B, Q, R, P)
    % Performs one step if riccati iteration.
    P = Q + A'*P*A - A'*P*B*((B'*P*B + R)\(B'*P*A));
    K = -(B'*P*B + R)\(B'*P*A);