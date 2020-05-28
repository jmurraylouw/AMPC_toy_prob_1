function Xi = sparsifyDynamics(Theta,dXdt,lambda,n)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

% compute Sparse regression: sequential least squares
Xi = Theta\dXdt;  % initial guess: Least-squares

% lambda is our sparsification knob.
for k=1:10
    small_indexes = (abs(Xi)<lambda);   % find small coefficients
    Xi(small_indexes)=0;                % set small coeffs to 0 (threshold)
    for index = 1:n                     % n is state dimension
        big_indexes = ~small_indexes(:,index);
        % Regress dynamics onto remaining terms to find sparse Xi
        Xi(big_indexes,index) = Theta(:,big_indexes)\dXdt(:,index); 
    end
end