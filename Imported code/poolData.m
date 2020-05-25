function Theta = poolData(x,polyorder)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

% Cannot handle heigher than 3rd order polynomial yet.

% x is matrix of states. Each column is a time series of a single state

n = size(x,2); % Number of states
m = size(x,1); % Number of data samples of each state
ind = 1;
% poly order 0
Theta(:,ind) = ones(m,1);
ind = ind+1;

% poly order 1
for i=1:n
    Theta(:,ind) = x(:,i);
    ind = ind+1;
end

if(polyorder>=2)    % poly order 2
    for i=1:n
        for j=i:n
            Theta(:,ind) = x(:,i).*x(:,j);
            ind = ind+1;
        end
    end
end

if(polyorder>=3)    % poly order 3
    for i=1:n
        for j=i:n
            for k=j:n
                Theta(:,ind) = x(:,i).*x(:,j).*x(:,k);
                ind = ind+1;
            end
        end
    end
end