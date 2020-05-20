function [Ar,Br,Cr,Dr,Ur,Sr,Vr] = ERA(YY,m,n,nin,nout,r)
% Code by Steve Brunton 
% Based on An Eigensystem Realization Algorithm
% for Modal Parameter Identification and Model Reduction
% by Juang and Pappa, 1985

for i=1:nout
     for j=1:nin
         Dr(i,j) = YY(i,j,1);
         Y(i,j,:) = YY(i,j,2:end);
     end
 end

% Yss = Y(1,1,end);
% Y = Y-Yss;
% Y(i,j,k)::
% i refers to i-th output
% j refers to j-th input
% k refers to k-th timestep

% nin,nout number of inputs and outputs
% m,n dimensions of Hankel matrix
% r, dimensions of reduced model

assert(length(Y(:,1,1))==nout);
assert(length(Y(1,:,1))==nin);
assert(length(Y(1,1,:))>=m+n);

for i=1:m
    for j=1:n
        for Q=1:nout
            for P=1:nin
                H(nout*i-nout+Q,nin*j-nin+P) = Y(Q,P,i+j-1);
                H2(nout*i-nout+Q,nin*j-nin+P) = Y(Q,P,i+j);
            end
        end
    end
end

[U,S,V] = svd(H,'econ');
Sr = S(1:r,1:r);
Ur = U(:,1:r);
Vr = V(:,1:r);
Ar = Sr^(-.5)*Ur'*H2*Vr*Sr^(-.5);
Br = Sr^(-.5)*Ur'*H(:,1:nin);
Cr = H(1:nout,:)*Vr*Sr^(-.5);
HSVs = diag(S);