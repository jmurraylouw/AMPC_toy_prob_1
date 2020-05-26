% This file will generate the necessary LHS library to identify the PDE.
% Please modify this file to code the library you want.
%
% Last Updated: 2019/06/19
% Coded By: K

function [Data,Sym_Struct]=LHS_Guess_Implicit_PDE_mKdV(U,Ut,Utt,Ux,Uxx,Uxxx)
%% First get the size of the u vector.
[Data_Length,~]=size(U);

%Also create the symbolic variable
syms u ut utt ux uxx uxxx

Data=[];
Index=1;

Org_Data=[U Ut Ux Uxx];
Org_Sym=[u ut ux uxx];

Data(:,Index)=Org_Data(:,2);
Sym_Struct{1,Index}=Org_Sym(1,2);
Index=Index+1;

for i=2:size(Org_Data,2)
    Data(:,Index)=Org_Data(:,1).*Org_Data(:,i);
    Sym_Struct{1,Index}=Org_Sym(1,1)*Org_Sym(1,i);
    Index=Index+1;
end






