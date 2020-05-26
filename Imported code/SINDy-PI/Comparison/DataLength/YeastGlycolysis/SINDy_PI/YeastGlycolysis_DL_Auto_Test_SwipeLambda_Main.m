%% This file is the main file of using the SINDy-PI method to
% infer the Yeast Glycolysis Model. We will figure out what is the minimum
% data length needed for SINDy-PI to accurately discover the six-th state
% of the Yeast Glycolysis Model. This file will be used to swipe through
% different sparse parameter.
%
% Date: 2019/06/12
% Coded By: K

%% Close all, clear all, clc
close all;clear all; clc;
set(0,'defaulttextInterpreter','latex')
addpath('Functions')
addpath('Datas')
%% Define some parameters
% Define whehter you have control, if you have it, please define it
n_control=0;u=0;

% Run the ODE files and gather the simulation data.
% (We use the same data for both the iSINDy method and SINDy-PI method for better comparision)
% (Baseline data for the data length comparison of state 1,2,3,4,5,6,7)
% The data is noise clean. It is generated using 900 different initial
% conditions with 51 points of each initial condition.
load('TrainingData.mat')

% Choose whether you want to display actual ODE or not
disp_actual_ode=1;

% If the ODEs you want to display is the actual underlyting dynamics of the
% system, please set actual as 1
actual=1;

% Define how many states we have in our example
n_state=7;

% Print the actual ODE we try to discover
digits(4)
Print_ODEs(@(t,y)YeastGlycolysis_ODE(t,y),n_state,n_control,disp_actual_ode,actual);

% Create symbolic states
dz=sym('dz',[n_state,1]);

% Now we first create the parameters of the function right hand side
Highest_Poly_Order_Guess=0;
Highest_Trig_Order_Guess=0;
Highest_U_Order_Guess=0;

% Then create the right hand side library parameters
Highest_Trig_Order=0;
Highest_U_Order=0;
Highest_dPoly_Order=1;

% Set the parameter normLib=1 to normalize the librry
normLib=1;

% Set how amny iterations you want for the sparse regression
N_iter=10;

% Set whether you want to display the ODE or not
disp=0;

% Set the library size
Highest_Poly_Order_Lib=[6;6;3;3;3;6;3];

% Determine how many percent of data you want.
percent_start=0.03;d_percent=0.005;percent_end=0.05;

% Determine the lambda
lam_start=1;lam_end=20;d_lambda=0.02;

% Define a cell matrix to store the result
Xi=cell(round(percent_end-percent_start)/d_percent+1,2,lam_end-lam_start+1);

% Determine which states you want to discover
Which_State=6;

% Determine the left hand side guess number
if Which_State==1 || Which_State==2 || Which_State==6
    LHS_Num=2;
else
    LHS_Num=1;
end

% Get the poly-order
Highest_Poly_Order=Highest_Poly_Order_Lib(Which_State);

% Determine which LHS you want to use, the first one or the second one
LHS_Pin=2;

% Create the new directory to save the result
FolderName=strcat('Result_DL_SINDy_Data_Length_Compare_State_',num2str(Which_State),'_LHS_Guess_',num2str(LHS_Pin),'_V2');
[fld_status, fld_msg, fld_msgID]=mkdir(FolderName);

tic
% Create the library data using all the data points
[SINDy_Data_Full,SINDy_Struct]=SINDyLib(xt,dxt(:,Which_State),Which_State,u,Highest_Poly_Order,Highest_Trig_Order,Highest_U_Order,Highest_dPoly_Order);
toc

% Create left hand side guess
[LHS_Data_Full,LHS_Sym]=GuessLib(xt,dxt(:,Which_State),Which_State,u,Highest_Poly_Order_Guess,Highest_Trig_Order_Guess,Highest_U_Order_Guess);

tic
% Create the right hand side, exclude the guess from SINDy library
[RHS_Data_Full,RHS_Struct]=ExcludeGuess(SINDy_Data_Full,SINDy_Struct,LHS_Sym{LHS_Pin});
toc

parlooooop=parpool(2)
%% Start!
for Total_Run=1:20
    fprintf('\n \n Get the result for the %i time...\n',Total_Run)
    
    % Set a pin to count which iteration it is now
    pinpin=0;
    
    % Define a cell matrix to store the value of discovery result
    Xi=cell(round(percent_end-percent_start)/d_percent+1,lam_end-lam_start+1);
    ODE_Guess=cell(round(percent_end-percent_start)/d_percent+1,lam_end-lam_start+1);
    ODE=cell(round(percent_end-percent_start)/d_percent+1,lam_end-lam_start+1);
    
    for percent=percent_start:d_percent:percent_end
        
        fprintf('\n Testing the percentage as %i ...\n',percent*100)
        
        pinpin=pinpin+1;
        
        % Define the new data length
        new_length=round(percent*length(xt));
        
        % Shuffel the original data
        Sequence=randperm(size(SINDy_Data_Full,1));
        Sequence_Trimed=Sequence(1:new_length);
        RHS_Data=RHS_Data_Full(Sequence_Trimed,:);
        LHS_Data=LHS_Data_Full(Sequence_Trimed,:);
        
        fprintf('\n \t Calculating the %i expression...\n',Which_State)
        
        % Print the left hand side that we are testing
        fprintf('\t Testing the left hand side as %s:\n',char(LHS_Sym{LHS_Pin}))
        
        % Set up dummy variables for parfor
        LHS_Data_Dum=LHS_Data(:,LHS_Pin);
        LHS_Sym_Dum=LHS_Sym{LHS_Pin};
        dz_dum=dz(Which_State);
        
        % Start parfor and sweep the different values of lambda
        parfor pin=lam_start:lam_end
            % Perform the sparse regression problem
            tic
            [Xi{pinpin,pin},ODE{pinpin,pin}]=sparsifyDynamics(RHS_Data,LHS_Data_Dum,LHS_Sym_Dum,d_lambda*pin,N_iter,RHS_Struct,disp,normLib);
            toc
            % Perform sybolic calculation and solve for dX
            digits(6)
            ODE_Guess{pinpin,pin}=vpa(solve(LHS_Sym_Dum==ODE{pinpin,pin},dz_dum));
        end
    end
        
    % Save the calculation result of current iteration
    cc=clock;
    ResultName=strcat(FolderName,'/DL_SINDY_Data_Length_Result_',num2str(Total_Run),'__','LHS',num2str(LHS_Pin),'_',num2str(cc(3)),'_',num2str(cc(4)),'_',num2str(cc(5)),'_',num2str(round(cc(6))),'_P2.mat');
    save(ResultName,'Xi','percent_start','d_percent','percent_end','lam_start','lam_end','d_lambda')
    
end



