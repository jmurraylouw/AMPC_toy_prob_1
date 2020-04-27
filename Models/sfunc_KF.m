function sfunc_KF(block)
%sfunc_DMDc_MW - Applies Dynamic Mode Decomposition with Control and a
%Moving Window of past input data to obtain the state space system matrixes

%%
%% The setup method is used to set up the basic attributes of the
%% S-function such as ports, parameters, etc. Do not add any other
%% calls to the main body of the function.
%%
setup(block);
% Get dimentions of state and input vectors

%endfunction

%% Function: setup ===================================================
%% Abstract:
%%   Set up the basic characteristics of the S-function block such as:
%%   - Input ports
%%   - Output ports
%%   - Dialog parameters
%%   - Options
%%
%%   Required         : Yes
%%   C MEX counterpart: mdlInitializeSizes
%%
function setup(block)
% Register parameters
% parameter 1 = F
% parameter 2 = G
% parameter 3 = H
% parameter 4 = D
% parameter 5 = Ts (sample time)
% parameter 6 = Q
% parameter 7 = R
% parameter 8 = x0
% parameter 9 = P0
% A,B,C,D,Ts,Q,R,x0,P0
block.NumDialogPrms     = 9;

% Read dialog parameters
size_G = size(block.DialogPrm(2).Data);
size_H = size(block.DialogPrm(3).Data);
nx = size_G(1);
nu = size_G(2);
ny = size_H(1);

% Register number of ports
block.NumInputPorts  = 2;
block.NumOutputPorts = 1;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override INPUT port properties
% u (input vector)
block.InputPort(1).Dimensions        = nu;
block.InputPort(1).DatatypeID        = 0;  % double
block.InputPort(1).Complexity        = 'Real';

% y (measurement vector)
block.InputPort(2).Dimensions        = ny;
block.InputPort(2).DatatypeID        = 0;  % double
block.InputPort(2).Complexity        = 'Real';


% Override OUTPUT port properties
% x_hat (estimated state vector)
block.OutputPort(1).DatatypeID       = 0; % double
block.OutputPort(1).Complexity       = 'Real';
block.OutputPort(1).SamplingMode     = 'Sample';
block.OutputPort(1).Dimensions       = [nx, 1];

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
block.SampleTimes = [block.DialogPrm(5).Data 0]; % Set sample time

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'CustomSimState',  < Has GetSimState and SetSimState methods
%    'DisallowSimState' < Error out when saving or restoring the model sim state
block.SimStateCompliance = 'DefaultSimState';

%% -----------------------------------------------------------------
%% The MATLAB S-function uses an internal registry for all
%% block methods. You should register all relevant methods
%% (optional and required) as illustrated below. You may choose
%% any suitable name for the methods and implement these methods
%% as local functions within the same file. See comments
%% provided for each function for more information.
%% -----------------------------------------------------------------

block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
block.RegBlockMethod('InitializeConditions', @InitializeConditions);
block.RegBlockMethod('Start', @Start);
block.RegBlockMethod('Outputs', @Outputs);     % Required
block.RegBlockMethod('Update', @Update);
block.RegBlockMethod('Derivatives', @Derivatives);
block.RegBlockMethod('Terminate', @Terminate); % Required
block.RegBlockMethod('SetInputPortSamplingMode', @SetInputPortSamplingMode);
  
%end setup

%%
%% PostPropagationSetup:
%%   Functionality    : Setup work areas and state variables. Can
%%                      also register run-time methods here
%%   Required         : No
%%   C MEX counterpart: mdlSetWorkWidths
%%
function DoPostPropSetup(block)

    block.NumDworks = 2;

    size_G = size(block.DialogPrm(2).Data);
    nx = size_G(1);
    
    block.Dwork(1).Name            = 'x_hat';
    block.Dwork(1).Dimensions      = nx;
    block.Dwork(1).DatatypeID      = 0;      % double
    block.Dwork(1).Complexity      = 'Real'; % real
    block.Dwork(1).UsedAsDiscState = true;

    % Needs to reshape because Dwork cannot store matrix    
    block.Dwork(2).Name            = 'P';
    block.Dwork(2).Dimensions      = nx*nx;
    block.Dwork(2).DatatypeID      = 0;      % double
    block.Dwork(2).Complexity      = 'Real'; % real
    block.Dwork(2).UsedAsDiscState = true;
  
%%
%% InitializeConditions:
%%   Functionality    : Called at the start of simulation and if it is 
%%                      present in an enabled subsystem configured to reset 
%%                      states, it will be called when the enabled subsystem
%%                      restarts execution to reset the states.
%%   Required         : No
%%   C MEX counterpart: mdlInitializeConditions
%%
function InitializeConditions(block)

%end InitializeConditions


%%
%% Start:
%%   Functionality    : Called once at start of model execution. If you
%%                      have states that should be initialized once, this 
%%                      is the place to do it.
%%   Required         : No
%%   C MEX counterpart: mdlStart
%%
function Start(block)
    F = block.DialogPrm(1).Data;   
    Q = block.DialogPrm(6).Data;
    size_G = size(block.DialogPrm(2).Data);
    nx = size_G(1);
    
    % Inititialise variables    
    x_hat = block.DialogPrm(8).Data;
    P = block.DialogPrm(9).Data;
    
    x_hat = F*x_hat; % Priori estimation / Extrapolate state
    P = F*P*F' + Q; % Extrapolate uncertainty

    % Assign to memory/Dwork
    block.Dwork(1).Data = x_hat;
    block.Dwork(2).Data = reshape(P, 1, nx*nx);
%end Start

%%
%% Outputs:
%%   Functionality    : Called to generate block outputs in
%%                      simulation step
%%   Required         : Yes
%%   C MEX counterpart: mdlOutputs
%%
function Outputs(block)
    % Dialog parameters
    F = block.DialogPrm(1).Data;
    G = block.DialogPrm(2).Data;
    H = block.DialogPrm(3).Data;
    Q = block.DialogPrm(6).Data;
    R = block.DialogPrm(7).Data;
    
    size_G = size(G);
    size_H = size(H);
    nx = size_G(1);
    nu = size_G(2);
    ny = size_H(1);

    % Input
    u       = block.InputPort(1).Data;
    y       = block.InputPort(2).Data;
    
    %Dwork data
    x_hat = block.Dwork(1).Data; % Priori estimation 
    P = reshape(block.Dwork(2).Data, nx, nx);
    
    % Update step
    K = (P*H')/(H*P*H' + R); % Compute Kalman gain (b*inv(A) -> b/A)
    x_hat = x_hat + K*(y - H*x_hat); % Posteriori / With measurement
    KH_term = (eye(nx) - K*H);
    P = KH_term*P*KH_term' + K*R*K'; % Update estimate uncertainty
   
    % Output
    block.OutputPort(1).Data = x_hat;
    
    % Extrapulate/Prediction for next time step
    x_hat = F*x_hat + G*u; % Priori estimation / Extrapolate state
    P = F*P*F' + Q; % Extrapolate uncertainty

    % Update Dwork memory   
    block.Dwork(1).Data = x_hat;
    block.Dwork(2).Data = reshape(P, 1, nx*nx);

%end Outputs

%%
%% Update:
%%   Functionality    : Called to update discrete states
%%                      during simulation step
%%   Required         : No
%%   C MEX counterpart: mdlUpdate
%%
function Update(block)


%end Update

%%
%% Derivatives:
%%   Functionality    : Called to update derivatives of
%%                      continuous states during simulation step
%%   Required         : No
%%   C MEX counterpart: mdlDerivatives
%%
function Derivatives(block)

%end Derivatives

function SetInputPortSamplingMode(block, port, mode)
    block.InputPort(port).SamplingMode = mode;
    
%end SetInputPortSamplingMode

%%
%% Terminate:
%%   Functionality    : Called at the end of simulation for cleanup
%%   Required         : Yes
%%   C MEX counterpart: mdlTerminate
%%


function Terminate(block)

%end Terminate

