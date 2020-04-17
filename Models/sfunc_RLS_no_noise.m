function sfunc_RLS(block)
%MSFUNTMPL_BASIC A Template for a Level-2 MATLAB S-Function
%   The MATLAB S-function is written as a MATLAB function with the
%   same name as the S-function. Replace 'msfuntmpl_basic' with the 
%   name of your S-function.

%   Copyright 2003-2018 The MathWorks, Inc.

%%
%% The setup method is used to set up the basic attributes of the
%% S-function such as ports, parameters, etc. Do not add any other
%% calls to the main body of the function.
%%
setup(block);

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

% Register number of ports
block.NumInputPorts  = 2;
block.NumOutputPorts = 3;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic();
block.SetPreCompOutPortInfoToDynamic();

% Override INPUT port properties
% f (Inputs to plant)
block.InputPort(1).Dimensions        = 1;
block.InputPort(1).DatatypeID        = 0;  % double
block.InputPort(1).Complexity        = 'Real';
block.InputPort(1).DirectFeedthrough = true;

% x (Actual outputs of plant/Desired outputs of model)
block.InputPort(2).Dimensions        = 1;
block.InputPort(2).DatatypeID        = 0;  % double
block.InputPort(2).Complexity        = 'Real';
block.InputPort(2).DirectFeedthrough = true;

% Override OUTPUT port properties
% y (Predicted output using previous model)
block.OutputPort(1).Dimensions       = 1;
block.OutputPort(1).DatatypeID       = 0; % double
block.OutputPort(1).Complexity       = 'Real';
block.OutputPort(1).SamplingMode     = 'Sample';

% eps (error of model / eps(n)  = x(n) - y(n))
block.OutputPort(2).Dimensions       = 1;
block.OutputPort(2).DatatypeID       = 0; % double
block.OutputPort(2).Complexity       = 'Real';
block.OutputPort(2).SamplingMode     = 'Sample';

% Theta (new model parameters)
block.OutputPort(3).Dimensions       = block.DialogPrm(1).Data; % May need to use set outputportdimentions
block.OutputPort(3).DatatypeID       = 0; % double
block.OutputPort(3).Complexity       = 'Real';
block.OutputPort(3).SamplingMode     = 'Sample';

% Register parameters
% parameter 1 = M (number of parameters in w)
% parameter 2 = lambda (forgetting factor)
% parameter 3 = sample time

block.NumDialogPrms     = 3;

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
block.SampleTimes = [block.DialogPrm(3).Data 0];

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
    M = block.DialogPrm(1).Data;   
    
    block.NumDworks = 3;
  
    % Vector of next phi
    block.Dwork(1).Name            = 'phi';
    block.Dwork(1).Dimensions      = M;
    block.Dwork(1).DatatypeID      = 0;      % double
    block.Dwork(1).Complexity      = 'Real'; % real
    block.Dwork(1).UsedAsDiscState = true;

    % Previous value of model parameters
    block.Dwork(2).Name            = 'Theta';
    block.Dwork(2).Dimensions      = M;
    block.Dwork(2).DatatypeID      = 0;      % double
    block.Dwork(2).Complexity      = 'Real'; % real
    block.Dwork(2).UsedAsDiscState = true;
    
    % Previous value of P, in reshaped vector form, P2
    % Needs to reshape because Dwork cannot store matrix
    block.Dwork(3).Name            = 'P2';
    block.Dwork(3).Dimensions      = M*M; % Matrix is reshaped as vector, Dwork can only be vecotr
    block.Dwork(3).DatatypeID      = 0;      % double
    block.Dwork(3).Complexity      = 'Real'; % real
    block.Dwork(3).UsedAsDiscState = true;


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
    M       = block.DialogPrm(1).Data;
    
    % Inititialise variables    
    phi     = zeros(M,1);
    Theta   = zeros(M,1); % Inititialise parameter vectors, Time steps downwards
    P       = 100*eye(M); % Initialise large P
    P2      = reshape(P, M*M, 1); % Form a vector to store in Dwork

    % Assign to memory/Dwork
    block.Dwork(1).Data = phi;
    block.Dwork(2).Data = Theta;
    block.Dwork(3).Data = P2;

%end Start

%%
%% Outputs:
%%   Functionality    : Called to generate block outputs in
%%                      simulation step
%%   Required         : Yes
%%   C MEX counterpart: mdlOutputs
%%
function Outputs(block)
    f       = block.InputPort(1).Data;
    x       = block.InputPort(2).Data;

    phi     = block.Dwork(1).Data;
    Theta   = block.Dwork(2).Data;
    P2    	= block.Dwork(3).Data;
    
    M       = block.DialogPrm(1).Data;
    lambda  = block.DialogPrm(2).Data;
    
    P       = reshape(P2, M, M); % Reshape vector P2 in matrix P
    % phi   = [f(n-1); x(n-1); x(n-2); eps(n-1); eps(n-2)]
    pi      = phi'*P; % Use M previous values of U as well
    gamma   = lambda + pi*phi;
    k       = pi'/gamma;
    y       = Theta'*phi;
    eps     = x - y;
    Theta   = Theta + k*eps; % Estimate next parameter vector

    P_prime = k*pi;
    P       = 1/lambda*(P - P_prime); % Update P for next use
    P2      = reshape(P, M*M, 1); % Form a vector to store in Dwork   
    phi     = [f; x; phi(2)];
    
    block.OutputPort(1).Data = y;
    block.OutputPort(2).Data = eps;
    block.OutputPort(3).Data = Theta;
    
    block.Dwork(1).Data = phi;
    block.Dwork(2).Data = Theta;
    block.Dwork(3).Data = P2;

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
    block.OutputPort(1).SamplingMode = mode;
    block.OutputPort(2).SamplingMode = mode;
    block.OutputPort(3).SamplingMode = mode;
%end SetInputPortSamplingMode

%%
%% Terminate:
%%   Functionality    : Called at the end of simulation for cleanup
%%   Required         : Yes
%%   C MEX counterpart: mdlTerminate
%%
function Terminate(block)

%end Terminate

