function sfunc_DMDc_MW(block)
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
% parameter 1 = Ts (sample time)
% parameter 2 = window (width of data memory window in data steps)
% parameter 3 = nx (length of state vector)
% parameter 4 = nu (length of input vector)
block.NumDialogPrms     = 4;

% Read dialog parameters
nx = block.DialogPrm(3).Data;
nu = block.DialogPrm(4).Data;

% Register number of ports
block.NumInputPorts  = 2;
block.NumOutputPorts = 3;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override INPUT port properties
% u (input vector)
block.InputPort(1).Dimensions        = nu;
block.InputPort(1).DatatypeID        = 0;  % double
block.InputPort(1).Complexity        = 'Real';

% x (state vector)
block.InputPort(2).Dimensions        = nx;
block.InputPort(2).DatatypeID        = 0;  % double
block.InputPort(2).Complexity        = 'Real';


% Override OUTPUT port properties
% A (Flattened system matrix)
block.OutputPort(1).DatatypeID       = 0; % double
block.OutputPort(1).Complexity       = 'Real';
block.OutputPort(1).SamplingMode     = 'Sample';
block.OutputPort(1).Dimensions       = [nx, nx];

% B (Input matrix)
block.OutputPort(2).DatatypeID       = 0; % double
block.OutputPort(2).Complexity       = 'Real';
block.OutputPort(2).SamplingMode     = 'Sample';
block.OutputPort(2).Dimensions       = [nx, nu];

% Mean Squared Error
block.OutputPort(3).Dimensions       = 1;
block.OutputPort(3).DatatypeID       = 0; % double
block.OutputPort(3).Complexity       = 'Real';
block.OutputPort(3).SamplingMode     = 'Sample';

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
block.SampleTimes = [block.DialogPrm(1).Data 0]; % Set sample time

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
    % w = Timestep width of window (T_window/Ts)    
    w = block.DialogPrm(2).Data; 
    
    % Get dimentions of state and input vectors
    nx = block.InputPort(2).Dimensions; % Length of state vector
    nu = block.InputPort(1).Dimensions; % Length of input vector
    
    block.NumDworks = 4;
  
    % Matrix of window of state vectors
    % Needs to reshape because Dwork cannot store matrix
    block.Dwork(1).Name            = 'X';
    block.Dwork(1).Dimensions      = w*nx;
    block.Dwork(1).DatatypeID      = 0;      % double
    block.Dwork(1).Complexity      = 'Real'; % real
    block.Dwork(1).UsedAsDiscState = true;

    % Matrix of window of input vectors
    % Needs to reshape because Dwork cannot store matrix
    block.Dwork(2).Name            = 'U';
    block.Dwork(2).Dimensions      = w*nu;
    block.Dwork(2).DatatypeID      = 0;      % double
    block.Dwork(2).Complexity      = 'Real'; % real
    block.Dwork(2).UsedAsDiscState = true;
    
    block.Dwork(3).Name            = 'A_prev';
    block.Dwork(3).Dimensions      = nx*nx;
    block.Dwork(3).DatatypeID      = 0;      % double
    block.Dwork(3).Complexity      = 'Real'; % real
    block.Dwork(3).UsedAsDiscState = true;

    block.Dwork(4).Name            = 'B_prev';
    block.Dwork(4).Dimensions      = nx*nu;
    block.Dwork(4).DatatypeID      = 0;      % double
    block.Dwork(4).Complexity      = 'Real'; % real
    block.Dwork(4).UsedAsDiscState = true;
    
    
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
    % w = Timestep width of window (T_window/Ts)    
    w = block.DialogPrm(2).Data;   
    
    % Get dimentions of state and input vectors
    nx = block.InputPort(2).Dimensions; % Length of state vector
    nu = block.InputPort(1).Dimensions; % Length of input vector
    
    % Inititialise variables    
    X = zeros(nx,w);
    U = zeros(nu,w);
    
    % Reshape to vector
    X_dwork = reshape(X, 1, w*nx);
    U_dwork = reshape(U, 1, w*nu);
    
    % Assign to memory/Dwork
    block.Dwork(1).Data = X_dwork;
    block.Dwork(2).Data = U_dwork;
%end Start

%%
%% Outputs:
%%   Functionality    : Called to generate block outputs in
%%                      simulation step
%%   Required         : Yes
%%   C MEX counterpart: mdlOutputs
%%
function Outputs(block)
    
    w = block.DialogPrm(2).Data; % Window width
    
    % Get dimentions of state and input vectors
    nx = block.InputPort(2).Dimensions; % Length of state vector
    nu = block.InputPort(1).Dimensions; % Length of input vector

    u       = block.InputPort(1).Data;
    x       = block.InputPort(2).Data;
    
    X_dwork = block.Dwork(1).Data;
    U_dwork = block.Dwork(2).Data;
    A_dwork = block.Dwork(3).Data;
    B_dwork = block.Dwork(4).Data;

    X = reshape(X_dwork, nx, w); % Reshape vector into matrix
    U = reshape(U_dwork, nu, w); % Reshape vector into matrix
    A = reshape(A_dwork, nx, nx);
    B = reshape(B_dwork, nx, nu);
    
    % Add input data to X2
    X2 = [X(:, 2:end), x];
    
    % Calculate Mean Squared Error
    X2_calc = A*X; % Calculate X2 according to A
    MSE = mean((X2 - X2_calc).^2, 'all'); % Mean Squared Error

    % Calculate A and B
    % Based on DMD control example video by Steve Brunton
    XU = [X; U];
    AB = X2*pinv(XU);
    A  = AB(:,1:nx);
    B  = AB(:,(nx+1):end);
       
    % Output
    block.OutputPort(1).Data = A;
    block.OutputPort(2).Data = B;
    block.OutputPort(3).Data = MSE;
    
    % Update Dwork memory 
    X_dwork = reshape(X2, 1, w*nx);
    U_dwork = reshape([U(:, 2:end), u], 1, w*nu);
    A_dwork = reshape(A, 1, nx*nx);
    B_dwork = reshape(B, 1, nu*nx);
    
    block.Dwork(1).Data = X_dwork;
    block.Dwork(2).Data = U_dwork;
    block.Dwork(3).Data = A_dwork;
    block.Dwork(4).Data = B_dwork;
   
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

