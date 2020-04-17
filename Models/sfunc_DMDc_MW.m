function sfunc_DMDc_MW(block)
%sfunc_DMDc_MW - Applies Dynamic Mode Decomposition with Control and a
%Moving Window of past input data to obtain the state space system matrixes

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
block.NumInputPorts  = 4;
block.NumOutputPorts = 3;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override INPUT port properties
% f (Force applied to model)
block.InputPort(1).Dimensions        = 1;
block.InputPort(1).DatatypeID        = 0;  % double
block.InputPort(1).Complexity        = 'Real';
block.InputPort(1).DirectFeedthrough = true;

% x (Plant output distance)
block.InputPort(2).Dimensions        = 1;
block.InputPort(2).DatatypeID        = 0;  % double
block.InputPort(2).Complexity        = 'Real';
block.InputPort(2).DirectFeedthrough = true;

% x_dot (Plant output velocity)
block.InputPort(3).Dimensions        = 1;
block.InputPort(3).DatatypeID        = 0;  % double
block.InputPort(3).Complexity        = 'Real';
block.InputPort(3).DirectFeedthrough = true;

% enable
block.InputPort(4).Dimensions        = 1;
block.InputPort(4).DatatypeID        = 0;  % double
block.InputPort(4).Complexity        = 'Real';
block.InputPort(4).DirectFeedthrough = true;


% Override OUTPUT port properties
% A (Flattened system matrix)
block.OutputPort(1).Dimensions       = [2, 2];
block.OutputPort(1).DatatypeID       = 0; % double
block.OutputPort(1).Complexity       = 'Real';
block.OutputPort(1).SamplingMode     = 'Sample';

% B (Input matrix)
block.OutputPort(2).Dimensions       = [2, 1];
block.OutputPort(2).DatatypeID       = 0; % double
block.OutputPort(2).Complexity       = 'Real';
block.OutputPort(2).SamplingMode     = 'Sample';

% A error (debug)
block.OutputPort(3).Dimensions       = [2, 2];
block.OutputPort(3).DatatypeID       = 0; % double
block.OutputPort(3).Complexity       = 'Real';
block.OutputPort(3).SamplingMode     = 'Sample';

% Register parameters
% parameter 1 = T_window (time width of data memory window)
% parameter 2 = Ts (sample time)
block.NumDialogPrms     = 2;

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
block.SampleTimes = [block.DialogPrm(2).Data 0]; % Set sample time

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
    w = block.DialogPrm(1).Data;    
    nx = 2; % Size of state vector 
    nu = 1; % Size of input vector
    
    block.NumDworks = 2;
  
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
    w = block.DialogPrm(1).Data;    
    nx = 2; % Size of state vector 
    nu = 1; % Size of input vector
    
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
    
    w = block.DialogPrm(1).Data; % Window width
    nx = 2; % Size of state vector 
    nu = 1; % Size of input vector

    f       = block.InputPort(1).Data;
    x       = block.InputPort(2).Data;
    x_dot   = block.InputPort(3).Data;
    
    X_dwork = block.Dwork(1).Data;
    U_dwork = block.Dwork(2).Data;

    X = reshape(X_dwork, nx, w); % Reshape vector into matrix
    U = reshape(U_dwork, nu, w); % Reshape vector into matrix
    
    % Add input data to X2
    X2 = [X(:, 2:end), [x_dot; x]];
   
    % Calculate A and B
    % Based on DMD control example video by Steve Brunton
    XU = [X; U];
    AB = X2*pinv(XU);
    A  = AB(:,1:2);
    B  = AB(:,end);
    A2=A;
    
    
    % If enable = 0, then output original model
    option = block.InputPort(4).Data;
    
    if option == 0
        % Model from mpc object
        A = [0.969169519504925  -0.246386829627017;  0.049277365925403  0.993808202467626];
        B = [0.049277365925403; 0.001238359506475];
        block.OutputPort(3).Data = A2-A ;
    end
    
    
    % Output
    block.OutputPort(1).Data = A;
    block.OutputPort(2).Data = B;
    
    % Update Dwork memory 
    X_dwork = reshape(X2, 1, w*nx);
    U_dwork = reshape([U(:, 2:end), f], 1, w*nu);
    
    block.Dwork(1).Data = X_dwork;
    block.Dwork(2).Data = U_dwork;
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
%     block.OutputPort(1).SamplingMode = mode;
%     block.OutputPort(2).SamplingMode = mode;
%     block.OutputPort(3).SamplingMode = mode;
    
%end SetInputPortSamplingMode

%%
%% Terminate:
%%   Functionality    : Called at the end of simulation for cleanup
%%   Required         : Yes
%%   C MEX counterpart: mdlTerminate
%%
function Terminate(block)

%end Terminate

