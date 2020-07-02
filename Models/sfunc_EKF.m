function sfunc_EKF(block)
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
% parameter 2 = f ( dx = f(x,u) Continuous state function handle)
% parameter 3 = g ( y = g(x,u) Continuous measurement function handle)
% parameter 4 = Q (Process Noise/Model Uncertainty)
% parameter 5 = R (Measurement Noise Uncertainty)
% parameter 6 = x0 (Estimate of initial state)
% parameter 7 = P0 (Initial estimate uncertainty)
% parameter 8 = u0 (Estimate of initial input)
% Ts,f,g,Q,R,x0,P0,u0
block.NumDialogPrms = 8;

% Read dialog parameters
Ts = block.DialogPrm(1).Data;
g = block.DialogPrm(3).Data;
x0 = block.DialogPrm(6).Data;
u0 = block.DialogPrm(8).Data;

y0 = g(x0,u0);
nx = length(x0);
nu = length(u0);
ny = length(y0);

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
block.SampleTimes = [Ts 0]; % Set sample time

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
%block.RegBlockMethod('jaccsd', @jaccsd);  
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
    
    % Read dialog parameters
    x0 = block.DialogPrm(6).Data;
    nx = length(x0);
    
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
    % Read dialog parameters
    Ts = block.DialogPrm(1).Data;
    f = block.DialogPrm(2).Data;
    Q = block.DialogPrm(4).Data;
    x0 = block.DialogPrm(6).Data;
    P0 = block.DialogPrm(7).Data;
    u0 = block.DialogPrm(8).Data;

    % Dimensions
    nx = length(x0);
    
    x_hat = x0;
    P = P0;
    u = u0;
    
    % Linearise by calculating Jacobians
    F = jaccsd(block,f,x_hat,u); % Calculate Jacobian of continuous system
   
    % Extrapolate
    % ??? Change this if f is dependent on time.
    x_hat = x_hat + f(x_hat,u)*Ts; % Numeric integration to extrapolate state

    Phi = eye(nx) + F*Ts + 1/2*(F*Ts)^2; % ??? where is this from? 2nd order Taylor expansion? (continuous to discrete)

    P = Phi*P*Phi' + Q; % Extrapolate uncertainty

%     % Enforce positive parameters constraint
%     min_param = 0.01; % Min allowable param value
%     if nx > 4 % If paramters are estimated
%         for i = 5:nx 
%             x_hat(i) = max(x_hat(i), min_param);
%         end
%     end
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
    % Read dialog parameters
    Ts = block.DialogPrm(1).Data;
    f = block.DialogPrm(2).Data;
    g = block.DialogPrm(3).Data;
    Q = block.DialogPrm(4).Data;
    R = block.DialogPrm(5).Data;
    x0 = block.DialogPrm(6).Data;
    u0 = block.DialogPrm(8).Data;
 
    nx = length(x0);
    
    % Input
    u = block.InputPort(1).Data;
    y = block.InputPort(2).Data;
    
    % Dwork data
    x_hat = block.Dwork(1).Data; % Priori estimation 
    P = reshape(block.Dwork(2).Data, nx, nx);
    
    % Update step
    H = jaccsd(block,g,x_hat,u); % Linearise measurement function
    K = (P*H')/(H*P*H' + R); % Compute Kalman gain (b*inv(A) -> b/A)
    x_hat = x_hat + K*(y - H*x_hat); % Posteriori / With measurement
    KH_term = (eye(nx) - K*H);
    P = KH_term*P*KH_term' + K*R*K'; % Update estimate uncertainty
   
    % Output
    block.OutputPort(1).Data = x_hat;
    
    % Extrapulate/Prediction for next time step
    x_hat = x_hat + f(x_hat,u)*Ts; % Numeric integration (extrapolate state)
    % ??? Needs to change to proper num integration if dependent on time
    
    F = jaccsd(block,f,x_hat,0); % Calculate Jacobian of continuous system
    Phi = eye(nx) + F*Ts + 0.5*(F*Ts)^2; % 2nd order Taylor expansion (continuous to discrete)
    P = Phi*P*Phi' + Q; % Extrapolate uncertainty

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

function J = jaccsd(block, f,x,u) % ??? Maybe should use simbolic diff for more exact
    % JACCSD Jacobian through complex step differentiation
    % By Yi Cao at Cranfield University, 02/01/2008
    % [z J] = jaccsd(f,x)
    % z = f(x)
    % J = f'(x)
    %
    f_x = f(x,u);
    n = numel(x);
    m = numel(f_x);
    J = zeros(m,n);
    h = n*eps;
    for k=1:n
        x1 = x;
        x1(k) = x1(k)+ h*1i;
        J(:,k) = imag(f(x1,u))/h;
    end

function Terminate(block)

%end Terminate

