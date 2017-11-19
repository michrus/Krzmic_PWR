function [sys,x0,str,ts,simStateCompliance] = reaktor_control(t,x,u,flag)
% SFUNTMPL General MATLAB S-Function Template
%   With MATLAB S-functions, you can define you own ordinary differential
%   equations (ODEs), discrete system equations, and/or just about
%   any type of algorithm to be used within a Simulink block diagram.
%
%   The general form of an MATLAB S-function syntax is:
%       [SYS,X0,STR,TS,SIMSTATECOMPLIANCE] = SFUNC(T,X,U,FLAG,P1,...,Pn)
%
%   What is returned by SFUNC at a given point in time, T, depends on the
%   value of the FLAG, the current state vector, X, and the current
%   input vector, U.
%
%   FLAG   RESULT             DESCRIPTION
%   -----  ------             --------------------------------------------
%   0      [SIZES,X0,STR,TS]  Initialization, return system sizes in SYS,
%                             initial state in X0, state ordering strings
%                             in STR, and sample times in TS.
%   1      DX                 Return continuous state derivatives in SYS.
%   2      DS                 Update discrete states SYS = X(n+1)
%   3      Y                  Return outputs in SYS.
%   4      TNEXT              Return next time hit for variable step sample
%                             time in SYS.
%   5                         Reserved for future (root finding).
%   9      []                 Termination, perform any cleanup SYS=[].
%
%
%   The state vectors, X and X0 consists of continuous states followed
%   by discrete states.
%
%   Optional parameters, P1,...,Pn can be provided to the S-function and
%   used during any FLAG operation.
%
%   When SFUNC is called with FLAG = 0, the following information
%   should be returned:
%
%      SYS(1) = Number of continuous states.
%      SYS(2) = Number of discrete states.
%      SYS(3) = Number of outputs.
%      SYS(4) = Number of inputs.
%               Any of the first four elements in SYS can be specified
%               as -1 indicating that they are dynamically sized. The
%               actual length for all other flags will be equal to the
%               length of the input, U.
%      SYS(5) = Reserved for root finding. Must be zero.
%      SYS(6) = Direct feedthrough flag (1=yes, 0=no). The s-function
%               has direct feedthrough if U is used during the FLAG=3
%               call. Setting this to 0 is akin to making a promise that
%               U will not be used during FLAG=3. If you break the promise
%               then unpredictable results will occur.
%      SYS(7) = Number of sample times. This is the number of rows in TS.
%
%
%      X0     = Initial state conditions or [] if no states.
%
%      STR    = State ordering strings which is generally specified as [].
%
%      TS     = An m-by-2 matrix containing the sample time
%               (period, offset) information. Where m = number of sample
%               times. The ordering of the sample times must be:
%
%               TS = [0      0,      : Continuous sample time.
%                     0      1,      : Continuous, but fixed in minor step
%                                      sample time.
%                     PERIOD OFFSET, : Discrete sample time where
%                                      PERIOD > 0 & OFFSET < PERIOD.
%                     -2     0];     : Variable step discrete sample time
%                                      where FLAG=4 is used to get time of
%                                      next hit.
%
%               There can be more than one sample time providing
%               they are ordered such that they are monotonically
%               increasing. Only the needed sample times should be
%               specified in TS. When specifying more than one
%               sample time, you must check for sample hits explicitly by
%               seeing if
%                  abs(round((T-OFFSET)/PERIOD) - (T-OFFSET)/PERIOD)
%               is within a specified tolerance, generally 1e-8. This
%               tolerance is dependent upon your model's sampling times
%               and simulation time.
%
%               You can also specify that the sample time of the S-function
%               is inherited from the driving block. For functions which
%               change during minor steps, this is done by
%               specifying SYS(7) = 1 and TS = [-1 0]. For functions which
%               are held during minor steps, this is done by specifying
%               SYS(7) = 1 and TS = [-1 1].
%
%      SIMSTATECOMPLIANCE = Specifices how to handle this block when saving and
%                           restoring the complete simulation state of the
%                           model. The allowed values are: 'DefaultSimState',
%                           'HasNoSimState' or 'DisallowSimState'. If this value
%                           is not speficified, then the block's compliance with
%                           simState feature is set to 'UknownSimState'.


%   Copyright 1990-2010 The MathWorks, Inc.

%
% The following outlines the general structure of an S-function.
%
switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes;

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1,
    sys=mdlDerivatives(t,x,u);

  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2,
    sys=mdlUpdate(t,x,u);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,x,u);

  %%%%%%%%%%%%%%%%%%%%%%%
  % GetTimeOfNextVarHit %
  %%%%%%%%%%%%%%%%%%%%%%%
  case 4,
    sys=mdlGetTimeOfNextVarHit(t,x,u);

  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
  case 9,
    sys=mdlTerminate(t,x,u);

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));

end

% end sfuntmpl

%
% =============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
 function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%
% Note that in this example, the values are hard coded.  This is not a
% recommended practice as the characteristics of the block are typically
% defined by the S-function parameters.
%
sizes = simsizes;

sizes.NumContStates  = 21;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 21;
sizes.NumInputs      = 2;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0 = [1.0 4625.71 1488.75 548.303 557.106 1506.38 565.909 575.711 1523.96 583.514 592.5 592.5 592.5 539.5 539.5 592.5 539.5 566 566 0.0 0.0];

% x0(1) = Pcor0
% x0(2) = C0
% x0(3) = Tf10
% x0(4) = Tmo10
% x0(5) = Tmo20
% x0(6) = Tf20
% x0(7) = Tmo30
% x0(8) = Tmo40
% x0(9) = Tf30
% x0(10) = Tmo50
% x0(11) = Tmo60
% x0(12) = Tup
% x0(13) = Thl
% x0(14) = Tlp
% x0(15) = Tcl

%x0(16) = Thlp0
%x0(17) = Tclp0
%x0(18) = Tset0
%x0(19) = Taves0
%x0(20) = AUXco0
%x0(21) = ROHex0


%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [0 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'DisallowSimState' < Error out when saving or restoring the model sim state
simStateCompliance = 'UnknownSimState';

% end mdlInitializeSizes

%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,u)
    
% ############################ Wejscia ####################################
 P = u(1);
 Tpo = u(2);


%% ##################################
%  #  STALE MODELU REAKTORA         #
%  ##################################
BETA1 = 0.000209
BETA2 = 0.001414
BETA3= 0.001309
BETA4 = 0.002727
BETA5 = 0.000925
BETA6 = 0.000314
BETAt = 0.00689
LANDA1 = 0.0125
LANDA2 = 0.0308
LANDA3 = 0.114
LANDA4 = 0.307
LANDA5 = 1.19
LANDA6 = 3.19
NGT = 17.9E-6
ALPHAf =-1.1E-05
ALPHAc = -2.0E-04
Mft = 222739
Cpf = 0.059
FR = 0.974

Vup = 1376
Vlp = 1791
Vhl = 1000
Vcl = 2000
Vmo = 540
DENSmo = 45.71
Cpmo = 1.39

Pcor10 = 3.436e9
Wpth =1.57559e8

Tf10 = 1488.75;
Tmo10 = 548.303;
Tmo20 = 557.106;
Tf20 = 1506.38;
Tmo30 = 565.909;
Tmo40 = 574.711;
Tf30 = 1523.96;
Tmo50 = 583.514;
Tmo60 = 592.5;

React=0.003375;
TOU1a1=10;
TOU1a2=5;
TOUrtd=4;
POWERi=691244.4199;
TOUset=30;
X8=1225.85;
K9=-0.035;
C1=1.2195;
Tfi=434.3;
CPfw=1.2181;
TOU1e=80;

% ###############        Koniec stalych reaktora     ###############

%% #########################################################
%  #  Stale parametry modelu reaktora - Obliczenia         #
%  #########################################################

Pcor1 = (Pcor10*3.41)/(3600.0*3.0)
Hfm = 200.0/3600.0
Afm = 59900.0/3.0
Wpt = Wpth/3600.0

LANDA = BETAt/( BETA1/LANDA1 + BETA2/LANDA2 + BETA3/LANDA3 + BETA4/LANDA4 + BETA5/LANDA5 + BETA6/LANDA6)

Mf=Mft/3.0
Mmo=Vmo*DENSmo/3.0
Mup=Vup*DENSmo
Mlp=Vlp*DENSmo
Mhl=Vhl*DENSmo
Mcl=Vcl*DENSmo

TOUmo=Mmo/(2.0*Wpt)
TOUup=Mup/Wpt
TOUlp=Mlp/Wpt
TOUhl=Mhl /Wpt
TOUcl=Mcl /Wpt
%############## Koniec stalych reaktora - obliczonych ##############

%% #########################################
%  #  Rownania rozniczkowe reaktora        #
%  #########################################

Pcor = x(1); 
C = x(2); 
Tf1 = x(3); 
Tmo1 = x(4);
Tmo2 = x(5);
Tf2 = x(6); 
Tmo3 = x(7);
Tmo4 = x(8);
Tf3 = x(9); 
Tmo5 = x(10);
Tmo6 = x(11);
Tup = x(12); 
Thl = x(13); 
Tlp = x(14); 
Tcl = x(15); 

Thlp = x(16);
Tclp = x(17);
Tset = x(18);
Taves = x(19);
AUXco = x(20);
ROHex = x(21);

ROH=ROHex+ALPHAf*((Tf1-Tf10)+(Tf2-Tf20)+(Tf3-Tf30))/3 + ALPHAc*((Tmo1-Tmo10)+(Tmo2-Tmo20)+(Tmo3-Tmo30)+(Tmo4-Tmo40)+(Tmo5-Tmo50)+(Tmo6-Tmo60))/6

DPcor=(ROH - BETAt)*Pcor/NGT + LANDA*C
xdot(1) = DPcor;

DC=BETAt*Pcor/NGT - LANDA*C
xdot(2) = DC;

DTf1=(FR*Pcor1*Pcor)/(Mf*Cpf)+(Hfm*Afm)*(Tmo1-Tf1)/(Mf*Cpf)
xdot(3) = DTf1;

DTmo1=((1-FR)*Pcor1*Pcor)/(Mmo*Cpmo)+(Hfm*Afm)*(Tf1-Tmo1)/(Mmo*Cpmo)+(Tlp-Tmo1)/TOUmo
xdot(4) = DTmo1;

DTmo2=((1-FR)*Pcor1*Pcor)/( Mmo* Cpmo )+( Hfm*Afm ) * ( Tf1 -Tmo1 )/( Mmo *Cpmo ) + ( Tmo1 - Tmo2 ) /TOUmo
xdot(5) = DTmo2;

DTf2=(FR*Pcor1 *Pcor ) / (Mf * Cpf ) + ( Hfm* Afm) * (Tmo3 - Tf2 ) / (Mf * Cpf)
xdot(6) = DTf2;

DTmo3 = ( ( 1-FR) * Pcor1 * Pcor ) / (Mmo *Cpmo ) + (Hfm*Afm ) * (Tf2-Tmo3 ) /(Mmo*Cpmo ) + (Tmo2 - Tmo3) /TOUmo
xdot(7) = DTmo3;

DTmo4 = ( ( 1 -FR) *Pcor1 * Pcor ) / (Mmo *Cpmo ) + ( Hfm*Afm ) * ( Tf2-Tmo3)/(Mmo *Cpmo) + (Tmo3 - Tmo4)/TOUmo
xdot(8) = DTmo4;

DTf3 = ( FR* Pcor1 * Pcor ) / (Mf * Cpf ) + ( Hfm*Afm) * (Tmo5 - Tf3 ) / (Mf *Cpf)
xdot(9) = DTf3;

DTmo5 = ( ( 1 -FR) *Pcor1 * Pcor) / (Mmo * Cpmo ) + ( Hfm*Afm ) * (Tf3-Tmo5 )/( Mmo *Cpmo ) + (Tmo4 - Tmo5 ) /TOUmo
xdot(10) = DTmo5;

DTmo6 = ( ( 1 -FR) * Pcor1 * Pcor ) / (Mmo *Cpmo ) + ( Hfm*Afm) * (Tf3-Tmo5 ) /( Mmo*Cpmo ) + (Tmo5 - Tmo6 ) /TOUmo
xdot(11) = DTmo6;

DTup=(Tmo6-Tup )/TOUup
xdot(12) = DTup;

DThl=(Tup-Thl)/TOUhl
xdot(13) = DThl;

DTlp = (Tcl - Tlp)/TOUlp
xdot(14) = DTlp;

DTcl = (Tpo - Tcl)/TOUcl
xdot(15) = DTcl;

Tave = (Tcl+Thl)/2.0;

DThlp = (Thl-Thlp)/TOUrtd;
xdot(16) = DThlp;

DTclp = (Tcl - Tclp)/TOUrtd;
xdot(17) = DTclp;

Tavep = (Thlp+Tclp)/2;
Hst = X8 + K9*P;
POWERs = P*C1*(Hst-CPfw*Tfi);

DTset = (0.19*(POWERs*100.0/POWERi)+547.0-Tset)/TOUset;
xdot(18) = DTset;

DTaves = AUXco;
xdot(19) = DTaves;

DAUXco = (Tavep-Taves-(TOU1a1+TOU1a2)*AUXco+TOU1e*(DTclp+DThlp)/2.0)/(TOU1a1*TOU1a2);
xdot(20) = DAUXco;

e = Tset - Taves;

if e >= 5
    DROHex=72*React*0.00689/60.0;
elseif e >= 3
    DROHex=React*0.00689*(32*(Tset-Taves)-88)/60.0;
elseif e >= 1
    DROHex=8*React*0.00689/60.0;
elseif abs(e) < 1
    DROHex = 0;
elseif e <= -1
    DROHex=-8*React*0.00689/60.0;
elseif e <= -3
    DROHex=React*0.00689*(32*(Tset-Taves)+88)/60.0;
else
    DROHex=-72*React*0.00689/60.0;
end

xdot(21) = DROHex;

%############  Koniec rownan rozniczkowych reaktora   ############
sys = xdot;

%%
% end mdlDerivatives

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,x,u)

sys = [];

% end mdlUpdate

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)

sys = x;
% end mdlOutputs

%
%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.  Note that the result is
% absolute time.  Note that this function is only used when you specify a
% variable discrete-time sample time [-2 0] in the sample time array in
% mdlInitializeSizes.
%=============================================================================
%
function sys=mdlGetTimeOfNextVarHit(t,x,u)

sampleTime = 1;    %  Example, set the next hit to be one second later.
sys = t + sampleTime;

% end mdlGetTimeOfNextVarHit

%
%=============================================================================
% mdlTerminate
% Perform any end of simulation tasks.
%=============================================================================
%
function sys=mdlTerminate(t,x,u)

sys = [];

% end mdlTerminate
