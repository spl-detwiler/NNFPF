% Ellis_Artificial_Removal
% 
% Solver of Ellis-fluid flow through a fracture with sand pack. 
% The code solves for the steady state velocity distribution, once the
% velocity distribution is solved, the shear rate everywhere is calculated.
% At every pixel-location where the mobilization condition is met (e.g., 
% phi<phi_mob or shearRate > shearRate_mob, the solids are artificially 
% mobilized by setting phi=0. The process is repeated until the number of 
% iterations (erodeITER) is met.
%
%
% INPUT			DESCRIPTION
% phiNAME           = concentration file name
% AP_LIST           = aperture file name
% ROI_1             = Region of interest [Row_start Row_end Col_start Col_end]
%                     If solving for whole field leave blank, e.g., ROI=[]
%
% hin               = head on the left boundary
% shearRate_mob     = shear rate threshold (for solid removal)
% phi_mob           = phi threshold (for solid removal)
%
% maxITER           = number of iterations for flow solver (h-field)
% erodeITER         = number of iterations to erode
%
%
%~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
%
% Copyright (c) 2018 Ricardo Medina and Russell Detwiler
% Contact:   University of California, Irvine, Irvine, CA
% E-mail: ricarm3@uci.edu  and  detwiler@uci.edu
% Citation: Medina, R., R.L. Detwiler, R. Prioul, W. Xu, and J.E. Elkhoury (2018), Settling and mobilization of sand-fiber proppants ina deformable fracture, Water Reources Research. 
%
%~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
%%
clear all; close all; clc;

%Set path to location of files:
path1 = '/Set/Local/Path/to/file/locations/'; %pwd


%~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% General Settings
%~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

% Concentration and Aperture Field (File names): Change the Name of both
% files
phiNAME = 'PHI_TestA_S_0.0kPa_EndofInj.fits'; % Name of $\phi$-field
AP_LIST = 'AP_TestA_S_0.0kPa_EndofInj.fits'; % Name of corresponding Aperture field
ROI =[];  % [RowSTART RowEND ColSTART ColEND]


%~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Simulation Settings
%~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
maxITER   =200; % number of iterations for flow solver (h-field)
erodeITER =100; % number of iterations to erode

% Mobilization Threshold
phi_mob = 0.55; % [-] mobilization concentration
shearRate_mob = 5;%[1/s] mobilization shear rate


hin=0.5;  % head at lhs (inlet)
hout=0; % head at rhs (outlet)



%~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Fluid Properties (ELLIS FLUID) - Constants
%~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
grav = 9.81; % [m/s2] : gravity
rho =  1180; % [kg/m3] : fluid density
gamma=rho*grav; % [Pa/m]=[kg/m^2/s^2] : specific weight of fluid

tau12 = 4.33;       % [Pa] : half-shear stress
mu0=2.4;            % [Pa.s] : zero-shear viscosity
alpha = 2.77;       % [-] : shear thinning index


%~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% GEOMETRY:
%~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Pixel size
dx=86e-6; % [m] : Assumes dy=dx

% FRACTURE APERTURE
ap = fitsread([path1 apNAME])*1e-6; %[m] convert aperture from microns to meters
ap = ap/2; %half-aperture

% Set ROI for aperture (if necessary)
if ~isempty(ROI)
    ap=ap(ROI(1):ROI(2), ROI(3):ROI(4));
end

% pad aperture array with one extra cell around entire domain
ap=padarray(ap,[1 1],'symmetric');

% size of domain (padded domain)
nx=size(ap,1); % # of rows
ny=size(ap,2); % # of columns

% set ap=0 in 'ghost' cells adjacent to no flow boundaries
ap(1,:)=0; ap(nx,:)=0;

% POROUS 'GEOMETRY'
a  = 150e-6; % Particle radius [m]
dp = a*2; % Particle diameter [m]

% PHI-FIELD
phi = fitsread([path1 phiNAME]);

% Set ROI for phi (if necessary)
if ~isempty(ROI)
    phi=phi(ROI(1):ROI(2), ROI(3):ROI(4));
end

% Set boundaries on phi:
phiMIN = 0;
phi(phi<=.005)=phiMIN; % for stability pusposes
phiMAX = 0.75;
phi(phi>=phiMAX)=phiMAX; % for stability pusposes


% Initialize vectors
b=zeros(nx-2,ny-2); % right-hand-side vector
h=zeros(nx-2,ny-2).*linspace( hin-(0.05*hin),hout+(0.05*hout),ny-2 ); % Init head field: linearly decreasing [used for calculating rhs]

for ITER=1:erodeITER %REMOVAL ITERATIONS
    
    % VOID FRACTION
    e = 1-phi;
    e = padarray(e,[1 1],'symmetric');
    
    % KOZENY-CARMAN - K(phi)
    k_CK = (dp^2/180) .* ( (e.^3)./( (1-e).^2 ));%Kozeny-Carman
    
    %~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    % SYSTEM PARAMETERS
    %~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    % FRACTURE CONSTANTS
    A = -gamma*ap.^2/3/mu0;% linear term
    C = 3/(alpha+2)*(gamma*ap./tau12).^(alpha-1); % non-linear term
    
    % SOLID MATRIX CONSTANTS (Effective stress formulation)
    D = -gamma.*k_CK/mu0; % linear term
    E = abs(( ((k_CK.*e).^.5)*gamma)./tau12 ).^(alpha-1);% non-linear term
    
    % GEOMETRIC AVERAGE (MIXED FRACTURE-MATRIX CONSTANTS)
    AD = A./(1+A./D);
    CE = C./(1+C./E);
    
    
    % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    % NON-NEWTONIAN FLOW SOLVER - PICARD SOLVER
    % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    
    % Initialize vectors
    %b=zeros(nx-2,ny-2); % right-hand-side vector
    % Initial head field: linearly decreasing (faster convergence)[used for calculating rhs]
    % h=zeros(nx-2,ny-2);%.*linspace(hin-0.05,hout+0.05,ny-2);
    
    ch_bound=1:nx-2;        % indices of constant head boundaries
    
    % Build coefficient matrix [A] and RHS vector {b}
    i=2:nx-1; j=2:ny-1;
    
    % off-diagonal terms
    % harmonic average of linear term at each face already set in ip,jp,in,jn
    ip = 2./(1./AD(i+1,j) + 1./AD(i,j));
    in = 2./(1./AD(i-1,j) + 1./AD(i,j));
    jp = 2./(1./AD(i,j+1) + 1./AD(i,j));
    jn = 2./(1./AD(i,j-1) + 1./AD(i,j));
    
    % harmonic average of non-linear term
    Bipa = 2./(1./CE(i+1,j) + 1./CE(i,j));
    Bina = 2./(1./CE(i-1,j) + 1./CE(i,j));
    Bjpa = 2./(1./CE(i,j+1) + 1./CE(i,j));
    Bjna = 2./(1./CE(i,j-1) + 1./CE(i,j));
    
    % Picard-Factor
    w=0.7;
    for k = 1:maxITER
        
        % Add boundary cells to head array and include BCs
        hwb=padarray(h,[1 1],'symmetric');
        hwb(:,1)=hin; hwb(:,ny)=hout;
        
        % head gradient across each face^(1-alpha)
        hipa = abs((hwb(i+1,j)-hwb(i,j))/dx).^(alpha-1);
        hina = abs((hwb(i-1,j)-hwb(i,j))/dx).^(alpha-1);
        hjpa = abs((hwb(i,j+1)-hwb(i,j))/dx).^(alpha-1);
        hjna = abs((hwb(i,j-1)-hwb(i,j))/dx).^(alpha-1);
        
        % off-diagonal terms for residual and Jacobian
        iip = ip.*(1+Bipa.*hipa);
        iin = in.*(1+Bina.*hina);
        jjp = jp.*(1+Bjpa.*hjpa);
        jjn = jn.*(1+Bjna.*hjna);
        
        % diagonal terms
        diag = -iip - iin - jjp - jjn;
        
        % modify right hand side for constant head B.C.s on left and right
        b(ch_bound,1)=-jjn(ch_bound,1).*hwb(ch_bound+1,1);
        b(ch_bound,ny-2)=-jjp(ch_bound,ny-2).*hwb(ch_bound+1,ny);
        
        % assemble equations
        B=[jjn(:) iin(:) diag(:) iip(:) jjp(:)];
        d=[(nx-2),1,0,-1,-(nx-2)];
        neq=(nx-2)*(ny-2);
        aellis=spdiags(B,d,neq,neq)';
        
        % solve for head field
        hnew = aellis\b(:);
        hnew = reshape(hnew,nx-2,ny-2);
        
        % Calculate and store residuals
        RES(k)= sum((hnew(:)-h(:)).^2);
        
        % UPDATE head field (Picard Iteration)
        h = w*hnew + (1-w)*h;
        
    end
    
    
    % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    % NON-NEWTONIAN VELOCITY CALCULATION
    % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    
    % Pad head field (need to calculate flux
    hh=padarray(h,[1 1],'symmetric');
    hh(:,1)=hin; hh(:,ny)=hout;
    
    % qx from left(L) and right(R)
    qxL = HMean(AD,i,j,i,j-1).*(hh(i,j-1)-hh(i,j))./dx.*( 1+HMean(CE,i,j,i,j-1).*abs((hh(i,j-1)-hh(i,j))./dx).^(alpha-1) );
    qxR = HMean(AD,i,j,i,j+1).*(hh(i,j)-hh(i,j+1))./dx.*( 1+HMean(CE,i,j,i,j+1).*abs((hh(i,j)-hh(i,j+1))./dx).^(alpha-1) );
    
    % qy from top(T) and bottom(B)
    qyT = HMean(AD,i,j,i-1,j).*(hh(i-1,j)-hh(i,j))./dx.*(1+HMean(CE,i,j,i-1,j).*abs((hh(i-1,j)-hh(i,j))./dx).^(alpha-1));
    qyB = HMean(AD,i,j,i+1,j).*(hh(i,j)-hh(i+1,j))./dx.*(1+HMean(CE,i,j,i+1,j).*abs((hh(i,j)-hh(i+1,j))./dx).^(alpha-1));
    
    % Average qx and qy components at each node
    qx = (qxL + qxR)/2;
    qy = (qyT + qyB)/2;
    
    % Velocity components
    ux = qx./e(i,j);
    uy = qy./e(i,j);
    
    % VELOCITY MAGNITUDE
    qmag = (qx.*qx + qy.*qy).^0.5;
    
    % Calculate shear-rate
    gx = (qxR-qxL)./e(i,j)./dx; % du/dx
    gy = (qyB-qyT)./e(i,j)./dx; % dv/dx
    g  = sqrt(gx.*gx + gy.*gy); % magnitude of shear rate
    
    
    % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    % ARTIFICIAL REMOVAL
    % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    
    % SET removal criteria:
    gLIM_rem  = shearRate_mob;        % [1/m] shear rate threshold
    phiLIM_rem= phi_mob;              % [-]   concentration theshold
    
    % FIND locations where solids meet removal criterion:
    clear idxREM
    idxREM = (g(:)>=gLIM_rem & phi(:)<=phiLIM_rem);
    phi0=phi;
    
    % Remove solids at those locations
    phi(idxREM)=phiMIN;
    
    % SAVE FILE (Every Eroded iteration
    SAVENAME = [phiNAME(5:18) sprintf('_hin_%1.2f_SR%.2f_Phi%.2f_IT%04g.mat',hin,gLIM_rem,phiLIM_rem,ITER)];
    save(SAVENAME,'phi0','phi','h','qx','qy','RES');
    
end

