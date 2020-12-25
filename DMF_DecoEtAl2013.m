% 
% This code makes a direct simulation of a network for which connectivity is
% given (DTI matrix). It compares the model correlation matrix to an
% empirical one (BOLD functionnal connectivity).
%
% Each node contains a reduced model à la Wang.
% Nodes are connected via a DTI matrix.
%
% Uses: BOLD.m  : Balloon-Windkessel hemodynamic model to convert the
% synaptic gating variables into BOLD signals.
%
% Outputs: 
%    
%       Cmat    : correlation coeff.  (MxMxnumWg)
%
%       where M: nb. of nodes, numWg: nb. of tested global couplings (wg)
%       
%       Coef    : Cmat vs. FC_emp fitting
%
%--------------------------------------------------------------------------

clear all;

% Load SC and FC data:
%---------------------

load Human_66.mat C Order anat_lbls
load Corbetta.mat Cb Pv

C=C(Order,Order);

FC_emp=Cb;

% C : anatomical connectivity.
% FC_emp : functional connectivity.


N=size(C,1);

%--------------

Isubdiag = find(tril(ones(N),-1)); % Values below the diagonal.
fc       = FC_emp(Isubdiag); % Vector of all FC values below the diagonal 

dtt   = 1e-3;  % Sampling rate of simulated neuronal activity (seconds)
ds   = 100;    % BOLD downsampling rate
f_c  = 0.25;   % Frequency cut-off to remove spurious BOLD


% simulation length and binsize:
dt=0.1;
tmax=200000; % WARNING: tmax should be larger than 50000  !
tspan=0:dt:tmax;

% Fixed parameters 
%(see documentation)
%----------------------
tao=100;
gamma=0.641;
sigma=0.001; %0.001;
JN11=0.2609;
I0=0.3;
w=0.9; % local excitatory recurrence



% transfert function:
%---------------------------------------
a=270;
b=108;
d=0.154;
H=@(x) (a*x-b)./(1-exp(-d*(a*x-b)));
%--------------------------------------



wgs=0:0.1:5; % global couplings
numWg=length(wgs);


Coef     =zeros(1,numWg);
Maxrate  =zeros(1,numWg);
neuro_act=zeros(tmax,N);


nTrials=10; % number of trials for each G

ix=1;

for we=wgs;

 display(sprintf('we=%g',we))
 
 
 MaxR=0;
 cb=zeros(length(Isubdiag),1);
  
for tr=1:nTrials
    
    display(tr)
 
y=0.001*ones(N,1);
nn=1;
j=0;
for i=2:1:length(tspan)
 x=w*JN11*y+I0+JN11*we*(C*y);
 r=feval(H,x);
 y=y+dt*(-y/tao+(1-y)*gamma.*r./1000)+sqrt(dt)*sigma*randn(N,1);
 j=j+1;
 if j==(1/dt)
  neuro_act(nn,:)=y';   % Down-sampling
  nn=nn+1;
  j=0;
 end
end

nn=nn-1;


% Friston BALLOON MODEL
%----------------------------------------------
T = nn*dtt; % Total time in seconds

B = BOLD(T,neuro_act(1:nn,1)'); % B=BOLD activity, bf=Foutrier transform, f=frequency range)
BOLD_act = zeros(length(B),N);
BOLD_act(:,1) = B;  

for nnew=2:N
    B = BOLD(T,neuro_act(1:nn,nnew));
    BOLD_act(:,nnew) = B;
end

% Downsampling and reordering removing the first 500ms
bds = BOLD_act(500:ds:end,:);


% BOLD correlation matrix = Simulated Functional Connectivity
% (fisher-transformed):
%-----------------------------------------------------------
Cb  = corrcoef(bds);
cb  = cb + atanh(Cb(Isubdiag))/nTrials; % Vector containing all the FC values below the diagonal 
        
clear BOLD BOLD_act bds

MaxR = MaxR + max(mean(neuro_act(end-10000:end,:)))/nTrials;

end
% FITTING:
%----------------------------
r_c    = corrcoef(cb,fc);
Coef(ix)=r_c(2);

Maxrate(ix)=MaxR;

ix=ix+1;

disp(num2str(Coef(ix-1)));


end;

set(figure,'Position',[100 130 600 400],'Color','w')

plot(wgs,Coef,'ko-','markerfacecolor','k')
xlabel('Coupling strength (G)','FontSize',9)
ylabel('Similarity','FontSize',9)


