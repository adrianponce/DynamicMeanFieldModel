%
%
% This prog. simulates the FIC model + associated BOLD time series,
% and calculates the fitting between the empirical FC and
% the model prediction for varying global couplings (G)
%
% For comparison with the empirical data, we considered the FC of simulated
% BOLD signals which are obtained by transforming the model synaptic activity 
% through a hemodynamic model. 
%
% The feedback inhition weights need to be calculated previously
% using Get_balanced_weights.m
%
% It uses: BOLD_HemModel.m
%
% see:
% Deco et al. (2014) J Neurosci.
% http://www.jneurosci.org/content/34/23/7886.long
% Ponce-Alvarez et al. (2015)
% http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.100
% 4445
%
% Adrian Ponce-Alvarez
%--------------------------------------------------------------------------


clear all;

% HERE --> Load C = the DTI matrix 
% N = size(C,1);
% HERE --> Load FC_emp = the emprirical FC

load Human_66.mat C FC_emp

N = size(C,1);



Isubdiag = find(tril(ones(N),-1)); % Indexes of all the values below the diagonal.
fc       = atanh(FC_emp(Isubdiag)); % Vector containing all the FC values below the diagonal, z-transform 
dsb   = 100;    % BOLD downsampling rate
dtt   = 1e-3;   % Sampling rate of simulated neuronal activity (seconds)
npairs = length(Isubdiag);


% Parameters:
%---------------
dt=.1;
tmax=100000;
tspan=0:dt:tmax;
ds=10; % downsampling stepsize
Tds=tmax/(ds*dt);
res=ds*dt/1000;
T = (Tds-1)*ds*dt/1000; % Total time in seconds (for bold model)


w=1.4; % local recurrence
tau_exc=100;
tau_inh=10;
gamma=0.641;
JN=0.15;
I0=0.382;
Jexte=1;
Jexti=0.7;
I_exc = I0*Jexte;
I_inh = I0*Jexti;
Io = [I_exc*ones(N,1); I_inh*ones(N,1)];
beta=0.001;  % additive ("finite size") noise

% number of stochastic realizations:
nTrials=5;


% transfer functions:
% transfert function: excitatory
%--------------------------------------------------------------------------
ae=310;
be=125;
de=0.16;
He=@(x) (ae*x-be)./(1-exp(-de*(ae*x-be)));

% transfert function: inhibitory
%--------------------------------------------------------------------------
ai=615;
bi=177;
di=0.087;
Hi=@(x) (ai*x-bi)./(1-exp(-di*(ai*x-bi)));


% load the Ji weights (previously calculated with Get_balanced_weights.m):
%--------------------------------------------------------------------------
load Bifurcations_BalancedModel_stochastic wes JI Se_init Si_init

numG=length(wes);

% initialization:
fittcorr = zeros(1,numG);
Steady_FR=zeros(N,numG);
neuro_act=zeros(Tds,N);
neuro_act2=zeros(Tds,N);
FC_z = zeros(npairs,numG);

for kk=1:numG 

    we=wes(kk);

    display(sprintf('Global coupling = %g',we))

 % feedback inhibition weights:   
 J=JI(:,kk);
 mu0=[Se_init(:,kk);Si_init(:,kk)];
 Ut=zeros(N,Tds);
 Rt=zeros(N,Tds);


       % Coupling matrix:
      %----------------------------------
      W11 = JN*we*C + w*JN*eye(N);
      W12 = diag(-J);
      W21 = JN*eye(N);
      W22 = -eye(N);
      Wmat = [W11 W12;W21 W22];
      

 cb = zeros(length(Isubdiag),1);

 for tr=1:nTrials
 
     display(sprintf('   trial:%g',tr))
     mu=mu0;

 
     % Warm-Up to reach stationarity:
     %--------------------------------

      for t=2:1000         
       u = Wmat*mu + Io;
       re = feval(He,u(1:N));
       re = gamma*re./1000;
       ri = feval(Hi,u(N+1:end));
       ri = ri./1000;       
       ke = -mu(1:N)/tau_exc+(1-mu(1:N)).*re;
       ki = -mu(N+1:end)/tau_inh + ri;     
       kei = [ke;ki];     
       mu = mu + dt*kei + sqrt(dt)*beta*randn(2*N,1);      
       mu(mu>1)=1;
       mu(mu<0)=0;
      end
      
      % Simulate dynamical model:
      % -------------------------
  
      nn=0;
      for t=2:length(tspan)          
       u = Wmat*mu + Io;
       re = feval(He,u(1:N));
       re = gamma*re./1000;
       ri = feval(Hi,u(N+1:end));
       ri = ri./1000;       
       ke = -mu(1:N)/tau_exc+(1-mu(1:N)).*re;
       ki = -mu(N+1:end)/tau_inh + ri;     
       kei = [ke;ki];     
       mu = mu + dt*kei + sqrt(dt)*beta*randn(2*N,1);      
       mu(mu>1)=1;
       mu(mu<0)=0;    
      
          if mod(t,ds)==0
          nn=nn+1;    
          Ut(:,nn)=u(1:N); %excitatory synaptic activity
          Rt(:,nn)=re;     %excitatory firing rate
          end
          
      end


%%%% BOLD Model
% Friston BALLOON MODEL
%--------------------------------------------------------------------------
B = BOLD_HemModel(T,Ut(1,:),res); % B=BOLD activity
Tnew=length(B);
BOLD_act = zeros(Tnew,N);
BOLD_act(:,1) = B;  
for i=2:N
    B = BOLD_HemModel(T,Ut(i,:),res);
    BOLD_act(:,i) = B;
end
% Downsampling
bds=BOLD_act(500:dsb:end,:);
% BOLD correlation matrix = Simulated Functional Connectivity
Cb  = corrcoef(bds);
cb  = cb + atanh(Cb(Isubdiag))/nTrials; % Vector containing all the FC values below the diagonal 
% Firing rate:
meanrate=mean(Rt,2)*1000/gamma;
Steady_FR(:,kk) = Steady_FR(:,kk) + meanrate/nTrials;


 end  %end loop over trials
 
 
Coef    = corrcoef(cb,fc);
fittcorr(kk)=Coef(2);
FC_z(:,kk)=cb;

end

figure
plot(wes,fittcorr);

figure
plot(wes,max(Steady_FR),'.');


