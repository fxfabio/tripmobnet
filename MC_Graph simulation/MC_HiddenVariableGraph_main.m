function [A,Fxy]=MC_HiddenVariableGraph_main(N,T,multi)
% Monte Carlo simulation of hidden variable mobility multi-graph model
% with multilink procedure or not in string 'multi'.
% @FVanni Gennaio 2023

%AD= 10; % average degree
%N=300; % number of nodes

%% INPUT parameters
alpha=1; % as in-visit exponent
beta=1; % as out-trip exponent
alphabeta=[alpha beta];

mu_x=2.001;
a=0.001;
b=N-1;
%%
 if nargin<3; multi='no', gx=0; 
 else
     if strcmpi(multi,'yes')==1 || strcmpi(multi,'y')==1 || strcmpi(multi,'si')==1
         gx=1;
     elseif strcmpi(multi,'no')==1 || strcmpi(multi,'n')==1 || strcmpi(multi,'not')==1
          gx=0;
     else disp('error; not valid string on multi')
     end
 end

%% generation of the network

% variable distribution for the locations
xmin=a;
xmax=b;
%F=MC_variable_distributions(N,'xmin',xmin,'truncated',mu,xmax);
%F=MC_variable_distributions(N,'xmin',xmin,'powerlaw',mu);
Fy=MC_variable_distributions(N,'exponential',1/1000);
%Fy=MC_variable_distributions(N,'xmin',xmin,'powerlaw',mu_x);
Fx=MC_variable_distributions(N,'xmin',xmin,'powerlaw',mu_x);
%F = unifrnd(xmin,xmax,[N,1]);
%F = lognrnd(0.01,0.001,[N,1]);
%F=MC_variable_distributions(n,'xmin',xmin,'powercutoff',mu,lambda);

Fxy=[Fx(:),Fy(:)];
% generation of the interbank Linking Function
linkF=MC_HiddenVariable_linkingP(Fxy,alphabeta);

%% for loops
n=length(Fx);

ad=T./n;

totlink=round(ad*n);
if totlink>=n*(n-1), warning('warning: multi-link procedure activated');end
t=1;
%linkFF=linkF-diag(diag(linkF));
linkFF=linkF; % self-loop included
L=linkFF(:);
Lcut=1-L;
lin_ind=zeros(size(L));

   tic  

while t<=totlink
    cL=[0; cumsum(L)];
    ra =  cL(end)* rand(1,1) ;
    %[k,~]=find( cL >= ra );
    %m=k(1);
    m=find(cL>=ra,1);
    %lin_ind(m-1)=1; %'cause cumsum added a 0
    lin_ind(m-1)= lin_ind(m-1)+1;% link multipli
    L(m-1)=gx*L(m-1); % una volta creato il link non si puo ripetere o si
    t=t+1;
end
    toc

A=(reshape(lin_ind,size(linkF)));
% so that:
% sum(A,1); % is the out-Degree distribution - Departures
% sum(A,2); % is the  in-Degree distribution - Visits