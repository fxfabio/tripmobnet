function L=MC_HiddenVariable_linkingP(Fd,alphabeta)
%       alfabeta=[alpha beta] - of the linking function
% this code allows self-loop (removed in other code in the evolution)
%%@F.Vanni2023 

n=length(Fd);
%% LINKING FUNCTION
alpha=alphabeta(1);
beta=alphabeta(2);
p=(1+alpha)*(1+beta) ; % random network connecting value
p=1;
A_sim=Fd(:); % force to be column vector 

Fx=Fd(:,1);
Fy=transpose(Fd(:,2));


linkF=(1./p).*(Fx.^alpha * Fy.^beta) ./ (max(Fx).^(alpha).*max(Fy).^(beta)) ;
c=1e3;
%linkF=(1./p).*(exp(Fx.*alpha./c)*exp(Fy.*beta./c))./( max(exp(Fx.*alpha./c).*max(Fy.*beta./c)) );
%C=400;
%TRE=double((A_sim+A_sim'-C)>=0);

L=linkF;

%L=TRE;
%L=1+0.*TRE;

%linkF=erf( (A_sim.^alpha * A_sim'.^beta)./(max(A_sim).^(alpha+beta)) );
%L=linkF;

%linkF=((A_sim.^alpha * A_sim'.^beta) ./ (max(A_sim).^(alpha+beta)) )./ (1+(A_sim.^alpha * A_sim'.^beta)./ (max(A_sim).^(alpha+beta))) ;

x=A_sim./max(A_sim);
y=A_sim'./max(A_sim);

xx=A_sim;
yy=A_sim';
c=0.0001;

%linkF=xx*yy./(1+c.*xx*yy) ;


%linkF= (x*y./(1+x*y));

%C=1;
%linkF=1./(1+exp(-(xx+yy-C)));
%L=linkF;

