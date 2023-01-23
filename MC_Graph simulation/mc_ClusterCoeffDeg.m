function [k_degree , C_cluster,C_clusterx]=mc_ClusterCoeffDeg(A,binnum,author)
% Cluster coefficient spectrum vs degree,for weighted and directed mobility network
% Input: 
%       A= weighted adjacency matrix
%       binum = number of bins of the degrees (resolution of the spectrum)
%       author = the type of cluster coefficient definitions, 3 string:
%            'Fagiolo' - In-CC as in Fagiolo (2007),doi:10.1103/PhysRevE.76.026107
%            'Clemente'- In-CC as in Clemente-Grassi (2018),doi:10.1016/j.chaos.2017.12.007
%            'Fardet'  - FanIn-CC as in Fardet-Levina (2021),doi:10.1103/PhysRevResearch.3.043124
%
%   Note:   All weights must be between 0 and 1.
% Output:
%       k_degree= binned strength (degree) grouped in binnum number of bins
%       C_cluster= Cluster coeeficient for nodes of a certain strength (normalized to 1) 
%  
%@F.Vanni2023 
if nargin<3
    author='Fagiolo';
end

W=transpose(A)./max(abs(A(:)));    % normalized scale by maximal weight
%W=A./max(sum(A,2));

% (CRin1)Fagiolo 2007 o (CRin2) Clemente 2017 o (CRin3)  Fardet 2021 
[Ctot,Cin,Kin,CRin1,CRin2,CRin3]=mc_clustering_coef(W); % Fagiolo 2007 o Clemente 2017

if  strcmpi(author,'Fagiolo'),CCin=CRin1;
elseif strcmpi(author,'Clemente'),CCin=CRin2;
elseif strcmpi(author,'Fardet'),CCin=CRin3;
else, CCin=Ctot; disp('total Cluster coefficient');
end

CCin(~isfinite(CCin))=NaN;
C=Cin;
%k_degree= sum(A,2)+sum(A,1)'; % total degree

k_degree= sum(A,2); % in degree
%k_degree= Kin;
max(k_degree)

x=unique(k_degree);

k_bin=linspace(min(x),max(x),binnum);

   cc=zeros(length(k_bin),1); 
   ccx=zeros(length(k_bin),1); 
   
for i=1:length(k_bin)-1
    %ix=(k_degree==x(i));
    ix=(k_degree>=k_bin(i) & k_degree<k_bin(i+1) );
    cc(i)=nansum(C(ix));
    ccx(i)=nansum(CCin(ix));
end

k_degree=k_bin;
C_cluster=cc';
C_cluster(C_cluster==0)=NaN;

C_clusterx=ccx';


C_clusterx(C_clusterx==0)=NaN; % excluding node with clustering zero 
C_clusterx(~isfinite(C_clusterx))=NaN;


function [C,Cin,Kin,CR_in1,CR_in2,CR_in3]=mc_clustering_coef(W)
%CLUSTERING_COEF_WD     Fan-in Clustering coefficient
%
%   The weighted clustering coefficient is the average "intensity"
%   (geometric mean) of all triangles associated with each node.
%
%   Input:      W,      weighted directed connection matrix
%                       (all weights must be between 0 and 1)
%
%   Output:     C,      clustering coefficient vector
%
%
%   2022: modification of https://rdrr.io/cran/DirectedClustering/src/R/ClustF.R
%   The weighted modification is as follows:
%   - The numerator: adjacency matrix is replaced with weights matrix ^ 1/3
%   - The denominator: no changes from the binary version


A=W~=0;                     %binary djacency matrix
K=sum(A+A',2);            	%total degree (in + out)
Kin=sum(A,2);               % in-degree

%% Original (Fagiolo) In-CC
Sin=(W.').^(1/3)*(W.^(1/3))^2; % a-la Fagiolo
cyc3in=diag(Sin);
KKin=Kin;
KKin(cyc3in==0)=inf;  

CYC3in=KKin.*(KKin-1);
Cin=cyc3in./(CYC3in);               %clustering coefficient

%% Total CC
S=W.^(1/3)+(W.').^(1/3);	%symmetrized weights matrix ^1/3
cyc3=diag(S^3)/2;           %number of 3-cycles (ie. directed triangles)
K(cyc3==0)=inf;             %if no 3-cycles exist, make C=0 (via K=inf)
CYC3=K.*(K-1)-2*diag(A^2);	%number of all possible 3-cycles
C=cyc3./CYC3;               % TOTAL clustering coefficient


%% other 3 types of In-CCs
%Sin=(W.').*W.^2;

%% Fagiolo
Ksin=sum(W,2);
W_hat = W.^(1/3);

triaIn_1=diag( (transpose(W_hat)*W_hat^2) );
Kin1=Kin;
Kin1(triaIn_1==0)=inf;
triaTot_1=(Kin1.*(Kin1 - 1));
CR_in1 = triaIn_1./triaTot_1; % Fagiolo


%% Clemente
triaIn_2=0.5*diag(transpose(W)*(A+A')*A);
Kin2=Kin;
Kin2(triaIn_2==0)=inf;
Ksin2=Ksin;
Ksin2(triaIn_2==0)=inf;

triaTot_2=(Ksin2.*(Kin2-1));
CR_in2 = triaIn_2./triaTot_2; % Clemente


%% Fardet
Ksin312=sum(W.^(1/2),2);
triaIn_3=diag( transpose(W.^(2/3))*((W.^(2/3))^2) );
Kin3=Kin;
Kin3(triaIn_3==0)=inf;
Ksin3=Ksin;
Ksin3(triaIn_3==0)=inf;
Ksin312(triaIn_3==0)=inf;

triaTot_3=(Ksin312.^2-Ksin3);
%triaTot_3=(Kin1.*(Kin1 - 1));

CR_in3 =triaIn_3./triaTot_3; % Fardet 2021





