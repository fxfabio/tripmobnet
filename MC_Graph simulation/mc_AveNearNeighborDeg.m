function [aveK,k_bin]=mc_AveNearNeighborDeg(Adj,bin,inout);
%gives a list of k-mean degree connectivity for the network for successive
%k=0,1,2....
% The average degree connectivity is also known as average nearest neighbor degree.
% The k-average degree connectivity is the average of the mean neighbor
% degrees of vertices of degree k
% 
%       Input: Adj - adjacencey matrix (also weighted and directed)
%                bin = subdivision of the degree interval on which find the
%                          average connectivity
%                inout = choose if "in", "out" or "total"
% it works also for weighted and directed networks
%@F.Vanni2023 

if nargin<3
    inout='total'
end


if strcmpi(inout,'in');Adj=Adj'; ave='in';
elseif strcmpi(inout,'out');Adj=Adj; ave='in';
else ave='total';
end

[deg,indeg,outdeg]=degrees(Adj);

switch ave
    case 'in'
        degr=max(indeg);
        k_bin=linspace(0,degr,bin);
        aveK=zeros(1,length(k_bin));
        k=indeg;
        knn=aveNeighborDegWW(Adj,'in');

        for i=1:length(k_bin)-1
           ix=(k>k_bin(i) & k<=k_bin(i+1) );
           aveK(i)=mean(knn(ix));
        end
        
        
    case 'out'
        degr=length(outdeg);
        k_bin=linspace(0,degr,bin);
        aveK=zeros(1,length(k_bin));
        k=outdeg;
        knn=aveNeighborDegWW(Adj,'out');

        for i=1:length(k_bin)-1
           ix=(k>k_bin(i) & k<=k_bin(i+1) );
           aveK(i)=mean(knn(ix));
        end
        
    case 'total'
        degr=max(deg);
        k_bin=linspace(0,degr,bin);
        aveK=zeros(1,length(k_bin));
        k=deg;
        knn=aveNeighborDegWW(Adj,'total');

        for i=1:length(k_bin)-1
           ix=(k>k_bin(i) & k<=k_bin(i+1) );
           aveK(i)=mean(knn(ix));
        end

end




function ave_n_deg=aveNeighborDegWW(adj,inout)
%  Average weighted and directed neighbor degree


ave_n_deg=zeros(1,length(adj));   % initialize output vector

A=double(adj>0); % adjacency matrix
W=adj; % weight matrix



if nargin<2
    inout='total';
end

switch lower(inout)
    case 'in'
        adj=adj;
        deg = sum(adj);
    case 'out'
        adj=transpose(adj);
       deg = sum(adj);
    case 'total'
        [deg,~,~]=degrees(adj);
end


for i=1:length(adj)  % across all nodes
  
  neigh=kneighbors(adj,i,1);  % neighbors of i, one link away
  peso=W(i,neigh);
  if isempty(neigh); ave_n_deg(i)=0; continue; end
  stren(i)=sum(peso);
  ave_n_deg(i)=sum(deg(neigh).*peso)/(stren(i));
end
