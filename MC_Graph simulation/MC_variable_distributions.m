function [Fd]=MC_variable_distributions(n,varargin)
% function that simulates a scale-free network with preferential attachment
% with fitness preferences as in WP_LEM.
%INPUT:
%         n - total number of nodes (size of the financial bank system)
%         mu - the power law coeff. for the fitness distribution (bank assets)
%         ad = average degree for directed networks ad=L/N
%         alfabeta=[alpha beta] - the power law parameters of the linking function
%         xmin - the min size of the size distribution (assets size)
%         xmax - the max size of the size distribution (assets size) (in truncated pareto)
%  varargin -> decide the distribution of the fitness variables (assets size of banks)
%
%OUTPUT
%        Fd = the size distribution from the fitness distribution
% Example:
%        fitness_density_MFa(n,'powerlaw',mu);
%        fitness_density_MFa(n,'xmin',xmin,'powerlaw',mu);
%        fitness_density_MFa(n,'xmin',xmin,'truncated',mu,xmax);
%        fitness_density_MFa(n,'xmin',xmin,'powercutoff',mu,lambda);
%        fitness_density_MFa(n,'exponential',lambda);
%        fitness_density_MFa(n,'stretched',lambda,bs);
% FabioVanni©2015



%default values 
xmin=5;
xmax=50;
mu=2;

%% ARGUMENTS selections
% parse command-line parameters; trap for bad input
if isempty(varargin), disp(' !!! ATTENTION !!! Using as power law fitness with mu=2'),
    type='PL';
    xmin=1;mu=2;
end
i=1; 
while i<length(varargin), 
  argom = 1; 
  if ischar(varargin{i}), 
    switch varargin{i},
        case 'xmin',            xmin = varargin{i+1}; i = i + 2;
        case 'powerlaw',     type = 'PL'; mu  = varargin{i+1}; i = i + 1;
        case 'truncated',    type = 'TC'; mu = varargin{i+1}; xmax = varargin{i+2};i = i + 2;
        case 'powercutoff',          type = 'PC'; mu  = varargin{i+1}; lambda = varargin{i+2}; i = i + 2;
        case 'exponential',     type = 'EX'; lambda = varargin{i+1}; i = i + 3;
        case 'stretched',       type = 'ST'; lambda = varargin{i+1}; bs = varargin{i+2}; i = i + 4;
        otherwise, argom=0;disp(['(scalefree) !!! INVALID INPUTS !!! , using default truncated pareto' ]);i=5 ; 
            n=250;xmin=5;xmax=50;mu=2; type='TC';
    end
  end
end

switch type
    case 'EX',% EXPONENTIAL (random graph)
        x = xmin - (1/lambda)*log(1-rand(n,1));
    case 'TC', % TRUNCATED Pareto
            alfa=mu-1;
            u=rand(n,1);
            l=xmin;
            h=xmax;
            if l<h,
            x=(-(u*h^alfa - u*l^alfa - h^alfa)./( h^alfa * l^alfa ) ).^(-1/alfa);
            else error('wrong ranges')
            end
    case 'ST', % stretched exponential truncation
        x = (xmin^bs - (1/lambda)*log(1-rand(n,1))).^(1/bs);
    case 'PC', % Exponential cutoff
        x = [];
        y = xmin - (1/lambda)*log(1-rand(10*n,1));
            while true
                y(rand(10*n,1)>=(y./xmin).^(-mu)) = [];
                x = [x; y];
                q = length(x)-n;
                if (q==0), break;
                end;
                if (q>0),
                    r = randperm(length(x));
                    x(r(1:q)) = [];
                    break;
                end;
                if (q<0),
                    y = xmin - (1/lambda)*log(1-rand(10*n,1));
                end;
            end; 
    case 'PL', % POWERLAW with xmin
         x = xmin*(1-rand(n,1)).^(-1/(mu-1));
    otherwise, x = xmin*(1-rand(n,1)).^(-1/(mu-1));
end;
%%
  
Fd=x; 



