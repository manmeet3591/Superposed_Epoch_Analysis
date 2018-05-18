function [SEA_means,SEA_prctiles]=SEA(T,ev,before,after)
%function [SEA_means,SEA_prctiles]=SEA(T,ev,before,after,varargin)

% Aim: Superposed Epoch Analysis (SEA) with confidence intervals estimated
% using block bootstrap, the procedure is freely inspired from Adams et al.
% (2003)
% 
% Input: T: time serie of any proxy (work both for dendrochonological
%          or paleoecological proxies) provided as a vector, note that the
%          time serie must be evenly spaced 
%        ev: list of key events from witch to compute SEA
%        before: Lag time units before key events for SEA  
%        after: Lag time units after key events for SEA 
%        nbboot: mumber of boostraps (default 10000)
%        b: block width for the block bootstrap procedure if unspecified b  
%           is computed such as b=2*tau where tau=1/log(r), r is the 
%           lag-one autocorrelation coefficient for the T serie 
%           (see Adams et al. 2003)
%        prctilesIN: percentiles associated with the block bootstrap
%                    by default prctilesIN=[97.5,2.5,95,5]
%        display: set to zero for no display
% 
% Output: SEA_means: Averaged zscores following the SEA procedure
%         SEA_prctiles: Associated percentiles
%[SEA_means,SEA_prctiles]=SEA(T,ev,20,20,'prctilesIN',[5,95],'nbboot',9999,'b',0)

% Example:
%{                     
t=1:500;
d=2.5*sin(2*pi*t/100)+1.5*sin(2*pi*t/41)+1*sin(2*pi*t/21);
T=(d+2*randn(size(d))).';
figure()
plot(t,T,'g')
ev=[50,93,131,175,214,257,297,337,381,428,470].';
hold on
%plot(ev,repmat(4,length(ev),1),'v')
%[SEA_means,SEA_prctiles]=SEA(T,ev,20,20,'prctilesIN',[5,95])
%}
%         
% Citation: Adams, J. B., et al. 2003. Nature 426:274-278.
% 
% 
% 
% Copyright: Olivier Blarquez
% If you found this function usefull you can cite: 
% Blarquez O., Carcaillet C. 2010. Fire, fuel composition and resilience 
% threshold in subalpine ecosystem. PLoS ONE 5(8) : e12480. 
% doi:10.1371/journal.pone.0012480
% blarquez@gmail.com
% Date: 14/02/2013   
    

%Set default values
%pnames = {'nbboot', 'b', 'display','prctilesIN'};
%dflts  = {10000, 0, 1, [97.5,2.5,95,5]};
%[nbboot,b,display,prctilesIN] = parseArgs(pnames,dflts,varargin{:});    
nbboot = 10000
b = 0
display = 0
prctilesIN = [97.5, 2.5, 95, 5]


%Calculates block width
if b==0
r=corrcoef(T(2:end),T(1:end-1));
b=round(2*(-1/log(r(2))));
end

% Arrange data
ev=ev+before;
T=([NaN(before,1);T;NaN(after,1)]);

% Reconstruct the key event matrix
%size(T)
for i=1:length(ev)
    T_matrix(:,i)=T(ev(i)-before:ev(i)+after);
end
T_matrix=[repmat(-999,b,length(ev));T_matrix;repmat(-999,b,length(ev))];

% Block bootstrap procedure
n=ceil((before+after+1+(b*2))/b);
for k=1:nbboot
    for i=1:size(T,2)   
            q=fix((size(T_matrix,1)-b))+1;
            y=zeros(b,q);
            for j=1:q,   
              y(:,j)=T_matrix((j-1)+1:(j-1)+b,i);
            end; 
        o=randsample(1:q,n*2,'true') ;
        yy{i}=reshape(y(:,o),n*2*b,1);
        yy{i}(yy{i}==-999)=[];
        y_n(:,i)=yy{i}(1:(before+after+1));
    end
  yy_m(:,k)=zscore(nanmean(y_n,2));
end

SEA_means=zscore(nanmean(T_matrix(b+1:size(T_matrix,1)-b,:),2));

for i=1:length(prctilesIN)
SEA_prctiles(:,i)=prctile(yy_m.',prctilesIN(i)).';
end

if display==1
figure()
    bar(-before:after,SEA_means,'facecolor','w', 'edgecolor','k')
for i=1:length(prctilesIN)
hold on
plot(-before:after,SEA_prctiles(:,i),'k-')
end
xlabel('Lag time from key event')
ylabel('Proxy zscores')
end
%source http://paleoecologie.umontreal.ca/public/code/SEA.m

