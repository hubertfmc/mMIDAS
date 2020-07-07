
function [be_values]=plot_b_weights(ind,k1,k2)
%Beta weights function.
if length(ind)>1
    u=(ind(1)-ind)/(ind(1)-ind(end));
    u(1)=u(1)+eps;
    u(end)=u(end)-eps;
end

if isscalar(ind)
    u=linspace(eps,1-eps,ind)';
end

% k1=(k1*(k1>0)+1e-8*~(k1>0))*(k1<300)+300*~(k1<300);% 0<k1<300
% k2=(k2*(k2>0)+1e-8*~(k2>0))*(k2<300)+300*~(k2<300);% 0<k1<300

be_values=u.^(k1-1).*(1-u).^(k2-1);%u.^(k1-1).*(1-u).^(k2-1)./beta(k1,k2);
be_values=be_values/sum(be_values);
end