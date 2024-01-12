function [fnext] = Yield_function(Qnext,P,tau1,c, alpha1) 
s=P*Qnext;
n=(s-alpha1)/sqrt((s-alpha1)'*(s-alpha1));
%     fnext=1.5.*(1-c.*dot(n',alpha1))*(s-alpha1)'*(s-alpha1)-tau1.^2
     fnext=sqrt(1.5*(1-c.*dot(n',alpha1)))*sqrt((s-alpha1)'*(s-alpha1))-tau1;
%      fnext=1.5.*(1-c*(Qbb'*P*(Qnext-Qbb))/(sqrt(Qnext-Qbb)'*P*(Qnext-Qbb)))*((Qnext-Qbb)'*P*(Qnext-Qbb))-tau1.^2;
end
