function [QQnext] = return_mapping(Qnow,Y,K,qnext,b,r,qnow,C)
epsilon=10e-8;
% C=inv(K);
  for n = 1:200
                          lamda_n(1,1)=0;
                          %phi=0.5*Q'*Y*Q+Q'*p-0.5*r^2
                          %Q=((C+lamda_n(1,n).*Y)^(-1))*(C*Qt(:,1)+q(:,i+1)-qt(:,1)-lamda_n(:,n).*p)
                          %∆lambda(n+1)=∆lambda(n)-Phi/(d(phi)/d(∆lambda)) 
                          %sigma(n)=sigma_0+∆sigma=(C+lambda(n)*Y)^(-1)(q_e+∆q-∆lambda*p)
                          %d(phi)/d(∆lambda)=-(Y*Q(n)+p)'(C+∆lambda*Y)^(-1)*(Y*Q(n)+p)
                          lamda_n(1,n+1)=lamda_n(1,n)-(0.5*transpose(((C+lamda_n(1,n).*Y)^(-1))*(C*Qnow(:,1)+qnext-qnow...
                              -lamda_n(:,n).*b))*Y*((C+lamda_n(1,n).*Y)^(-1))*(C*Qnow(:,1)+qnext-qnow-lamda_n(:,n).*b)...
                              +transpose(((C+lamda_n(1,n).*Y)^(-1))*(C*Qnow(:,1)+qnext-qnow-lamda_n(:,n).*b))*b-0.5.*r.^2)./...
                                  (-transpose(Y*((C+lamda_n(1,n).*Y)^(-1))*(C*Qnow(:,1)+qnext-qnow-lamda_n(:,n).*b)+b)*...
                                  (C+lamda_n(1,n).*Y)^(-1)*(((C+lamda_n(1,n).*Y)^(-1))*Y*(C*Qnow(:,1)+qnext-qnow-lamda_n(:,n).*b)+b));                     
                           f_val(1,n)=(0.5*transpose(((C+lamda_n(1,n+1).*Y)^(-1))*(C*Qnow(:,1)+qnext-qnow-lamda_n(1,n+1).*b))...
                               *Y*((C+lamda_n(1,n+1).*Y)^(-1))*(C*Qnow(:,1)+qnext-qnow-lamda_n(1,n+1).*b)+transpose(((C+lamda_n(1,n+1).*Y)^(-1))...
                               *(C*Qnow(:,1)+qnext-qnow-lamda_n(1,n+1).*b))*b-0.5.*r.^2);
%                           lamfa_roop(i+1)=n;
                           if abs(f_val(1,n)) < epsilon  & abs(f_val(1,n-1)) > epsilon                             
                               QQnext=((C+lamda_n(1,n+1).*Y)^(-1))*(C*Qnow(:,1)+qnext-qnow-lamda_n(1,n+1).*b);
                           lambda=lamda_n(1,n+1);
                    
                     gggg=n;

                               break;

                           end
%                            if n==200
%                            [Xnext] = return_free(Qnow,uu,D,W,S,dt,Phi,g,PhiStar,r);
%                            QQnext=Xnext(1:6);
%                            
% kkkk=1
% 
% break;
                           end
                           
                       
end