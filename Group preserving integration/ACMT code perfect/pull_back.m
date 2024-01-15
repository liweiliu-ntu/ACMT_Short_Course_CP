function [ts,nextQ,nextq,nextqe,nextdq,nextdQ,nextyf,nextqp] = pull_back (iflag,dQ,q,qe,ke,dt,dq,Qy,Q,ee,ii)
 syms ts
f=ke*q(ii-1,:)+ke*dq*ts;
 f1=norm(f)-Qy;
 ts1=solve(f1,ts);
 ts2=double(ts1);
 for jj=1:length(ts2)
     if ts2(jj)>0
         ts=ts2(jj);
     end
 end
 ts;
%  [nextQ,nextq,nextqe,nextdq,nextdQ,nextyf,nextqp]=off_phase(Q,ke,dt,Qy,q,qe,ee,ii);
nextq=q(ii,:);
nextdq=(q(ii-1,:)-q(ii-2,:))/dt;
nextdQ=nextdq*ke;
nextqe=nextq;
nextqp=round(nextq-nextqe,8);
nextQ=(ke*nextdq*ts+Q);
nextyf=norm(nextQ)-Qy;

% tolrance=1e-9;
% t_start=0.0;
% tau=0.01;
% tempp=-Qy;
% iit=0;
% %
% while abs(tempp) > tolrance
%     iit=iit+1;
%     tau=tau/10;
%     for jj=2:10
%         if iflag==[0]
%             [nextQ,nextq,nextdq,nextdQ,nextQb,nextQa,nextyf,nextqp]=off_phase(q,ke,kp,t_start+tau*(jj-1),Qy,Q,Qb,ee,ii);
%         end
%         tempp=nextyf;
%         if tempp > tolrance
%             t_start=t_start+tau*(jj-1);
%             break      
%         end
%     end
%     if (iit > 10000)
%         fprintf('pull-back warning,');
%         fprintf('\t tempp=%8.5e,',tempp);
%         fprintf('\t t_start=%8.5e,',t_start); 
%         break
%     end
% end
% ts=t_start+tau*jj;
         