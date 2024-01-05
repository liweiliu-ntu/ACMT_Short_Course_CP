function [ts,nextq,nextdq,nextdQ,nextQb,nextQa,nextyf,nextqp,nextQ] = pull_back (iflag,dQ,q,ke,kp,dt,Qy,Q,Qb,ee,ii)
 syms ts
 f=Q(ii-1,:)+dQ*ts;
 f1=norm(f)-Qy;
 ts1=solve(f1,ts);
 ts=double(ts1);
 [nextQ,nextq,nextdq,nextdQ,nextQb,nextQa,nextyf,nextqp]=off_phase(q,ke,kp,ts,Qy,Q,Qb,ee,ii);

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
         