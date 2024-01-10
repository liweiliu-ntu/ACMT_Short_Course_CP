clc
clear all
close all
syms c11 c22 c33 c12 c23 c13 t11 t22 t33 t12 t13 t23 Qy d11 d22 d33  d12 d13 d23 x xx yy zz
load("Propor strain loading.csv");

q=transpose(Propor_strain_loading(:,2:7));
t=transpose(Propor_strain_loading(:,1));


% material constant
Es=200e3;
v=0.3;
dt=0.01;

% Hardening parameter
kappa1=10000;
kappa2=0.005;

% Kinematic parameter
a1=50000;
a2=0.01;

%distortional model parameter
c=0.0008;

% stiffness modulus matreix
P= [ 2/3 -1/3 -1/3 0      0      0;
     -1/3 2/3 -1/3 0      0      0;
     -1/3 -1/3 2/3  0     0      0;
     0     0     0     1   0      0;
     0     0     0     0      1  0;
     0     0     0     0      0      1];
P1=P/((1+v)*(1-2*v));
K=Es*P1;

%intial stress
Q(:,1)=[0;100;0;0;0;0];

%intial back stress
Qb(:,1)=[0;0;0;0;0;0];

alpha(:,1)=[0;0;0;0;0;0];
%yield function
tau(1)=150;
err=10e-9;


Q = zeros(size(q));
ff=[];
Iflag=[];


Iflag(1)=[0];
for i = 1:length(q)-1
    q_dot(:,i)=(q(:,i+1)-q(:,i))/dt;
    Q_dot(:,i)=K*(q(:,i+1)-q(:,i))/dt;
end


for i = 1:length(q)-1

    
     if Iflag==[0]  

      Qnext = Q(:,i) + K*(q(:,i+1)-q(:,i));

      alpha1=alpha(:,i);

      tau1=tau(i);

      [fnext] = Yield_function(Qnext,P,tau1,c, alpha1) ;

          if fnext>err 
          %elastic to plastic  
            Iflag=[1];

            iflag(i)=1;

            Q_now=Q(:,i);q_now=q(:,i); Qb_now=Qb(:,i);tton=t(i);Q_dot1=Q_dot(:,i);q_dot1=q_dot(:,i);

           taunow=tau(i);
           %time_on
           [t_on,Qton,qton,d_ton] = ton(Q_dot1,Q_now,Qb_now,c, q_now,q_dot1,taunow,tton,P);
           
            Qnow=Qton;
            
            qnow=qton;

            Qbnow=Qb(:,i);

            qnext=q(:,i+1);

            tau1=tau(i);

            alpha2=alpha(:,i);

           [ QQnext, alphanext, taunext,fnext] = flowRule(Qnow, Qbnow, tau1, Q_dot1, c, alpha2,kappa1, kappa2, a1, a2, err,dt,K,P);

            tau(i+1)= taunext;

            alpha(:,i+1)=alphanext;

            Qb(:,i+1)=alphanext;

            Q(:,i+1)=QQnext;

            ff(i)=fnext;

          else  
          %elastic to elastic
          Iflag=[0];

           iflag(i)=2;

           tau(i+1)=tau(i);

           alpha(:,i+1)=alpha(:,i);

           Qb(:,i+1)=Qb(:,i);

           Q(:,i+1)=Qnext;

           ff(i)=fnext;     
       
          end

     elseif Iflag==[1] 

       Qnext = Q(:,i) + K*(q(:,i+1)-q(:,i));   

       alpha1=alpha(:,i);

       tau1=tau(i);

      [fnext] = Yield_function(Qnext,P,tau1,c, alpha1) ;


           if fnext<=0

         %plastic to elastic

          Iflag=[0];

          iflag(i)=3;

          tau(i+1)=tau(i);

          Qb(:,i+1)=Qb(:,i);

          alpha(:,i+1)=alpha(:,i)

          Q(:,i+1)=Qnext;

          ff(i)=fnext;

         else
  
          %plastic to plastic

           Iflag=[1];

            iflag(i)=4;

            Qnow=Q(:,i);

            qnow=q(:,i);

            qnext=q(:,i+1);
              
            Qbnow=Qb(:,i);
      
            tau1=tau(i);

            alpha2=alpha(:,i);

           [ QQnext, alphanext, taunext,fnext] = flowRule(Qnow, Qbnow, tau1, Q_dot1, c, alpha2,kappa1, kappa2, a1, a2, err,dt,K,P);

            tau(i+1)= taunext;

            alpha(:,i+1)=alphanext;

            Qb(:,i+1)=alphanext;
            
            Q(:,i+1)=QQnext;

            ff(i)=fnext;

            
     end
      
     end


end



 figure(1)
 plot(q(1,:),Q(1,:))
 xlabel('strain q1');ylabel('Stress Q1 (MPa)')

%  figure(2)
% plot3(Q(1,:),Q(2,:),Q(3,:))
% grid on
%  xlabel('Stress Q1 (Mpa)');ylabel('Stress Q2 (Mpa)');zlabel('Stress Q3 (Mpa)')


figure(3)
  plot(t(1,10:length(ff)),ff(1,10:length(ff)))
