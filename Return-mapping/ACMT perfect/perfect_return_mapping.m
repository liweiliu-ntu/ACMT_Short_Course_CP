clc
clear all
close all
syms c11 c22 c33 c12 c23 c13 t11 t22 t33 t12 t13 t23 Qy d11 d22 d33  d12 d13 d23 x xx yy zz
load("Propor strain loading.csv");

% material constant
Es=200e3;
v=0.3;
tau_y=sqrt(0.5)*150;

% operator matrix
P=zeros(6,6);
Pd=P;
Pd(1:3,1:3)=-1/3*ones(3,3);
Pd=Pd+eye(6);
P(1:3,1:3)=v*ones(3,3)-(2*v-1)*eye(3);
P(4:6,4:6)=(1-2*v)/2*eye(3);
P=P/((1+v)*(1-2*v));

% stiffness modulus matreix
K=Es*P;


q=transpose(Propor_strain_loading(:,2:7));
t=transpose(Propor_strain_loading(:,1));
Q = zeros(size(q));
ff=[];
dt=0.01;
%intial stress
Q(:,1)=[100;0;0;0;0;0];
%yield function
Qy=1000;
err=10e-8;

Iflag(1)=[0];
for i = 1:length(q)-1
    q_dot(:,i)=(q(:,i+1)-q(:,i))/dt;
    Q_dot(:,i)=K*(q(:,i+1)-q(:,i))/dt;
end
tStart = cputime;
for i = 1:length(q)-1
     ff(i)=sqrt(Q(:,i)'*Q(:,i))-Qy;

    Qnext = Q(:,i) + K*(q(:,i+1)-q(:,i));
     
     if Iflag==[0]  

      [fnext] = Yield_function(Qnext,Qy);
          if fnext>err 
              Iflag=[1];
              iflag(i)=1;
           %elastic to plastic

           Q_off=Q(:,i); q_off=q(:,i);tton=t(i);Q_dot1=Q_dot(:,i);q_dot1=q_dot(:,i);
           %time_on
           [t_on,Qt,qt,d_ton] = ton(Q_dot1,Q_off,q_dot1,tton, q_off,Qy);

            Qnow=Qt;
            qnow=qt;
            qnext=q(:,i+1);

            [QQnext] = return_mapping(Qnow,K,q_dot1,Qy,dt);

             Q(:,i+1)=QQnext;
             Qnext=QQnext;
             [fnext] = Yield_function(Qnext,Qy);
             f(i)=fnext;
          
          else  
iflag(i)=2;
          %elastic to elastic
          Iflag=[0];

          Q(:,i+1)=Qnext;
           [fnext] = Yield_function(Qnext,Qy);
          f(i)=fnext;
       
          end

     elseif Iflag==[1] 
          q_dot1=q_dot(:,i);
          Qnow=Q(:,i);
          [Sg] = Straining_condition(q_dot1,Qnow);
           if Sg<=0
         %plastic to elastic
iflag(i)=3;
          Q(:,i+1)=Qnext;
           [fnext] = Yield_function(Qnext,Qy);
          f(i)=fnext;
          Iflag=[0];
         else
iflag(i)=4;  
          %plastic to plastic

             Qnow=Q(:,i);
             qnow=q(:,i);
             qnext=q(:,i+1);
             [QQnext] = return_mapping(Qnow,K,q_dot1,Qy,dt);
              Q(:,i+1)=QQnext;
              Qnext=QQnext;
             [fnext] = Yield_function(Qnext,Qy);
             f(i)=fnext;
              Iflag=[1];
     end
   
     

      
         
     end
end

tEnd = cputime - tStart


figure(1)
 plot(q(1,:),Q(1,:))
 xlabel('Strain q1 ');ylabel('Stress Q1 (Mpa)')
 
figure(2)
  plot(t(1,68:length(f)),f(1,68:length(f)))

%  figure(3)
% plot3(Q(1,:),Q(2,:),Q(3,:))
% grid on
%  xlabel('Stress Q1 (Mpa)');ylabel('Stress Q2 (Mpa)');zlabel('Stress Q3 (Mpa)')

