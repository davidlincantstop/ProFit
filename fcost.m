% cost function for the non-linear fitting function fmincon 
% is a sum of fit residual squares plus the sum of squares of penalty
% factors for regularisation 

function [f Bt_r sr]=fcost(x0)
% x0=[delta1,phi0,theta1,theta2,delta2,sigmae,sigmag]
delta1=x0(1);
phi0=x0(2);
if length(x0)==9 
     delta2=x0(3:5);
     sigmae=x0(6:8);
     sigmag=x0(9);
     M=3;
elseif length(x0)==11
     theta1=x0(3);
     theta2=x0(4);
     delta2=x0(5:7);
     sigmae=x0(8:10);
     sigmag=x0(11);
     M=3;
else
      theta1=x0(3);
      theta2=x0(4);
      delta2=x0(5:17);
      sigmae=x0(18:30);
      sigmag=x0(31);
      M=13;
end
if M==3
load basis3;
for i=1:3
basis(:,:,i)=basis3(:,:,i)';
end


else 
    load basis13;
    for i=1:13
       basis(:,:,i)=basis13(:,:,i)';
    end
end
load sig2;
sig=sig2';
rho_sigmae=8;
rho_delta2=1;
sigma_etyp=0;
T_theta1=17e-3; T_theta2=100e-3;
t1=linspace(20e-3,140e-3,25);
t2=linspace(0,511*0.25e-3,512);
if length(x0)==9 
    upsilon=zeros(1,512);
else
    
    upsilon=theta1*exp(-t2/T_theta1)+theta2*exp(-t2/T_theta2);
end
s=zeros(size(sig)); 
% tau=zeros(size(sig));

% f=linspace(-2e3,2e3,512);

b=zeros(25,512,M);
for k=1:25
   for l=1:512
    
       s(k,l)=sig(k,l)*exp(1i*(delta1*t1(k)-phi0+upsilon(l)));
%        tau(k,l)=exp(1i*pi*t1(k)*f(l));
       
         for m=1:M
             if M==3
             b(k,l,m)=basis(k,l,m)*exp(1i*delta2(m)*t2(l)-pi*sigmae(m)*(t1(k)+t2(l))-(sigmag*pi*t2(l))^2/(2*log(2)));
             else
             b(k,l,m)=basis(k,l,m)*exp(1i*delta2(m)*t2(l)-pi*sigmae(m)*(t1(k)+t2(l))-(sigmag*pi*t2(l))^2/(2*log(2))); 
             end
         end
   end
end


  bb=zeros(128,512,M);
  
for m=1:M
   
    bb(:,:,m)=fliplr(fftshift(fft(fftshift(fft(b(:,:,m),512,2),2),128,1),1));
   
end
s=fftshift(fft(fftshift(fft(s,512,2),2),128,1),1);
K=128;L=512;


Bt_r=zeros(K*L,M);
for m=1:M
    for k=1:K
        for l=1:L
            
            if ~(((k>=59)&&(k<=70))&&((l>=123)&&(l<=230)))
                bb(k,l,m)=0;
                if m==1
                s(k,l)=0;
                end
            end
            
            
        end
    end
end
for m=1:M
    Bt_r(:,m)=real(reshape(bb(:,:,m)',K*L,1));
end


sr=real(reshape(s',K*L,1));

c=pinv(Bt_r)*reshape(sr,K*L,1);
f=0;
for k=59:70
    for l=123:230
        f=f+(sr((k-1)*512+l)-Bt_r((k-1)*512+l,:)*c)^2;
    end
end
   f=f+rho_sigmae/(size(b,3)-1)*sum((sigmae-sigma_etyp).^2)+rho_delta2*var(delta2);
