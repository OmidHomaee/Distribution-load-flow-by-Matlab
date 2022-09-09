%                      POWER FLOW IN DISTRIBUTION NETWORKS
%**************************************************************************
%function [DLF_result]=distributionloadflow(network_data_result)

clc
clear all
model=SelectModel()

Sbase=200;
n=size(model,1);   %please enter bus numbers:
L=[model(:,3)-1 1.6*(model(:,6)+i*model(:,7))/Sbase]; 

Vbase=12660;
Zbase=(Vbase^2)/(Sbase*1000);

%**************************************************************************
[sal,sol]=size(L);
Busno=L(:,1);
phaz=L(:,2);
v=0;
s=0;
%***********************************************************************
%enter the line impedance like below:
%   [sd   re    z  ]
%               ...                                                                                                                                                                                                       ]  
zl=[model(:,2:1:3)-1 (model(:,4)+i*model(:,5))/Zbase];
%**************************************************************************
v=0;
s=0;
sl=0;
time=1;
for k0=2:1:(time+1)
v1=v;
tt=1.;
t1=tt;
v=repmat(t1,1,n);
[sa,so]=size(zl);
b=zl(:,1);
bb=zl(:,2);
I=zeros(1,n);
J=zeros(1,n);
h=0;
nn=0 ;
vv=0;
mm=1;
%for pp=1:20
while mm >0.0001
%**************************************************************************
%                   Nodal current calculation
nn=nn+1;
 if ~mod(nn,10)
    fprintf('%d\n',nn)
    fprintf('%d\n',mm)
  end
I=zeros(1,n); 
for k1=1:n
    i=(-1)^0.5;
     for k2=1:sal
         if Busno(k2)==k1
            I(1,k1)=I(1,k1)+((conj((L(k2,k0))/(v(1,k1)))));
         end
     end
          J(1,k1)=I(1,k1);
end
%**************************************************************************
%              Backward sweep—section current calculation
i=0;
for i=n:-1:1
    b1=b';
    indx=find(b1==i);
    [sai soi]=size(indx);
    if soi>0
    for l=1:soi
        k1=indx(1,l);
       J(:,i)=J(:,i)+J(:,k1);
    end
    end
end
%**************************************************************************
%               Forward sweep—nodal voltage calculation
for i=0:1:n
if i==0
   indx=find(b1==i);
   [sai soi]=size(indx);
    for l=1:soi
     k1=indx(1,l);
     v(1,k1)=t1-[zl(k1,3)]*J(1,k1);
    
    end
else
    indx=find(b1==i);
   [sai soi]=size(indx);
   for l=1:soi
     k1=indx(1,l);
     v(1,k1)=v(1,i)-[zl(k1,3)]*J(:,k1);
   
   end
end
end
%**************************************************************************
%                      Convergence Criterion
dv=v-vv;
mm=abs(max(max(dv)));
vv=v;
end
%**************************************************************************
voltage(k0-1,:)=abs(v);
%**************************************************************************
%                        power loses calculation
    sl(k0-1)=sum(zl(:,3).*((abs(J(:))).^2));
   % sl(k0-1)=sum(zl(:,3).*((J(:)).^2));

 ActiveLoss(1,k0-1) = real(sl(k0-1))*(Sbase);
 ReactiveLoss(1,k0-1) = imag(sl(k0-1))*(Sbase);


SHVnetwork(k0-1)=t1*conj(J(1))*(Sbase);
end


DLF_result.SHVnetwork=SHVnetwork;
DLF_result.voltage=voltage;
DLF_result.ConvergenceIteration=nn;
DLF_result.ActiveLoss=ActiveLoss;
DLF_result.ReactiveLoss=ReactiveLoss;

Active_Losses=ActiveLoss
Reactive_Loss=ReactiveLoss
Reactive_HV_network=SHVnetwork