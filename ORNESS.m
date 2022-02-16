%ORNESS Fusion method
% in this fusion code 4 sattelite-based streamflow datasets (TRMM-CMORPH-PERSIANN-MSWEP) are fused with
% ground-based streamflow simulations
clc
clear
tic
load('USGS0706930520022009','USGS0706930520022009');
USGS0706930520022009=USGS0706930520022009*0.0283168466 ;
USGS0706930520022009mm=(USGS0706930520022009*3600/(2216.179*10^6))*1000;
load('cmorph','ADA');
cmorph=ADA*1;
load('ground','ADA');
ground=ADA*1;
load('mswep','ADA');
mswep=ADA*1;
load('persiann','ADA');
persiann=ADA*1;
load('trmmv7','ADA');
trmmv7=ADA*1;
%%
%calculate average daily precipitation 
for j=1:10
lb=zeros(5,1);
ub=[1;1;1;1;1];
aeq=[1 1 1 1 1;1 0.75 0.5 0.25 0];
beq=[1;j/10];
options=gaoptimset('PopulationSize',200,'Generations',500,'PlotFcns',@gaplotbestf);
%,'StallGenLimit',500)
[x,fval(j)] = ga(@ornesscostfunction,5,[],[],aeq,beq,lb,ub,[]);
ss(j,:)=x;
for i=1:1000
optstpr(i,1)=x(1)*persiann(i)+x(2)*ground(i)+x(3)*cmorph(i)+x(4)*trmmv7(i)+x(5)*mswep(i);
end
qave=mean(USGS0706930520022009mm(1:1000));
for i=2:1000
    dife(i)=(USGS0706930520022009mm(i)-qave)^2;
end
for i=2:1000
    se(i)=((USGS0706930520022009mm(i))-(optstpr(i)))^2;
end
nse(j)=1-(sum(se(250:1000),'omitnan')/sum(dife(250:1000),'omitnan'))

end
toc