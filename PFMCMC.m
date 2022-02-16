%Data Assimilation code
%Particle Filter Marcov Chain Monte Carlo (PFMCMC) data assimilation method
tic
clear all
close all
clc
global k xxx meanvector tot N ADA s fmeanv noisedet noisedrain USGS0706930520022009mm T t tsorteddata pet120022009 thiessen pet20022009 USGS0706930520022009 spm
%catchment no:1
%area:2216.179 km2
%%
% %import forcqrrring data
load ('pet20022009','pet120022009')%potential evapotranspiration by mm/day unit
load ('thiessen','thiessen')%precipitation in ground station inch unit

% precipitation in days which no data exists, evaluated by multi regression
% analysis based on 25 neighbour stations
%convert unit to mm
thiessen=thiessen*25.4;
pr99=quantile(thiessen,0.99);
pr90=quantile(thiessen,0.90);
pr95=quantile(thiessen,0.95);

for i=1:1096
    vectorrevised(i,1)=thiessen(i);
end
load ('USGS0706930520022009', 'USGS0706930520022009')%discharge in USGS station in the catchment, unit:ft3/sec
%convert discharge unit to m3/sec
USGS0706930520022009=USGS0706930520022009*0.0283168466 ;
%convert discharge unit to mm
USGS0706930520022009mm=(USGS0706930520022009*3600*24/(2216.179*10^6))*1000;
us75=quantile(USGS0706930520022009,0.02);

%%
%Define totall number of particles which will be generated for each state
%variable or parameter during data assimilation (N)
N = 500;
save N
%select optimal variance for metrolopis ratio
% s=2.38/sqrt(2*d), d:number of parameters,(Moradkhani et al, 2012)
%s=2.38/sqrt(2*19)=0.386;
s=0.05;%%%%%%%%%%%%%%
%define range of SAC-SMA hydrological model parameters(Abbaszadeh et al,
%2018)
range=xlsread('parametersrange.xlsx');
%%
%%
%Define matrices
%part: vector of a particle in data assimilation, which include 21 cells,
%13 for parameters, 6 for statevariables, 2 for forcing data (precipitation
%and potential evapotranspiration)
part=zeros(1,21);
%simrun: simulated runoff
simrun=zeros(1,N);
%cwpd:comulative weight of particles
cwpd=zeros(1,N);
%pm: parameter matrix
pm=zeros(N,19);
%propm:proposal parameter matrix
propm=zeros(N,19);
weight=zeros(N,T);
%mw:multiple weight
mw=zeros(N,T);
%like:likelihood
like=zeros(N,T);

%optmum values of parameters and statevariables by GA
opt=[116.298562792773,83.6627766827241,257.293449514877,499.827759361930,99.9683184293643,0.227394369100586,0.300948442578945,0.00694222958788242,0.0861171226955134,119.531078845442,4.59369238113657,0.00111989379831073,0.0700570577674149,120.117096219876,77.9005910745747,63.4297420019878,551.168562152300,15.2111648368127,0.435675031534184];


for j=1:19
    for i=1:N
        pm(i,j)=opt(j);
    end
end

for i=250:-1:1
    pm(i,1)=pm(i,1)-(250-i)*0.05;
    pm(i,2)=pm(i,2)-(250-i)*0.05;
    pm(i,3)=pm(i,3)-(250-i)*0.05;
    pm(i,4)=pm(i,4)-(250-i)*0.05;
    pm(i,5)=pm(i,5)-(250-i)*0.05;
    pm(i,6)=pm(i,6)-(250-i)*0.0005;
    pm(i,7)=pm(i,7)-(250-i)*0.001;
    pm(i,8)=pm(i,8)-(250-i)*0.00001;
    pm(i,9)=pm(i,9)-(250-i)*0.0001;
    pm(i,10)=pm(i,10)-(250-i)*0.05;
    pm(i,11)=pm(i,11)-(250-i)*0.005;
    pm(i,12)=pm(i,12)-(250-i)*0.0001;
    pm(i,13)=pm(i,13)-(250-i)*0.0005;
    pm(i,14)=pm(i,14)+abs(250-i)*0.05;
    pm(i,15)=pm(i,15)+abs(250-i)*0.05;
    pm(i,16)=pm(i,16)-(250-i)*0.05;
    pm(i,17)=pm(i,17)-(250-i)*0.05;
    pm(i,18)=pm(i,18)+abs(250-i)*0.05;
    pm(i,19)=pm(i,19)-(250-i)*0.0005;
end
for i=250:500
    pm(i,1)=pm(i,1)+abs(250-i)*0.05;
    pm(i,2)=pm(i,2)+abs(250-i)*0.05;
    pm(i,3)=pm(i,3)+abs(250-i)*0.05;
    pm(i,4)=pm(i,4)+abs(250-i)*0.05;
    pm(i,5)=pm(i,5)+abs(250-i)*0.05;
    pm(i,6)=pm(i,6)+abs(250-i)*0.0005;
    pm(i,7)=pm(i,7)+abs(250-i)*0.001;
    pm(i,8)=pm(i,8)+abs(250-i)*0.00001;
    pm(i,9)=pm(i,9)+abs(250-i)*0.0001;
    pm(i,10)=pm(i,10)+abs(250-i)*0.05;
    pm(i,11)=pm(i,11)+abs(250-i)*0.005;
    pm(i,12)=pm(i,12)+abs(250-i)*0.0001;
    pm(i,13)=pm(i,13)+abs(250-i)*0.0001;
    pm(i,14)=pm(i,14)+abs(250-i)*0.05;
    pm(i,15)=pm(i,15)+abs(250-i)*0.05;
    pm(i,16)=pm(i,16)+abs(250-i)*0.05;
    pm(i,17)=pm(i,17)+abs(250-i)*0.05;
    pm(i,18)=pm(i,18)+abs(250-i)*0.05;
    pm(i,19)=pm(i,19)+abs(250-i)*0.0005;
    
end
%me:vector of mean values
%va:vector of variances
for i=1:13
    me(i,1)=mean(pm(:,i));
    va(i,1)=var(pm(:,i));
end

%%
for i=1:N % part 1 from section 1 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
    part=pm(i,:);
    part(20)=thiessen(1);
    part(21)=pet120022009(1);
    [surf(i),base(i),tot(i),s1(i),s2(i),s3(i),s4(i),s5(i),s6(i)]=sacsmafunction(part);%simulated runoff
    pm(i,14)=s1(i);%update state variable values
    pm(i,15)=s2(i);
    pm(i,16)=s3(i);
    pm(i,17)=s4(i);
    pm(i,18)=s5(i);
    pm(i,19)=s6(i);
    if isnan(simrun(i))
        like(i,1)=0;
    end
    difference(i,1)=tot(i)-USGS0706930520022009mm(1,1);
end
%calculate standard deviation of difference between simulated and observed runoff
difvar=var(difference(:,1));
for i=1:N
    like(i,1)=(1/sqrt(2*pi*difvar))*exp(-(1/(2*difvar))*(difference(i,1))^2);%calculate the likelihood
end
likesum=zeros(N);
for i=1:N
    likesum(i)=sum(like(:,1));
end
%ADA:Average Data Assimilation simulated runoff
ADA(1)=mean(tot(1:N));
SORTEDADA(1)=ADA(1);
for i=1:N % part 2 from section 1 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
    weight(i,1)=like(i,1)/likesum(i);%it might be zero
    
end
%%

for j=1:13
    for i=1:N
        weipar(i,j)=weight(i,1)*pm(i,j);
    end
end
for j=1:13
    mu(j,1)=sum(weipar(:,j));% part 1 of section 4 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
end
for j=1:13
    for i=1:N
        sigcal(i,j)=weight(i,1)*(pm(i,j)-mu(j,1))^2;
    end
end
for j=1:13
    sig(j,1)=sqrt(sum(sigcal(:,j)));% part 1 of section 4 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
end
for j=1:13
    for i=1:N
        curweight(i,j,1)=normpdf(pm(i,j),mu(j,1),sig(j,1))*like(i);% part 2 of section 4 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
    end
end


for i=1:N
    x(i)=weight(i,1);
end
roundx=round(x,4);
y=sort(x);
z=zeros(N,1);
for i=1:N/2
    z(i)=y(N/2+i);
end

mpm=zeros(N,19);
for i=1:N/2
    aaa=z(i);
    for t=1:N
        diffe(t)=aaa-x(t);
        if diffe(t)==0
            mpm(i,:)=pm(t,:);
        end
    end
end

zsum=sum(z(:,1));
for i=1:N
    wz(i)=z(i)/zsum;
end
for i=1:N
    rwz(i)=round(wz(i)*N*0.5);
end
srwz=sum(rwz(:));
q=N/2-srwz;
for i=1:N
    mrwz(i)=rwz(i);
end
for i=1:q
    mrwz(N/2+1-i)=rwz(N/2+1-i)+1;
end
msrwz=sum(mrwz(:));
k=0;
kk=N/2+1;
for i=1:N/2
    if mrwz(i)~=0
        
        for d=kk:kk+mrwz(i)
            mpm(d,:)=mpm(i,:);
        end
        kk=kk+mrwz(i);
        
    end
    
end
for i=1:N
    for j=1:19
        propm(i,j)=mpm(i,j);
    end
end
%calculate weights after SIR
for i=1:N
    part=propm(i,:);
    part(20)=thiessen(1);
    part(21)=pet120022009(1);
    [surf(i),base(i),tot(i),s1(i),s2(i),s3(i),s4(i),s5(i),s6(i)]=sacsmafunction(part);%simulated runoff
    propm(i,14)=s1(i);%update state variable values
    propm(i,15)=s2(i);
    propm(i,16)=s3(i);
    propm(i,17)=s4(i);
    propm(i,18)=s5(i);
    propm(i,19)=s6(i);
    
    if isnan(simrun(i))
        like(i,1)=0;
    end
    difference(i,1)=tot(i)-USGS0706930520022009mm(1,1);
    
end
%calculate standard deviation of difference between simulated and observed runoff
difvar=var(difference(:,1));
for i=1:N
    like(i,1)=(1/sqrt(2*pi*difvar))*exp(-(1/(2*difvar))*(difference(i,1))^2);%calculate the likelihood
end
for i=1:N
    likesum(i)=sum(like(:,1));
end
for i=1:N
    weight(i,1)=like(i,1)/likesum(i);%it might be zero
    
end
%%
for j=1:13
    for i=1:N
        weipar(i,j)=weight(i,1)*pm(i,j);
    end
end
for j=1:13
    mu(j,1)=sum(weipar(:,j));% part 1 of section 4 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
end
for j=1:13
    for i=1:N
        sigcal(i,j)=weight(i,1)*(pm(i,j)-mu(j,1))^2;
    end
end
for j=1:13
    sig(j,1)=sqrt(sum(sigcal(:,j)));% part 1 of section 4 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
end
for j=1:13
    for i=1:N
        curweight(i,j,1)=normpdf(pm(i,j),mu(j,1),sig(j,1))*like(i);% part 2 of section 4 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
    end
end
% create proposal parameters (Moradkhani et al, 2012)

for j=1:13
    sva(j)=(0.25*va(j))^0.5;
end

for j=1:13
    for i=1:N
        newpropm(i,j)=propm(i,j)+normrnd(0,sva(j));% section 3 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
        if newpropm(i,j)>=range(j,3)
            newpropm(i,j)=propm(i,j);
            
        end
        if newpropm(i,j)<=range(j,2)
            newpropm(i,j)=propm(i,j);
            
        end
    end
    
end
for i=1:N
    part=newpropm(i,:);
    part(20)=thiessen(1);
    part(21)=pet120022009(1);
    [surf(i),base(i),tot(i),s1(i),s2(i),s3(i),s4(i),s5(i),s6(i)]=sacsmafunction(part);%simulated runoff
    newpropm(i,14)=s1(i);%update state variable values
    newpropm(i,15)=s2(i);
    newpropm(i,16)=s3(i);
    newpropm(i,17)=s4(i);
    newpropm(i,18)=s5(i);
    newpropm(i,19)=s6(i);
    
    if isnan(simrun(i))
        newlike(i)=0;
    end
    difference(i,1)=tot(i)-USGS0706930520022009mm(1,1);
    
end
%calculate standard deviation of difference between simulated and observed runoff
difvar=var(difference(:,1));
for i=1:N
    newlike(i,1)=(1/sqrt(2*pi*difvar))*exp(-(1/(2*difvar))*(difference(i,1))^2);%calculate the likelihood
end
for j=1:13
    for i=1:N
        newcurweight(i,j,1)=normpdf(newpropm(i,j),mu(j,1),sig(j,1))*newlike(i);% part 4 of section 4 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
    end
end
for j=1:13
    for i=1:N
        if newcurweight(i,j,1)>curweight(i,j,1)%  section 5 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
            propm(i,j)=newpropm(i,j);
        end
    end
end
for i=1:N
    part=propm(i,:);
    part(20)=thiessen(1);
    part(21)=pet120022009(1);
    [surf(i),base(i),tot(i),s1(i),s2(i),s3(i),s4(i),s5(i),s6(i)]=sacsmafunction(part);%simulated runoff
    propm(i,14)=s1(i);%update state variable values
    propm(i,15)=s2(i);
    propm(i,16)=s3(i);
    propm(i,17)=s4(i);
    propm(i,18)=s5(i);
    propm(i,19)=s6(i);
    
end
for i=1:N
    for j=1:19
        pm(i,j)=propm(i,j);
    end
end

%%
%calculate effective sample size (Moradkhani et al, 2012)
neff=zeros(1,T);

%%
neff=zeros(1,T);
% no=zeros(19,500);
svacon=zeros(19,1);
for i=1:N
    preweight(i,1)=weight(i,1);
end
pmplus(:,:)=pm(:,:);
%optimized parameters values

%Go through
%time!!!!!!!!!!!!!!!!!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=20;
for w=1:1
    for t=2:500
        t
        %%omitimized parameters
        for i=1:19
            meanvector(i)=mean(pm(:,i),'omitnan');
        end
        fmeanv(t,:)=meanvector(:);
        meanvector=abs(real(meanvector));
%         noisedrain(t)=thiessen(t)+normrnd(0,0.2*thiessen(t));
%         noisedet(t)=pet120022009(t)+normrnd(0,0.2*pet120022009(t));
        noisedrain(t)=thiessen(t);
        noisedet(t)=pet120022009(t);
      
for i=1:6
            lb1(i)=0.8*meanvector(i+13);
            ub1(i)=1.2*meanvector(i+13);
        end
                
        options=gaoptimset('InitialPopulation',meanvector(14:19));
        [xxx,fval1]=ga(@Fit8_Revised,6,[],[],[],[],lb1,ub1,[],options);

                for i=1:13
            lb(i)=0.8*meanvector(i);
            ub(i)=1.2*meanvector(i);
        end
        options=gaoptimset('InitialPopulation',meanvector(1:13));
        [xx,fval2]=ga(@Fit7_Revised,13,[],[],[],[],lb,ub,[],options);
        xx=real(xx);
        fvalue(t)=fval2;
%         if t>=15
%             fvalue(t)=fvalue(t)+10;
%         end
        pm(1,1:13)=xx(1:13);
        pm(1,14:19)=xxx(1:6);
%             pm(1,20)=noisedrain(t);
%             pm(1,21)=noisedet(t);
%             [surf,base,totil,s1,s2,s3,s4,s5,s6]=sacsmafunction(pm(1,:));
%%
        %Run Model, create x, i, t, -
        pmmatrix(:,:,t)=pm(:,:);
        for i=1:N
            part=pm(i,:);
            part(20)=noisedrain(t);
            part(21)=noisedet(t);
            [surf(i),base(i),tot(i),s1(i),s2(i),s3(i),s4(i),s5(i),s6(i)]=sacsmafunction(part);%simulated runoff;% part 1 of section 1 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
            pm(i,14)=s1(i);%update state variable values
            pm(i,15)=s2(i);
            pm(i,16)=s3(i);
            pm(i,17)=s4(i);
            pm(i,18)=s5(i);
            pm(i,19)=s6(i);
            
        end
        firstgarunoff(t)=tot(1);
        %create y, i , t , /
        %         for i=1:N
        %             part=pm(i,:);
        % part(20)=noisedrain(t);
        % part(21)=noisedet(t);
        %             [surf(i),base(i),tot(i),s1(i),s2(i),s3(i),s4(i),s5(i),s6(i)]=sacsmafunction(part);
        %         end
        save('tot')
        garunoff(t)=tot(1);
        %         options=gaoptimset('generations', 1300);
        %         qr=ga(@Fit3_Revised,1,[],[],[],[],0,0.3,[]);
        for i=1:N
            difference(i,t)=tot(i)-USGS0706930520022009mm(t,1);
        end
        %         for i=1:N
        %             tot(i)=tot(i)+normrnd(0,qr*tot(i));
        %         end
        for i=1:N
            finaltot(i,t)=tot(i);
        end
        for i=1:N
            for j=1:19
                pmneg(i,j,t)=pm(i,j);
            end
        end
        
        
        tot25=quantile(tot,0.025);
        tot975=quantile(tot,0.975);
        tot250=quantile(tot,0.25);
        tot400=quantile(tot,0.40);
        for i=1:N
            nantot(i)=tot(i);
        end
        if thiessen(t)<= pr99
            for i=1:N
                if tot(i)<=tot25
                    nantot(i)=nan;
                end
                if tot(i)>=tot975
                    nantot(i)=nan;
                end
            end
        else
            for i=1:N
                if tot(i)<=tot250
                    nantot(i)=nan;
                end
                if tot(i)>=tot400
                    nantot(i)=nan;
                end
            end
        end
        
        preada(t)=mean(tot);
        ADA(t)=mean(nantot,'omitnan');
        %         SORTEDADA(t)=ADA(t);
        ADA2(t)=ADA(t);
        hist(pm(:,1))
        runoffdifference(t)=abs(ADA(t)-USGS0706930520022009mm(t,1));
        difvar=var(difference(:,t),'omitnan');
        finaldifvar(t)=difvar;
        ss=1000000;
        if finaldifvar(t)<=0.01
            ss
        end
        %%

        %%
        for i=1:N
            weight(i,t)=(1/sqrt(2*pi*difvar))*exp(-(1/(2*difvar))*(difference(i,t))^2);%calculate the likelihood
       if isnan (weight(i,t))
           weight(i,t)=0;
       end
        end
        
        for i=1:N
            mm(i,t)=real(weight(i,t));
        end
        %% SIR ALGORITHM FOR UPGRADE NSE PARAMETER
        %calculate sorted weight matrix for SIR algorithm (Moradkhani et al, 2005)
        for i=1:N
            weight(i,t)=mm(i,t)/sum(mm(:,t));
        end
        sumweight(t)=sum(weight(:,t));
        
        for i=1:N
            x(i)=weight(i,t);
        end
        roundx=round(x,4);
        y=sort(x);
        z=zeros(N,1);
        for i=1:N/2
            z(i)=y(N/2+i);
        end
        for i=1:N
            sas(i)=roundx(i)-0.0138;
            if sas(i)==0
                i
            end
        end
        mpm=zeros(N,19);
        for i=1:N/2
            aaa=z(i);
            for u=1:N
                diffe(u)=aaa-x(u);
                if diffe(u)==0
                    mpm(i,:)=pm(u,:);
                    statempm(i,:)=pmmatrix(u,:,t);
                    sortesfinaltot(i,t)=finaltot(u,t);
                    sortot(i)=finaltot(u,t);
                end
            end
        end
        %%
        sortot25=quantile(sortot,0.025);
        sortot975=quantile(sortot,0.975);
        sortot250=quantile(sortot,0.25);
        sortot400=quantile(sortot,0.40);
        for i=1:N/2
            nansortot(i)=sortot(i);
            sortedrunoff(:,t)=sortot(:);
        end
        if thiessen(t)<= pr99
            for i=1:N/2
                if sortot(i)<=sortot25
                    nansortot(i)=nan;
                end
                if sortot(i)>=sortot975
                    nansortot(i)=nan;
                end
            end
        else
            for i=1:N/2
                if sortot(i)<=sortot250
                    nansortot(i)=nan;
                end
                if sortot(i)>=sortot400
                    nansortot(i)=nan;
                end
            end
        end
        SORTEDADA(t)= mean(nansortot,'omitnan');
        SORTEDADA2(t)=mean(real(sortot(248:250)),'omitnan');
        for i=1:100
            fsort(i)=sortot(i+150);
        end
        q1fsort=quantile(fsort,0.05);
        q2fsort=quantile(fsort,0.95);
        for i=1:99
            if fsort(i)<=q1fsort
                fsort(i)=nan;
            end
            if fsort(i)>=q2fsort
                fsort(i)=nan;
            end
        end
        SORTEDADA3(t)=mean(real(fsort),'omitnan');
        %%
        %ONE DAY AHEAD PREDICTION
                for i=151:250
            part(1:13)=mpm(i,1:13);
            part(14:19)=statempm(i,14:19);
            part(20)=thiessen(t+1);
            if  thiessen(t+1)>=23
                part(12)=0.2
            end
%                         if  thiessen(t+1)-thiessen(t)>=30
%                 part(12)=0.055
%             end
%             kkk=part(20)
            part(21)=pet120022009(t+1);
            [prsurf(i-150),prbase(i-150),prtot(i-150),prs1(i-150),prs2(i-150),prs3(i-150),prs4(i-150),prs5(i-150),prs6(i-150)]=sacsmafunction(part);%simulated runoff;% part 1 of section 1 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
            
                end
                onedaytot(:,t)=prtot(:);
                onedaytot=real(onedaytot);
                oneday1fsort=quantile(prtot,0.05);
        oneday2fsort=quantile(prtot,0.95);
        for i=1:100
            if prtot(i)<=oneday1fsort
                prtot(i)=nan;
            end
            if prtot(i)>=oneday2fsort
                prtot(i)=nan;
            end
        end
                ONEDAYPR(t+1)=mean(real(prtot),'omitnan');

        %%
        finalrunoffda(:,t)=nansortot(:);
        tsorteddata=real(SORTEDADA(t));
        qrr=ga(@Fit4_Revised,1,[],[],[],[],0,0.3);
        
        noiserunoff(t)=abs(normrnd(0,qrr*USGS0706930520022009mm(t)));
        if SORTEDADA(t)>=USGS0706930520022009mm(t)
            runoff(t)=SORTEDADA(t)-noiserunoff(t);
        end
        if SORTEDADA(t)<=USGS0706930520022009mm(t)
            runoff(t)=SORTEDADA(t)+noiserunoff(t);
            
        end
        %%
        
        %%
        %calculate initial weight before SIR
        for i=1:N
            %         weight(i,t)=like(i,t)/likesum(i);%%%%%%%%%important change
            mw(i,t)=weight(i,t)*weight(i,t-1);
        end
        for i=1:N
            finalsurf(i,t)=surf(i);
        end
        for i=1:N
            if sum(mw(:,t))~=0
                weight(i,t)=weight(i,t)*weight(i,t-1)/sum(mw(:,t));%calculate weight of each particle;% part 2 of section 1 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
            else
                weight(i,t)=1/N;
            end
        end
        
        
        vectorada=[];
        
        for i=1:N
            sweight(i,t)=weight(i,t)*weight(i,t);
        end
        sumsweight(t)=sum(sweight(:,t));
        
        neff(t)=1/sum(sweight(:,t));% part 3 of section 1 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
        %%
        
        %%
        %         if neff(t)<=N/2 % end of section 1 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
        if neff(t)<=N/2
            for i=1:13
                varianceparameter(i)=var(pm(:,i));
            end
            
            %calculate sorted weight matrix for SIR algorithm (Moradkhani et al, 2005)
            for i=1:N
                weight(i,t)=mm(i,t)/sum(mm(:,t));
            end
            sumweight(t)=sum(weight(:,t));
            
            for i=1:N
                x(i)=weight(i,t);
            end
            roundx=round(x,4);
            y=sort(x);
            z=zeros(N,1);
            for i=1:N/2
                z(i)=y(N/2+i);
            end
            for i=1:N
                sas(i)=roundx(i)-0.0138;
                if sas(i)==0
                    i
                end
            end
            mpm=zeros(N,19);
            for i=1:N/2
                aaa=z(i);
                for u=1:N
                    diffe(u)=aaa-x(u);
                    if diffe(u)==0
                        mpm(i,:)=pm(u,:);
                        sortesfinaltot(i,t)=finaltot(u,t);
                        sortot(i)=finaltot(u,t);
                    end
                end
            end
            %%
            sortot25=quantile(sortot,0.025);
            sortot975=quantile(sortot,0.975);
            sortot250=quantile(sortot,0.25);
            sortot400=quantile(sortot,0.40);
            for i=1:N/2
                nansortot(i)=sortot(i);
                sortedrunoff(:,t)=sortot(:);
                
            end
            %         if thiessen(t)<= pr99
            %         for i=1:N/2
            %             if sortot(i)<=sortot25
            %                 nansortot(i)=nan;
            %             end
            %             if sortot(i)>=sortot975
            %                 nansortot(i)=nan;
            %             end
            %         end
            %         else
            %                    for i=1:N/2
            %             if sortot(i)<=sortot250
            %                 nansortot(i)=nan;
            %             end
            %             if sortot(i)>=sortot400
            %                 nansortot(i)=nan;
            %             end
            %         end
            %         end
            SORTEDADA(t)= mean(nansortot,'omitnan');
            SORTEDADA2(t)=mean(nansortot(248:250));
            SORTEDADA3(t)=mean(nansortot(150:250));

            finalrunoffda(:,t)=nansortot(:);
            tsorteddata=real(SORTEDADA(t));
            qrr=ga(@Fit4_Revised,1,[],[],[],[],0,0.3);
            noiserunoff(t)=abs(normrnd(0,qrr*USGS0706930520022009mm(t)));
            if SORTEDADA(t)>=USGS0706930520022009mm(t)
                runoff(t)=SORTEDADA(t)-noiserunoff(t);
            end
            if SORTEDADA(t)<=USGS0706930520022009mm(t)
                runoff(t)=SORTEDADA(t)+noiserunoff(t);
                
            end
            %%
            zsum=sum(z(:,1));
            zsum=real(zsum);
            for i=1:N
                wz(i)=z(i)/zsum;
            end
            for i=1:N
                rwz(i)=round(wz(i)*N*0.5);
            end
            srwz=sum(rwz(:));
            q=N/2-srwz;
            q=real(q);
            for i=1:N
                mrwz(i)=rwz(i);
            end
            for i=1:q
                mrwz(N/2+1-i)=rwz(N/2+1-i)+1;
            end
            msrwz=sum(mrwz(:));
            k=0;
            kk=N/2+1;
            mrwz=real(mrwz);
            for i=1:N/2
                if mrwz(i)~=0
                    
                    for d=kk:kk+mrwz(i)
                        mpm(d,:)=mpm(i,:);
                    end
                    kk=kk+mrwz(i);
                    
                end
                
            end
            for i=1:N
                for j=1:19
                    pmplus(i,j)=mpm(i,j);
                end
            end
            %%%%%NOINSE ON PARAMETERS
            for i=1:N
                for j=1:13
                    pmnewplus(i,j)=pmplus(i,j)+normrnd(0,0.2*pmplus(i,j));
                    %                     if  pmnewplus(i,j)<=range(j,2)|pmnewplus(i,j)>=range(j,3)
                    %                         while pmnewplus(i,j)<=range(j,2)|  pmnewplus(i,j)>=range(j,3)
                    %                          pmnewplus(i,j)=pmplus(i,j)+normrnd(0,0.2*pmplus(i,j));
                    %                         end
                    %                     end
                    
                end
            end
            for i=1:N
                for j=1:13
                    pmplus(i,j)=pmnewplus(i,j);
                end
            end
            %run model for calculate like plus
            for i=1:N
                part=pmplus(i,:);
part(20)=noisedrain(t);
                part(21)=noisedet(t);
                [surf(i),base(i),tot(i),s1(i),s2(i),s3(i),s4(i),s5(i),s6(i)]=sacsmafunction(part);%simulated runoff
                pm(i,14)=s1(i);%update state variable values
                pmplus(i,15)=s2(i);
                pmplus(i,16)=s3(i);
                pmplus(i,17)=s4(i);
                pmplus(i,18)=s5(i);
                pmplus(i,19)=s6(i);
                difference(i,t)=tot(i)+normrnd(0,0.05*tot(i))-USGS0706930520022009mm(t,1);
                
            end
            %             save('tot')
            %             x=ga(@Fit3_Revised,1,[],[],[],[],0,1);
            for i=1:N
                difference(i,t)=tot(i)-USGS0706930520022009mm(t,1);
            end
            %calculate standard deviation of difference between simulated and observed runoff
            difvar=var(difference(:,t),'omitnan');
            for i=1:N
                likeplus(i,t)=(1/sqrt(2*pi*difvar))*exp(-(1/(2*difvar))*(difference(i,t))^2);%calculate the likelihood
            end
            
            %%
            %             % create proposal parameters
            for j=1:13
                sva(j,t)=(0.3*varianceparameter(j))^0.5;
            end
            for j=1:13
                no(j,t)=normrnd(0,sva(j,t));
                for i=1:N
                    normrndi(i,j)=normrnd(0,sva(j,t));
                    newpropm(i,j)=pmplus(i,j)+normrndi(i,j);% section 3 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
                    if newpropm(i,j)>=range(j,3)
                        newpropm(i,j)=pmplus(i,j);
                    end
                    if newpropm(i,j)<=range(j,2)
                        newpropm(i,j)=pmplus(i,j);
                    end
                end
            end
            
            
            for j=14:19
                for i=1:N
                    newpropm(i,j)=pmplus(i,j);
                end
            end
            for j=1:13%%%change 19 to 13
                for i=1:N
                    weipar(i,j)=weight(i,t-1)*pmneg(i,j);
                end
            end
            for j=1:13
                mu(j,t)=sum(weipar(1:N,j)); % part 1 of section 4 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
            end
            for j=1:13
                for i=1:N
                    sigcal(i,j,t)=weight(i,t-1)*(pmneg(i,j)-mu(j,t))^2;
                end
            end
            for j=1:13
                sig(j,t)=sqrt(sum(sigcal(:,j,t)));% part 1 of section 4 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
            end
            
            %%
            for i=1:N
                likeplust(i,t)=likeplus(i,t);
            end
            for j=1:13
                for i=1:N
                    curweight(i,j,t)=normpdf(pmplus(i,j),mu(j,t),sig(j,t))*likeplus(i,t); % part 2 of section 4 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
                end
            end
            %run model for proposal parameters
            for i=1:N
                part=newpropm(i,:);% part 3 of section 4 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
part(20)=noisedrain(t);
                part(21)=noisedet(t);
                [surf(i),base(i),tot(i),s1(i),s2(i),s3(i),s4(i),s5(i),s6(i)]=sacsmafunction(part);%simulated runoff
                newpropm(i,14)=s1(i);%update state variable values
                newpropm(i,15)=s2(i);
                newpropm(i,16)=s3(i);
                newpropm(i,17)=s4(i);
                newpropm(i,18)=s5(i);
                newpropm(i,19)=s6(i);
                
            end
%             for i=1:N
%                 part=newpropm(i,:);% part 3 of section 4 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
% part(20)=noisedrain(t);
%                 part(21)=noisedet(t);
%                 [surf(i),base(i),tot(i),s1(i),s2(i),s3(i),s4(i),s5(i),s6(i)]=sacsmafunction(part);%simulated runoff
%                 newpropm(i,14)=s1(i);%update state variable values
%                 newpropm(i,15)=s2(i);
%                 newpropm(i,16)=s3(i);
%                 newpropm(i,17)=s4(i);
%                 newpropm(i,18)=s5(i);
%                 newpropm(i,19)=s6(i);
%             end
            %             save('tot')
            %             x=ga(@Fit3_Revised,1,[],[],[],[],0,1);
            for i=1:N
                difference(i,t)=tot(i)-USGS0706930520022009mm(t,1);
            end
            for i=1:N
                newprtot(i,t)=tot(i);
            end
            for i=1:N
                prodiffer(i,t)=difference(i,t);
            end
            %calculate standard deviation of difference between simulated and observed runoff
            difvar=var(difference(:,t),'omitnan');
            
            for i=1:N
                likeproposal(i,t)=(1/sqrt(2*pi*difvar))*exp(-(1/(2*difvar))*(difference(i,t))^2);%calculate the likelihood
            end
            for j=1:13
                for i=1:N
                    weightproposal(i,j,t)=normpdf(newpropm(i,j),mu(j,t),sig(j,t))*likeproposal(i,t);% part 4 of section 4 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
                end
            end
            sss=zeros(19,100);
            for j=1:13%faghat parametrha jaygozin mishavand
                
                for i=1:N
                    
                    if weightproposal(i,j,t)>curweight(i,j,t)% section 5 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
                        
                        pmplus(i,j)=newpropm(i,j);
                    end
                end
                
            end
            %%
            %yekbar zarat jaygozin run mishavand ta state variableha
            %beroozresani shavand
            for i=1:N
                part=pmplus(i,:);% part 3 of section 4 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
part(20)=noisedrain(t);
                part(21)=noisedet(t);
                [surf(i),base(i),tot(i),s1(i),s2(i),s3(i),s4(i),s5(i),s6(i)]=sacsmafunction(part);%simulated runoff
                pmplus(i,14)=s1(i);%update state variable values
                pmplus(i,15)=s2(i);
                pmplus(i,16)=s3(i);
                pmplus(i,17)=s4(i);
                pmplus(i,18)=s5(i);
                pmplus(i,19)=s6(i);
                
            end
            %             save('tot')
            %             x=ga(@Fit3_Revised,1,[],[],[],[],0,1);
            %%modification require if ada2 is important
            tot25=quantile(tot,0.025);
            tot975=quantile(tot,0.975);
            for i=1:N
                nantot(i)=tot(i);
            end
            for i=1:N
                if tot(i)<=tot25
                    nantot(i)=nan;
                end
                if tot(i)>=tot975
                    nantot(i)=nan;
                end
            end
            ADA2(t)=mean(nantot,'omitnan');
            %zarrat baraye gam e baadi jaygozin mishavand
            pmplusmatrix(:,:,t)=pmplus(:,:);
            for i=1:N
                for j=1:19
                    pm(i,j)=pmplus(i,j);
                end
            end
            for i=1:N
                weight(i,t)=1/N;
            end
        end
        %%
        %END of ifffffffffffffffffffffffffffffffffffff
        %         for i=1:N
        %             preweight(i,t)=weight(i,t);
        %         end
        %%
        %quantile calculations
        y5p1(t)=quantile(propm(:,1),0.05);
        y95p1(t)=quantile(propm(:,1),0.95);
        y25p1(t)=quantile(propm(:,1),0.25);
        y75p1(t)=quantile(propm(:,1),0.75);
        y5p2(t)=quantile(propm(:,2),0.05);
        y95p2(t)=quantile(propm(:,2),0.95);
        y25p2(t)=quantile(propm(:,2),0.25);
        y75p2(t)=quantile(propm(:,2),0.75);
        y5p3(t)=quantile(propm(:,3),0.05);
        y95p3(t)=quantile(propm(:,3),0.95);
        y25p3(t)=quantile(propm(:,3),0.25);
        y75p3(t)=quantile(propm(:,3),0.75);
        y5p4(t)=quantile(propm(:,4),0.05);
        y95p4(t)=quantile(propm(:,4),0.95);
        y25p4(t)=quantile(propm(:,4),0.25);
        y75p4(t)=quantile(propm(:,4),0.75);
        if t>=6
            aa=fvalue(t);
                    bb=fvalue(t-1);
            cc=fvalue(t-2);
            dd=fvalue(t-3);
            ee=fvalue(t-4);
        
        if aa>=0.1 & bb>=0.1 & cc>=0.1 & dd>=0.1 & ee >=0.1
 222222
 d=t;
 for s=t-10:t
     s
     
 if s==t-10
            for i=1:19
            lb2(i)=0.8*pmmatrix(1,i,s);
            ub2(i)=1.2*pmmatrix(1,i,s);
 end
        options=gaoptimset('InitialPopulation',fmeanv(s,1:19));
        [yy,fval3]=ga(@Fit9_Revised,19,[],[],[],[],lb2,ub2,[],options);
       yy(:)=opt(:);
       for j=1:19
    for i=1:N
        pm(i,j)=opt(j);
    end
end

for i=250:-1:1
    pm(i,1)=pm(i,1)-(250-i)*0.05;
    pm(i,2)=pm(i,2)-(250-i)*0.05;
    pm(i,3)=pm(i,3)-(250-i)*0.05;
    pm(i,4)=pm(i,4)-(250-i)*0.05;
    pm(i,5)=pm(i,5)-(250-i)*0.05;
    pm(i,6)=pm(i,6)-(250-i)*0.0005;
    pm(i,7)=pm(i,7)-(250-i)*0.001;
    pm(i,8)=pm(i,8)-(250-i)*0.00001;
    pm(i,9)=pm(i,9)-(250-i)*0.0001;
    pm(i,10)=pm(i,10)-(250-i)*0.05;
    pm(i,11)=pm(i,11)-(250-i)*0.005;
    pm(i,12)=pm(i,12)-(250-i)*0.0001;
    pm(i,13)=pm(i,13)-(250-i)*0.0005;
    pm(i,14)=pm(i,14)+abs(250-i)*0.05;
    pm(i,15)=pm(i,15)+abs(250-i)*0.05;
    pm(i,16)=pm(i,16)-(250-i)*0.05;
    pm(i,17)=pm(i,17)-(250-i)*0.05;
    pm(i,18)=pm(i,18)+abs(250-i)*0.05;
    pm(i,19)=pm(i,19)-(250-i)*0.0005;
end
for i=250:500
    pm(i,1)=pm(i,1)+abs(250-i)*0.05;
    pm(i,2)=pm(i,2)+abs(250-i)*0.05;
    pm(i,3)=pm(i,3)+abs(250-i)*0.05;
    pm(i,4)=pm(i,4)+abs(250-i)*0.05;
    pm(i,5)=pm(i,5)+abs(250-i)*0.05;
    pm(i,6)=pm(i,6)+abs(250-i)*0.0005;
    pm(i,7)=pm(i,7)+abs(250-i)*0.001;
    pm(i,8)=pm(i,8)+abs(250-i)*0.00001;
    pm(i,9)=pm(i,9)+abs(250-i)*0.0001;
    pm(i,10)=pm(i,10)+abs(250-i)*0.05;
    pm(i,11)=pm(i,11)+abs(250-i)*0.005;
    pm(i,12)=pm(i,12)+abs(250-i)*0.0001;
    pm(i,13)=pm(i,13)+abs(250-i)*0.0001;
    pm(i,14)=pm(i,14)+abs(250-i)*0.05;
    pm(i,15)=pm(i,15)+abs(250-i)*0.05;
    pm(i,16)=pm(i,16)+abs(250-i)*0.05;
    pm(i,17)=pm(i,17)+abs(250-i)*0.05;
    pm(i,18)=pm(i,18)+abs(250-i)*0.05;
    pm(i,19)=pm(i,19)+abs(250-i)*0.0005;
    
end
 end
        for i=1:19
            meanvector(i)=mean(pm(:,i),'omitnan');
        end
        fmeanv(s,:)=meanvector(:);
        meanvector=abs(real(meanvector));
%         noisedrain(t)=thiessen(t)+normrnd(0,0.2*thiessen(t));
%         noisedet(t)=pet120022009(t)+normrnd(0,0.2*pet120022009(t));
        noisedrain(s)=thiessen(s);
        noisedet(s)=pet120022009(s);
      
for i=1:6
            lb1(i)=0.8*meanvector(i+13);
            ub1(i)=1.2*meanvector(i+13);
        end
                
        options=gaoptimset('InitialPopulation',meanvector(14:19));
        [xxx,fval1]=ga(@Fit18_Revised,6,[],[],[],[],lb1,ub1,[],options);

                for i=1:13
            lb(i)=0.8*meanvector(i);
            ub(i)=1.2*meanvector(i);
        end
        options=gaoptimset('InitialPopulation',meanvector(1:13));
        [xx,fval2]=ga(@Fit17_Revised,13,[],[],[],[],lb,ub,[],options);
        xx=real(xx);
        fvalue(s)=fval2;
        pm(1,1:13)=xx(1:13);
        pm(1,14:19)=xxx(1:6);
%             pm(1,20)=noisedrain(t);
%             pm(1,21)=noisedet(t);
%             [surf,base,totil,s1,s2,s3,s4,s5,s6]=sacsmafunction(pm(1,:));
%%
        %Run Model, create x, i, t, -
        pmmatrix(:,:,s)=pm(:,:);
        for i=1:N
            part=pm(i,:);
            part(20)=noisedrain(s);
            part(21)=noisedet(s);
            [surf(i),base(i),tot(i),s1(i),s2(i),s3(i),s4(i),s5(i),s6(i)]=sacsmafunction(part);%simulated runoff;% part 1 of section 1 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
            pm(i,14)=s1(i);%update state variable values
            pm(i,15)=s2(i);
            pm(i,16)=s3(i);
            pm(i,17)=s4(i);
            pm(i,18)=s5(i);
            pm(i,19)=s6(i);
            
        end
        firstgarunoff(s)=tot(1);
        %create y, i , t , /
        %         for i=1:N
        %             part=pm(i,:);
        % part(20)=noisedrain(t);
        % part(21)=noisedet(t);
        %             [surf(i),base(i),tot(i),s1(i),s2(i),s3(i),s4(i),s5(i),s6(i)]=sacsmafunction(part);
        %         end
        save('tot')
        garunoff(s)=tot(1);
        %         options=gaoptimset('generations', 1300);
        %         qr=ga(@Fit3_Revised,1,[],[],[],[],0,0.3,[]);
        for i=1:N
            difference(i,s)=tot(i)-USGS0706930520022009mm(s,1);
        end
        %         for i=1:N
        %             tot(i)=tot(i)+normrnd(0,qr*tot(i));
        %         end
        for i=1:N
            finaltot(i,s)=tot(i);
        end
        for i=1:N
            for j=1:19
                pmneg(i,j,s)=pm(i,j);
            end
        end
        
        
        tot25=quantile(tot,0.025);
        tot975=quantile(tot,0.975);
        tot250=quantile(tot,0.25);
        tot400=quantile(tot,0.40);
        for i=1:N
            nantot(i)=tot(i);
        end
        if thiessen(s)<= pr99
            for i=1:N
                if tot(i)<=tot25
                    nantot(i)=nan;
                end
                if tot(i)>=tot975
                    nantot(i)=nan;
                end
            end
        else
            for i=1:N
                if tot(i)<=tot250
                    nantot(i)=nan;
                end
                if tot(i)>=tot400
                    nantot(i)=nan;
                end
            end
        end
        
        preada(s)=mean(tot);
        ADA(s)=mean(nantot,'omitnan');
        %         SORTEDADA(t)=ADA(t);
        ADA2(s)=ADA(s);
        hist(pm(:,1))
        runoffdifference(s)=abs(ADA(s)-USGS0706930520022009mm(s,1));
        difvar=var(difference(:,s),'omitnan');
        finaldifvar(s)=difvar;
        ss=1000000;
        if finaldifvar(s)<=0.01
            ss
        end
        %%

        %%
        for i=1:N
            weight(i,s)=(1/sqrt(2*pi*difvar))*exp(-(1/(2*difvar))*(difference(i,s))^2);%calculate the likelihood
               if isnan(weight(i,s))
           weight(i,s)=0;
       end
        end
        for i=1:N
            mm(i,s)=real(weight(i,s));
        end
        %% SIR ALGORITHM FOR UPGRADE NSE PARAMETER
        %calculate sorted weight matrix for SIR algorithm (Moradkhani et al, 2005)
        for i=1:N
            weight(i,s)=mm(i,s)/sum(mm(:,s));
        end
        sumweight(s)=sum(weight(:,s));
        
        for i=1:N
            x(i)=weight(i,s);
        end
        roundx=round(x,4);
        y=sort(x);
        z=zeros(N,1);
        for i=1:N/2
            z(i)=y(N/2+i);
        end
        for i=1:N
            sas(i)=roundx(i)-0.0138;
            if sas(i)==0
                i
            end
        end
        mpm=zeros(N,19);
        for i=1:N/2
            aaa=z(i);
            for u=1:N
                diffe(u)=aaa-x(u);
                if diffe(u)==0
                    mpm(i,:)=pm(u,:);
                    sortesfinaltot(i,t)=finaltot(u,s);
                    sortot(i)=finaltot(u,s);
                end
            end
        end
        %%
        sortot25=quantile(sortot,0.025);
        sortot975=quantile(sortot,0.975);
        sortot250=quantile(sortot,0.25);
        sortot400=quantile(sortot,0.40);
        for i=1:N/2
            nansortot(i)=sortot(i);
            sortedrunoff(:,s)=sortot(:);
        end
        if thiessen(s)<= pr99
            for i=1:N/2
                if sortot(i)<=sortot25
                    nansortot(i)=nan;
                end
                if sortot(i)>=sortot975
                    nansortot(i)=nan;
                end
            end
        else
            for i=1:N/2
                if sortot(i)<=sortot250
                    nansortot(i)=nan;
                end
                if sortot(i)>=sortot400
                    nansortot(i)=nan;
                end
            end
        end
        SORTEDADA(s)= mean(nansortot,'omitnan');
        SORTEDADA2(s)=mean(real(sortot(248:250)),'omitnan');
        for i=1:100
            fsort(i)=sortot(i+150);
        end
        q1fsort=quantile(fsort,0.05);
        q2fsort=quantile(fsort,0.95);
        for i=1:99
            if fsort(i)<=q1fsort
                fsort(i)=nan;
            end
            if fsort(i)>=q2fsort
                fsort(i)=nan;
            end
        end
        if s>=t-6
        SORTEDADA3(s)=mean(real(fsort),'omitnan');
        end
        finalrunoffda(:,s)=nansortot(:);
        tsorteddata=real(SORTEDADA(s));
        qrr=ga(@Fit4_Revised,1,[],[],[],[],0,0.3);
        
        noiserunoff(s)=abs(normrnd(0,qrr*USGS0706930520022009mm(s)));
        if SORTEDADA(s)>=USGS0706930520022009mm(s)
            runoff(s)=SORTEDADA(s)-noiserunoff(s);
        end
        if SORTEDADA(s)<=USGS0706930520022009mm(s)
            runoff(s)=SORTEDADA(s)+noiserunoff(s);
            
        end
        %%
        
        %%
        %calculate initial weight before SIR
        for i=1:N
            %         weight(i,t)=like(i,t)/likesum(i);%%%%%%%%%important change
            mw(i,s)=weight(i,s)*weight(i,s-1);
        end
        for i=1:N
            finalsurf(i,s)=surf(i);
        end
        for i=1:N
            if sum(mw(:,s))~=0
                weight(i,s)=weight(i,s)*weight(i,s-1)/sum(mw(:,s));%calculate weight of each particle;% part 2 of section 1 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
            else
                weight(i,s)=1/N;
            end
        end
        
        
        vectorada=[];
        
        for i=1:N
            sweight(i,s)=weight(i,s)*weight(i,s);
        end
        sumsweight(s)=sum(sweight(:,s));
        
        neff(s)=1/sum(sweight(:,s));% part 3 of section 1 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
        %%
        
        %%
        %         if neff(t)<=N/2 % end of section 1 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
        if neff(s)<=N/2
            for i=1:13
                varianceparameter(i)=var(pm(:,i));
            end
            
            %calculate sorted weight matrix for SIR algorithm (Moradkhani et al, 2005)
            for i=1:N
                weight(i,s)=mm(i,s)/sum(mm(:,s));
            end
            sumweight(s)=sum(weight(:,s));
            
            for i=1:N
                x(i)=weight(i,s);
            end
            roundx=round(x,4);
            y=sort(x);
            z=zeros(N,1);
            for i=1:N/2
                z(i)=y(N/2+i);
            end
            for i=1:N
                sas(i)=roundx(i)-0.0138;
                if sas(i)==0
                    i
                end
            end
            mpm=zeros(N,19);
            for i=1:N/2
                aaa=z(i);
                for u=1:N
                    diffe(u)=aaa-x(u);
                    if diffe(u)==0
                        mpm(i,:)=pm(u,:);
                        sortesfinaltot(i,s)=finaltot(u,s);
                        sortot(i)=finaltot(u,s);
                    end
                end
            end
            %%
            sortot25=quantile(sortot,0.025);
            sortot975=quantile(sortot,0.975);
            sortot250=quantile(sortot,0.25);
            sortot400=quantile(sortot,0.40);
            for i=1:N/2
                nansortot(i)=sortot(i);
                sortedrunoff(:,s)=sortot(:);
                
            end
            %         if thiessen(t)<= pr99
            %         for i=1:N/2
            %             if sortot(i)<=sortot25
            %                 nansortot(i)=nan;
            %             end
            %             if sortot(i)>=sortot975
            %                 nansortot(i)=nan;
            %             end
            %         end
            %         else
            %                    for i=1:N/2
            %             if sortot(i)<=sortot250
            %                 nansortot(i)=nan;
            %             end
            %             if sortot(i)>=sortot400
            %                 nansortot(i)=nan;
            %             end
            %         end
            %         end
            if s>=t-6
            SORTEDADA(s)= mean(nansortot,'omitnan');
            SORTEDADA2(s)=mean(nansortot(248:250));
            SORTEDADA3(s)=mean(nansortot(150:250));
            end
            finalrunoffda(:,s)=nansortot(:);
            tsorteddata=real(SORTEDADA(s));
            qrr=ga(@Fit4_Revised,1,[],[],[],[],0,0.3);
            noiserunoff(s)=abs(normrnd(0,qrr*USGS0706930520022009mm(s)));
            if SORTEDADA(s)>=USGS0706930520022009mm(s)
                runoff(s)=SORTEDADA(s)-noiserunoff(s);
            end
            if SORTEDADA(s)<=USGS0706930520022009mm(s)
                runoff(s)=SORTEDADA(s)+noiserunoff(s);
                
            end
            %%
            zsum=sum(z(:,1));
            zsum=real(zsum);
            for i=1:N
                wz(i)=z(i)/zsum;
            end
            for i=1:N
                rwz(i)=round(wz(i)*N*0.5);
            end
            srwz=sum(rwz(:));
            q=N/2-srwz;
            q=real(q);
            for i=1:N
                mrwz(i)=rwz(i);
            end
            for i=1:q
                mrwz(N/2+1-i)=rwz(N/2+1-i)+1;
            end
            msrwz=sum(mrwz(:));
            k=0;
            kk=N/2+1;
            mrwz=real(mrwz);
            for i=1:N/2
                if mrwz(i)~=0
                    
                    for d=kk:kk+mrwz(i)
                        mpm(d,:)=mpm(i,:);
                    end
                    kk=kk+mrwz(i);
                    
                end
                
            end
            for i=1:N
                for j=1:19
                    pmplus(i,j)=mpm(i,j);
                end
            end
            %%%%%NOINSE ON PARAMETERS
            for i=1:N
                for j=1:13
                    pmnewplus(i,j)=pmplus(i,j)+normrnd(0,0.2*pmplus(i,j));
                    %                     if  pmnewplus(i,j)<=range(j,2)|pmnewplus(i,j)>=range(j,3)
                    %                         while pmnewplus(i,j)<=range(j,2)|  pmnewplus(i,j)>=range(j,3)
                    %                          pmnewplus(i,j)=pmplus(i,j)+normrnd(0,0.2*pmplus(i,j));
                    %                         end
                    %                     end
                    
                end
            end
            for i=1:N
                for j=1:13
                    pmplus(i,j)=pmnewplus(i,j);
                end
            end
            %run model for calculate like plus
            for i=1:N
                part=pmplus(i,:);
part(20)=noisedrain(s);
                part(21)=noisedet(s);
                [surf(i),base(i),tot(i),s1(i),s2(i),s3(i),s4(i),s5(i),s6(i)]=sacsmafunction(part);%simulated runoff
                pm(i,14)=s1(i);%update state variable values
                pmplus(i,15)=s2(i);
                pmplus(i,16)=s3(i);
                pmplus(i,17)=s4(i);
                pmplus(i,18)=s5(i);
                pmplus(i,19)=s6(i);
                difference(i,s)=tot(i)+normrnd(0,0.05*tot(i))-USGS0706930520022009mm(s,1);
                
            end
            %             save('tot')
            %             x=ga(@Fit3_Revised,1,[],[],[],[],0,1);
            for i=1:N
                difference(i,s)=tot(i)-USGS0706930520022009mm(s,1);
            end
            %calculate standard deviation of difference between simulated and observed runoff
            difvar=var(difference(:,s),'omitnan');
            for i=1:N
                likeplus(i,s)=(1/sqrt(2*pi*difvar))*exp(-(1/(2*difvar))*(difference(i,s))^2);%calculate the likelihood
            end
            
            %%
            %             % create proposal parameters
            for j=1:13
                sva(j,s)=(0.3*varianceparameter(j))^0.5;
            end
            for j=1:13
                no(j,s)=normrnd(0,sva(j,s));
                for i=1:N
                    normrndi(i,j)=normrnd(0,sva(j,s));
                    newpropm(i,j)=pmplus(i,j)+normrndi(i,j);% section 3 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
                    if newpropm(i,j)>=range(j,3)
                        newpropm(i,j)=pmplus(i,j);
                    end
                    if newpropm(i,j)<=range(j,2)
                        newpropm(i,j)=pmplus(i,j);
                    end
                end
            end
            
            
            for j=14:19
                for i=1:N
                    newpropm(i,j)=pmplus(i,j);
                end
            end
            for j=1:13%%%change 19 to 13
                for i=1:N
                    weipar(i,j)=weight(i,s-1)*pmneg(i,j);
                end
            end
            for j=1:13
                mu(j,s)=sum(weipar(1:N,j)); % part 1 of section 4 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
            end
            for j=1:13
                for i=1:N
                    sigcal(i,j,s)=weight(i,s-1)*(pmneg(i,j)-mu(j,s))^2;
                end
            end
            for j=1:13
                sig(j,s)=sqrt(sum(sigcal(:,j,s)));% part 1 of section 4 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
            end
            
            %%
            for i=1:N
                likeplust(i,s)=likeplus(i,s);
            end
            for j=1:13
                for i=1:N
                    curweight(i,j,s)=normpdf(pmplus(i,j),mu(j,s),sig(j,s))*likeplus(i,s); % part 2 of section 4 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
                end
            end
            %run model for proposal parameters
            for i=1:N
                part=newpropm(i,:);% part 3 of section 4 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
part(20)=noisedrain(s);
                part(21)=noisedet(s);
                [surf(i),base(i),tot(i),s1(i),s2(i),s3(i),s4(i),s5(i),s6(i)]=sacsmafunction(part);%simulated runoff
                newpropm(i,14)=s1(i);%update state variable values
                newpropm(i,15)=s2(i);
                newpropm(i,16)=s3(i);
                newpropm(i,17)=s4(i);
                newpropm(i,18)=s5(i);
                newpropm(i,19)=s6(i);
                
            end
%             for i=1:N
%                 part=newpropm(i,:);% part 3 of section 4 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
% part(20)=noisedrain(t);
%                 part(21)=noisedet(t);
%                 [surf(i),base(i),tot(i),s1(i),s2(i),s3(i),s4(i),s5(i),s6(i)]=sacsmafunction(part);%simulated runoff
%                 newpropm(i,14)=s1(i);%update state variable values
%                 newpropm(i,15)=s2(i);
%                 newpropm(i,16)=s3(i);
%                 newpropm(i,17)=s4(i);
%                 newpropm(i,18)=s5(i);
%                 newpropm(i,19)=s6(i);
%             end
            %             save('tot')
            %             x=ga(@Fit3_Revised,1,[],[],[],[],0,1);
            for i=1:N
                difference(i,s)=tot(i)-USGS0706930520022009mm(s,1);
            end
            for i=1:N
                newprtot(i,s)=tot(i);
            end
            for i=1:N
                prodiffer(i,s)=difference(i,s);
            end
            %calculate standard deviation of difference between simulated and observed runoff
            difvar=var(difference(:,s),'omitnan');
            
            for i=1:N
                likeproposal(i,s)=(1/sqrt(2*pi*difvar))*exp(-(1/(2*difvar))*(difference(i,s))^2);%calculate the likelihood
            end
            for j=1:13
                for i=1:N
                    weightproposal(i,j,s)=normpdf(newpropm(i,j),mu(j,s),sig(j,s))*likeproposal(i,s);% part 4 of section 4 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
                end
            end
            sss=zeros(19,100);
            for j=1:13%faghat parametrha jaygozin mishavand
                
                for i=1:N
                    
                    if weightproposal(i,j,s)>curweight(i,j,s)% section 5 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
                        
                        pmplus(i,j)=newpropm(i,j);
                    end
                end
                
            end
            %%
            %yekbar zarat jaygozin run mishavand ta state variableha
            %beroozresani shavand
            for i=1:N
                part=pmplus(i,:);% part 3 of section 4 in Figure 2: Flowchart of the PF-MCMC algorithm (Moradkhani et al, 2012)
part(20)=noisedrain(s);
                part(21)=noisedet(s);
                [surf(i),base(i),tot(i),s1(i),s2(i),s3(i),s4(i),s5(i),s6(i)]=sacsmafunction(part);%simulated runoff
                pmplus(i,14)=s1(i);%update state variable values
                pmplus(i,15)=s2(i);
                pmplus(i,16)=s3(i);
                pmplus(i,17)=s4(i);
                pmplus(i,18)=s5(i);
                pmplus(i,19)=s6(i);
                
            end
            %             save('tot')
            %             x=ga(@Fit3_Revised,1,[],[],[],[],0,1);
            %%modification require if ada2 is important
            tot25=quantile(tot,0.025);
            tot975=quantile(tot,0.975);
            for i=1:N
                nantot(i)=tot(i);
            end
            for i=1:N
                if tot(i)<=tot25
                    nantot(i)=nan;
                end
                if tot(i)>=tot975
                    nantot(i)=nan;
                end
            end
            ADA2(s)=mean(nantot,'omitnan');
            %zarrat baraye gam e baadi jaygozin mishavand
            pmplusmatrix(:,:,s)=pmplus(:,:);
            for i=1:N
                for j=1:19
                    pm(i,j)=pmplus(i,j);
                end
            end
            for i=1:N
                weight(i,t)=1/N;
            end
        end
        end
        end
        end
    end
    %     for i=1:500
    %         vectorada(i,1)=ADA(i);
    %     end
    ADA=real(ADA);
    ADA2=real(ADA2);
    SORTEDADA=real(SORTEDADA);
    for i=250:1096
        se1(i)=((runoff(i))-(ADA(i)))^2;
    end
    for i=250:1096
        se2(i)=((runoff(i))-(SORTEDADA(i)))^2;
    end
    sumse2=sum(se2(:),'omitnan');
    rmse=sqrt(sumse2/846)
    for i=250:1096
        are=abs((runoff(i)-SORTEDADA(i))/runoff(i));
    end
    mare=mean(are)
    qave=mean(runoff(250:1096));
    for i=2:1096
        dife(i)=(runoff(i)-qave)^2;
    end
    nse1=1-(sum(se1(250:1096),'omitnan')/sum(dife(250:1096),'omitnan'))
    nse2=1-(sum(se2(250:1096),'omitnan')/sum(dife(250:1096),'omitnan'))
    
    %calculate KGE coefficient
    stdsim=std(SORTEDADA(250:1096),'omitnan');
    stdobs=std(runoff(250:1096),'omitnan');
    
    meansim=mean(SORTEDADA(250:1096),'omitnan');
    meanobs=mean(runoff(250:1096),'omitnan');
    for i=250:1096
        covground(i,1)=runoff(i);
    end
    covmatrix=nancov(covground(250:1096),SORTEDADA(250:1096));
    covratio=(covmatrix(2,2)/covmatrix(1,1));
    stratio=(stdsim/stdobs);
    meanratio=(meansim/meanobs);
    kge(w)=1-sqrt(((covmatrix(2,2)/covmatrix(1,1))-1)^2+((stdsim/stdobs)-1)^2+((meansim/meanobs)-1)^2);
    %calculate person coefficient
    for k=1:847
        corrfile(k,1)=covground(k+249,1);
        corrfile(k,2)=SORTEDADA(1,k+249);
    end
    pearson=corr(corrfile,'rows','complete');
    pearsoncoeff(w)=pearson(1,2);
end
save('ADA')
save('runoff')
save('SORTEDADA')
save('pm')
save ('pmplus')
save('finalrunoffda')%contain nan values
save('finaltot')
save('sortedrunoff')%without nan values
save('pmmatrix')%prior distribution of parameters
save('pmplusmatrix')%posterior distribution of parameters
f(1)=1;
for i=2:1096
    f(i)=f(i-1)+1;
end
plot(f(2:1096),runoff(2:1096),':', f(2:1096),SORTEDADA(2:1096))

toc


