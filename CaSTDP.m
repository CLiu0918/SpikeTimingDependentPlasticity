
% Author Chang Liu
% Project:Calcium-dependent storage of long-term potential induced by 
% spike-timing-dependent plasticity

% Variables:
% EPSP -- the onset of presynaptic spikes
% BPAP -- back-propogated action potential that is injected into the
% postsynaptic neuron soma for stimulating postsynaptic spikes
% preSpike -- presynaptic spikes generated using Possion spike generator.
%          -- with the firing rate =1 Hz.
% Casp -- Ca2+ dynamics induced by single presynaptic spikes stimulation 

% Note: The arrival time of the postsynaptic spike is controlled by the onset of
% BPAP.

clc
clear all
%% Single Presynaptic spike stimulation 
%generate of presynaptic spikes
fr=0.001;   % 1 Hz =0.001 Hz/ms
tSim = 1000; % 10000 ms   10s
nTrials = 1;
[preSpike,tVec]=possionSpikeGen(fr,tSim,nTrials);
spikepos=find(preSpike==1);

%%
EPSP=ones(1,length(tVec))*(-65);
tauep1=50;  tauep2=5; s=1; norm=1;
for count=1:length(tVec)-1
    
    for num=1:count
    
        if preSpike(num) ==1
        
           dv = s/norm *(exp(-(tVec(count)-tVec(num))/tauep1) -...
                            exp(-(tVec(count)-tVec(num))/tauep2))  ;
           EPSP(count)=EPSP(count)+dv;
        end  
    end 
    
 end

% single presynaptic stimulation:
postVol=EPSP;
% NMDA channel current and Ca2+ dynamics
Mg =1; Vr=130;
BV= 1./ (1+exp(-0.062*postVol)*Mg/3.57);
HV= BV.*(postVol-Vr);

P0=0.5; Gnmda = -1/500; 
tauf=50; taus=200; 
If=0.5; Is=0.5;
tauCa=50;   Ca=0;

for count=1:length(tVec)
    if tVec(count)<=spikepos
        Inmda(count)=0;
    else
        Inmda(count) = P0* Gnmda* (If*exp(-tVec(count)/tauf)+...
                        Is*exp(-tVec(count)/taus)).*HV(count);
    end
    dCa = Inmda(count)-(1/tauCa)*Ca(count);
    Ca(count+1)=Ca(count)+dCa;
    
end 

Casp=Ca;

%% BPAP post-pre
 Ibsf=0.75; Ibss=0.25;
 taubsf=3; taubfs=25;
 BPAP = 100* (Ibsf*exp(-tVec/taubsf) + Ibss* exp(-tVec/taubfs));
 BPAP = [zeros(1,spikepos-60) BPAP]; BPAP=BPAP(1:length(tVec));
postVol=-65;
postSpike=(zeros(1,length(tVec)));
for count=1:length(tVec)
    
    postVol(count)= BPAP(count)+EPSP(count);
   
    if postVol (count)>=-40 
        postSpike(count)=100;
    end
    if postVol(count)>30
        postVol(count)=-65;
    end
end
%% NMDA channel current and Ca2+ dynamics
Mg =1; Vr=130;
BV= 1./ (1+exp(-0.062*postVol)*Mg/3.57);
HVn10= BV.*(postVol-Vr);

P0=0.5; Gnmda = -1/500; 
tauf=50; taus=200; 
If=0.5; Is=0.5;
tauCa=50;   Ca=0;
for count=1:length(tVec)
    if tVec(count)<=spikepos
        Inmda(count)=0;
    else
        Inmda(count) = P0* Gnmda* (If*exp(-tVec(count)/tauf)+...
                        Is*exp(-tVec(count)/taus)).*HVn10(count);
    end
    dCa = Inmda(count)-(1/tauCa)*Ca(count);
    Ca(count+1)=Ca(count)+dCa;
    
end  
Can60=Ca;
%% pre-post spike stimulation

Ibsf=0.75; Ibss=0.25;
taubsf=3; taubfs=25;
BPAP = 100* (Ibsf*exp(-tVec/taubsf) + Ibss* exp(-tVec/taubfs));
BPAP = [zeros(1,spikepos+60) BPAP]; BPAP=BPAP(1:length(tVec));

% postsynaptic spikes
postVol=-65;
postSpike=(zeros(1,length(tVec)));
for count=1:length(tVec)
    
    postVol(count)= BPAP(count)+EPSP(count);
    if postVol (count)>=-40 
        postSpike(count)=100; % PostSpike is in the same scale as preSpike  
    end
    if postVol(count)>30
        postVol(count)=-65;
    end
end

% NMDA channel current and Ca2+ dynamics
Mg =1; Vr=130;
BV= 1./ (1+exp(-0.062*postVol)*Mg/3.57);
HV10= BV.*(postVol-Vr);

P0=0.5; Gnmda = -1/500; 
tauf=50; taus=200; 
If=0.5; Is=0.5;
tauCa=50;   Ca=0;
for count=1:length(tVec)
    if tVec(count)<=spikepos
        Inmda(count)=0;
    else
        Inmda(count) = P0* Gnmda* (If*exp(-tVec(count)/tauf)+...
                        Is*exp(-tVec(count)/taus)).*HV10(count);
    end
    dCa = Inmda(count)-(1/tauCa)*Ca(count);
    Ca(count+1)=Ca(count)+dCa;
    
end  
Ca60=Ca;

%% figures
figure(1)
[Plothandle,spikePos1] = plotRaster(preSpike,tVec,'red');
xlabel('\bf time, ms');
title('\bf Pre-spike');
set(gca,'fontsize',18);

figure(2)
semilogx(tVec,EPSP,'r','linewidth',2);
%semilogx(tVec,postVol,'g','linewidth',2);
legend('\bf EPSP');
title('\bf Single Presyn stimualtion');
xlabel('\bf Time ms');
ylabel('\bf postsynaptic voltage mV');
set(gca,'fontsize',14);

figure(3)
plot(Casp,'k','linewidth',2);
xlabel('\bf Time ms');
ylabel('\bf [Ca^{2+}] \muM');
set(gca,'fontsize',14);

%% possion spike generation function
%function [spikeMat,tVec]=possionSpikeGen(fr,tSim,nTrials)
%dt=1; %1ms
%nbins = floor(tSim/dt);
%spikeMat = rand(nTrials,nbins)<fr*dt;
%tVec=0:dt:tSim;
%end

%% PlotRaster
%function [Plothandle,spikePos] = plotRaster(spikeMat,tVec,c)
%hold all
%for Ntrail = 1:size(spikeMat,1)
 %   spikePos = tVec(spikeMat(Ntrail,:));
  %  for spikeCount =1:length(spikePos)
   %     Plothandle =  plot([spikePos(spikeCount) spikePos(spikeCount)],...
    %        [Ntrail-0.6 Ntrail+0.2] , 'color',c);
    %end
    %ylim([0 size(spikeMat,1)+0.5]);
%end
%end






        
      
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    