function [xt_recov_IMFsLowToHigh,xt_recov_IMFsHighToLow] = FBDM(xt,Fs,alfa)
    threshold=0;
    N=length(xt);
    Xk=FB1(xt,alfa);     
    nTotalHormonics=N;
        xt_recov_IMFsLowToHigh=zeros(1,N)';
        IFOfSignalIF=zeros(1,N)'; 
        xt_AnalyticFIBF=zeros(1,N)';
        init=1;
        p=1;     
        while(init<=N)          
            [kk,xt_recov_FIBF,xt_recov_AnalyticFIBF,IFOfSignal]=getIMFsScanAllLowToHigh(Xk,Fs,init,nTotalHormonics,threshold,alfa);
            init=kk+1;
            xt_recov_IMFsLowToHigh(:,p)=xt_recov_FIBF;
            IFOfSignalIF(:,p)=IFOfSignal;
            xt_AnalyticFIBF(:,p)=xt_recov_AnalyticFIBF;            
            p=p+1;  
        end         
   
        %sp_PlotTF(xt_AnalyticFIBF,IFOfSignalIF,t,Fs,Fs/2);
        
        %title('FBDM: TFR of FBIBFs (LTH-FS)','FontSize',32,'FontName','Times');
        xt_recov_IMFsHighToLow=zeros(1,N)'; 
        IFOfSignalIF=zeros(1,N)'; 
        xt_AnalyticFIBF=zeros(1,N)';
        mm=0;        
        final=N;        
        while(final>=1)
            mm=mm+1; 
            [kk,xt_recov_FIBF,xt_recov_AnalyticFIBF,IFOfSignal]=getIMFsScanAllHighToLow(Xk,Fs,final,threshold,alfa);
            final=kk-1;
            xt_recov_IMFsHighToLow(:,mm)=xt_recov_FIBF;
            xt_AnalyticFIBF(:,mm)=xt_recov_AnalyticFIBF;
            IFOfSignalIF(:,mm)=IFOfSignal;            
        end
    
        %sp_PlotTF(xt_AnalyticFIBF,IFOfSignalIF,t,Fs,Fs/2);   
        %title('FBDM: TFR of FBIBFs (HTL-FS)','FontSize',32,'FontName','Times');        
%% Calculate Energy Leakge
%     FMF_energy=0;
%     for np=1:length(xt_recov_IMFsHighToLow(1,:))         
%          FMF_energy=FMF_energy+sum(xt_recov_IMFsHighToLow(:,np).*xt_recov_IMFsHighToLow(:,np));         
%     end   
%     SignalEnergy=sum(xt(1:end-1).*xt(1:end-1)); 
%     EnergyLekage=SignalEnergy-FMF_energy;          
  %%%
  function [kk,xt_recov_FIBF,xt_recov_AnalyticFIBF,IFOfSignal]=getIMFsScanAllLowToHigh(Xk,Fs,init,final,threshold,alfa)
    xt_recov_FIBF=0; 
    xt_recov_AnalyticFIBF=0;    
    N=length(Xk);
    EndPoints=3;  
    Analytic_zt=0;
    n=1:length(Xk);
    col=1;
    PositiveIMF=1;
    for kk=init:final
         Analytic_zt=Analytic_zt+hilbert((Xk(kk)).*besselj(0,alfa(kk)/N*n)); 
         IMFs_phase=unwrap(angle(Analytic_zt));
         tmp1=(IMFs_phase);
         tmp=((tmp1(3:end)-tmp1(1:end-2))/2)*(Fs/(2*pi)); 
         tmp=[tmp(1) tmp tmp(end)];
                
         if((min(tmp(EndPoints:end-(EndPoints-1)))<threshold) && (max(xt_recov_FIBF)~=0) && (PositiveIMF==1) && (kk<=final)) % check -ve IF
             PositiveIMF=0;             
             recordIMF(:,col)=xt_recov_FIBF;
             recordAnalyticFIBFs(:,col)=xt_recov_AnalyticFIBF;
             recordIMF_IF(:,col)=xt_IF; % IF
             recordIMFkk(col)=kk-1;
             col=col+1;             
         end
         
         if(min(tmp(EndPoints:end-(EndPoints-1)))>=threshold) % check +ve IF, First value may be negative so 2:end
             PositiveIMF=1;             
         end         
         xt_recov_FIBF=real(Analytic_zt); % *2 as taking only half of component 
         xt_recov_AnalyticFIBF=Analytic_zt;
         xt_IF=tmp; % instantaneous Freq
         
         % First Value
         if(kk==init)
             recordIMF(:,col)=xt_recov_FIBF;
             recordAnalyticFIBFs(:,col)=xt_recov_AnalyticFIBF;
             recordIMF_IF(:,col)=xt_IF; % IF
             recordIMFkk(col)=kk;
             col=col+1;      
         end
       
         % for last value
         if((kk==final) && (PositiveIMF==1))             
             recordIMF(:,col)=xt_recov_FIBF;
             recordAnalyticFIBFs(:,col)=xt_recov_AnalyticFIBF;
             recordIMF_IF(:,col)=xt_IF; % IF
             recordIMFkk(col)=kk;             
         end         
    end 
    getHighestValue=length(recordIMFkk);
    xt_recov_FIBF=recordIMF(:,getHighestValue); % return last value IMF
    xt_recov_AnalyticFIBF=recordAnalyticFIBFs(:,getHighestValue);
    IFOfSignal=recordIMF_IF(:,getHighestValue);
    kk=recordIMFkk(end); % return last value
    
    
    function [kk,xt_recov_FIBF,xt_recov_AnalyticFIBF,IFOfSignal]=getIMFsScanAllHighToLow(Xk,Fs,final,threshold,alfa)
    xt_recov_FIBF=0; 
    xt_recov_AnalyticFIBF=0;    
    N=length(Xk);
    EndPoints=3;   
    Analytic_zt=0;
    n=1:length(Xk);
    col=1;
    PositiveIMF=1;
    for kk=final:-1:1
        Analytic_zt=Analytic_zt+hilbert((Xk(kk)).*besselj(0,alfa(kk)/N*n));        
        IMFs_phase=unwrap(angle(Analytic_zt));
        tmp1=(IMFs_phase);
        tmp=((tmp1(3:end)-tmp1(1:end-2))/2)*(Fs/(2*pi));
        tmp=[tmp(1) tmp tmp(end)];
                
         if( (min(tmp(EndPoints:end-(EndPoints-1)))<threshold) && (max(xt_recov_FIBF)~=0) && (PositiveIMF==1) && (kk>=1) ) % check -ve IF
             PositiveIMF=0;             
             recordIMF(:,col)=xt_recov_FIBF;
             recordAnalyticFIBFs(:,col)=xt_recov_AnalyticFIBF;
             recordIMF_IF(:,col)=xt_IF; % IF
             recordIMFkk(col)=kk+1;
             col=col+1;             
         end

         if(min(tmp(EndPoints:end-(EndPoints-1)))>=threshold) % check +ve IF, First value & last value may be negative so 2:end
             PositiveIMF=1;             
         end         
         xt_recov_FIBF=real(Analytic_zt); % *2 as taking only half of component
         xt_recov_AnalyticFIBF=Analytic_zt;
         xt_IF=tmp; % IF
         % First Value
         if(kk==final)
             recordIMF(:,col)=xt_recov_FIBF;
             recordAnalyticFIBFs(:,col)=xt_recov_AnalyticFIBF;
             recordIMF_IF(:,col)=xt_IF;
             recordIMFkk(col)=kk;
             col=col+1;      
         end
         
         %% for last value
         if((kk==2) && (PositiveIMF==1))             
             recordIMF(:,col)=xt_recov_FIBF;
             recordAnalyticFIBFs(:,col)=xt_recov_AnalyticFIBF;
             recordIMF_IF(:,col)=xt_IF; % IF
             recordIMFkk(col)=kk;             
         end         
    end 
    getHighestValue=length(recordIMFkk);
    xt_recov_FIBF=recordIMF(:,getHighestValue); % return last value IMF
    xt_recov_AnalyticFIBF=recordAnalyticFIBFs(:,getHighestValue);
    IFOfSignal=recordIMF_IF(:,getHighestValue);
    kk=recordIMFkk(end);
    
%     function sp_PlotTF(xt_recov_FIBFs,xt_recov_IMFs_IF,t,Fs,fw1)        
%     xt_recov_FIBFs_abs=abs(xt_recov_FIBFs);
%     %% To plot Use code of RCADA    
%     %[nt,tscale,fscale]=PlotTF_FFT(freq,amp,t0,t1,fres,tres,fw0,fw1,tw0,tw1);
%     [nt,tscale,fscale]=PlotTF_FFT(xt_recov_IMFs_IF(:,1:end),xt_recov_FIBFs_abs(:,1:end),t,Fs,fw1); % magnitude value 
%     %[nt,tscale,fscale]=PlotTF_FFT(xt_recov_IMFs_IF(:,1:end),xt_recov_IMFs_AnalyticAbs(:,1:end),t,Fs,fw1); % magnitude value
%     q=fspecial('gaussian',7,0.6);  
%     nsu=filter2(q,nt);
%     nsu=filter2(q,nsu);
%     if 1
%     norme = nsu/max(abs(nsu(:)));   
%     for k=1:length(norme(:,1))
%     E1(k)=sum(norme(k,1:length(norme(k,:))).^3);
%     E2(k)=sum(norme(k,1:length(norme(k,:))));
%     end
%     RE=(1/(1-3))*log2(sum(E1)/sum(E2))
%     end
%     figure;
%     %subplot(5,1,1);
%     imagesc(tscale,fscale,nsu.^.5); 
%     axis xy;
%     xlabel('Time (second)')
%     ylabel('Frequency (Hz)')
%     colorbar;
%     %title('FMD based Time Frequncy Plot');
%     %title('FBDM: Time-Frequency plot of FBIBFs','FontSize',32,'FontName','Times'); 
%     %% Use end code of RCADA
%     FntSize=32;
%     set(gca,'FontSize',FntSize,'FontName','Times')
%     h1 = get(gca, 'xlabel');
%     set(h1,'FontSize',FntSize,'FontName','Times')
%     h1 = get(gca, 'ylabel');
%     set(h1,'FontSize',FntSize,'FontName','Times') 
%     
%     function [nt,tscale,fscale]=PlotTF_FFT(freq,amp,t,Fs,fw1)
%     
%     t0=0;t1=t(end);    
%     multfactor=4;
%     if(length(freq(:,1))>=100*multfactor) 
%         fres=100*multfactor; tres=100*multfactor;
%     else
%         fres=length(freq(:,1));    tres=fres;
%     end
%     
%     %fres=268;    tres=fres; % for earthquake data
%     
%     %fw0=0; fw1=max(max(freq));
%     fw0=min(min(freq)); fw1=max(max(Fs/2));
%     if(fw0<0)          fw0=0;    end
%     %fw1=Fs/2; % max frequency in plot
%     tw0=t0;     tw1=t1;
%     %  4.call nspplote.m to plot the figure of time-frequency spectrum 
%     %----- Get the values to plot the spectrum
%     lscale=0;
%     [nt,tscale,fscale] = nspplote(freq,amp,t0,t1,fres,tres,fw0,fw1,tw0,tw1,lscale); 
%    