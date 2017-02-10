    % CHECKING AND VISUALIZATION
% -----------------------------------------
%    DATA PIPELINE
% -----------------------------------------
%
% (1) RAW_PROFILES     (2)  OVERSAMPLED (cubic/linear interp)
%    Data(i).I   -->          curveOverSampled             -->
%    Data(i).I2
%
%  (3) GAUSS_FIT (2-4 gaussians) (4)  Gauss_FIT_SAMPLED 
%        curvefitA3   -->                  sI2
%






% Fill the Intensity for dots that passed the checks (in Data)
clear all;
close all;

FV = 10;    % Field of view used to validate spots detected
TF = 36;     % Number of time points to check
th = 0:pi/50:2*pi;
r = 5;       % radius of marker in global image
PrePostFrames = 10;   % Number of frames recorded before/after A3 event was first and last detected
WithWeights = 1;    % Introduces a weight as a function of the amplitud of the curve


% --------------------------
% Asking input from User:
fprintf('Checking level on the events: (each level includes previous level check)\n');
fprintf('1 = Events with a minimum intensity level\n');
fprintf('2 = Events with good correlation A3-Dyn\n');
fprintf('3 = Events with gaussian profile shape for A3 channel\n');
nThresholds=input('Indicate Event check level: ');

fprintf('Do you want to generate the subfigures where event was detected?\n');
PlotDots=input('Indicate (Yes=1, No=0): ');

fprintf('Do you want to store the .fig Matlab plots of each profile?\n');
SaveFig=input('Indicate (Yes=1, No=0): ');


fprintf('Do you want to fit to RAW data (1 Gaussian) or Model adjusted data (several Gaussians)?\n');
RawData=input('Indicate (RAW (SC) =1, MODEL=0): ');



% Get filename of file produced by postuTrack
[spots_name,spots_folder] = uigetfile('*.mat','Select the A3 .mat Events file');

% If dots are to be ploted we need to open the .tiff stack
if (PlotDots)
    [image_name,image_folder] = uigetfile({'*.tif';'*.tiff'},'Select the A3 .tif stack file');
    
    % Extract File properties from TIFF STACK
    File_props = imfinfo([image_folder image_name]);

    % Load original Image stack
    fprintf('   Removing Background from TIFF File [%s]\n',[image_folder image_name]);
    se = strel('disk',12);
    for ( i=1:length(File_props) )
        S(i).ima = imtophat(imread([image_folder image_name],'Index',i),se);
        Io(i).ima = imread([image_folder image_name],'Index',i);
    end;
end;

% File with Statistics for each event
stats_filename = [spots_folder spots_name(1:end-4) '_stats.txt'];
stats_filenamexls = [spots_folder spots_name(1:end-4) '_stats.xls'];

stats_f = fopen(stats_filename,'w');

% Load data produced by postuTrack
load([spots_folder spots_name]);
isDyn = isfield(Data,'DI2');
if (isDyn)
   fprintf(' > Two channels detected\n');
   fprintf(stats_f,'Data for 2 channels:\n');
else
  fprintf(' > One channels detected\n');
  fprintf(stats_f,'Data for 1 channel:\n');
end

BlackI = zeros(FV*2+1,FV*2+1);
dLen = length(Data);
aux = [Data(:).Thres];
nSel = sum(aux==nThresholds);
selected = find(aux == nThresholds);

% # POINT, MaxA3, MaxDyn, timeA3, timeDyn, WidthA3, WidthDyn,   
AllStats = zeros(nSel,33);     % Matrix to store all stats in Excel

fprintf('Total number of events DETECTED is: %g\n',dLen);
fprintf('Total number of events SELECTED is: %g\n',nSel);

fprintf(stats_f,'Total number of events SELECTED is: %g\n',nSel);

tgap = zeros(1,nSel);
GSelA3 = zeros(1,nSel);
GSelDyn = zeros(1,nSel);
%options = fitoptions('gauss2'); % Set Robust fitting




fprintf('Producing files for Manual Validation\n')
    for s = 1:nSel  
        i = selected(s);  % Select the events from the list index
        len = length(Data(i).fr);
        x_i = round(Data(i).x(1));
        y_i = round(Data(i).y(1)); 
        fprintf('Event [%.0d/%.0d]\n',s,nSel);
        ini = Data(i).fr(1); % First frame in which it was detected
        fin = Data(i).fr(end);
        
        % Data to file
        fprintf(stats_f,'POINT %g: ',i); 
        AllStats(s,1) = i; % Event number

        tgap(i) = 0; % Initialization
        Mcor = 0;
        options.Weights = [];

        % Check for repeated time points (it happens some times)
        Rep = len-length(unique(Data(i).fr'));
        
        % Pre-fit to smooth the input data
        t = Data(i).fr(1+Rep):0.25:Data(i).fr(end); % All the range including gaps OVERSAMPLED
        
        
        if (RawData == 0)
            Data(i).I2 = Data(i).I2 - min(Data(i).I2); % Removing the base intensity
            curve = fit(Data(i).fr(1+Rep:end)',Data(i).I2(1+Rep:end)','cubicinterp');
            curveOverSampled = curve(t);  
        else
            curve = fit(Data(i).fr(1+Rep:end)',Data(i).I2(1+Rep:end)','linearinterp');
            curveOverSampled = curve(t);
        end
        
        if (isDyn)
            
            if (RawData == 0)
                Data(i).DI2 = Data(i).DI2 - min(Data(i).DI2); % Removing the base intensity
                curveD = fit(Data(i).fr(1+Rep:end)',Data(i).DI2(1+Rep:end)','cubicinterp');
                curveOverDSampled = curveD(t);
            else
                curveD = fit(Data(i).fr(1+Rep:end)',Data(i).DI2(1+Rep:end)','linearinterp');
                curveOverDSampled = curveD(t);
            end
        end
        % -------------------------
        % Checking THE BEST FIT based on GoF (comparing 2,3,4 gaussians)
        % Gaussian GOF
        if (WithWeights)
            options.Weights = Data(i).I2/max(Data(i).I2);
        end
        
        
        % Fitting 2 GAUSSIANS to Curve oversampled with cubic/linear interp
        [curvefitA3_2,gofA3,outputA3] = fit(t',curveOverSampled,'gauss2');
        gofA3 = gofA3.rsquare;
        sI2 = curvefitA3_2(t);
        if (abs(max(sI2)/max(Data(i).I2))<2)
            % Then there is NOT overfit -> KEEP
            GSelA3(s) = 2; % Choose gauss2
            curvefitA3 = curvefitA3_2;
        end
        
        if (RawData == 0)
            if (len > 9)  % At least 9 points to fit 3 Gaussians
                % Fitting 3 GAUSSIANS to Curve oversampled with cubic/linear interp
                [curvefitA3_3,gofA3Aux,outputA3] = fit(t',curveOverSampled,'gauss3');
                aux = curvefitA3_3(t);
                if (abs(max(aux)/max(Data(i).I2))<2) && (gofA3Aux.rsquare>gofA3)
                    % Then there is NOT overfit -> KEEP
                    GSelA3(s) = 3; % Choose gauss3
                    curvefitA3 = curvefitA3_3;
                    gofA3 = gofA3Aux.rsquare;
                    sI2 = aux;
                end
            end
            if (len > 12)  % At least 12 points to fit 3 Gaussians
                % Fitting 4 GAUSSIANS to Curve oversampled with cubic/linear interp
                [curvefitA3_4,gofA3Aux,outputA3] = fit(t',curveOverSampled,'gauss4');
                aux = curvefitA3_4(t);
                if (abs(max(aux)/max(Data(i).I2))<2) && (gofA3Aux.rsquare>gofA3)
                    % Then there is NOT overfit -> KEEP
                    GSelA3(s) = 4; % Choose gauss4
                    curvefitA3 = curvefitA3_4;
                    gofA3 = gofA3Aux.rsquare;
                    sI2 = aux;
                end
            end    
        end
        % -------------------------
        % Getting some A3 statistics 
        % sI2 = Best possible Gaussian fit SAMPLED
        MaxPeakA3 = max(findpeaks(sI2));
        MinPeakA3 = min(sI2);
        if isempty(MaxPeakA3)
            [MaxPeakA3, PeakA3]= max(sI2);
            fprintf(stats_f,'MaxA3[%1.2f] ',MaxPeakA3);
            AllStats(s,2) = MaxPeakA3; % Maximum value of A3
        else
            [PeakA3, ~]= find(sI2==MaxPeakA3);
            fprintf(stats_f,'MaxA3[%1.2f] ',MaxPeakA3);
            AllStats(s,2) = MaxPeakA3;
        end;
        
        
        AllStats(s,8) = GSelA3(s);
        AllStats(s,4) = t(PeakA3);
        fprintf(stats_f,'tA3[%1.2f] ', t(PeakA3));
        % If there is DYN Channel, fit the curves
        
        
        if (isDyn)        
            % -------------------------
            % Checking THE BEST FIT based on GoF
            % Gaussian GOF
            if (WithWeights)
                options.Weights = Data(i).DI2/max(Data(i).DI2);
            end

            [curvefitDyn_2,gofDyn,outputDyn] = fit(t',curveOverDSampled,'gauss2');
            gofDyn = gofDyn.rsquare;
            sDI2 = curvefitDyn_2(t);
            if (abs(max(sDI2)/max(Data(i).I2))<2)
                % Then there is NOT overfit -> KEEP
                GSelDyn(s) = 2; % Choose gauss2
                curvefitDyn = curvefitDyn_2;
            end

            if (RawData == 0)
                if (len > 9)  % At least 9 points to fit 3 Gaussians
                    [curvefitDyn_3,gofDynAux,outputDyn] = fit(t',curveOverDSampled,'gauss3');
                    aux = curvefitDyn_3(t);
                    if (abs(max(aux)/max(Data(i).I2))<2) && (gofDynAux.rsquare>gofDyn)
                        % Then there is an overfit -> Discard
                        GSelDyn(s) = 3; % Choose gauss3
                        curvefitDyn = curvefitDyn_3;
                        gofDyn = gofDynAux.rsquare;
                        sDI2 = aux;
                    end
                end
                if (len > 12)  % At least 12 points to fit 3 Gaussians
                    [curvefitDyn_4,gofDynAux,outputDyn] = fit(t',curveOverDSampled,'gauss4');
                    aux = curvefitDyn_4(t);
                    if (abs(max(aux)/max(Data(i).I2))<2) && (gofDynAux.rsquare>gofDyn)
                        % Then there is an overfit -> Discard
                        GSelDyn(s) = 4; % Choose gauss4
                        curvefitDyn = curvefitDyn_4;
                        gofDyn = gofDynAux.rsquare;
                        sDI2 = aux;
                    end
                end 
            end
           % -------------------------
           % Getting some Dyn statistics                                
            MaxPeakDyn = max(findpeaks(sDI2));
            MinPeakDyn = min(sDI2);
            
            if isempty(MaxPeakDyn)
                [MaxPeakDyn, PeakDyn]= max(sDI2);
                AllStats(s,3) = MaxPeakDyn;
                fprintf(stats_f,'MaxDyn[%1.2f] ',MaxPeakDyn);
            else
                [PeakDyn, ~]= find(sDI2==MaxPeakDyn);
                AllStats(s,3) = MaxPeakDyn;
                fprintf(stats_f,'MaxDyn[%1.2f] ',MaxPeakDyn);
            end;
            tgap(i) =  t(PeakA3)-t(PeakDyn);
            fprintf(stats_f,'tDyn[%1.2f] ',t(PeakDyn));            
            AllStats(s,5) = t(PeakDyn);
        end;
        
        
        
        % Find the Gaussian that corresponds to the Max peak
        % Measure the distance from the Max Peak to each Guassian
        
        
        PeakDist = len*ones(1,4); % Set to maximum distance in the sequence
        WA3 = zeros(1,4);
        AA3 = zeros(1,4);
        tA3 = zeros(1,4);
        SelGaussA3 = zeros(4,length(t));
        
        PeakDist(1) = abs(t(PeakA3) - curvefitA3.b1);
        WA3(1) = curvefitA3.c1;
        AA3(1) = curvefitA3.a1;
        tA3(1) = curvefitA3.b1;
        SelGaussA3(1,:) = curvefitA3.a1 * exp(-((t-curvefitA3.b1)/curvefitA3.c1).^2);
            
        PeakDist(2) = abs(t(PeakA3) - curvefitA3.b2);
        WA3(2) = curvefitA3.c2;
        AA3(2) = curvefitA3.a2;
        tA3(2) = curvefitA3.b2;
        SelGaussA3(2,:) = curvefitA3.a2 * exp(-((t-curvefitA3.b2)/curvefitA3.c2).^2);
        
        if (GSelA3(s) == 3) && (curvefitA3.a3 > 0)
            PeakDist(3) = abs(t(PeakA3) - curvefitA3.b3);
            WA3(3) = curvefitA3.c3;
            AA3(3) = curvefitA3.a3;
            tA3(3) = curvefitA3.b3;
            SelGaussA3(3,:) = curvefitA3.a3 * exp(-((t-curvefitA3.b3)/curvefitA3.c3).^2);
        end
        
        if (GSelA3(s) == 4) && (curvefitA3.a4 > 0)
            PeakDist(4) = abs(t(PeakA3) - curvefitA3.b4);
            WA3(4) = curvefitA3.c4;
            AA3(4) = curvefitA3.a4;
            tA3(4) = curvefitA3.b4;
            SelGaussA3(4,:) = curvefitA3.a4 * exp(-((t-curvefitA3.b4)/curvefitA3.c4).^2);
        end
        
        [~, IndPeakA3] = min(PeakDist);
        [~, IndSortedA3] = sort(AA3,'descend');
        
        AllStats(s,10) = WA3(IndSortedA3(1));
        AllStats(s,11) = AA3(IndSortedA3(1));
        AllStats(s,12) = tA3(IndSortedA3(1));
        
        AllStats(s,13) = WA3(IndSortedA3(2));
        AllStats(s,14) = AA3(IndSortedA3(2));
        AllStats(s,15) = tA3(IndSortedA3(2));
        
        %if (GSelA3(s) == 3)
            AllStats(s,16) = WA3(IndSortedA3(3));
            AllStats(s,17) = AA3(IndSortedA3(3));
            AllStats(s,18) = tA3(IndSortedA3(3));
        %end;
        
        %if (GSelA3(s) == 4)
            AllStats(s,19) = WA3(IndSortedA3(4));
            AllStats(s,20) = AA3(IndSortedA3(4));
            AllStats(s,21) = tA3(IndSortedA3(4));
       % end;
       fprintf(stats_f,'WA3[%1.2f] ',WA3(IndPeakA3));
       AllStats(s,6) = WA3(IndPeakA3);
       AllStats(s,8) = GSelA3(s);
       fprintf(stats_f,'GA3[%g] ',GSelA3(s));

       if (isDyn)
           WDyn = zeros(1,4);
           ADyn = zeros(1,4);
           tDyn = zeros(1,4);
           SelGaussDyn = zeros(4,length(t));
        
            PeakDist = len*ones(1,4); % Set to maximum distance in the sequence
            WDyn = zeros(1,4);
            
            PeakDist(1) = abs(t(PeakDyn) - curvefitDyn.b1);
            WDyn(1) = curvefitDyn.c1;
            ADyn(1) = curvefitDyn.a1;
            tDyn(1) = curvefitDyn.b1;
            %SelGaussDyn(1,:) = curvefitDyn.a1 * exp(-((t-curvefitDyn.b1)/curvefitDyn.c1).^2);
            
            PeakDist(2) = abs(t(PeakDyn) - curvefitDyn.b2);
            WDyn(2) = curvefitDyn.c2;
            ADyn(2) = curvefitDyn.a2;
            tDyn(2) = curvefitDyn.b2;
            %SelGaussDyn(2,:) = curvefitDyn.a2 * exp(-((t-curvefitDyn.b2)/curvefitDyn.c2).^2);
            
            if (GSelDyn(s) == 3) && (curvefitDyn.a3 > 0)
                PeakDist(3) = abs(t(PeakDyn) - curvefitDyn.b3);
                WDyn(3) = curvefitDyn.c3;
                ADyn(3) = curvefitDyn.a3;
                tDyn(3) = curvefitDyn.b3;
                %SelGaussDyn(3,:) = curvefitDyn.a3 * exp(-((t-curvefitDyn.b3)/curvefitDyn.c2).^3);
            end
            if (GSelDyn(s) == 4) && (curvefitDyn.a4 > 0)
                PeakDist(4) = abs(t(PeakDyn) - curvefitDyn.b4);
                WDyn(4) = curvefitDyn.c4;
                ADyn(4) = curvefitDyn.a4;
                tDyn(4) = curvefitDyn.b4;
                %SelGaussDyn(4,:) = curvefitDyn.a4 * exp(-((t-curvefitDyn.b4)/curvefitDyn.c4).^2);
            end
            [~, IndPeak] = min(PeakDist);
            [~, IndSortedDyn] = sort(ADyn,'descend');
            
            
           fprintf(stats_f,'WDyn[%1.2f] ',WDyn(IndPeak));
           AllStats(s,7) = WDyn(IndPeak);

            AllStats(s,22) = WDyn(IndSortedDyn(1));
            AllStats(s,23) = ADyn(IndSortedDyn(1));
            AllStats(s,24) = tDyn(IndSortedDyn(1));
        
            AllStats(s,25) = WDyn(IndSortedDyn(2));
            AllStats(s,26) = ADyn(IndSortedDyn(2));
            AllStats(s,27) = tDyn(IndSortedDyn(2));
        
            AllStats(s,28) = WDyn(IndSortedDyn(3));
            AllStats(s,29) = ADyn(IndSortedDyn(3));
            AllStats(s,30) = tDyn(IndSortedDyn(3));
            
            AllStats(s,31) = WDyn(IndSortedDyn(4));
            AllStats(s,32) = ADyn(IndSortedDyn(4));
            AllStats(s,33) = tDyn(IndSortedDyn(4));
           
           AllStats(s,9) = GSelDyn(s);
           fprintf(stats_f,'GDyn[%g] ',GSelDyn(s));
           fprintf(stats_f,'\n');
       end
        % Show the intensity profile (using Gaussian fitted)
        h1=figure(1); 
        
        %plot(Data(i).fr, sI2,'r'); % Intensity of A3 Channel
        h=figure(1);
       
        plot(Data(i).fr,Data(i).I2,'r--o');
        hold on;
        
        plot(curvefitA3,'r'); % Intensity of A3 Channel
        %plot(curve,'k'); % Intensity of A3 Channel
        
        % Commented: SHOW ALL GAUSSIANS FITTED and STORE them
        plot(t,SelGaussA3(IndSortedA3(2),:),'b--'); % Intensity of A3 Channel
        if (GSelA3(s) == 3)
            plot(t,SelGaussA3(IndSortedA3(3),:),'y--'); % Intensity of A3 Channel
        end
        if (GSelA3(s) == 4)
            plot(t,SelGaussA3(IndSortedA3(4),:),'c--'); % Intensity of A3 Channel
        end
        % Plot again the one with Highest Peak in Black color
        plot(t,SelGaussA3(IndPeakA3,:),'k--'); % Intensity of A3 Channel
         
        y1=get(gca,'ylim');
        if (isDyn)
            plot(curvefitDyn,'g'); % Intensity of Dyn Channel
            plot(Data(i).fr,Data(i).DI2,'g--o');
            y1=get(gca,'ylim');
            line([t(PeakDyn) t(PeakDyn)],y1,'Color',[0 1 0]);
            
            % Calculating Correlation between curves
            aux = corrcoef(sI2,sDI2);
            Mcor = aux(1,2); 
        end
        
        
        line([t(PeakA3) t(PeakA3)],y1,'Color',[1 0 0]);
        title(sprintf('Temporal gap %g',tgap(i)));
        xlabel('Frames');
        

        hold off;
        
        %Calculating lag-delay using cross-correlation maximum
        %[acor,lag]  = xcorr(sI2,sDI2);
        %[cor,I] = max(abs(acor));
        %tgap(i) = lag(I);
        

        % Saving plots in files
        if ((GSelA3(s)>0) && (GSelDyn(s)>0) || ~isDyn)
            if (tgap(i) > 0) 
                saveas(h1,[spots_folder sprintf('DynA3(%1.2f)(%1.2f)(%1.2f)(G%g%g)_prof_%g.jpg',abs(tgap(i)),Mcor,gofA3,GSelA3(s),GSelDyn(s),i)],'jpg');
                if SaveFig
                    saveas(h1,[spots_folder sprintf('DynA3(%1.2f)(%1.2f)(%1.2f)(G%g%g)_prof_%g.fig',abs(tgap(i)),Mcor,gofA3,GSelA3(s),GSelDyn(s),i)],'fig');
                end
            else
                saveas(h1,[spots_folder sprintf('A3Dyn(%1.2f)(%1.2f)(%1.2f)(G%g%g)_prof_%g.jpg',abs(tgap(i)),Mcor,gofA3,GSelA3(s),GSelDyn(s),i)],'jpg');
                if SaveFig
                    saveas(h1,[spots_folder sprintf('A3Dyn(%1.2f)(%1.2f)(%1.2f)(G%g%g)_prof_%g.fig',abs(tgap(i)),Mcor,gofA3,GSelA3(s),GSelDyn(s),i)],'fig');
                end
            end;
        end;
        
        
        % Code to sva the spots detected in each frame of the stack
        if (PlotDots)
            h2=figure(2);
            % Checking max and min image dimensions for the subfigures Field of View
            up_c = y_i-FV;
            down_c = y_i+FV;
            left_c = x_i-FV;
            right_c = x_i+FV;

            if (up_c <= 0) 
                up_c = 1;
            end
            if (down_c > File_props(1).Height)
                down_c = File_props(1).Height;
            end;
            if (left_c <=0)
                left_c = 1;
            end;
            if right_c > File_props(1).Width
                right_c = File_props(1).Width;
            end;

        % Plotting subfigures around the spot

            k = 1;
            if (ini>2) % plot 2 frames before detection
                subplot(6,6,1);
                imshow(Io(ini-2).ima(up_c:down_c,left_c:right_c));
                title(['fr=' num2str(ini-2) ],'Color','b');

                subplot(6,6,2);            
                imshow(Io(ini-1).ima(up_c:down_c,left_c:right_c));
                title(['fr=' num2str(ini-1) ],'Color','b');
                k = 3;
            end;

            for (j=1:TF-k+1)  
              if ((k <= TF) && (j<=len))  
                subplot(6,6,k);
                imshow(Io(Data(i).fr(j)).ima(up_c:down_c,left_c:right_c));
                title(['fr=' num2str(Data(i).fr(j))],'Color','r');
                k = k + 1;
              end;
              if (j>len)
                subplot(6,6,k);
                imshow(BlackI);
                title('No spot');
                k = k + 1;
              end;
            end;
            saveas(h2,[spots_folder sprintf('event%.0d',i)],'jpg');
        end;
    end; % end for
    
    
    fclose(stats_f);
    xlswrite(stats_filenamexls,AllStats);
    
    