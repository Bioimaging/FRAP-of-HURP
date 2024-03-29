%   This script is extracting all the intensity values inside each *.mat
%   generated by the bleaching control experiments. These *.mat files have
%   to be in the same folder. 
% 
%   The intensity of the 1st channel is corrected for the bleach, and can
%   be swapped to have the most intense peak at the end.
%
%   A running average filter along the distance is applied for each time
%   point of each experiments. And the position of the maximum HURP signal
%   for each half-spindle is extracted.
%
%   For the Kymograph: the data are normalized according to the value at
%   the first time point.
%
%   For the FRAP efficiencies: the data are 100*(I-min(I))/(max(I)-min(I)) 
%   for each half-spindle.
%
%   For the bleaching efficiencies: the data are 100*I/I(t=0) for each
%   half-spindle.
%
%
%   Nicolas Liaudet
%   Bioimaging Core Facility - UNIGE
%   https://www.unige.ch/medecine/bioimaging/en/bioimaging-core-facility/
% 
%   v1.0 17-Jul-2018 NL


clc
clear
close all

path =  uigetdir(pwd,'Select your FRAP data');
fnames = dir([path filesep '*.mat']);
fnames = {fnames.name};

mch1 = cell(1,length(fnames));
mch2 = cell(1,length(fnames));
d    = cell(1,length(fnames));
t    = cell(1,length(fnames));
for idxF = 1:length(fnames)
    load([path filesep fnames{idxF}])
    mch1(idxF) = {cat(1,DATA.rel_mch1{:})};
    mch2(idxF) = {cat(1,DATA.rel_mch2{:})};            
    d(idxF) = {DATA.rel_d{1}};
    t(idxF) = {DATA.t};
end


dl = cellfun(@(x) length(x), d);
d = cat(1,d{:});
d = unique(d(:));
if length(d)~=dl
   errordlg('The distance data are inconsistent within all the files...',...
       'Distance problem')
end


MaxTframe = max(cellfun(@(x) size(x,1), mch1));
T = unique(cat(2,t{:}));



MCH1 = nan(MaxTframe,length(d),length(fnames));
MCH2 = nan(MaxTframe,length(d),length(fnames));

path1 =  uigetdir(path,'Select the folder where your BLEACHING MODEL is');
load([path1 filesep 'bleach_model.mat'])
bleac_fit = fitobject(T);
bleac_fit = bleac_fit-bleac_fit(1);
for idxF = 1:length(fnames)
    idxT = ismember(t{idxF},T);
    MCH1(idxT,:,idxF) = mch1{idxF}(idxT,:)-bleac_fit;
    MCH2(idxT,:,idxF) = mch2{idxF}(idxT,:);
end

mch1 = MCH1;
mch2 = MCH2;



choice = 'No';
choice = questdlg('Should I flip the distances to have Ch1 max at the end?', ...
	'Distance flipping', ...
	'Yes','No','No');

switch choice
    case 'No'
    case 'Yes'
        for idxF = 1:length(fnames)
            [~,idx] = max(mch1(:,:,idxF),[],2);
            idx = idx<length(d)/2;            
            mch1(idx,:,idxF)= mch1(idx,[length(d):-1:1],idxF);
            mch2(idx,:,idxF)= mch2(idx,[length(d):-1:1],idxF);
        end
end

%Running average along the distance, for each time point
%get the maximum on the left and the one the right (ch1)
Max_LR_idx = zeros(size(mch1,3),2);
for idxF = 1:size(mch1,3)
    for idxT = 1:size(mch1,1)
        tmp1 = mch1(idxT,:,idxF);
        tmp1(isnan(tmp1)) = 0;    
        tmp2 = mch2(idxT,:,idxF);
        tmp2(isnan(tmp2)) = 0;

        windowSize = 3; 
        b = (1/windowSize)*ones(1,windowSize);
        a = 1;
        tmp1 = filter(b,a,tmp1);
        tmp2 = filter(b,a,tmp2);
        
        
        if idxT == 1
            Sz = length(tmp1);
            SzH = round(Sz/2);        
            [~,Max_LR_idx(idxF,1)] = max(tmp1(1:SzH));
            [~,Max_LR_idx(idxF,2)] = max(tmp1(SzH+1:end));
            Max_LR_idx(idxF,2) = Max_LR_idx(idxF,2)+SzH;
        end

    end
end




%normalization according to maximum of the first row, i.e. the first time point
mch1n = mch1./repmat(max(mch1(1,:,:)),[size(mch2,1) size(mch2,2) 1]);

MCH1m = mean(mch1n,3,'omitnan');
 

Mont_mch1 = reshape(mch1,[size(mch1,1),size(mch1,2),1,size(mch1,3)]);
Mont_mch2 = reshape(mch2,[size(mch2,1),size(mch2,2),1,size(mch2,3)]);
Mont_mch1n = reshape(mch1n,[size(mch1n,1),size(mch1n,2),1,size(mch1n,3)]);

close all
%% figure 1 Kymograph
%-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
figure('PaperPosition',[0.5 0.5 49.5 20.5],'Color',[1 1 1],...
            'PaperSize',[50 21],'units','normalized','outerposition',[0 0 1 1],'Name','Kymograph')

%__________________________________________________________________________
ax(1) = subplot(3,2,1);

im(1) = montage(Mont_mch1,'Size',[1 size(mch1,3)],'DisplayRange',[min(Mont_mch1(:)) max(Mont_mch1(:))]);
for idxC = 1:size(mch1,3)-1
    line( idxC*size(mch1,2)*ones(1,2)+0.5, [1-0.5 size(mch1,1)+0.5],'Color',[0 1 0],...
        'LineWidth',2,'Parent',ax(1))        
end
for idxC = 1:size(mch1,3)   
    line( (idxC-1)*size(mch1,2)*ones(1,2)+Max_LR_idx(idxC,:), ones(2,1),'Color',[0 1 0],...
        'LineWidth',2,'Parent',ax(1),'LineStyle','none','Marker','o','MarkerEdgeColor',[0 0 1])     
end

title('Ch1')
ylabel('t(s)')
xlabel('Cell')
set(ax(1),'Visible','on','DataAspectRatioMode','auto',...
    'XTick',[0:1:size(mch1,3)-1]*size(mch1,2)+size(mch1,2)/2,'XTickLabel',num2str([1:size(mch1,3)]'),...
    'YTick',[1:5:length(T)],'YTickLabel',num2str(T([1:5:length(T)])'))
colormap(ax(1),hot)
set(im(1),'CDataMapping','Scaled')
colorbar
 
%__________________________________________________________________________ 
ax(2) = subplot(3,2,3);

im(2) = montage(Mont_mch2,'Size',[1 size(mch2,3)],'DisplayRange',[min(Mont_mch2(:)) max(Mont_mch2(:))]);
for idxC = 1:size(mch2,3)-1
    line( idxC*size(mch2,2)*ones(1,2)+0.5, [1-0.5 size(mch2,1)+0.5],'Color',[0 1 0],...
        'LineWidth',2,'Parent',ax(2))
end

title('Ch2')
ylabel('t(s)')
xlabel('Cell')
set(ax(2),'Visible','on','DataAspectRatioMode','auto',...
    'XTick',[0:1:size(mch2,3)-1]*size(mch2,2)+size(mch2,2)/2,'XTickLabel',num2str([1:size(mch2,3)]'),...
    'YTick',[1:5:length(T)],'YTickLabel',num2str(T([1:5:length(T)])'))
colormap(ax(2),hot)
set(im(2),'CDataMapping','Scaled')
colorbar

%__________________________________________________________________________
ax(3) = subplot(3,2,5);

im(3) = montage(Mont_mch1n,'Size',[1 size(mch1n,3)],'DisplayRange',[min(Mont_mch1n(:)) max(Mont_mch1n(:))]);
for idxC = 1:size(mch1n,3)-1
    line( idxC*size(mch1n,2)*ones(1,2)+0.5, [1-0.5 size(mch1n,1)+0.5],'Color',[0 1 0],...
        'LineWidth',2,'Parent',ax(3))
end

title('Ch1 t_1')
ylabel('t(s)')
xlabel('Cell')
set(ax(3),'Visible','on','DataAspectRatioMode','auto',...
    'XTick',[0:1:size(mch1n,3)-1]*size(mch1n,2)+size(mch1n,2)/2,'XTickLabel',num2str([1:size(mch1n,3)]'),...
    'YTick',[1:5:length(T)],'YTickLabel',num2str(T([1:5:length(T)])'),...
    'CLim',[min(Mont_mch1n(:)) max(Mont_mch1n(:))])
colormap(ax(3),hot)
set(im(3),'CDataMapping','Scaled')
colorbar

%__________________________________________________________________________
ax(4) = subplot(3,2,[2 4 6]);%THIS ONE

im(4) = imshow(MCH1m,hot(512));

title('Ch1 t_1')
ylabel('t(s)')
xlabel('d (%)')
set(ax(4),'Visible','on','DataAspectRatioMode','auto',...
    'XTick',[1:5:length(d)],'XTickLabel',num2str(d([1:5:length(d)]')),...
    'YTick',[1:5:length(T)],'YTickLabel',num2str(T([1:5:length(T)])'))
set(im(4),'CDataMapping','Scaled')
set(ax(4),'CLim',[0.1 0.9]);
hc = colorbar(ax(4));
hc.Limits = [0 1];

saveas(gcf,[path filesep 'FRAP.fig'])

 
%% Computation of the bleaching and FRAP efficieny per half-spindles
 
 Bleach_eff = zeros(size(mch1,1),2,size(mch1,3));
 FRAP_eff   = zeros(size(mch1,1),2,size(mch1,3));
 for idxC = 1:size(mch1,3)
     tmp = mch1(:,:,idxC);
     Bleach_eff(:,1,idxC) = 100*tmp(:,Max_LR_idx(idxC,1))/tmp(1,Max_LR_idx(idxC,1));
     Bleach_eff(:,2,idxC) = 100*tmp(:,Max_LR_idx(idxC,2))/tmp(1,Max_LR_idx(idxC,2));
     FRAP_eff(:,1,idxC) = 100*(tmp(:,Max_LR_idx(idxC,1))-(min(tmp(:,Max_LR_idx(idxC,1)))))/ ( (max(tmp(:,Max_LR_idx(idxC,1)))) -(min(tmp(:,Max_LR_idx(idxC,1))))) ;
     FRAP_eff(:,2,idxC) = 100*(tmp(:,Max_LR_idx(idxC,2))-(min(tmp(:,Max_LR_idx(idxC,2)))))/ ( (max(tmp(:,Max_LR_idx(idxC,2)))) -(min(tmp(:,Max_LR_idx(idxC,2))))) ;
 end
 
 
 
cmap = jet(size(mch1,3));
  
  
%% figure 2  all traces concatenated
%-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_

figure('Name','all traces concatenated')

hold on
for idxT=1:size(mch1,1)
    tmp  = mch1(idxT,:,:);
    plot([1:size(mch1,2)*size(mch1,3)] ,...
        reshape(tmp,[size(mch1,2)*size(mch1,3),1]))
end
hold off
set(gca,'Color',[0.9 0.9 0.9],'Box','on')

%% figure 3  FRAP efficiencies
close all
%-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
cmap = jet(size(mch1,3));
figure('Name','FRAP efficiencies')


med_l = median(FRAP_eff(:,1,:),3)';
p25_l   = prctile(FRAP_eff(:,1,:),25,3)';
p75_l   = prctile(FRAP_eff(:,1,:),75,3)';

med_r = median(FRAP_eff(:,2,:),3)';
p25_r   = prctile(FRAP_eff(:,2,:),25,3)';
p75_r   = prctile(FRAP_eff(:,2,:),75,3)';

if med_r(3)-med_r(end) <0 % increase model i.e. FRAP
    subplot(2,3,1)
    
    for idxC=1:size(mch1,3)
        line(T,FRAP_eff(:,1,idxC),'LineStyle','-','Color',cmap(idxC,:))
    end
    
    title('FRAP left half-spindle')
    xlabel('t(s)')
    ylabel('\epsilon_{FRAP}(%)')
    set(gca,'Color',[0.9 0.9 0.9],'Box','on')
    
    %__________________________________________________________________________
    subplot(2,3,2)
    
    for idxC=1:size(mch1,3)
        line(T,FRAP_eff(:,2,idxC),'LineStyle','--','Color',cmap(idxC,:))
    end
    
    title('FRAP right half-spindle')
    xlabel('t(s)')
    ylabel('\epsilon_{FRAP}(%)')
    set(gca,'Color',[0.9 0.9 0.9],'Box','on')
    
    %__________________________________________________________________________
    subplot(2,3,3)
    
    for idxC=1:size(mch1,3)
        line(T,FRAP_eff(:,1,idxC),'LineStyle','-','Color',cmap(idxC,:))
        line(T,FRAP_eff(:,2,idxC),'LineStyle','--','Color',cmap(idxC,:))
    end
    
    title('FRAP')
    xlabel('t(s)')
    ylabel('\epsilon_{FRAP}(%)')
    set(gca,'Color',[0.9 0.9 0.9],'Box','on')
    
    %__________________________________________________________________________
    ax = subplot(2,3,4);%THIS ONE
    
    frapmodel = fittype('A*(1-exp(-lambda*t))','independent','t');
    RecoveryEffect = fit( T(2:end)'-T(2), med_l(2:end)', frapmodel,'StartPoint',[40 0.08]);
    
    patch([T, fliplr(T)]',[p75_l, fliplr(p25_l)]',1,...
        'FaceColor',[0 0 1 ],'EdgeColor','none',...
        'FaceAlpha',0.5,'Parent',ax);
    line(T,med_l,'LineStyle','-','Color',[0 0 1],'LineWidth',2);
    
    tt = [T(2):0.1:T(end)]-T(2);
    hl = line(tt+T(2),RecoveryEffect.A*(1-exp(-RecoveryEffect.lambda*tt)),...
        'Color',[0 0 0],'LineWidth',1);
    
    title(['FRAP median left half-spindle, 5\lambda^{-1}=' num2str(5/RecoveryEffect.lambda,'%03.2f') 's, t_{0.5}=' num2str(log(2)/RecoveryEffect.lambda,'%03.2f') ])
    xlabel('t(s)')
    ylabel('\epsilon_{FRAP}(%)')
    set(gca,'Color',[1 1 1],'Box','on')
    hl =legend(hl,['$\textrm{' num2str(RecoveryEffect.A,'%03.2f') '}\cdot',...
        '\left(\textrm{1}-\textrm{e}^{\textrm{-' num2str(RecoveryEffect.lambda,'%03.2f') 't}}\right)$']);
    hl.Interpreter='latex';
    hl.FontSize = 12;
    
    
    
    %
    %__________________________________________________________________________
    ax = subplot(2,3,5);%THIS ONE
    
    
    
    frapmodel = fittype('A*(1-exp(-lambda*t))','independent','t');
    RecoveryEffect = fit( T(2:end)'-T(2), med_r(2:end)', frapmodel,'StartPoint',[40 0.08]);
    tt = [T(2):0.1:T(end)]-T(2);
    
    patch([T, fliplr(T)]',[p75_r, fliplr(p25_r)]',1,...
        'FaceColor',[1 0 0],'EdgeColor','none',...
        'FaceAlpha',0.5,'Parent',ax);
    line(T,med_r,'LineStyle','-','Color',[1 0 0],'LineWidth',2)
    
    
    
    hl  = line(tt+T(2),RecoveryEffect.A*(1-exp(-RecoveryEffect.lambda*tt)),...
        'Color',[0 0 0],'LineWidth',1);
    
    title(['FRAP median right half-spindle, 5\lambda^{-1}=' num2str(5/RecoveryEffect.lambda,'%03.2f') 's, t_{0.5}=' num2str(log(2)/RecoveryEffect.lambda,'%03.2f') ])
    xlabel('t(s)')
    ylabel('\epsilon_{FRAP}(%)')
    set(gca,'Color',[1 1 1],'Box','on')
    
    hl =legend(hl,['$\textrm{' num2str(RecoveryEffect.A,'%03.2f') '}\cdot',...
        '\left(\textrm{1}-\textrm{e}^{\textrm{-' num2str(RecoveryEffect.lambda,'%03.2f') 't}}\right)$']);
    
    hl.Interpreter='latex';
    hl.FontSize = 12;
    
    
    %__________________________________________________________________________
    ax = subplot(2,3,6);%THIS ONE
    
    line(T,med_l,'LineStyle','-','Color',[0 0 1],'LineWidth',2)
    line(T,med_r,'LineStyle','-','Color',[1 0 0],'LineWidth',2)
    patch([T, fliplr(T)]',[p75_l, fliplr(p25_l)]',1,...
        'FaceColor',[0 0 1],'EdgeColor','none',...
        'FaceAlpha',0.5,'Parent',ax);
    
    
    
    title('FRAP medians')
    xlabel('t(s)')
    ylabel('\epsilon_{FRAP}(%)')
    set(gca,'Color',[1 1 1],'Box','on')
    
    
    
else
     subplot(2,1,1)
    
    for idxC=1:size(mch1,3)
        line(T,FRAP_eff(:,1,idxC),'LineStyle','-','Color',cmap(idxC,:))
    end
    
    title('FRAP left half-spindle')
    xlabel('t(s)')
    ylabel('\epsilon_{FRAP}(%)')
    set(gca,'Color',[0.9 0.9 0.9],'Box','on')
    
   
    
    %__________________________________________________________________________
    ax = subplot(2,1,2);%THIS ONE
    
    frapmodel = fittype('A*(1-exp(-lambda*t))','independent','t');
    RecoveryEffect = fit( T(2:end)'-T(2), med_l(2:end)', frapmodel,'StartPoint',[40 0.08]);
    
    patch([T, fliplr(T)]',[p75_l, fliplr(p25_l)]',1,...
        'FaceColor',[0 0 1 ],'EdgeColor','none',...
        'FaceAlpha',0.5,'Parent',ax);
    line(T,med_l,'LineStyle','-','Color',[0 0 1],'LineWidth',2);
    
    tt = [T(2):0.1:T(end)]-T(2);
    hl = line(tt+T(2),RecoveryEffect.A*(1-exp(-RecoveryEffect.lambda*tt)),...
        'Color',[0 0 0],'LineWidth',1);
    
    title(['FRAP median left half-spindle, 5\lambda^{-1}=' num2str(5/RecoveryEffect.lambda,'%03.2f') 's, t_{0.5}=' num2str(log(2)/RecoveryEffect.lambda,'%03.2f') ])
    xlabel('t(s)')
    ylabel('\epsilon_{FRAP}(%)')
    set(gca,'Color',[1 1 1],'Box','on')
    hl =legend(hl,['$\textrm{' num2str(RecoveryEffect.A,'%03.2f') '}\cdot',...
        '\left(\textrm{1}-\textrm{e}^{\textrm{-' num2str(RecoveryEffect.lambda,'%03.2f') 't}}\right)$']);
    hl.Interpreter='latex';
    hl.FontSize = 12;
    
    
    
   
    
    
    
end


saveas(gcf,[path filesep 'FRAP graphs.fig'])
 

%% figure 4  Bleaching efficiencies
%-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_

figure('Name','Bleaching efficiencies')

%__________________________________________________________________________
subplot(2,3,1)

cmap = lines(size(mch1,3));
for idxC=1:size(mch1,3)
    line(T,Bleach_eff(:,1,idxC),'LineStyle','-','Color',cmap(idxC,:))
end

title('Bleaching left half-spindle')
xlabel('t(s)')
ylabel('\epsilon_{Bleach}(%)')
set(gca,'Color',[0.9 0.9 0.9],'Box','on')

%__________________________________________________________________________
subplot(2,3,2)

cmap = lines(size(mch1,3));
for idxC=1:size(mch1,3)    
    line(T,Bleach_eff(:,2,idxC),'LineStyle','--','Color',cmap(idxC,:))
end
title('Bleaching right half-spindle')
xlabel('t(s)')
ylabel('\epsilon_{Bleach}(%)')
set(gca,'Color',[0.9 0.9 0.9],'Box','on')

%__________________________________________________________________________ 
subplot(2,3,3)

cmap = lines(size(mch1,3));
for idxC=1:size(mch1,3)    
    line(T,Bleach_eff(:,1,idxC),'LineStyle','-','Color',cmap(idxC,:))
    line(T,Bleach_eff(:,2,idxC),'LineStyle','--','Color',cmap(idxC,:))
end

title('Bleaching')
xlabel('t(s)')
ylabel('\epsilon_{Bleach}(%)')
set(gca,'Color',[0.9 0.9 0.9],'Box','on')

%__________________________________________________________________________ 
ax = subplot(2,3,4);%THIS ONE

med_l = median(Bleach_eff(:,1,:),3)';
p25_l   = prctile(Bleach_eff(:,1,:),25,3)';
p75_l   = prctile(Bleach_eff(:,1,:),75,3)';

patch([T, fliplr(T)]',[p75_l, fliplr(p25_l)]',1,...
    'FaceColor',[0 0 1 ],'EdgeColor','none',...
    'FaceAlpha',0.5,'Parent',ax);
line(T,med_l,'LineStyle','-','Color',[0 0 1],'LineWidth',2)

title('Bleaching median left half-spindle')
xlabel('t(s)')
ylabel('\epsilon_{Bleach}(%)')
set(gca,'Color',[1 1 1],'Box','on')

%__________________________________________________________________________ 
ax = subplot(2,3,5);%THIS ONE

med_r = median(Bleach_eff(:,2,:),3)';
p25_r   = prctile(Bleach_eff(:,2,:),25,3)';
p75_r   = prctile(Bleach_eff(:,2,:),75,3)';

patch([T, fliplr(T)]',[p75_r, fliplr(p25_r)]',1,...
    'FaceColor',[1 0 0 ],'EdgeColor','none',...
    'FaceAlpha',0.5,'Parent',ax);

line(T,med_r,'LineStyle','-','Color',[1 0 0],'LineWidth',2)
title('Bleaching median right half-spindle')
xlabel('t(s)')
ylabel('\epsilon_{Bleach}(%)')
set(gca,'Color',[1 1 1],'Box','on')

%__________________________________________________________________________ 
ax = subplot(2,3,6);%THIS ONE

line(T,med_l,'LineStyle','-','Color',[0 0 1],'LineWidth',2)
line(T,med_r,'LineStyle','-','Color',[1 0 0],'LineWidth',2)
patch([T, fliplr(T)]',[p75_l, fliplr(p25_l)]',1,...
    'FaceColor',[0 0 1],'EdgeColor','none',...
    'FaceAlpha',0.5,'Parent',ax);
patch([T, fliplr(T)]',[p75_r, fliplr(p25_r)]',1,...
    'FaceColor',[1 0 0],'EdgeColor','none',...
    'FaceAlpha',0.5,'Parent',ax);

xlabel('t(s)')
ylabel('\epsilon_{Bleach}(%)')
set(gca,'Color',[1 1 1],'Box','on')
title('Bleaching median')

saveas(gcf,[path filesep 'Bleaching graphs.fig'])

%% figure 3
%-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
% close all
figure('Name','Initial-Final values')
 
%__________________________________________________________________________
ax(1) =subplot(2,2,1);

cmap = lines(2);
h = gobjects(2,size(mch1,3));
for idxC =1:size(mch1,3)
    tmp  = mch1(1,:,idxC);
    h(1,idxC) = line(d,tmp,'Color',cmap(1,:));
    tmp  = mch1(end,:,idxC);
    h(2,idxC) = line(d,tmp,'Color',cmap(2,:));
end

legend(h(:,end),'First frame','Last frame')
xlabel('distance (%)')
ylabel('mch1 intensity(-)')
title('Individual line profiles')

%__________________________________________________________________________
ax(2) = subplot(2,2,3);%THIS ONE

ave_mch1_T0 = mean(squeeze(mch1(1,:,:)),2);
std_mch1_T0 = std(squeeze(mch1(1,:,:)),[],2);
ave_mch1_TF = mean(squeeze(mch1(end,:,:)),2);
std_mch1_TF = std(squeeze(mch1(end,:,:)),[],2);

hl(1) = patch([d; flipud(d)],[ave_mch1_T0+std_mch1_T0; flipud(ave_mch1_T0-std_mch1_T0)],1,...
    'FaceColor',cmap(1,:),'EdgeColor','none',...
    'FaceAlpha',0.3,'Parent',ax(2));
hl(2) = patch([d; flipud(d)],[ave_mch1_TF+std_mch1_TF; flipud(ave_mch1_TF-std_mch1_TF)],1,...
    'FaceColor',cmap(2,:),'EdgeColor','none',...
    'FaceAlpha',0.3,'Parent',ax(2));
line(d,ave_mch1_T0,'Color',cmap(1,:))
line(d,ave_mch1_TF,'Color',cmap(2,:))


xlabel('distance (%)')
ylabel('mch1 intensity(-)')
title('Mean \pm standard deviation line profile')
legend(hl,'First frame','Last frame')

%-------------------------------
StartPL = zeros(size(mch1,3),1);
StartPR = zeros(size(mch1,3),1);
StopPL = zeros(size(mch1,3),1);
StopPR = zeros(size(mch1,3),1);
for idxC =1:size(mch1,3)
    StartPL(idxC) = mch1(1,Max_LR_idx(idxC,1),idxC);
    StartPR(idxC) = mch1(1,Max_LR_idx(idxC,2),idxC);
    StopPL(idxC) = mch1(end,Max_LR_idx(idxC,1),idxC);
    StopPR(idxC) = mch1(end,Max_LR_idx(idxC,2),idxC);
end

%__________________________________________________________________________
ax(3) = subplot(2,2,2);%THIS ONE

boxplot([StartPL;StartPR],[repmat({'Left half-spindle'},[length(StartPL) 1]);repmat({'Right half-spindle'},[length(StartPR) 1])] )
xL = 0.2*(rand(size(mch1,3),1)-0.5)+1;
xR = 0.2*(rand(size(mch1,3),1)-0.5)+2;
line(xL, StartPL,...
    'LineStyle','none','Marker','o','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor',[0.2 0.2 0.2])
line(xR, StartPR,...
    'LineStyle','none','Marker','o','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor',[0.2 0.2 0.2])
for idxP = 1:length(StartPL)
    line([xL(idxP) xR(idxP)],[StartPL(idxP) StartPR(idxP)],'LineStyle','--','Color',[0 0 0])        
end

[h,p] = ttest2(StartPL,StartPR);
line([1 2], ones(2,1)*max([StartPL;StartPR])*1.1)
% text(1.5,max([StartPL;StartPR])*1.15, ['p=' num2str(p)],'HorizontalAlignment','center')
ylabel('mch1 intensity(-)')
title(['First frame, paired t-test: p=' num2str(p)])
%__________________________________________________________________________
ax(4) = subplot(2,2,4);%THIS ONE

boxplot([StopPL;StopPR],[repmat({'Left half-spindle'},[length(StopPL) 1]);repmat({'Right half-spindle'},[length(StopPR) 1])] )
xL = 0.2*(rand(size(mch1,3),1)-0.5)+1;
xR = 0.2*(rand(size(mch1,3),1)-0.5)+2;
line(xL, StopPL,...
    'LineStyle','none','Marker','o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor',[0.2 0.2 0.2])
line(xR, StopPR,...
    'LineStyle','none','Marker','o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor',[0.2 0.2 0.2])
for idxP = 1:length(StopPL)
    line([xL(idxP) xR(idxP)],[StopPL(idxP) StopPR(idxP)],'LineStyle','--','Color',[0 0 0])        
end
[h,p] = ttest2(StopPL,StopPR);
line([1 2], ones(2,1)*max([StopPL;StopPR])*1.1)
% text(1.5,max([StopPL;StopPR])*1.15, ['p=' num2str(p)],'HorizontalAlignment','center')
ylabel('mch1 intensity(-)')
title(['Last frame, paired t-test: p=' num2str(p)])

saveas(gcf,[path filesep 'FRAP ratios.fig'])



