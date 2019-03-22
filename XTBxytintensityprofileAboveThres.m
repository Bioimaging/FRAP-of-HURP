% Extract intensity profile between two spots in 2D+t
%
%    <CustomTools>
%      <Menu name="Bioimaging XT">
%       <Submenu name="2D+t">
%        <Item name="I-profile above threshold between 2 spots" icon="Matlab" tooltip="Intensity profile  above threshold between 2 spots in 2D+t">
%          <Command>MatlabXT::XTBxytintensityprofileAboveThres(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%    </CustomTools>
%
% Copyright (c) Mar 2017, Bioimaging Core Facility
% v1. Nicolas Liaudet




function [ output_args ] = XTBxytintensityprofileAboveThres(aImarisApplicationID)
%XTNSEG Summary of this function goes here
%   Extract intensity profile between two spots in 2D+t.
%    <CustomTools>
%      <Menu>
%        <Item name="I-profile above threshold between 2 spots" icon="Matlab" tooltip="Intensity profile above threshold between 2 spots in 2D+t">
%          <Command>MatlabXT::XTBxytintensityprofileAboveThres(%i)</Command>
%        </Item>
%      </Menu>
%    </CustomTools>



%% connect to Imaris interface
disp('Imaris connection...')
tic
if ~isa(aImarisApplicationID, 'Imaris.IApplicationPrxHelper')
    javaaddpath ImarisLib.jar
    vImarisLib = ImarisLib;
    if ischar(aImarisApplicationID)
        aImarisApplicationID = round(str2double(aImarisApplicationID));
    end
    vImarisApplication = vImarisLib.GetApplication(aImarisApplicationID);
else
    vImarisApplication = aImarisApplicationID;
end
vDataSet = vImarisApplication.GetDataSet;
 %vImarisDataSet = vImarisApplication.GetDataSet.Clone;
% vImarisDataSet = vImarisApplication.GetDataSet;
toc



%% Get DATA geometry
DataSizeXYZTC = [vDataSet.GetSizeX, vDataSet.GetSizeY,...
    vDataSet.GetSizeZ, vDataSet.GetSizeT,...
    vDataSet.GetSizeC];

DataPosXYZ = [vDataSet.GetExtendMinX, vDataSet.GetExtendMaxX;...
              vDataSet.GetExtendMinY, vDataSet.GetExtendMaxY;...
              vDataSet.GetExtendMinZ, vDataSet.GetExtendMaxZ];
          
DataResXYZ = diff(DataPosXYZ,1,2)./DataSizeXYZTC(1:3)';

%% get the spots into tracks

vSpots = vImarisApplication.GetFactory.ToSpots(vImarisApplication.GetSurpassSelection);
if isempty(vSpots)
    errordlg('Please select a "Spots" object!','Selection Error');        
    return
end

XYZ        = vSpots.GetPositionsXYZ;
Tidx       = vSpots.GetIndicesT;
TrackIds   = vSpots.GetTrackIds;
TrackEdges = vSpots.GetTrackEdges;



if rem(length(Tidx),2) ~= 0
     errordlg('Some spots are not always there!','Dataset Error');        
    return
end

tracks = unique(TrackIds);
if length(tracks) ~= 2
     errordlg('There is not 2 continuous tracks!','Dataset Error');        
    return
end

%Sort and connect the spots according to their tracks
SpotAidx = ismember(TrackIds,tracks(1));
SpotAidx = TrackEdges(SpotAidx,:);
SpotAidx = unique(SpotAidx)+1;

SpotBidx = ismember(TrackIds,tracks(2));
SpotBidx = TrackEdges(SpotBidx,:);
SpotBidx = unique(SpotBidx)+1;

%Reassign the A or B labelling according to the highligthed spot in the current scene 
SelectedSpots = vSpots.GetSelectedIndices;

% any(SpotAidx == SelectedSpots)

if ~isempty(SelectedSpots)
    if (any(SpotAidx == SelectedSpots))
        
        tmp1 = SpotAidx;
        tmp2 = SpotBidx;
        
        SpotBidx = tmp1;
        SpotAidx = tmp2;
    end
end



TidxA = Tidx(SpotAidx);
TidxB = Tidx(SpotBidx);

[~,pTA] = sort(TidxA);
[~,pTB] = sort(TidxB);

SPOTA_XYZ_um  = XYZ(SpotAidx,:);
SPOTB_XYZ_um  = XYZ(SpotBidx,:);

SPOTA_XYZ_um = SPOTA_XYZ_um(pTA,:);
SPOTB_XYZ_um = SPOTB_XYZ_um(pTB,:);

SPOTA_XYZ = SPOTA_XYZ_um-repmat(DataPosXYZ(:,1)',size(SPOTA_XYZ_um,1),1);
SPOTB_XYZ = SPOTB_XYZ_um-repmat(DataPosXYZ(:,1)',size(SPOTB_XYZ_um,1),1);
SPOTA_XYZ = SPOTA_XYZ./repmat(DataResXYZ',size(SPOTA_XYZ_um,1),1); 
SPOTB_XYZ = SPOTB_XYZ./repmat(DataResXYZ',size(SPOTB_XYZ_um,1),1);
SPOTA_XYZ = round(SPOTA_XYZ);
SPOTB_XYZ = round(SPOTB_XYZ);
% 
% SPOTB_XYZ_um  = XYZ(SpotBidx,:);

t = cell(1,length(TidxA));
for idxt = 1:length(TidxA)
    t(idxt) = vSpots.GetTimePoint(TidxA(idxt));
end
t = t(pTA);



%%

% absolute
abs_d = cell(1,length(TidxA));%distances in absolute units
abs_mch1 = cell(1,length(TidxA));%average intensity at a given distance in Ch1
abs_sch1 = cell(1,length(TidxA));%std intensity at a given distance in Ch1

abs_mch2 = cell(1,length(TidxA));%average intensity at a given distance in Ch2
abs_sch2 = cell(1,length(TidxA));%std intensity at a given distance in Ch2

% relative
rel_d = cell(1,length(TidxA));%distances in relative units 100% is centrosome-centrosome
rel_mch1 = cell(1,length(TidxA));%average intensity at a given distance in Ch1
rel_sch1 = cell(1,length(TidxA));%std intensity at a given distance in Ch1

rel_mch2 = cell(1,length(TidxA));%average intensity at a given distance in Ch2
rel_sch2 = cell(1,length(TidxA));%std intensity at a given distance in Ch2


prompt = {'Enter the radius (um):'};
dlg_title  = 'Integration cylinder';
num_lines  = 1;
defaultans = {'1'};
rcyl = inputdlg(prompt,dlg_title,num_lines,defaultans);
rcyl = str2double(rcyl{:});




DataSizeXYZTC(5) = DataSizeXYZTC(5)+2;
      
vDataSet.SetSizeC(DataSizeXYZTC(5));
  
        
        
vDataSet.SetChannelName(DataSizeXYZTC(5)-2,['Link'])
vDataSet.SetChannelRange(DataSizeXYZTC(5)-2,0,1)

vDataSet.SetChannelName(DataSizeXYZTC(5)-1,['Geodesic distance'])
vDataSet.SetChannelRange(DataSizeXYZTC(5)-1,0,max((sum((SPOTB_XYZ-SPOTA_XYZ).^2,2)).^0.5))




hwbar = waitbar(0, ['Processing...']);


for idxT=1:length(t)
    hwbar = waitbar(idxT/length(t), hwbar,...
        ['Processing z-stack ',...
        num2str(idxT) '/' num2str(length(t))]);
    
    if strcmp(vDataSet.GetType,'eTypeUInt8')
        tmp = zeros([DataSizeXYZTC(1:3)],'uint8');
    elseif strcmp(vDataSet.GetType,'eTypeUInt16')
        tmp = zeros([DataSizeXYZTC(1:3)],'uint16');
    elseif strcmp(vDataSet.GetType,'eTypeFloat')
        tmp = zeros([DataSizeXYZTC(1:3)],'single');
    end
    
    a_xyz = SPOTA_XYZ(idxT,:);
    b_xyz = SPOTB_XYZ(idxT,:);
   
    %make the cylinder axis    
    [link_x,link_y,link_z] = bresenham_line3d(a_xyz,b_xyz, 0);               
    link_ind = sub2ind(DataSizeXYZTC([1 2 3]),link_x,link_y,link_z);
    tmp(link_ind) = 1;        
    vDataSet.SetDataVolumeAs1DArrayShorts(tmp(:),DataSizeXYZTC(5)-2 , idxT-1 )
    
    
    %keep voxels closer than xxx um away from the axis
    D = bwdist(logical(tmp));    
    D = uint16(D<= rcyl/DataResXYZ(1)) ;            
    idx_D = find(D~=0);    
    [Dx,Dy,Dz] = ind2sub(DataSizeXYZTC([2 1 3]),idx_D);
    Dxyz = [Dx,Dy,Dz];
    
    
    %kill voxel behind A    
    n_AB = b_xyz-a_xyz;
    D_a_xyz = Dxyz-repmat(a_xyz,[size(Dxyz,1) 1 ]);    
    Pscal_abDa = dot(repmat(n_AB,[size(D_a_xyz,1) 1 ]),D_a_xyz,2);
    idx_rmv_A = Pscal_abDa<0;    
    D(idx_D(idx_rmv_A)) = 0;


    %kill voxel behind B
    n_BA = -n_AB;
    D_b_xyz = Dxyz-repmat(b_xyz,[size(Dxyz,1) 1 ]);    
    Pscal_baDb = dot(repmat(n_BA,[size(D_b_xyz,1) 1 ]),D_b_xyz,2);
    idx_rmv_B = Pscal_baDb<0;
    D(idx_D(idx_rmv_B)) = 0;
    
    
   
    
    %find "A face"
    Bc = strel(ones(3,3,3));
    gD = D-imerode(D,Bc);
    
    tmp = false(size(D));
    dm = floor(rcyl/DataResXYZ(1));
    
    tmp(a_xyz(1),a_xyz(2),a_xyz(3)) = true;
    
    sw = (2*dm-1)/2;
    ses2 = ceil(2*dm/2);
    [y,x,z] = meshgrid(-sw:sw, -sw:sw, -sw:sw);
    m = sqrt(x.^2 + y.^2 + z.^2);
    bc = (m <= m(ses2, ses2, 2*dm));
    Bc = strel('arbitrary', bc);    
    
    tmp = imdilate(tmp,Bc);
        
    msk = gD~=0;
    msk = msk&tmp;
    
    
    %distance along the cylinder
    D = bwdistgeodesic(D~=0,msk,'quasi-euclidean');
    
    
    vDataSet.SetDataVolumeAs1DArrayShorts(uint16(D(:)),DataSizeXYZTC(5)-1 , idxT-1 )
    
    %group voxels according to the distance with ds um steps
    D = D*DataResXYZ(1);%in um...
    ds = 0.15;
    abs_edges = [0:ds:floor(max(D(:)))];
    
    %group voxels according to percentage of the max ditance
%     D = D*DataResXYZ(1);%in um...
    dperc = 0.02;
        
    rel_edges = linspace(0,max(D(:)),round(1/dperc+1));
    

    %Get channels 1 and 2 value for the current time point    
    if strcmp(vDataSet.GetType,'eTypeUInt8')
        Ch1 = zeros([DataSizeXYZTC(1:3)],'uint8');
        Ch1(:) = typecast(vDataSet.GetDataVolumeAs1DArrayShorts(...
            0,idxT-1), 'uint8');
        Ch2 = zeros([DataSizeXYZTC(1:3)],'uint8');
        Ch2(:) = typecast(vDataSet.GetDataVolumeAs1DArrayShorts(...
            1,idxT-1), 'uint8');
    elseif strcmp(vDataSet.GetType,'eTypeUInt16')
        Ch1 = zeros([DataSizeXYZTC(1:3)],'uint16');
        Ch1(:) = typecast(vDataSet.GetDataVolumeAs1DArrayShorts(...
            0,idxT-1), 'uint16');
        Ch2 = zeros([DataSizeXYZTC(1:3)],'uint16');
        Ch2(:) = typecast(vDataSet.GetDataVolumeAs1DArrayShorts(...
            1,idxT-1), 'uint16');
    elseif strcmp(vDataSet.GetType,'eTypeFloat')
        Ch1 = zeros([DataSizeXYZTC(1:3)],'single');
        Ch1(:) = typecast(vDataSet.GetDataVolumeAs1DArrayShorts(...
            0,idxT-1), 'single');
        Ch2 = zeros([DataSizeXYZTC(1:3)],'single');
        Ch2(:) = typecast(vDataSet.GetDataVolumeAs1DArrayShorts(...
            1,idxT-1), 'single');
    end
    
    
    %get average and std intensities at a given distance
    m1 = zeros(1,length(abs_edges)-1);
    s1 = zeros(1,length(abs_edges)-1);
    m2 = zeros(1,length(abs_edges)-1);
    s2 = zeros(1,length(abs_edges)-1);    
    
    for didx = 1:length(abs_edges)-1   
        
        Binidx = D(:)>=abs_edges(didx)&D(:)<abs_edges(didx+1);
        valch1 = Ch1(Binidx);
        valch2 = Ch2(Binidx);
        
        valch1(valch1<1000) = [];
        
        
        m1(didx) = mean(valch1);
        s1(didx) = std(single(valch1));
        m2(didx) = mean(valch2);
        s2(didx) = std(single(valch2));
        
    end
    
    abs_d(1,idxT)={abs_edges(1:end-1)+ds/2};
    abs_mch1(1,idxT) = {m1};
    abs_sch1(1,idxT) = {s1};
    abs_mch2(1,idxT) = {m2};
    abs_sch2(1,idxT) = {s2};

    %get average and std intensities at a given percentage distance
    m1 = zeros(1,length(rel_edges)-1);
    s1 = zeros(1,length(rel_edges)-1);
    m2 = zeros(1,length(rel_edges)-1);
    s2 = zeros(1,length(rel_edges)-1);    
    
    for didx = 1:length(rel_edges)-1           
        Binidx = D(:)>=rel_edges(didx)&D(:)<rel_edges(didx+1);
        valch1 = Ch1(Binidx);
        valch2 = Ch2(Binidx);
                
        valch1(valch1<1000) = [];
        
        m1(didx) = mean(valch1);
        s1(didx) = std(single(valch1));
        m2(didx) = mean(valch2);
        s2(didx) = std(single(valch2));
        
    end
    
    rel_d(1,idxT)={100*[dperc/2:dperc:1-dperc/2]};
    rel_mch1(1,idxT) = {m1};
    rel_sch1(1,idxT) = {s1};
    rel_mch2(1,idxT) = {m2};
    rel_sch2(1,idxT) = {s2};

    
%     unique() 
   
end




close(hwbar)



dt = diff(datetime(t,'Format','yyyy-MM-dd HH:mm:ss.SSSSSSSSS'));
dt.Format = 'm';
dt = round(mean(seconds(dt)))
nt = [0:1:length(t)-1]*dt;




FilePathName = char(vImarisApplication.GetCurrentFileName);

FilePathName = FilePathName(1:end-3);


%absolute version
FileParam = ['d absolute, r=' num2str(rcyl) 'um'];
isAbs = true;
makefigure(SPOTA_XYZ_um,SPOTB_XYZ_um,abs_d,abs_mch1,abs_mch2,nt,FilePathName,FileParam,isAbs)


%relative version
FileParam = ['d relative, r=' num2str(rcyl) 'um'];
isAbs = false;
makefigure(SPOTA_XYZ_um,SPOTB_XYZ_um,rel_d,rel_mch1,rel_mch2,nt,FilePathName,FileParam,isAbs)


DATA.FileName = [FilePathName '.ims'];
DATA.rcyl  = rcyl;

DATA.ds    = ds;
DATA.dperc = dperc;


DATA.SPOTA_XYZ_um = SPOTA_XYZ_um;
DATA.SPOTB_XYZ_um = SPOTB_XYZ_um;
DATA.t = nt;

DATA.rel_d    = rel_d;
DATA.rel_mch1 = rel_mch1;
DATA.rel_mch2 = rel_mch2;

DATA.abs_d    = abs_d;
DATA.abs_mch1 = abs_mch1;
DATA.abs_mch2 = abs_mch2;



save([FilePathName 'mat'],'DATA')



    function makefigure(SPOTA_XYZ_um,SPOTB_XYZ_um,d,mch1,mch2,nt,FilePathName,FileParam,isAbs)
        
        
        
        
        % non normalized distance
        figure('PaperPosition',[0.635 0.635 28.431 19.731],'Color',[1 1 1],...
            'PaperSize',[29.7 21],'units','normalized','outerposition',[0 0 1 1],'Name',FileParam)
        %
        %
        % 'PaperOrientation','Portrait','PaperType','A4','PaperPosition',[0.3 0.63 29.1 19.73],...
        %     'units','normalized','outerposition',[0 0 1 1],'Name','d absolute')
        ax= zeros(1,9);
        for idx = 1:9
            ax(idx) = subplot(3,3,idx);
        end
        
        delete(ax(7))
        
        %--------------------------------------------------------------------------
        hp(1) = line(SPOTA_XYZ_um(:,1),SPOTA_XYZ_um(:,2),SPOTA_XYZ_um(:,3),'Parent',ax(1),...
            'Color',[1 0 0],'Marker','o','MarkerEdgeColor',[0.2 0.2 0.2],'MarkerFaceColor',[1 0 0]);
        hp(2) = line(SPOTB_XYZ_um(:,1),SPOTB_XYZ_um(:,2),SPOTB_XYZ_um(:,3),'Parent',ax(1),...
            'Color',[0 1 0],'Marker','o','MarkerEdgeColor',[0.2 0.2 0.2],'MarkerFaceColor',[0 1 0]);
        xlabel(ax(1),'x(\mum)')
        ylabel(ax(1),'y(\mum)')
        zlabel(ax(1),'z(\mum)')
        title(ax(1),'Centrosomes motion')
        legend(ax(1),hp,{'Centrosome A','Centrosome B'},'Location','Best')
        set(ax(1),'View',[-37.5 30],'Box','on','BoxStyle','full',...
            'XGrid','on','YGrid','on','ZGrid','on')
        
        
        %--------------------------------------------------------------------------
        line(nt,(sum((SPOTA_XYZ_um-SPOTB_XYZ_um).^2,2)).^0.5,'Parent',ax(4),...
            'Color',[0 0 1],'Marker','o','MarkerEdgeColor',[0.2 0.2 0.2],'MarkerFaceColor',[0 0 1]);
        set(ax(4),'XGrid','on','YGrid','on')
        xlabel(ax(4),'t(s)')
        ylabel(ax(4),'d_{AB}(\mum)')
        title(ax(4),'Centrosomes distances')
        
        
        %--------------------------------------------------------------------------
        cmap = parula(length(nt));
        for idxT=1:length(nt)
            line(d{idxT},mch1{idxT},'Color',cmap(idxT,:),'Parent',ax(2))
            line(d{idxT},mch2{idxT},'Color',cmap(idxT,:),'Parent',ax(3))
        end
        colormap(ax(2),cmap)
        set(ax(2),'CLim',[nt(1) nt(end)],'XGrid','on','YGrid','on')
        colormap(ax(3),cmap)
        set(ax(3),'CLim',[nt(1) nt(end)],'XGrid','on','YGrid','on')
        
        if isAbs
            xlabel(ax(2),'d(\mum)')
        else
            xlabel(ax(2),'d(%)')
        end
        ylabel(ax(2),'Ch1(-)')
        title(ax(2),'Ch1 vs. d @ t')
        hcb(1) = colorbar(ax(2));
        hcb(1).Label.String = 't(s)';
        
        
        if isAbs
            xlabel(ax(3),'d(\mum)')
        else
            xlabel(ax(3),'d(%)')
        end
        ylabel(ax(3),'Ch2(-)')
        title(ax(3),'Ch2 vs. d @ t')
        hcb(2) = colorbar(ax(3));
        hcb(2).Label.String = 't(s)';
        
        
        
        %--------------------------------------------------------------------------
        dtot = unique(cat(2,d{:}));
        Z_ch1 = nan(length(d),length(dtot));
        Z_ch2 = nan(length(d),length(dtot));
        for idxT=1:length(d)
            [~,Locb] = ismember(d{idxT},dtot);
            Z_ch1(idxT,Locb) = mch1{idxT};
            Z_ch2(idxT,Locb) = mch2{idxT};
        end
        
        
        waterfall(ax(5),dtot,nt,Z_ch1)
        colormap(ax(5),hot)
        set(ax(5),'ZLim',get(ax(2),'YLim'),'YLim',[nt(1) nt(end)],...
            'YDir','reverse','CLim',[min(Z_ch1(:)) max(Z_ch1(:))],...
            'Box','on','BoxStyle','full',...
            'XGrid','on','YGrid','on','ZGrid','on')
        if isAbs
            xlabel(ax(5),'d(\mum)')
        else
            xlabel(ax(5),'d(%)')
        end
        ylabel(ax(5),'t(s)')
        zlabel(ax(5),'Ch1(-)')
        title(ax(5),'Ch1 vs. d @ t')
        hcb(3) = colorbar(ax(5));
        hcb(3).Label.String = 'Ch1(-)';
        
        waterfall(ax(6),dtot,nt,Z_ch2)
        colormap(ax(6),hot)
        set(ax(6),'ZLim',get(ax(3),'YLim'),'YLim',[nt(1) nt(end)],...
            'YDir','reverse','CLim',[min(Z_ch2(:)) max(Z_ch2(:))],...
            'Box','on','BoxStyle','full',...
            'XGrid','on','YGrid','on','ZGrid','on')        
        if isAbs
            xlabel(ax(6),'d(\mum)')
        else
            xlabel(ax(6),'d(%)')
        end
        ylabel(ax(6),'t(s)')
        zlabel(ax(6),'Ch2(-)')
        title(ax(6),'Ch2 vs. d @ t')
        hcb(4) = colorbar(ax(6));
        hcb(4).Label.String = 'Ch2(-)';
        
        
        %--------------------------------------------------------------------------
        colormap(ax(8),hot)
        imshow(Z_ch1,[],'Parent',ax(8))
        
        if isAbs
            xlabel(ax(8),'d(\mum)')
        else
            xlabel(ax(8),'d(%)')
        end
        
        ylabel(ax(8),'t(s)')
        title(ax(8),'Ch1 kymograph')
        set(ax(8),'Visible','on','CLim',[min(Z_ch1(:)) max(Z_ch1(:))])
        hcb(5) = colorbar(ax(8));
        hcb(5).Label.String = 'Ch1(-)';
        
        
        set( ax(8),'DataAspectRatioMode','auto');
        
        txt_nt = num2cell(nt);
        txt_nt = cellfun(@(x) num2str(x),txt_nt,'UniformOutput',false);
        txt_dtot = num2cell(dtot);
        txt_dtot = cellfun(@(x) num2str(x),txt_dtot,'UniformOutput',false);
        
        tmpYTickidx = [1:2:length(nt)];% get(ax(8),'YTick');
        set( ax(8),'YTick',tmpYTickidx,'YTickLabel',txt_nt(tmpYTickidx));
        
        tmpXTickidx = [1:6:length(dtot)];%get(ax(8),'XTick');
        set( ax(8),'XTick',tmpXTickidx,'XTickLabel',txt_dtot(tmpXTickidx));
        
        
        colormap(ax(9),hot)
        imshow(Z_ch2,[],'Parent',ax(9))
        if isAbs
            xlabel(ax(9),'d(\mum)')
        else
            xlabel(ax(9),'d(%)')
        end
        ylabel(ax(9),'t(s)')
        title(ax(9),'Ch2 kymograph')
        set(ax(9),'Visible','on','CLim',[min(Z_ch2(:)) max(Z_ch2(:))])
        hcb(6) = colorbar(ax(9));
        hcb(6).Label.String = 'Ch2(-)';
        
        set( ax(9),'DataAspectRatioMode','auto');
        
        
        set( ax(9),'YTick',tmpYTickidx,'YTickLabel',txt_nt(tmpYTickidx));
        set( ax(9),'XTick',tmpXTickidx,'XTickLabel',txt_dtot(tmpXTickidx));
        
        
        
        
        print('-dpng','-r300',[FilePathName,FileParam],'-loose')
        
    end





end