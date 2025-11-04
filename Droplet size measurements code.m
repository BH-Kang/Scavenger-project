clear all;close all
tic
set(0,'DefaultFigureVisible','on')

MakeJPGFile = 1;

%%%%%%%%%%%%%%%%%%%%% CHANGE %%%%%%%%%%%%%%%%%%%%%%%%
%% designate a file path
path_in = 'XXX'; 

%% microscope parameter
objective = 10;
magnification = 1;
binsize = 2; % 2048x2048 -> 1; 1024x1024 -> 2; 512x512 -> 4
scale = 10; %um scale bar length you want to make 

thr = 0.5; % Brightness Threshold  
thrarea_min = 100; % Size threshold, small
thrarea_max = inf; % Size threshold, large

nbins = 10; % number of range for histogram
%%%%%%%%%%%%%%%%%%%%% CHANGE %%%%%%%%%%%%%%%%%%%%%%%%

path_out = fullfile(path_in,'output_ObjectMeasure');
%% calculate conversion from pixels to microns for M2
pixelsize = 2440*(10/objective)*(binsize/4)*(1/magnification); %% nm/pixel
pixelsize = pixelsize*1e-3; %%% micron/pixel


if ~isfolder(fullfile(path_out))
    mkdir(fullfile(path_out));
end

%% look for a file with .seq extension, your experiment
cd(path_in)
fileListing = dir('**/*.tif*');
numberOfFiles = length(fileListing);

xo = []; yo = []; A = []; C = [];
Area_avg_um2 = zeros(numberOfFiles,1);
Cir_avg = zeros(numberOfFiles,1);
Area_std_um2 = zeros(numberOfFiles,1);
Diameter_avg_um=zeros(numberOfFiles,1);
Diameter_std_um=zeros(numberOfFiles,1); 
fileName_list = cell(numberOfFiles,1);
for k = 1 : numberOfFiles
    path_tif = fullfile(fileListing(k).folder, fileListing(k).name);
    fprintf('Processing image %d of %d : %s...\n', k, numberOfFiles, path_tif);
    
    %% read tiff file
    % fluorphore and bright field image
    warning off
    tiff_info = imfinfo(path_tif);
    tiff_stack = Tiff(path_tif,'r');
    nframe = length(tiff_info);
    warning on
    
    %% find the file name
    filesepLoc = strfind(fileListing(k).folder, filesep);
    path_folder = fileListing(k).folder(1:filesepLoc(end));
    FileName = fileListing(k).folder(filesepLoc(end)+1:end);
%     FileName_upper = fileListing(k).folder(filesepLoc(end-1)+1:filesepLoc(end)-1);

    %% allocate path out file
    jpg_out = fullfile(path_out,strcat(FileName, '_label.jpg'));
    fig_out = fullfile(path_out,strcat(FileName,'_histogram.jpg'));
    avi_out = fullfile(path_out,strcat(FileName, '_label'));
    hist_out = fullfile(path_out,strcat(FileName, '_Diam_Area_List.xlsx'));
    % delete old versions of save files and open new video writers
    if exist(jpg_out)
        delete(jpg_out)
    end
    
    % get frame height and width
    frameHeight = tiff_info.Height;
    frameWidth = tiff_info.Width;
    
    % set scale bar text
    xl_text = [0.8*frameWidth 0.8*frameWidth-scale/pixelsize];
    yl_text = 0.95*frameHeight;
    
    % scale bar position
    leftToRight = round(0.8*frameWidth + [-1 1]*scale/pixelsize/2);
    TopToBottom = round(0.99*frameHeight + [-1 1]*frameHeight/170);

    for i_frame = 1:nframe
        % Use setDirectory method to access the image in the stack
        warning off
        tiff_stack.setDirectory(i_frame);
        warning on
    
        % Read image data
        frame = tiff_stack.read();
        frame = imadjust(frame);
    
        frame = double(frame)/(2^(tiff_info.BitDepth)-1);%max(frame(:));
        frame(frame>1) = 1;
        
        %%% find particle aggregation
        %%% Binary image of the frame
        Bim = frame > thr;
        % Bim = ~Bim;
        Bim = bwareafilt(Bim,[100 inf]);
        Bim_stats = regionprops(Bim,'centroid','Area',"Circularity");
        % for perfect circle the circularity value is 1

        N = length(Bim_stats); Area_um2 = zeros(N,1); Cir = zeros(N,1);xc = zeros(N,1);yc = zeros(N,1);
        id = 1;
        for i = 1:N
            %save particle centers
            xc(id,1) = (Bim_stats(i).Centroid(:,1));
            yc(id,1) = (Bim_stats(i).Centroid(:,2));
            Area_um2(id,1) = (Bim_stats(i).Area)*(pixelsize^2); % um2
            Cir(id,1) = Bim_stats(i).Circularity;
            id = id + 1;
        end

        % Delete particles at each edge
        Diameter_um = sqrt(4*Area_um2/pi); % um
        Diameter_avg_um(k) = mean(Diameter_um,'all'); % um
        Diameter_std_um(k) = std(Diameter_um,0,'all'); % um
        Area_avg_um2(k) = mean(Area_um2,'all'); % um2
        Area_std_um2(k) = std(Area_um2,0,'all'); % um2
        Cir_avg(k) = mean(Cir, 'all');
        fileName_list{k} = FileName;

        fig1 = figure;
        axe1 = axes();
        set(fig1,'Visible','off');

        % label scale bar
        frame_scalebar = frame;
        frame_scalebar(TopToBottom(1):TopToBottom(2), leftToRight(1):leftToRight(2)) = 2^16-1;
%         imshow(frame_scalebar); 
        imshow(imoverlay(frame_scalebar, Bim, 'r')); hold on

        % label scale bar 
        text(axe1,mean(xl_text),yl_text,strcat(num2str(round(scale,1)),' µm'),'Color','w','FontSize',frameWidth/100)

%         imshow(frame)
%         figure;imshow(Bimline)
%         imshowpair(frame,imoverlay(frame, Bim, 'r'),'montage')
%         set(fig1,'Visible','on');

        plot(axe1,xc,yc, 'or','MarkerSize',4,'MarkerFaceColor','k');
        Area_str = num2str(Area_um2,'%.1f');
        N_id = 1:N;
        N_id = N_id';
        str = strcat('#',num2str(N_id),'=',Area_str, 'μm^2');
        text(axe1,xc,yc,str,'Color','w','FontSize',10);

        frame_label = frame2im(getframe(gca));
    %     IM=getframe(gcf);
        IM=getframe(gca); % this is only for frame not including colorbar

        if MakeJPGFile ~= 0
            % save jpg file
            imwrite(frame_label, jpg_out);
        end

        % histogram
        fig2 = figure();
                set(fig2,'Visible','off');
                axe2 = axes();
                fig2.Units = "inches";
                axe2.LineWidth = 1.5;
                axe2.FontSize  = 16;
                axe2.NextPlot  = "add";
                axe2.Box       = "on";
                axe2.XLabel.String = "Diameter/ μm";
                axe2.XAxis.FontWeight = 'bold';
                axe2.XLabel.FontSize = 24;
                axe2.XLabel.FontWeight = 'normal';
                axe2.YLabel.String = "Number of samples";
                axe2.YAxis.FontWeight = 'bold';
                axe2.YLabel.FontSize = 24;
                axe2.YLabel.FontWeight = 'normal';
                axe2.XMinorTick = "on";
                axe2.YMinorTick = "on";
                axe2.TickLength = [0.03 0.05];
                legend('off')
                axis square
        h = histogram(Diameter_um, nbins);
        %%%%%%%%%%%%%%%%%%%%%%%%%%% save figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        exportgraphics(fig2,fig_out,'Resolution',300)
        % savefig(fig2,fig_out);
       
        % save data
        % save_temp = [N_id Diameter_um Area_um2];
        % writematrix(save_temp, hist_out,'Delimiter','tab')
        save_temp = table(N_id, Diameter_um, Area_um2);
        writetable(save_temp,hist_out)  

        clf(fig2)
    end
    
    clf
    close all hidden
    close all
    
    toc
end
T = table(fileName_list,Area_avg_um2,Area_std_um2,Diameter_avg_um,Diameter_std_um,Cir_avg); 
path_save = fullfile(path_out,'fileName_AreaAvg_AreaSTD_DiamAvg_DiamSTD_Cir.xlsx');
writetable(T,path_save)  

% writematrix(T, path_save,'Delimiter','tab')

