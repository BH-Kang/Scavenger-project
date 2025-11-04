% IntensityMeasure BW

close all;
clc, clear;

addpath(genpath(cd));

% initial setting
dir = 'C:\Users\XXX'; % Name of folder that contains gel image file
filename = 'YYY'; % Name of gel image file (*.tif)
bit = 8;  % bit
inversion = 1; % 0: no inversion, 1: inversion

row = 100;
col = 100;
frame = ones(col,row);
bound = zeros(col,row);
bound(1,:) = 1;
bound(:,1) = 1;
bound(end,:) = 1;
bound(:,end) = 1;
x = 0;
y = 0;

while(1)
    temp = inputdlg({'Directory','Filename','Bit (8 or 16)','Inversion (0 or 1)'},' ',1,{sprintf('%s',dir),sprintf('%s',filename),sprintf('%d',bit),sprintf('%d',inversion)});
    if size(temp,1)==0
        break;
    end
    dir = temp{1};
    filename = temp{2};
    bit = str2double(temp{3});
    inversion = str2double(temp{4});
    
    scale = 1;
    IMG = im2double(imread(sprintf('%s/%s.tif',dir,filename)));
    if inversion~=0
        IMG = 1-IMG;
    end
    
    U = 1;
    D = size(IMG,1);
    L = 1;
    R = size(IMG,2);
    img = IMG(:,:,1);
    
    while(1)
        if x~=0 && y~=0
            temp(y+1:y+col,x+1:x+row) =  temp(y+1:y+col,x+1:x+row)+bound;
            imshow(temp), title(sprintf('x %d',scale));
        else
            temp = img*scale;
            imshow(temp), title(sprintf('x %d',scale));
        end
        [x,y,button] = ginput(1);
        x = round(x);
        y = round(y);
        if button==1            
            fprintf('Filename:   %s\nPosition:    %d %d\nIntensity:    %.2f\n\n',filename,x,y,sum(sum(img(y+1:y+col,x+1:x+row).*(2^bit-1))));
        else
            x = 0;
            y = 0;
            if button==3        
                temp = inputdlg({'Frame Row', 'Frame Column', 'U (Min : 1)', sprintf('D (Max : %d)',size(IMG,1)), 'L (Min : 1)', sprintf('R (Max : %d)',size(IMG,2))},' ',1,{sprintf('%d',row),sprintf('%d',col),sprintf('%d',U),sprintf('%d',D),sprintf('%d',L),sprintf('%d',R)});
                if size(temp,1)==0
                    break;
                end
                row = str2double(temp{1});
                col = str2double(temp{2});
                U = str2double(temp{3});
                D = str2double(temp{4});
                L = str2double(temp{5});
                R = str2double(temp{6});
                img = IMG(U:D,L:R,1);
                frame = ones(col,row);
                bound = zeros(col,row);
                bound(1,:) = 1;
                bound(:,1) = 1;
                bound(end,:) = 1;
                bound(:,end) = 1;
            elseif button==30   % ก่: Contrast increase
                scale = scale+1;
            elseif button==31   % ก้: Contrast decrease
                scale = scale-1;
            else
                break;
            end
        end
    end
end