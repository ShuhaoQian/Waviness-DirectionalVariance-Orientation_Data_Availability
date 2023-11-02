 %% Here to calculate the waviness and directional variance of collagen in a 3D context.
% This is the main program

%% Here to load the 3D image stack of fiber-like structure.
clear;
sz1 = 512;
sz2 = 512;
sz3 = 61; % modify 1
% Define sz1, sz2, and sz3 according to the size of the 3D image.
imgpath = 'D:\Matlab code for Data Availability\Examples\'; % modify 2 
% Define the path of the folder where the images are saved.
fibstack = zeros(sz1,sz2,sz3);
for i = 1:sz3
    fibstack(:,:,i) = imread(strcat(imgpath,num2str(i),'.tif')); % modify 3 
    % 'imgstack' is the 3D image stack used for waviness analysis.
    % Define the image format.
end

%% Here to create the mask selecting fiber regions.
threshint = 10;  % modify 4 
% 'threshint' is the raw background intensity acquired by averaging intensity of several regions identified as background
maska3d = fibstack > threshint; % maska3d is a raw binary mask based on the raw background intensity

meanint = mean(fibstack(maska3d)); % calculate the mean intensity within the raw binary mask
threshinta = 0.45*meanint; % modify 5 
% 'threshinta' is an inmproved intensity threshold, and the factor '0.45' can be adjusted according to different sample
maska3da = fibstack > threshinta; % 'maska3da' is an improved binary mask

% Here to create the finalmask
finalmaska = maska3d.*maska3da; % 'finalmaska' is the binary mask acquired based on the above two masks

finalmask = zeros(sz1,sz2,sz3);
for j = 1:sz3
    Maskf1 = finalmaska(:,:,j); 
    Maskf2 = miprmdebrisn(Maskf1,15); % 'miprmdebrisn.m' is a function which is able to remove very tiny structures
    finalmask(:,:,j) = Maskf2;
end
finalmask = logical(finalmask); % 'finalmask' is the final binary mask selecting the fiber-only region

%% Here to calculate the voxel-wise 3D orientation.
armdx = 1; % modify 6
% the width of window in both 'x' and 'y' dimensions is both '2*armdx+1'.
% 'armdx' is chosen based on the diamater of fibers. Typically '2*armdx+1'
% is 2 to 3 times the diameter of fiber so as to provide optimal accuracy,
% as published in Biomed. Opt. Express 2015 (6), 2294-2310
armdz = 1; % modify 7
% '2*armdz+1' is the width of the window in 'z' dimension. We would expect
% the actual width of the window in 'x', 'y' dimensions is equal to the one
% in 'z' dimension. Therefore, 'armdz' is determined by the sampling 
% frequency in the 'xy' dimension and 'z' dimension  
filton = 1; % modify 8
% 'filton' is a binary variable that tells program whether to filter the
% data prior to analysis. 1: filtering; 0: no filtering
para = 1; % modify 9
% 'para' is the ratio of sampling frequency between 'xy' dimension and 'z'
% dimension
[aSDE,pSDER,bSDER,gSDER] = calcfibang3D(fibstack,armdx,armdz,filton,para); 
% 'calcfibang3D.m' is the function which calculates the voxel-wise 3D
% orientation of a 3D image. 'aSDE' is the calculated theta stack, 'pSDER'
% is the calculated phi stack, 'bSDER' is the calculated beta stack, and
% 'gSDER' is the calculated gamma stack

%% Here to calculate the voxel-wise fiber waviness.
wavarmdx = 1;% modify 10
wavarmdy = 1;% modify 11
wavarmdz = 1;% modify 12
% '2*wavarmdx+1', '2*wavarmdy+1' and '2*wavarmdz+1' are width of the window
% in 'x', 'y' and 'z' dimensions respectively for the calculation of the waviness. 
%'wavarmdx' and 'wavarmdx' are chosen based on the diamater of fibers. 
% Typically '2*wavarmdx/wavarmdy+1' is 5 to 6 times the diameter of fiber. 
% We would expect the actual width of the window in 'x', 'y' dimensions is
% equal to the one in 'z' dimension. Therefore, 'wavarmdz' is determined by
% the sampling frequency in the 'xy' dimension and 'z' dimension.  
doublefm = double(finalmask);
wav_theta = wavwin_cal3D(wavarmdy,wavarmdx,wavarmdz,aSDE,doublefm);
wav_beta = wavwin_cal3D(wavarmdy,wavarmdx,wavarmdz,bSDER,doublefm);
wav_gama = wavwin_cal3D(wavarmdy,wavarmdx,wavarmdz,gSDER,doublefm);
wavmatr3D = 1/3*wav_theta + 1/3*wav_beta + 1/3*wav_gama;
% 'wavwin_cal3D.m' is the function which calculates the voxel-wise
% waviness of a 3D fiber stack. Width of window of three dimensions, 
% orientation matrix, and fiber-background binary matrix are the input 
% parameters for the waviness calculation function. Waviness of three
% orienation angles are calculated respectively, and 'wavmatr3D' is the 
% final voxel-wise waviness matrix for the fiber stack.
nanind = isnan(wavmatr3D);
wavmatr3D(nanind) = 0;
finalmask(nanind) = 0;
Meanwav =mean(wavmatr3D(finalmask));
% 'Meanwav' is the mean waviness of all the fiber pixels selected  by 
% 'finalmask' in the image.

%% Here to calculate the voxel-wise 3D directional variance based on a localized window
bwhole = sqrt(1./(tan(2*bSDER*pi/180).^2)+1./(tan(2*gSDER*pi/180).^2));
Cwhole = (bwhole./sqrt(1+bwhole.^2)).*cos(2*aSDE*pi/180);
Swhole = (bwhole./sqrt(1+bwhole.^2)).*sin(2*aSDE*pi/180);
Zwhole = zeros(sz1,sz2,sz3);
% 'bwhole', 'Cwhole', 'Swhole' and 'Zwhole' are all assisting variables
% used to acquire the directional variance results.  
for m = 1:sz1
   for n = 1:sz2
     for o = 1:sz3
        if bSDER(m,n,o) <= 90   
            Zwhole(m,n,o) = 1/sqrt(1+bwhole(m,n,o)^2);
        else
            Zwhole(m,n,o) = -1/sqrt(1+bwhole(m,n,o)^2);
        end
     end
   end
end

Cwholemean = zeros(sz1,sz2,sz3);
Swholemean = zeros(sz1,sz2,sz3);
Zwholemean = zeros(sz1,sz2,sz3);
ll = 0;

armdxx = 1; % modify 13
armdzz = 1; % modify 14
% 'armdxx' and 'armdzz' are defining the size of the window used to
% acquired the 3D directional variance (xy: 2*armdxx+1; z: 2*armdzz+1).
% Typically, these two parameters are defined the same as 'armdx' and
% 'armdz', so that the window used to calculate the voxel-wise orientation 
% and the one used to acquire localized directional variance have the same 
% size. However, they can be set to other values according to the purpuse 
% of analysis 

h = waitbar(0,'Please Wait');
ijk = 0;
nijk = (2*armdxx+1)*(2*armdxx+1)*(2*armdzz+1);
tic

for  i = -armdxx:armdxx
    for j = -armdxx:armdxx
        for k = -armdzz:armdzz
            ijk = ijk+1;
            waitbar(ijk/nijk)
            Cim = circshift(Cwhole,[i j k]);
            Sim = circshift(Swhole,[i j k]);
            Zim = circshift(Zwhole,[i j k]);
                        
            Cwholemean = Cwholemean+Cim;
            Swholemean = Swholemean+Sim;
            Zwholemean = Zwholemean+Zim;
            ll = ll+1;
        end
    end
end
    
toc
close(h)
Cwholemean = Cwholemean/ll;
Swholemean = Swholemean/ll;
Zwholemean = Zwholemean/ll;
DVmatr = 1 - sqrt(Cwholemean.^2+Swholemean.^2+Zwholemean.^2);
DVmatr(isnan(DVmatr)) = 0;
% 'DVmatr' is the voxel-wise 3D directional variance stack, with every voxel
% filled with an directional variance value showing the fiber organization
% within the localized window region surrounding it

%% Here to do post-processing.
% For post-processing, we prepare the 'pretty' images of orientation and
% waviness. In these 'pretty' images, the raw intensity image 
% is used to provide the contrast of fiber features, and the orientation or
% waviness maps are labeled by different colors to show the orientation 
% or waviness information

% first,prepare the pretty orientation stack
fibstackre = fibstack;
fibstackre = fibstackre/max(max(max(fibstackre)));

uplim = 180;
botlim = 0;
bright = 0.99; 
dark = 0.01;

prettytheta = zeros(sz1,sz2,3,sz3);
for mm = 1:sz3
    shgima = fibstackre(:,:,mm);
    thetaima = aSDE(:,:,mm);
    prettyima = prettymap(thetaima,shgima,'none',hsv(64),uplim,botlim,bright,dark);
    % 'prettymap.m' is the function used to create 'pretty' images.
    % 'thetaima' is the theta orientation map, 'shgima' is the raw SHG
    % image. 'hsv(64)' designates the color scheme. 'uplim' and 'botlim' 
    % are the upper and bottom limits of the orientation index. The range 
    % of both theta and phi is from 0 to 180. 'bright' and 'dark' are used 
    % to enhance the contrast of the image
    prettytheta(:,:,:,mm) = prettyima;
end
prettyphi = zeros(sz1,sz2,3,sz3);
for mm = 1:sz3
    shgima = fibstackre(:,:,mm);
    phiima = pSDER(:,:,mm);
    prettyima = prettymap(phiima,shgima,'none',hsv(64),uplim,botlim,bright,dark);
    % Here 'phiima' is the phi orientation map 
    prettyphi(:,:,:,mm) = prettyima;
end
prettygama = zeros(sz1,sz2,3,sz3);
for mm = 1:sz3
    shgima = fibstackre(:,:,mm);
    phiima = gSDER(:,:,mm);
    prettyima = prettymap(phiima,shgima,'none',hsv(64),uplim,botlim,bright,dark);
    % Here 'phiima' is the gama orientation map 
    prettygama(:,:,:,mm) = prettyima;
end
prettybeta = zeros(sz1,sz2,3,sz3);
for mm = 1:sz3
    shgima = fibstackre(:,:,mm);
    phiima = bSDER(:,:,mm);
    prettyima = prettymap(phiima,shgima,'none',hsv(64),uplim,botlim,bright,dark);
    % Here 'phiima' is the beta orientation map 
    prettybeta(:,:,:,mm) = prettyima;
end


% Second, prepare the pretty waviness and directional variance stack.
uplim = 1;
botlim = 0;
prettywav = zeros(sz1,sz2,3,sz3);
for mm = 1:sz3
    shgima = fibstackre(:,:,mm);
    wavima = wavmatr3D(:,:,mm);
    prettyima = prettymap(wavima,shgima,'none',jet(64),uplim,botlim,bright,dark);
    % Here 'prettywav' is the waviness map.
    prettywav(:,:,:,mm) = prettyima;
end

prettyDV = zeros(sz1,sz2,3,sz3);
for mm = 1:sz3
    shgima = fibstackre(:,:,mm);
    DVima = DVmatr(:,:,mm);
    prettyima = prettymap(DVima,shgima,'none',jet(64),uplim,botlim,bright,dark);
    % Here 'prettyDV' is the directional variance map.
    prettyDV(:,:,:,mm) = prettyima;
end