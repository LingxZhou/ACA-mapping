%% Here to calculate the multiparametrics of a 3D image with adaptive compression analysis
% This is the main program
clear;

%%--data loading--
sz1 = 612; sz2 = 612; sz3 = 192; % modify 1
% define sz, sz2, and sz3 according to the size of the 3D image

I = zeros(sz1,sz2,sz3); Im = zeros(sz1,sz2,sz3); 
for i = 1:sz3
    I(:,:,i) = imread(['C:\Users\97916\Desktop\ACA mapping\Example\',num2str(i),'.tif']);
    Im(:,:,i) = imread(['C:\Users\97916\Desktop\ACA mapping\U-Net_Pytorch\output\',num2str(i),'.png']); 
    % 'I' is the intensity image, and 'Im' is the segmentation image output by U-net.
    % modify 2，modify the image reading path accordingly
end
Im = im2double(Im); I = im2double(I);

%%--Parameter settings--
armdx = 12; armdz = 12;  % modify 3,
% define caculation window. the width of window in both 'x' and 'y' dimensions is 
% both '2*armdx+1'.the width of window in 'z' dimensions is '2*armdz+1'. 'armdx' 
% and 'armdz' is chosen based on the diamater of vessels. Typically '2*armdx(z)+1'
% is 2 to 3 times the diameter of vessel so as to provide optimal accuracy.
ske_layer = 16; % modify 4
% This determines the number of image layers used for fast image skeletonization.
% The image will be divided into N*ske_layers for skeletonization. We recommand 
% that ske_layer be no more than 32.
loc = 0.1; % modify 5
% This determines the thickness distributin position used as compression factor.
% The value range is 0~1.
filton = 1; % modify 6
% 'filton' is a binary variable that tells program whether to filter the
% data prior to analysis. 1: filtering; 0: no filtering
para = 1; % modify 7
% 'para' is the ratio of sampling frequency between 'xy' dimension and 'z'
% dimension
threshint = 0; 

%%--thickness calculation and compression factor(CF) determination--
[thickness,CF] = fast_thick_3D(Im,ske_layer,loc);

%%--multiparametric analysis--
Im = imresize3(Im,1/CF,'nearest');
[rsz1,rsz2,rsz3]=size(Im);
armdx=ceil(armdx/CF); armdz=ceil(armdz/CF);
% Here implement resampling of data and calculation windows

% theta and phi orientation calculation
maska3d = Im > threshint;
meanint = mean(Im(maska3d)); threshinta = 0.45*meanint; maska3da = Im > threshinta;
finalmaska = maska3d.*maska3da;
finalmask = zeros(rsz1,rsz2,rsz3);
for j=1:rsz3
    Maskf1= finalmaska(:,:,j);
    Maskf2 = miprmdebrisn(Maskf1,15);
    finalmask(:,:,j) = Maskf2;
end
finalmask = logical(finalmask);
[aSDE,pSDER,bSDER,gSDER] = calcfibang3D_opt(Im,armdx,armdz,filton,para);
% 'aSDE' is the calculated theta stack, 'pSDER' is the calculated phi stack, 
% 'bSDER' is the calculated beta stack, and 'gSDER' is the calculated gamma stack.


% directioanl variance calculation
[Vmatr,Vvaluelocal,Vvalueentire] = calcfibvarspeed_nf(armdx,armdz,aSDE,bSDER,gSDER,finalmask);
% 'Vmatr' is the voxel-wise 3D directional variance stack, with every voxel
% filled with a directional variance value showing the vessel organization
% within the localized window region surrounding it.
% 'Vvaluelocal' is a mean value of the directional variance of this 3D 
% stack, with directional variance acquired by the localized window. Notice
% that only the regions identified as vessels are considered and contribute
% to the mean value.
% 'Vvalueentire' is the 3D directional variance value acquired from the
% orientation information of all the vessels within the entire 3D stack

%waviness calculation
doublefm = double(finalmask);
wav_theta = wavwin_cal3D(armdx,armdx,armdz,aSDE,doublefm);
wav_beta = wavwin_cal3D(armdx,armdx,armdz,bSDER,doublefm);
wav_gama = wavwin_cal3D(armdx,armdx,armdz,gSDER,doublefm);
wavmatr3D = 1/3*wav_theta + 1/3*wav_beta + 1/3*wav_gama;
% 'wavmatr3D' is the voxel-wise 3D waviness stack, with every voxel
% filled with a waviness value showing the vessel organization
% within the localized window region surrounding it

nanind = isnan(wavmatr3D);
wavmatr3D(nanind) = 0;
finalmaskw = finalmask;
finalmaskw(nanind) = 0;
Meanwav = mean(wavmatr3D(finalmaskw));
% 'Meanwav' is a mean value of the waviness of this 3D stack,
% with waviness acquired by the localized window. Notice
% that only the regions identified as vessels are considered
% and contribute to the mean value.

% local coverage calculation
doublefm = double(finalmask);
armdx = armdx*2; armdz = armdz*2;
winnum = (2*armdx+1)*(2*armdx+1)*(2*armdz+1);
win_th = 15; % modify 8
% 'win_th' is the threshold window size where the computational time for LC calculation 
% is roughly equal in both the spatial and frequency domains.
LCmatr3D = LCcalc3D_SorF(armdx,armdx,armdz,doublefm,winnum,win_th);
% 'LCmatr3D' is the voxel-wise 3D local coverage stack, with every voxel
% filled with a local coverage value showing the vessel organization
% within the localized window region surrounding it

nanind = isnan(LCmatr3D);
LCmatr3D(nanind) = 0;
finalmaskl = finalmask;
finalmaskl(nanind) = 0;
Meanlc = mean(LCmatr3D(finalmaskw));
% 'Meanlc' is a mean value of the local coverage of this 3D 
% stack, with local coverage acquired by the localized window.
% Notice that only the regions identified as vessels are 
% considered and contribute to the mean value

%%--post-processing--
% For post-processing, we prepare the 'pretty' images of theta orientation,
% phi orientation,directional variance, waviness and local coverage.
% In these 'pretty' images, the intensity image 
% is used to provide the contrast of vessel features, and characteristic
% maps are labeled by different colors to show the information
I = I/max(max(max(I)));

bright = 0.99; 
dark = 0.01;

% prepare the pretty theta orientation stack
uplim = 180;
botlim = 0;
prettyima = zeros(sz1,sz2,3,sz3);
aSDE = imresize3(aSDE,[sz1 sz2 sz3],'nearest');
for mm = 1:sz3
    shgima = I(:,:,mm);
    thetaima = aSDE(:,:,mm);
    prettyima(:,:,:,mm) = prettymap(thetaima,shgima,'none',hsv(64),uplim,botlim,bright,dark);
    % Here 'thetaima' is the phi orientation map. For theta orientation,
    % the range is from 0 to 180. 
    imwrite(double(prettyima(:,:,:,mm)),['C:\Users\97916\Desktop\ACA mapping\Results\theta\z_',num2str(mm),'.tif']); % modify 9
end

% prepare the pretty phi orientation stack
uplim = 180;
botlim = 0;
pSDER = imresize3(pSDER,[sz1 sz2 sz3],'nearest');
for mm = 1:sz3
    shgima = I(:,:,mm);
    phiima = pSDER(:,:,mm);
    prettyima(:,:,:,mm) = prettymap(phiima,shgima,'none',hsv(64),uplim,botlim,bright,dark);
    % Here 'phiima' is the phi orientation map. For phi orientation,
    % the range is from 0 to 180. 
    imwrite(double(prettyima(:,:,:,mm)),['C:\Users\97916\Desktop\ACA mapping\Results\phi\z_',num2str(mm),'.tif']); % modify 10
end

% prepare the pretty directional variance stack
uplim = 1; % modify 11
botlim = 0; % modify 12
Vmatr = imresize3(Vmatr,[sz1 sz2 sz3],'nearest');
for mm = 1:sz3
    shgima = I(:,:,mm);
    varima = Vmatr(:,:,mm);
    prettyima(:,:,:,mm) = prettymap(varima,shgima,'none',jet(64),uplim,botlim,bright,dark);
    % Here 'varima' is the directional variance map. For directional
    % variance, the range is from 0 to 1. Therefore, 'uplim' and 'botlim'
    % are modified accordingly.
    imwrite(double(prettyima(:,:,:,mm)),['C:\Users\97916\Desktop\ACA mapping\Results\variance\z_',num2str(mm),'.tif']); % modify 13
end

% prepare the pretty waviness stack
uplim = 0.7; % modify 14
botlim = 0; % modify 15
wavmatr3D = imresize3(wavmatr3D,[sz1 sz2 sz3],'nearest');
for mm = 1:sz3
    shgima = I(:,:,mm);
    wavima = wavmatr3D(:,:,mm);
    prettyima(:,:,:,mm) = prettymap(wavima,shgima,'none',jet(64),uplim,botlim,bright,dark);
    % Here 'wavima' is the waviness map. For waviness is 
    % usually lower than directional variance, therefore, 
    % 'uplim' and 'botlim' are modified accordingly.
    imwrite(double(prettyima(:,:,:,mm)),['C:\Users\97916\Desktop\ACA mapping\Results\waviness\z_',num2str(mm),'.tif']); % modify 16
end

% prepare the pretty local coverage stack
uplim = 0.8; % modify 17
botlim = 0; % modify 18
LCmatr3D = imresize3(LCmatr3D,[sz1 sz2 sz3],'nearest');
for mm = 1:sz3
    shgima = I(:,:,mm);
    LCima = LCmatr3D(:,:,mm);
    prettyima(:,:,:,mm)= prettymap(LCima,shgima,'none',jet(64),uplim,botlim,bright,dark); 
    % Here 'LCima' is the local coverage map. For local
    % coverage, the range is from 0 to 1. Therefore, 'uplim' and 'botlim'
    % are modified accordingly
    imwrite(double(prettyima(:,:,:,mm)),['C:\Users\97916\Desktop\ACA mapping\Results\LC\z_',num2str(mm),'.tif']); % modify 19
end