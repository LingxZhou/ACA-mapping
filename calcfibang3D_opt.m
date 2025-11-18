function [aSDE,pSDER,bSDER,gSDER] = calcfibang3D_opt(shgim,armdx,armdz,filton,para)
% This function calculates the voxel-wise orientation. The output and input
% indices are explained in the main program
im = shgim;

im = double(im);
sz = size(im,1);
sz2 = size(im,2);
sz3 = size(im,3);

% Here starts the determination of theta
terimc = zeros(sz,sz2,sz3 + 2*armdz); 
terimc(:,:,1:armdz) = im(:,:,sz3-armdz+1 : sz3); 
terimc(:,:,(armdz + 1):(armdz + sz3)) = im; 
terimc(:,:,(armdz + sz3 + 1):(armdz + sz3 + armdz)) = im(:,:,1 : armdz); 

imc = zeros(sz,sz2,sz3);

for i = 1:sz3
imc(:,:,i) = mean(terimc(:,:,i:i+armdz*2),3); 
end
clear terimc

meanc = zeros(sz,sz2,sz3);
meanc = calcfibangspeed(imc,armdx,filton);
clear imc

% Here starts the determination of beta
terimj = zeros( sz + armdz*2 , sz2,sz3); 
terimj(1 : armdz,:,:) = im(sz-armdz+1:sz,:,:); 
terimj(armdz + 1 : armdz + sz,:,:) = im;
terimj(armdz + sz + 1: armdz + sz + armdz,:,:) = im(1:armdz,:,:); 

imj = zeros(sz,sz2,sz3);
for i = 1:sz
    imj(i,:,:) = mean(terimj(i:i+armdz*2,:,:),1);
end
clear terimj


meanj = zeros(sz,sz2,sz3);
imj = permute(imj,[2,3,1]);
meanj = calcfibangspeed(imj,armdx,filton);
meanj = permute(meanj,[3,1,2]);
clear imj
meanj = (pi/2-meanj).*(meanj<=pi/2)+(3*pi/2-meanj).*(meanj>pi/2);

% Here starts the determination of gamma

terimg = zeros(sz,sz2+armdz*2,sz3); 
terimg(:,1:armdz,:) = im(:,sz2 - armdz + 1:sz2,:); 
terimg(:,armdz + 1:armdz+sz2,:) = im; 
terimg(:,armdz + sz2 + 1: armdz + sz2 + armdz,:) = im(:,1:armdz,:);

img = zeros(sz,sz2,sz3);
for i = 1:sz2
    img(:,i,:) = mean(terimg(:,(i):(i+armdz*2),:),2); 
end
clear terimg

meang = zeros(sz,sz2,sz3);
img = permute(img,[1,3,2]);
meang = calcfibangspeed(img,armdx,filton);
meang = permute(meang,[1,3,2]);
clear img
meang = (pi/2+meang).*(meang<=pi/2)+(meang-(pi/2)).*(meang>pi/2);

% Here to acquire the orientation of phi
meanp = atan(sqrt(1./(tan(meanj).^2)+1./(tan(meang).^2)));
meanplim = meanp;
meanp = meanp.*(meanj <= pi/2) + (pi-meanp).*(meanj > pi/2);

% Here to acquire the real phi orientation regardless of Z resolution
meanpreal = (pi/2)-(atan(para*tan((pi/2)-meanplim)));
meanpreal = meanpreal.*(meanj <= pi/2) + (pi-meanpreal).*(meanj > pi/2);

% Here to calculate the real beta angle regardless of Z resolution
meanjreal = atan(para*tan(meanj));
meanjreal = meanjreal.*(meanjreal>=0)+(meanjreal+pi).*(meanjreal<0);

% Here to calculate the real gamma angle regardless of Z resolution
meangreal = atan(para*tan(meang));
meangreal = meangreal.*(meangreal>=0)+(meangreal+pi).*(meangreal<0);

% Here prepares the final output
aSDE = meanc*180/pi;
pSDER = meanpreal*180/pi;
bSDER = meanjreal*180/pi;
gSDER = meangreal*180/pi;

end



  
  



