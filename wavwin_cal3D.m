function [wavmatr1] = wavwin_cal3D(xwin,ywin,zwin,orimatr,finalmask)
% This is a subroutine used to calculate waviness based on theta & phi orientation.
[sz1,sz2,sz3] = size(orimatr);
wavmatr = zeros(sz1,sz2,sz3,'double','gpuArray');
oripad = zeros(sz1,sz2,sz3+2*zwin,'double','gpuArray');
fmpad =  zeros(sz1,sz2,sz3+2*zwin,'double','gpuArray');
oridiff = zeros(sz1+2*xwin, sz2+2*ywin, sz3+2*zwin,'double','gpuArray');
fmsum = zeros(sz1+2*xwin, sz2+2*ywin, sz3+2*zwin,'double','gpuArray');

oripadshift = zeros(1+2*xwin, 1+2*ywin, 1+2*zwin,'double','gpuArray');
fmshift = zeros(1+2*xwin, 1+2*ywin, 1+2*zwin,'double','gpuArray');
oridiff_temp = zeros(1+2*xwin, 1+2*ywin, 1+2*zwin,'double','gpuArray');

for k = 1:zwin
    oripad(:,:,k) = orimatr(:,:,1);
    oripad(:,:,zwin+sz3+k) = orimatr(:,:,sz3);
    fmpad(:,:,k) = finalmask(:,:,1);
    fmpad(:,:,zwin+sz3+k) = finalmask(:,:,sz3);
end
oripad(:,:,zwin+1 : zwin+sz3) = orimatr;
fmpad(:,:,zwin+1 : zwin+sz3) = finalmask;
oripad = padarray(oripad,[xwin,ywin],'replicate');
fmpad = padarray(fmpad,[xwin,ywin],'replicate');

for k = -zwin:zwin
    for i = -xwin:xwin
        for j = -ywin:ywin
            if i~=0 || j~=0 || k~=0
                oripadshift = circshift(oripad, [i,j,k]);
                fmshift = circshift(fmpad, [i,j,k]);
                oridiff_temp = abs(oripad - oripadshift);
                oridiff_temp = (180-oridiff_temp).*(oridiff_temp >= 90) + (oridiff_temp < 90).*oridiff_temp;
                oridiff = oridiff + oridiff_temp .* fmshift;
                fmsum = fmsum + fmshift;
            end
        end
    end
end
clear oripadshift; clear fmshift; clear oridiff_temp; clear oripad; clear fmpad;

wavmatr = oridiff(xwin+1 : sz1+xwin, ywin+1 : sz2+ywin, zwin+1 : sz3+zwin);
fmsum = fmsum(xwin+1 : sz1+xwin, ywin+1 : sz2+ywin, zwin+1 : sz3+zwin);
wavmatr = (wavmatr / 90) ./ fmsum;
wavmatr(isnan(wavmatr)) = 0;
wavmatr1 = gather(wavmatr);
clear wavmatr; clear oridiff; clear fmsum;
end