function [LCmatr3D] = LCcalc3D_SorF(ywin,xwin,zwin,finalmask,winnum,win_th)
% This is a subroutine used for adaptive LC calculation. 
% Based on the relationship between the compressed image size
% 'ywin', 'xwin', 'zwin' and the threshold window size 'win_th',
% it adaptively selects  the faster LC calculation method (spatial 
% domain calculation or frequency domain calculation).

[sz1,sz2,sz3] = size(finalmask);
LCmatr3D = zeros(sz1,sz2,sz3);
fmpad = zeros(sz1,sz2,sz3+2*zwin);
for k = 1:zwin
    fmpad(:,:,k) = finalmask(:,:,1);
    fmpad(:,:,zwin+sz3+k) = finalmask(:,:,sz3);
end
fmpad(:,:,zwin+1 : zwin+sz3) = finalmask;
fmpad = padarray(fmpad,[ywin,xwin],'replicate');
F=ones(2*ywin+1,2*xwin+1,2*zwin+1)/winnum;

if ywin*xwin*zwin < win_th^3
    LC = convn(fmpad,F,'same');
    LCmatr3D = LC(ywin+1:ywin+sz1,xwin+1:xwin+sz2,zwin+1:zwin+sz3);
else
    FF = fftn(F,size(fmpad));
    fmpad =circshift(fmpad,[-xwin,-ywin,-zwin]);
    fmpadF = fftn(fmpad);
    LCF = fmpadF.*FF;
    LC = real(ifftn(LCF));
    LCmatr3D = LC(ywin+1:ywin+sz1,xwin+1:xwin+sz2,zwin+1:zwin+sz3);
end
end