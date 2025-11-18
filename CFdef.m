function CF = CFdef(D,loc)
% This is a subroutine used to determine the image compression factor(CF)
% under a specific information loss based on the image voxel thickness 
% distribution and the input position 'loc'.

S = sum(ceil(0.25*pi*D.^2));  L = length(D);
D1 = find(D==1); L1 = length(D1); N1 = ceil(0.25*pi*L1);
D2 = find(D==2); L2 = length(D2); N2 = ceil(0.25*pi*4*L2);
D3 = find(D==3); L3 = length(D3); N3 = ceil(0.25*pi*9*L3);
D4 = find(D==4); L4 = length(D4); N4 = ceil(0.25*pi*16*L4);
if L==0
    CF = Nan;
elseif ceil(loc*S)<=N1
    CF = 1;
elseif ceil(loc*S)<=(N1 + N2)
    CF = 2;
elseif ceil(loc*S)<=(N1 + N2 + N3)
    CF = 3;
elseif ceil(loc*S)<=(N1 + N2 + N3 + N4)
    CF = 4;
else
    CF = 5;
end
end