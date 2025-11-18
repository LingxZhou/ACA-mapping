function [Vmatr, Vvaluelocal, Vvalueentire]= calcfibvarspeed_nf(armdxx, armdzz, aSDE, bSDER, gSDER,finalmask)
% This is a subroutine used to calculate directional variance based on theta & phi orientation.
[sz, sz2, sz3] = size(aSDE);

bwhole = sqrt(1./(tan(2*bSDER*pi/180).^2)+1./(tan(2*gSDER*pi/180).^2));
bwhole(bwhole==Inf) = 0;
Cwhole = (bwhole./sqrt(1+bwhole.^2)).*cos(2*aSDE*pi/180);
Swhole = (bwhole./sqrt(1+bwhole.^2)).*sin(2*aSDE*pi/180);
Zwhole = zeros(sz,sz2,sz3);

bwhole(isnan(bwhole)) = 0;
Cwhole(isnan(Cwhole)) = 0;
Swhole(isnan(Swhole)) = 0;

mask = (bSDER <= 90);
Zwhole(mask) = (1./sqrt(1+(bwhole(mask)).^2));
Zwhole(logical(1-mask)) = (-1./sqrt(1+(bwhole(logical(1-mask))).^2));
Zwhole(isnan(Zwhole)) = 0;

Cwholemean = imboxfilt3(Cwhole,[2*armdxx+1, 2*armdxx+1, 2*armdzz+1],'padding','circular');
Swholemean = imboxfilt3(Swhole,[2*armdxx+1, 2*armdxx+1, 2*armdzz+1],'padding','circular');
Zwholemean = imboxfilt3(Zwhole,[2*armdxx+1, 2*armdxx+1, 2*armdzz+1],'padding','circular');

Rwholemean = sqrt(Cwholemean.^2+Swholemean.^2+Zwholemean.^2);
Vmatr = 1 - Rwholemean;
Vmatr(isnan(Vmatr)) = 0;
Vvaluelocal = mean(Vmatr(finalmask));

Cvalue = mean(Cwhole(finalmask));
Svalue = mean(Swhole(finalmask));
Zvalue = mean(Zwhole(finalmask));
Rvalue = sqrt(Cvalue^2+Svalue^2+Zvalue^2);
Vvalueentire = 1- Rvalue;

end