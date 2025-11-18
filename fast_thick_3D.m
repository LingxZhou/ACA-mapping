function [thickness,CF] = fast_thick_3D(Im,ske_layer,loc)
[sz1,sz2,sz3] = size(Im);

N = ceil(sz3/ske_layer);


%%--Calculation of edge voxel set and skeleton voxel set--
edg = zeros(sz1,sz2,sz3); ske = zeros(sz1,sz2,sz3);
% Initialize edge voxel set (edg) and skeleton voxel set (ske)
lp=fspecial('laplacian');
for k = 1:N-1
    I_t = 0;
    % Initialize the intensity overlay image
    for i = (k-1)*ske_layer+1:k*ske_layer
        edg(:,:,i) = imfilter(Im(:,:,i),lp);  %modify 1 
        % Choose an appropriate method to calculate the edge point set 
        % based on the data quality. For images with low signal-to-noise ratio, 
        % modify 1 can be changed to edg(:,:,i) = edge(Im(:,:,1),'canny').
        I_t = Im(:,:,i)+I_t;
    end
    I_t = bwmorph(I_t,'thin','inf');
    for i = (k-1)*ske_layer+1:k*ske_layer
        ske(:,:,i) = Im(:,:,i).*double(I_t);
    end
end

I_t = 0;
for i = (N-1)*ske_layer+1:sz3
    edg(:,:,i) = imfilter(Im(:,:,i),lp);  % modify 2
    % Choose an appropriate method to calculate the edge point set 
    % based on the data quality. Keep same with modify 1.
    I_t = Im(:,:,i)+I_t;
end
I_t = bwmorph(I_t,'thin','inf');
for i = (N-1)*ske_layer+1:sz3
    ske(:,:,i) = Im(:,:,i).*double(I_t);
end

M = 5; %modify 3
ske = imbinarize(ske); edg = imbinarize(edg);
ske = bwskel(ske, 'MinBranchLength', M);
% Remove erroneous skeletonization at intersections by pruning.
% 'M' is adjusted according to the thickness of the vessels.

%%--thickness calculation--

S = find(ske==1); b = ceil(S/sz1/sz2); n = ceil((S-(b-1)*sz1*sz2)/sz1); m = S-(b-1)*sz1*sz2-(n-1)*sz1; SS = [m,n,b];
E = find(edg==1); h = ceil(E/sz1/sz2); l = ceil((E-(h-1)*sz1*sz2)/sz1); g = E-(h-1)*sz1*sz2-(l-1)*sz1; ES = [g,l,h];
[nearest_indices,D] = knnsearch(ES,SS,'K',2);
% Determine the two closest voxels from the edge voxel set to each voxel 
% in the skeleton voxel set.

ske = double(ske);
for i = 1:length(m)
    ske(m(i),n(i),b(i))=D(i,1)+D(i,2)-1;
    % Here to calculate the thickness for each voxel in the skeleton voxel set.
end

%%--modification in axial direction of skeleton voxel set--
thickness = zeros(sz1,sz2,sz3);
for i = 1:sz1
    for j = 1:sz2
        n = 1;    
        for k = 1:sz3           
            if ske(i,j,k)==0||k<n
                continue
            end
            if k < sz3
                for n = k+1:sz3
                    if ske(i,j,n)==0
                        break;
                    end
                end
                a = ske(i,j,k:n-1);
                ind = find(a==max(a));
                thickness(i,j,ind+k-1) = ske(i,j,ind+k-1);
            else
                if ske(i,j,k-1) >= ske(i,j,k)
                    thickness (i,j,k) = 0;
                else
                    thickness(i,j,k) = ske(i,j,k); 
                end
            end
        end
    end
end
thickness = ceil(thickness);

%%--determination of compression factor--
thickness_d = reshape(thickness,1,[]);
thickness_d = sort(thickness_d);
NZ = find(thickness_d~=0);
D = thickness_d(min(NZ):length(thickness_d));
CF = CFdef(D,loc);

end