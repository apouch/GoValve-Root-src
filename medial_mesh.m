function [] = medial_root_mesh(fnvtk,fnmedout,fnbndout)

% apply rigid transform to orient the outflow tract along the vertical axis
[pts_rotxyz,T_init] = reorient_pts(fnvtk);

% obtain medial and boundary meshes from the reoriented point cloud
[~,mbnd,mmed] = edge_sampling(pts_rotxyz);

% apply inverse rigid transform to medial and boundary points so it's
% aligned with the original segmentation
T_init_inv = inv(T_init);
mmed.points = tform_apply(mmed.points,T_init_inv);
mbnd.points = tform_apply(mbnd.points,T_init_inv);

% write medial and boundary meshes
vtk_polydata_write(fnmedout,mmed);
vtk_polydata_write(fnbndout,mbnd);

end

% -------------------------------------------------------------------------

function [pts_resamp,mbnd,mmed] = edge_sampling(pts_rotxyz)

% number of boundary nodes per slice
nsamp = 48; 

% angles of rotation about z-axis
gamma = 0 : 10 : 350;
nrot = length(gamma);

% indices for 2D medial model
bsh = ((nsamp-2)/2 - 1)/2 + 1;
c1 = 1;
e1 = c1 + bsh;
c2 = e1 + bsh;
e2 = c2 + bsh;
nqs = bsh + 1;

% medial linkages
mi_ref_ord = horzcat([1:nsamp]',[e1:-1:2 1:c2 (c2-1):-1:e1+1]');

% initialize resampled bounday points
pts_resamp = zeros(nsamp,3,nrot);

% rotate point cloud about z-axis
for i = 1 : nrot

    % rotate points around vertical axis
    rg = deg2rad(gamma(i));
    Tz = rotz(rg);
    Tz_inv = rotz(-rg);
    pts_rotz = tform_apply(pts_rotxyz,Tz);
    
    % find a slab of points close to y = 0 along positive x-axis
    slice_ind = find(abs(pts_rotz(:,2)) < 1);
    pts_slab = pts_rotz(slice_ind,:);
    pts_slab_pos = pts_slab(find(pts_slab(:,1) > 0),:);

    % shift vertically so that all coordinates are positive 
    buff_vertical = 10;
    pts_slab_pos(:,3) = pts_slab_pos(:,3) + buff_vertical;
    pts_slab_pos = 10*pts_slab_pos;

    % find x and z extent and add a padding buffer
    buff_pad = 10;
    xrange = [1 ceil(max(pts_slab_pos(:,1)) + buff_pad)];
    zrange = [1 ceil(max(pts_slab_pos(:,3)) + buff_pad)];

    % rasterize the points
    img = zeros(zrange(2),xrange(2));
    ind = sub2ind(size(img),round(pts_slab_pos(:,3)),round(pts_slab_pos(:,1)));
    img(ind) = 1;
    
    % morphologically close, compute boundary contour
    %     figure(1)
    %     set(gca,'ydir','reverse');
    img_cl = imclose(img,strel('disk',50));
    h = fspecial('gaussian',[20 20],7); % [3 3], 0.5
    img_cl_sm = imfilter(img_cl,h);
    %     imagesc(img_cl_sm); axis image; colormap gray; hold on;
    contour(img_cl_sm, [0.5 0.5], 'r'); 

    C = contourc(double(img_cl_sm), [0.5, 0.5]);
    ncp = C(2,1);
    cuv = C(:,2:ncp+1);
    
    % check ordering of points
    if ~ispolycw(cuv(1,:)',cuv(2,:)')
        [cuv1cw,cuv2cw] = poly2cw(cuv(1,:)',cuv(2,:)');
        cuv(1,:) = cuv1cw';
        cuv(2,:) = cuv2cw';
    end

    % find vertical min and max and reorder point indices
    emax = find(cuv(2,:) == max(cuv(2,:)));
    emax = emax(1);
    emin = find(cuv(2,:) == min(cuv(2,:)));
    emin = emin(1);

    cuv_tcs = cuv';
    cuv_tcs(end,:) = [];
    cuv_tcs = circshift(cuv_tcs,ncp-emax,1)';
    cuv_tcs(:,end+1) = cuv_tcs(:,1);
    
    emax_cuv_tcs = 1;
    if emax > emin
        emin_cuv_tcs = emin + (ncp-emax);
    else
        emin_cuv_tcs = emin - emax + 1;
    end

    % determine c1,c2 indices on opposite sides of root near the vertical
    % mean
    mean_pt = mean(cuv_tcs,2);
    d = (cuv_tcs(2,:) - repmat(mean_pt(2),1,ncp));
    dsqrt = sqrt(d.*d);
    c1_cuv_tcs = find(dsqrt(emax_cuv_tcs:emin_cuv_tcs) == min(dsqrt(emax_cuv_tcs:emin_cuv_tcs)));
    c2_cuv_tcs = find(dsqrt(emin_cuv_tcs:end) == min(dsqrt(emin_cuv_tcs:end))) + emin_cuv_tcs - 1;

    % spline c1 -> emin
    x1 = c1_cuv_tcs:emin_cuv_tcs;
    f1 = spline(x1,cuv_tcs(:,x1));
    a1 = linspace(x1(1),x1(end),nqs);
    fap1 = ppval(f1,a1);
    fap1 = fap1';

    % spline emin -> c2
    x2 = emin_cuv_tcs:c2_cuv_tcs;
    f2 = spline(x2,cuv_tcs(:,x2));
    a2 = linspace(x2(1),x2(end),nqs);
    fap2 = ppval(f2,a2);
    fap2 = fap2';

    % spline c2 -> emax
    x3 = c2_cuv_tcs:ncp+1;
    f3 = spline(x3,[cuv_tcs(:,c2_cuv_tcs:ncp),cuv_tcs(:,1)]);
    a3 = linspace(x3(1),x3(end),nqs);
    fap3 = ppval(f3,a3);
    fap3 = fap3';

    % spline emax -> c1
    x4 = emax_cuv_tcs:c1_cuv_tcs;
    f4 = spline(x4,cuv_tcs(:,x4));
    a4 = linspace(x4(1),x4(end),nqs);
    fap4 = ppval(f4,a4);
    fap4 = fap4';

    % resampled 2D contour with all four quadrants
    fap = [fap1; fap2(2:end,:); fap3(2:end,:); fap4(2:end-1,:)];
    
    %     plot(fap(:,1),fap(:,2),'LineWidth',2.5); hold on;
    %     plot(fap(e1,1),fap(e1,2),'ro');
    %     plot(fap(e2,1),fap(e2,2),'mo');
    %     plot(fap(c1,1),fap(c1,2),'ko');
    %     plot(fap(c2,1),fap(c2,2),'go');
    %     hold off;
    
    % resampled 3D contour, rotated back 
    pts_slice = [fap(:,1) zeros(length(fap),1) fap(:,2)];
    pts_zrot_inv = tform_apply(pts_slice,Tz_inv);
    pts_resamp(:,:,i) = pts_zrot_inv/10;
    pts_resamp(:,3,i) = pts_resamp(:,3,i) - buff_vertical;
    
    %     figure(4)
    %     plot3(pts_resamp(:,1,i),pts_resamp(:,2,i),pts_resamp(:,3,i),'r','LineWidth',1.5);
    %     hold on
    
end

% for j = 1 : nsamp
%     pts = squeeze(pts_resamp(j,:,:));
%     pts = pts';
%     plot3(pts(:,1),pts(:,2),pts(:,3),'b','LineWidth',1.5);
% end

% assign an index to each boundary point
pts_resamp(:,4,:) = 0;
for i = 1:nrot
    pts_resamp(:,4,i) = repmat((i-1)*nsamp,nsamp,1) + [1:nsamp]';
end

% initialize boundary mesh data
pts_resamp_list = zeros(nsamp*nrot,4);
medind = zeros(nsamp*nrot,1);
triangles_bnd = [];
bnd_side = [];

% build triangulated mesh from boundary points, slice by slice
for s = 1:nrot
    
    if s == 1
        i = nrot;
    else
        i = s-1;
    end
    j = s;
    
    rstart = (s-1)*nsamp+1;
    pts_resamp_list(rstart:rstart+nsamp-1,:) = pts_resamp(:,:,s);
    medind(rstart:rstart+nsamp-1,1) = mi_ref_ord(:,2)+(s-1)*((nsamp-2)/2+2);
    
    for n = 1 : nsamp
        
        switch n
            case e1
                if mod(e1,2) == 1
                    t = [pts_resamp(n,4,i) pts_resamp(n-1,4,j) pts_resamp(n,4,j) ;
                         pts_resamp(n,4,i) pts_resamp(n,4,j) pts_resamp(n+1,4,j)];
                else
                    t = [pts_resamp(n,4,i) pts_resamp(n-1,4,i) pts_resamp(n,4,j);
                         pts_resamp(n,4,i) pts_resamp(n,4,j) pts_resamp(n+1,4,i)];
                end
                bside = [1;0];
            case c1
                t = [pts_resamp(n,4,i) pts_resamp(nsamp,4,j) pts_resamp(n,4,j) ;
                     pts_resamp(n,4,i) pts_resamp(n,4,j) pts_resamp(n+1,4,j) ];
                bside = [1;1];
            case e2
                if mod(e2,2) == 1
                    t = [pts_resamp(n,4,i) pts_resamp(n,4,j) pts_resamp(n+1,4,j) ;
                         pts_resamp(n,4,i) pts_resamp(n-1,4,j) pts_resamp(n,4,j) ];
                else
                    t = [pts_resamp(n,4,i) pts_resamp(n,4,j) pts_resamp(n+1,4,i) ;
                         pts_resamp(n,4,i) pts_resamp(n-1,4,i) pts_resamp(n,4,j) ];
                end
                bside = [1;0];
            case nsamp
                t = [pts_resamp(n,4,i) pts_resamp(n-1,4,i) pts_resamp(n,4,j);
                     pts_resamp(n,4,i) pts_resamp(n,4,j) pts_resamp(1,4,i)];
                bside = [1;1];
            otherwise
                if mod(n,2) == 1
                    t = [pts_resamp(n,4,i) pts_resamp(n-1,4,j) pts_resamp(n,4,j);
                         pts_resamp(n,4,i) pts_resamp(n,4,j) pts_resamp(n+1,4,j)];
                else
                    t = [pts_resamp(n,4,i) pts_resamp(n-1,4,i) pts_resamp(n,4,j);
                         pts_resamp(n,4,i) pts_resamp(n,4,j) pts_resamp(n+1,4,i)];
                end  

                if n < c2 
                    if n < e1
                        bside = [1;1];
                    else 
                        bside = [0;0];
                    end
                else
                    if n < e2
                        bside = [0;0];
                    else
                        bside = [1;1];
                    end
                end
                    
        end    
        
        triangles_bnd = [triangles_bnd; t];
        bnd_side = [bnd_side; bside];
    end
    
end

pts_resamp_list = pts_resamp_list(:,1:3);
% trimesh(double(triangles_bnd),pts_resamp_list(:,1),pts_resamp_list(:,2),pts_resamp_list(:,3));

% create vtk polydata for boundary mesh
mbnd = struct;
mbnd.hdr.name = 'File written by itkPolyDataMeshIO';
mbnd.hdr.type = 'ASCII';
mbnd.hdr.dst = 'POLYDATA';
mbnd.points = pts_resamp_list;
for i = 1 : length(triangles_bnd)
    mbnd.cells.polygons{i} = triangles_bnd(i,:)';
end
mbnd = vtk_add_point_data(mbnd,'MedialIndex',medind-1);
mbnd = vtk_add_cell_data(mbnd,'Label',ones(length(triangles_bnd),1));
mbnd = vtk_add_cell_data(mbnd,'Side',bnd_side);

% approximate the medial points, thickness, and other data that will be
% stored in the medial and boundary meshes
pts_cent = mean(pts_resamp_list,1);
pts_med = 0*pts_resamp_list;
pts_side = zeros(length(pts_med),1);
stj_flag = zeros(length(pts_med),1);
lvo_flag = zeros(length(pts_med),1);
med_thickness = zeros(length(pts_med),1);
bnd_thickness = zeros(length(pts_resamp_list),1);
for i = 1 : length(medind)
    xmi = medind(i);
    xb = find(medind == xmi);
    if length(xb) > 1
        
        mp = (pts_resamp_list(xb(1),:)+pts_resamp_list(xb(2),:))/2;
        pts_med(xb(1),:) = mp;
        pts_med(xb(2),:) = mp;
        
        [prow1,~] = find(triangles_bnd == xb(1));
        [prow2,~] = find(triangles_bnd == xb(2));
        pts_side(xb(1)) = bnd_side(prow1(1));
        pts_side(xb(2)) = bnd_side(prow2(1));
        
        med_thickness(xb(1)) = pdist([pts_resamp_list(xb(1),:);
                                    pts_resamp_list(xb(2),:)]);
        med_thickness(xb(2)) = med_thickness(xb(1));
        
        bnd_thickness(xb(1)) = med_thickness(xb(1));
        bnd_thickness(xb(2)) = med_thickness(xb(1));
        
    else
        
        pts_med(xb,:) = pts_resamp_list(xb,:);
        pts_side(xb) = 3;
        med_thickness(xb) = 0.3;
        bnd_thickness(xb) = 0.3;
        
        if pts_med(xb,3) > pts_cent(3)
            stj_flag(xb) = 1;
        else 
            lvo_flag(xb) = 1;
        end
        
    end
end

mbnd = vtk_add_point_data(mbnd,'Thickness',bnd_thickness);
mbnd = vtk_add_point_data(mbnd,'Radius',bnd_thickness/2);

% obtain triangulation of medial mesh
pts_side0_ind = find(pts_side == 1);
pts_side1_ind = find(pts_side ~= 1);
triangles_med = triangles_bnd;
for i = 1 : length(triangles_med)
    for j = 1 : 3
        
        node = triangles_med(i,j);
        if sum(ismember(pts_side0_ind,node))
            triangles_med(i,:) = 0;
            break;
        else
            triangles_med(i,j) = find(pts_side1_ind == node);
        end
  
    end
end
triangles_med(sum(triangles_med,2) == 0,:) = [];
pts_med = pts_med(pts_side1_ind,:);

% create vtk polydata for the medial mesh
mmed = struct;
mmed.hdr.name = 'File written by itkPolyDataMeshIO';
mmed.hdr.type = 'ASCII';
mmed.hdr.dst = 'POLYDATA';
mmed.points = pts_med;
for i = 1 : length(triangles_med)
    mmed.cells.polygons{i} = triangles_med(i,:)';
end
mmed = vtk_add_point_data(mmed,'Label',ones(length(pts_med),4));
mmed = vtk_add_point_data(mmed,'Thickness',med_thickness(pts_side1_ind));
mmed = vtk_add_point_data(mmed,'Radius',med_thickness(pts_side1_ind)/2);
mmed = vtk_add_point_data(mmed,'STJ',stj_flag(pts_side1_ind));
mmed = vtk_add_point_data(mmed,'LVO',lvo_flag(pts_side1_ind));

end

% -------------------------------------------------------------------------

function [pts_rotxyz,T] = reorient_pts(fnvtk)

% read medial model and point labels
m = vtk_polydata_read(fnvtk);
labels = m.point_data(1).data;

% assumes: STJ is 1, LVO is 2, IAS is 6
stj_ind = find(labels == 1);
lvo_ind = find(labels == 2);
septum_ind = find(labels == 6);

% outflow tract vector
outflow_vec = [mean(m.points(lvo_ind,:)); mean(m.points(stj_ind,:))];

% translate point clouds to origin
Tt = trans(-outflow_vec(1,:));
pts_origin = tform_apply(m.points,Tt);
outflow_vec_trans = tform_apply(outflow_vec,Tt);

% rotate eigenvector associated with lowest variance to align it with the
% positive z-axis

[Tx,Ty,~,~] = rotxy_apply(outflow_vec_trans(2,:)');

% now rotate points
pts_rotxy = tform_apply(tform_apply(pts_origin,Tx),Ty);

% next rotate the model around the z-axis
m_new = mean(pts_rotxy(septum_ind,:));
f_new = [100 0 0];

x0 = 0; 
f = @(x)rotz_apply(x,f_new,m_new);
x = fminunc(f,x0);

gamma = x;
Tz = rotz(gamma);

% finally rotate the template
pts_rotxyz = tform_apply(pts_rotxy,Tz);

% composite transform
T = Tz * Ty * Tx * Tt;

end

% -------------------------------------------------------------------------

function pts_new = tform_apply(pts,T)

pts = [pts'; ones(1,size(pts,1))];
pts_new = T*pts;
pts_new = pts_new(1:3,:)';

end

% -------------------------------------------------------------------------

function Tx = rotx(alpha)

Tx = [1     0           0           0;
      0     cos(alpha)  -sin(alpha) 0;
      0     sin(alpha)  cos(alpha)  0;
      0     0           0           1];

end

% -------------------------------------------------------------------------

function Ty = roty(beta)

Ty = [cos(beta)    0    sin(beta)   0;
      0            1    0           0;
      -sin(beta)   0    cos(beta)   0;
      0            0    0           1];
  
end

% -------------------------------------------------------------------------

function Tz = rotz(gamma)

Tz = [cos(gamma)    -sin(gamma) 0   0;
      sin(gamma)    cos(gamma)  0   0;
      0             0           1   0;
      0             0           0   1];

end

% -------------------------------------------------------------------------

function Tt = trans(x)

Tt = [1 0   0   x(1);
      0 1   0   x(2);
      0 0   1   x(3);
      0 0   0   1 ];
end

% -------------------------------------------------------------------------

function [Tx,Ty,alpha,beta] = rotxy_apply(v)

% Angle of rotation about x-axis
alpha = acos(v(3)/sqrt(v(2)^2 + v(3)^2));
if v(2) < 0
    alpha = -alpha;
end

Tx = rotx(alpha);
 
v_rx = Tx*[v; 1];

% Angle of rotation about y-axis
beta = acos(v_rx(3)/sqrt(v_rx(1)^2 + v_rx(3)^2));
if v_rx(1) > 0
    beta = -beta;
end

Ty = roty(beta);
  
v_rxy = Ty * v_rx;

% Check that rotation is correct
if (abs(v_rxy(1)) < 1E-6) && (abs(v_rxy(2)) < 1E-6) && ((abs(v_rxy(3))-1) < 1E-6)
    disp('Third eigenvector is correctly aligned');
else
    disp('Check third eigenvector');
end

end

% -------------------------------------------------------------------------
function D = rotz_apply(x,f,m)

% Rotation angle and transform
gamma = x;
Tz = rotz(gamma);

% Transformed points
m_new = tform_apply(m,Tz);

D = (f - m_new).*(f - m_new);
D = sum(D(:));

end
