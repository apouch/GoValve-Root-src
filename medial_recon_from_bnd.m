function [] = med_thickness2mesh(outdir,tag,fnbnd_ref,fnmed_ref,framenums,fref)

mbnd_ref = vtk_polydata_read(fnbnd_ref);
bnd_medind = vtk_get_point_data(mbnd_ref,'MedialIndex')+1;

mmed_ref = vtk_polydata_read(fnmed_ref);
med_label = vtk_get_point_data(mmed_ref,'Label');
med_stj = vtk_get_point_data(mmed_ref,'STJ');
med_vaj = vtk_get_point_data(mmed_ref,'VAJ');
nmpts = size(mmed_ref.points,1);

for k = 1 : length(framenums)
    
    fnum = framenums(k);
    
    if fnum ~= fref
    
        fnbnd = [outdir '/' tag '_bnd-prop_f' int2str(fnum) '.vtk'];
        mbnd = vtk_polydata_read(fnbnd);

        fnmed_out = [outdir '/' tag '_med-prop-recon_f' int2str(fnum) '.vtk'];

        pts_med = zeros(nmpts,3);
        med_thickness = zeros(nmpts,1);
        for i = 1 : nmpts
            xb = find(bnd_medind == i);
            if length(xb) > 1
                mp = (mbnd.points(xb(1),:) + mbnd.points(xb(2),:))/2;
                pts_med(i,:) = mp;
                med_thickness(i) = pdist([mbnd.points(xb(1),:);
                                                    mbnd.points(xb(2),:)]);
            else
                pts_med(i,:) = mbnd.points(xb,:);
                med_thickness(i) = 0.3;
            end
    %         display(['Med Pt ' int2str(i) ': Med Ind ' int2str(xb')]);
        end

        mmed = struct;
        mmed.hdr.name = 'File written by itkPolyDataMeshIO';
        mmed.hdr.type = 'ASCII';
        mmed.hdr.dst = 'POLYDATA';
        mmed.points = pts_med;
        mmed.cells = mmed_ref.cells;
        mmed = vtk_add_point_data(mmed,'Label',med_label);
        mmed = vtk_add_point_data(mmed,'Thickness',med_thickness);
        mmed = vtk_add_point_data(mmed,'Radius',med_thickness/2);
        mmed = vtk_add_point_data(mmed,'STJ',med_stj);
        mmed = vtk_add_point_data(mmed,'VAJ',med_vaj);
        vtk_polydata_write(fnmed_out,mmed);
    end

end
