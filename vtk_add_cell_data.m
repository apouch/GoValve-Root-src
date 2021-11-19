function p1 = vtk_add_cell_data(p, name, data, overwrite)
% Add a cell attribute array to a vtk mesh
% Usage:
%   p1 = vtk_add_cell_data(p, name, data, overwrite)
% Parameters
%   p         VTK mesh struct (from vtk_polydata_read)
%   name      Name of the new array (string)
%   data      An Nxk matrix of values to add
%   overwrite boolean, whether to overwrite existing array of the same name

arr.name = name;
arr.type = 'field';

% Count all the cells
cell_types = fieldnames(p.cells);
total_cells = 0;
for i = 1:length(cell_types)
    pc = p.cells.(cell_types{i});
    n_cell = length(pc); 
    total_cells = total_cells + n_cell;
end

if size(data, 1) == total_cells
    arr.data = data;
elseif size(data, 2) == total_cells
    arr.data = data';
else
    error('Data size does not match cell array size');
end

p1 = p;
if ~isfield(p1, 'cell_data')
    p1.cell_data(1) = arr;
else
    pos = strmatch(name, {p1.cell_data.name}, 'exact');
    if ~isempty(pos)
       if overwrite
          p1.cell_data(pos(1)) = arr;
       else
          error('This array already exists');
       end
    else         
       p1.cell_data(length(p1.cell_data)+1) = arr;
    end
end
