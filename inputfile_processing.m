% ======================================================================
%> @file inputfile_processing
%> @brief Parses an input file containing data re
%> @details
%>
%> @author Mafalda Gonçalves
%> @date April 2023
% ======================================================================

function [input_data] = inputfile_processing(inputfile_name)

% Open the input file for reading
fid = fopen(inputfile_name, 'r');

% Initialize the structure
input_data = struct();

% Loop through each line in the file
while ~feof(fid)
    % Read the current line
    line = fgetl(fid);
    
    % Skip empty lines or comments (lines starting with '#')
    if isempty(line) || line(1) == '#'
        continue
    elseif line(1) == '*'
        tag = line;
        continue
    else
        key_value = strsplit(line, ',');
    end
    
    % Store the key and value pairs in the structure
    if strcmp(tag, '**Domain')
        input_data.domainy = str2double(key_value{1});
        input_data.domainx = str2double(key_value{2});
    elseif strcmp(tag, '**Material')
        input_data.E0 = str2double(key_value{1});
        input_data.poisson = str2double(key_value{2});
    elseif strcmp(tag, '**Penal')
        input_data.penal = str2double(key_value{1});
    elseif strcmp(tag, '**Springs')
        input_data.spring_x = str2double(key_value{1});
        input_data.spring_y = str2double(key_value{2});
    elseif strcmp(tag, '**Vol')
        input_data.vol_frac = str2double(key_value{1});
    elseif strcmp(tag, '**Constraints')
        input_data.constraints = key_value{1};
    elseif strcmp(tag, '**Radius')
        input_data.radius = str2double(key_value{1});
    elseif strcmp(tag, '**BCs')
        input_data.ndis = str2double(key_value{1});
    elseif strcmp(tag, '**Fin')
        input_data.fin = str2double(key_value{1});
    elseif strcmp(tag, '**Analysis_type')
        input_data.analysis = key_value{1};
    elseif strcmp(tag, '**Material_behavior')
        input_data.matbehav = key_value{1};
    elseif strcmp(tag, '**Integration_type')
        input_data.integration = key_value{1};
    elseif strcmp(tag, '**Filename')
        input_data.filename = key_value{1};
    elseif strcmp(tag, '**Domain Type')
        input_data.dtype = key_value{1};
    end
end

% Close the input file
fclose(fid);

end

