function copy_data_files_from_dir(target_path, file_paths)
    mkdir(target_path);
    if isfolder(file_paths)
        filedirs = caraslab_lsdir(file_paths);
        filedirs = {filedirs.name};
    else
        filedirs = dir(file_paths);
    end
    for file_idx=1:length(filedirs)
        if ~isempty(filedirs)
            if isfolder(file_paths)
                cur_filedir = fullfile(file_paths, filedirs{file_idx});
                copyfile(cur_filedir, fullfile(target_path, filedirs{file_idx}));
            else
                cur_filedir = filedirs(file_idx);
                copyfile(fullfile(cur_filedir.folder, cur_filedir.name), fullfile(target_path, cur_filedir.name));
            end
        end
    end
end