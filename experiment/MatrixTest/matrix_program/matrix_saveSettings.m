function matrix_saveSettings(h, filename)    
    s = matrix_getSettings(h);
    save(filename, '-struct', 's');
end