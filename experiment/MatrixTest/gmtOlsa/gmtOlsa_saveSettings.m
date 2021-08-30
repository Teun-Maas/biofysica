function gmtOlsa_saveSettings(h, filename)    
    s = gmtOlsa_getSettings(h);
    save(filename, '-struct', 's');
end