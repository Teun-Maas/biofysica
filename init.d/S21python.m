function S21python
    global BIOFYSICA_PY_SYS
    global biofpy
    [v,e,loaded] = pyversion;
    if loaded && str2num(v)>=3
        BIOFYSICA_PY_SYS=py.eval('__import__(''sys'')',struct);
        PYPATH=[biofysica_root '/utilities/python'];
        BIOFYSICA_PY_SYS.path.append(PYPATH);
        biofpy=py.eval('__import__(''biofpy'')',struct);
    else
        warning('biofysica toolbox needs python version >= 3, you have %s',v);
    end

end
