function S21python
    global BIOFYSICA_PY_SYS
    global BIOFYSICA_PY_EXE
    global biofpy
    [v,e,loaded] = pyversion;
    if loaded && str2num(v)>=3
        BIOFYSICA_PY_EXE=e;
        BIOFYSICA_PY_SYS=py.eval('__import__(''sys'')',struct);
        PYPATH=[biofysica_root '/utilities/python'];
        BIOFYSICA_PY_SYS.path.append(PYPATH);
        biofpy=py.eval('__import__(''biofpy'')',struct);
    else
        if loaded
            warning('biofysica toolbox needs python version >= 3, you have %s',v);
        else
            warning('could not load your python %s interpreter ''%s''',v,e);
        end
    end

end
