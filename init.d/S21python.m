function S21python
    global BIOFYSICA_ROOT
    global BIOFYSICA_PY_SYS
    [v,e,loaded] = pyversion;
    if loaded && str2num(v)>=3
        BIOFYSICA_PY_SYS=py.eval('__import__(''sys'')',struct);
        PYPATH=[BIOFYSICA_ROOT '/utilities/python'];
        BIOFYSICA_PY_SYS.path.append(PYPATH);
    end

end
