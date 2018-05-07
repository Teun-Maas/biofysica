function h5demo01
    f=h5file('spherecfg.h5');
    info=f.info();
    disp(info);
    listGroup(f);   
end

function listGroup(g)
    listAttr(g);
    listDataset(g);
    gi=h5g_iterator(g);
    while gi.hasnext()
        g=gi.next();
        info=g.info();
        disp('---Group---');
        disp(info);
        listGroup(g);
    end
    
end

function listAttr(obj)
    ai=h5a_iterator(obj);
    while ai.hasnext()
        a=ai.next();
        disp('---Attribute---');
        info=a.info();
        disp(info);
    end
    
end

function listDataset(g)
    di=h5d_iterator(g);
    while di.hasnext()
        d=di.next();
        info=d.info();
        disp('---Dataset---');
        disp(info);
        listAttr(d);
    end
end
