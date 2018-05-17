#!/usr/bin/env python3
import pickle, numpy, scipy.io, os, sys, msgpack

src = sys.argv[1]
dst = sys.argv[2]

with open(src, 'rb') as f:
    pupil_data = msgpack.unpack(f, encoding='utf-8')

scipy.io.savemat(dst, mdict={'pupil_data': pupil_data})
