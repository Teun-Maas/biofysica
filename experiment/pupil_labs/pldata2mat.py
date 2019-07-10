#!/usr/bin/env python3
'''
(*)~---------------------------------------------------------------------------
Copyright (c) 2019 Gunter Windau

Code fragments from Pupil Labs .../pupil_src/shared_modules/file_methods.py
Copyright (C) 2012-2018 Pupil Labs
Distributed under the terms of the GNU
Lesser General Public License (LGPL v3.0).
See COPYING and COPYING.LESSER for license details.
---------------------------------------------------------------------------~(*)
'''
import collections
import os
import sys
import scipy.io
import msgpack
import numpy as np

PLData = collections.namedtuple('PLData', ['data', 'timestamps', 'topics'])

class _Empty(object):
    def purge_cache(self):
        pass

# an Immutable dict for dics nested inside this dict.
class _FrozenDict(dict):
    def __setitem__(self, key, value):
        raise NotImplementedError('Invalid operation')

    def clear(self):
        raise NotImplementedError()

    def update(self, *args, **kwargs):
        raise NotImplementedError()


class Serialized_Dict(object):
    __slots__ = ['_ser_data', '_data']
    cache_len = 100
    _cache_ref = [_Empty()] * cache_len
    MSGPACK_EXT_CODE = 13

    def __init__(self, python_dict=None, msgpack_bytes=None):
        if type(python_dict) is dict:
            self._ser_data = msgpack.packb(python_dict, use_bin_type=True,
                                           default=self.packing_hook)
        elif type(msgpack_bytes) is bytes:
            self._ser_data = msgpack_bytes
        else:
            raise ValueError("Neither mapping nor payload is supplied or wrong format.")
        self._data = None

    def _deser(self):
        if not self._data:
            self._data = msgpack.unpackb(self._ser_data, raw=False, use_list=False,
                                         object_hook=self.unpacking_object_hook,
                                         ext_hook=self.unpacking_ext_hook)
            self._cache_ref.pop(0).purge_cache()
            self._cache_ref.append(self)

    @classmethod
    def unpacking_object_hook(self,obj):
        if type(obj) is dict:
            return _FrozenDict(obj)

    @classmethod
    def packing_hook(self, obj):
        if isinstance(obj, self):
            return msgpack.ExtType(self.MSGPACK_EXT_CODE, obj.serialized)
        raise TypeError("can't serialize {}({})".format(type(obj), repr(obj)))

    @classmethod
    def unpacking_ext_hook(self, code, data):
        if code == self.MSGPACK_EXT_CODE:
            return self(msgpack_bytes=data)
        return msgpack.ExtType(code, data)

    def purge_cache(self):
        self._data = None

    @property
    def serialized(self):
        return self._ser_data

    def __setitem__(self, key, item):
        raise NotImplementedError()

    def __getitem__(self, key):
        self._deser()
        return self._data[key]

    def __repr__(self):
        self._deser()
        return 'Serialized_Dict({})'.format(repr(self._data))

    @property
    def len(self):
        '''Replacement implementation for __len__
        If __len__ is defined numpy will recognize this as nested structure and
        start deserializing everything instead of using this object as it is.
        '''
        self._deser()
        return len(self._data)

    def __delitem__(self, key):
        raise NotImplementedError()

    def get(self,key,default):
        try:
            return self[key]
        except KeyError:
            return default

    def clear(self):
        raise NotImplementedError()

    def copy(self):
        self._deser()
        return self._data.copy()

    def has_key(self, k):
        self._deser()
        return k in self._data

    def update(self, *args, **kwargs):
        raise NotImplementedError()

    def keys(self):
        self._deser()
        return self._data.keys()

    def values(self):
        self._deser()
        return self._data.values()

    def items(self):
        self._deser()
        return self._data.items()

    def pop(self, *args):
        raise NotImplementedError()

    def __cmp__(self, dict_):
        self._deser()
        return self._data.__cmp__(dict_)

    def __contains__(self, item):
        self._deser()
        return item in self._data

    def __iter__(self):
        self._deser()
        return iter(self._data)



def load_pldata_file(directory, topic):
    ts_file = os.path.join(directory, topic + '_timestamps.npy')
    msgpack_file = os.path.join(directory, topic + '.pldata')
    try:
        data = collections.deque()
        topics = collections.deque()
        data_ts = np.load(ts_file)
        with open(msgpack_file, "rb") as fh:
            for topic, payload in msgpack.Unpacker(fh, raw=False, use_list=False):
                data.append(Serialized_Dict(msgpack_bytes=payload))
                topics.append(topic)
    except FileNotFoundError:
        data = []
        data_ts = []
        topics = []

    return PLData(data, data_ts, topics)


datadir = sys.argv[1]
topic = sys.argv[2]
matfile = sys.argv[3]

topic_data=load_pldata_file(datadir, topic)
scipy.io.savemat(matfile, mdict={topic : topic_data})
