import io
import urllib.request
import h5py
from collections import OrderedDict

URL = 'https://ftp.cngb.org/pub/stomics/STT0000036/supp/SpatialTranscriptome.h5ad'

class RemoteFile(io.RawIOBase):
    def __init__(self):
        super().__init__()
        self.pos = 0
        self.cache = OrderedDict()
        self.bytes_read = 0
        with urllib.request.urlopen(urllib.request.Request(URL, headers={'Range': 'bytes=0-0'}), timeout=25) as r:
            assert r.status == 206
            self.size = int(r.headers['Content-Range'].split('/')[-1])
    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos
    def seek(self, offset, whence=0):
        self.pos = offset if whence == 0 else self.pos + offset if whence == 1 else self.size + offset
        return self.pos
    def read(self, n=-1):
        if n < 0: n = self.size - self.pos
        end = min(self.pos + n, self.size)
        blocks = []
        while self.pos < end:
            start = self.pos // 65536 * 65536
            if start not in self.cache:
                stop = min(start + 65536, self.size) - 1
                with urllib.request.urlopen(urllib.request.Request(URL, headers={'Range': f'bytes={start}-{stop}'}), timeout=40) as r:
                    assert r.status == 206
                    assert r.headers['Content-Range'].startswith(f'bytes {start}-')
                    b = r.read()
                assert len(b) == stop-start+1
                self.cache[start] = b
                self.bytes_read += len(b)
                if len(self.cache) > 512: self.cache.popitem(last=False)
            b = self.cache[start]
            take = min(end-self.pos, len(b)-(self.pos-start))
            blocks.append(b[self.pos-start:self.pos-start+take])
            self.pos += take
        return b''.join(blocks)
    def readinto(self, b):
        d = self.read(len(b))
        b[:len(d)] = d
        return len(d)

def walk(group, prefix, depth=0):
    for key in group:
        obj = group[key]
        name = prefix + '/' + key
        print(name, 'GROUP' if isinstance(obj,h5py.Group) else str(obj.shape)+' '+str(obj.dtype), flush=True)
        if isinstance(obj,h5py.Group) and depth < 5: walk(obj,name,depth+1)

if __name__ == '__main__':
    remote = RemoteFile()
    print('SIZE',remote.size,flush=True)
    with h5py.File(remote,'r') as h:
        print('ROOT',list(h),flush=True)
        for k in ['uns','obsm','obs']: print(k,list(h[k]),flush=True)
        walk(h['uns'],'uns')
    print('RANGE_BYTES',remote.bytes_read,flush=True)
