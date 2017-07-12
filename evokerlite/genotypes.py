from numpy import array
from struct import unpack
from math import ceil

class Genotypes(object):
    
    def __init__(self, file_path, n_samples, n_variants):
        self.file_path = file_path
        self.n_samples = n_samples
        self.n_variants = n_variants
        self.HEADER_MAGIC = (108, 27, 1)
        self.check_file()
        
        
    def check_file(self):
        with open(self.file_path, 'rb') as f:
            _bytes = f.read(3)
            magic = unpack('3B', _bytes)
            assert magic == self.HEADER_MAGIC

    def get_offset(self, variant_index):
        offset = len(self.HEADER_MAGIC)
        offset += ceil(self.n_samples / 4.) * variant_index
        return int(offset)
    
    def get_genotypes(self, variant_index):
        offset = self.get_offset(variant_index)
        n_bytes = int(ceil(self.n_samples / 4.))
        with open(self.file_path, 'rb') as f:
            f.seek(offset)
            _bytes = f.read(n_bytes)
        _bytes = unpack('B'*n_bytes, _bytes)
        genotypes = []
        for byte in _bytes:
            byte = '{:08b}'.format(byte)
            genotypes.extend([byte[i:i+2] for i in range(6, -2, -2)])
        return array(genotypes[:self.n_samples])
