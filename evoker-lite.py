%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
from math import ceil
from struct import unpack
import re

class Variants(object):
    
    def __init__(self, file_path):
        self.file_path = file_path
        self.load_variants()

    def load_variants(self):
        variants = []
        with open(self.file_path) as f:
            for line in f:
                tokens = line.split()
                chrom, name, pos, coord, A1, A2 = tokens
                variant = Variant(chrom, name, pos, coord, A1, A2)
                variants.append(variant)
        self.variants = variants
        
    def get_variants(self):
        return self.variants
        
    def get_n_variants(self):
        return len(self.variants)
    
    def get_index(self, variant_name):
        for index, variant in enumerate(self.variants):
            if variant.name == variant_name:
                return index
    
    def get_A1(self, variant_index):
        return self.variants[variant_index].get_A1()
    
    def get_A2(self, variant_index):
            return self.variants[variant_index].get_A2()
    
    def get_name(self, variant_index):
        return self.variants[variant_index].get_name()
    
    def __str__(self):
        return '\n'.join([str(variant) for variant in self.variants])


class Variant(object):
    
    def __init__(self, chrom, name, pos, coord, A1, A2):
        self.chrom = chrom
        self.name = name
        self.pos = pos
        self.coord = coord
        self.A1 = A1
        self.A2 = A2
        
    def get_A1(self):
        return self.A1
    
    def get_A2(self):
        return self.A2
    
    def get_name(self):
        return self.name
    
    def __str__(self):
        return '{}:{} {}'.format(self.chrom, self.pos, self.name)

class Samples(object):
    
    def __init__(self, file_path):
        self.file_path = file_path
        self.load_samples()

    def load_samples(self):
        samples = []
        with open(self.file_path) as f:
            for line in f:
                tokens = line.split()
                sample = tokens[0]
                samples.append(sample)
        self.samples = samples
        
    def get_samples(self):
        return self.samples
        
    def get_n_samples(self):
        return len(self.samples)


class Genotypes(object):
    
    def __init__(self, file_path, n_samples, n_variants):
        self.file_path = file_path
        self.n_samples = n_samples
        self.n_variants = n_variants
        self.HEADER_MAGIC = (108, 27, 1)
        self.check_file()
        
        
    def check_file(self):
        with open(self.file_path) as f:
            bytes = f.read(3)
            magic = unpack('3B', bytes)
            assert magic == self.HEADER_MAGIC

    def get_offset(self, variant_index):
        offset = len(self.HEADER_MAGIC)
        offset += ceil(self.n_samples / 4.) * variant_index
        return int(offset)
    
    def get_genotypes(self, variant_index):
        offset = self.get_offset(variant_index)
        n_bytes = int(ceil(self.n_samples / 4.))
        with open(self.file_path) as f:
            f.seek(offset)
            bytes = f.read(n_bytes)
        bytes = unpack('B'*n_bytes, bytes)
        genotypes = []
        for byte in bytes:
            byte = '{:08b}'.format(byte)
            genotypes.extend([byte[i:i+2] for i in range(6, -2, -2)])
        return np.array(genotypes[:self.n_samples])

    
class TextIntensity(object):
    
    def __init__(self, file_path, variants, samples):
        self.file_path = file_path
        self.variants = variants
        self.samples = samples
        self.map_sample_ids()
        
    def map_sample_ids(self):
        with open(self.file_path) as f:
            header = f.readline()
        tokens = header.split()
        header_samples = tokens[3::2]
        header_samples = [x[:-1] for x in header_samples]

        # Check if these variants exist in fam file
#         s = set(header_samples)
#         t = set(self.samples.get_samples())
#         intersection = s & t
#         if len(intersection) < min(len(s), len(t)):
        mapping = []
        for fam_sample in self.samples.get_samples():
            for index, header_sample in enumerate(header_samples):
                if fam_sample in header_sample:
                    mapping.append(index)
                    break
            else:
                raise Exception('Sample from fam: {} not found!'.format(fam_sample))
        self.mapping = np.array(mapping, dtype=int)
    
    def get_intensities(self, variant_index):
        variant_name = self.variants.get_name(variant_index)
        p = re.compile(variant_name)
        with open(self.file_path) as f:
            for line in f:
                if p.match(line):
                    data = line.split()
                    break
        xy = np.array(data[3:], dtype='float')
        xy = xy.reshape((-1, 2))
        xy = xy[self.mapping]
        return xy


class BinaryIntensity(object):
    
    def __init__(self, file_path, samples):
        self.file_path = file_path
        self.samples = samples
        self.HEADER_MAGIC = (26, 49)

        self.open_file()
        self.check_file()
        
    def open_file(self):
        self.f = open(self.file_path)
        
    def check_file(self):
        self.f.seek(0)
        magic = self.f.read(2)
        magic = unpack('2B', magic)
        assert  magic == self.HEADER_MAGIC

    def get_offset(self, variant_index):
        offset = len(self.HEADER_MAGIC)
        offset += 8 * variant_index * self.samples.get_n_samples()
        return offset
    
    def get_n_bytes(self):
        return 8 * self.samples.get_n_samples()
    
    def get_intensities(self, variant_index):
        offset = self.get_offset(variant_index)
        self.f.seek(offset)
        n_bytes = self.get_n_bytes()
        bytes = self.f.read(n_bytes)
        xy = unpack('f'*(n_bytes/4), bytes)
        xy = np.array(xy, dtype='float')
        xy = xy.reshape((-1, 2))
        return xy

    
class EvokerLite(object):
    
    def __init__(self, bfile_path, bnt_path=None, int_path= None, fam_path=None, bim_path=None, bed_path=None, exclude_list_path=None):
        self.samples = Samples(fam_path or (bfile_path + '.fam'))
        self.variants = Variants(bim_path or (bfile_path + '.bim'))
        self.genotypes = Genotypes(bed_path or (bfile_path + '.bed'),
                                   self.samples.get_n_samples(),
                                   self.variants.get_n_variants(),
                                  )
        
        if bnt_path:
            self.intensities = BinaryIntensity(bnt_path, self.samples)
        elif int_path:
            self.intensities = TextIntensity(int_path, self.variants, self.samples)
        else:
            raise Exception('Must provide a binary or intensity file')

        if exclude_list_path:
            exclude_list = []
            with open(self.exclude_list_path) as f:
                for line in f:
                    sample = line.strip()
                    exclude_list.append(sample)
            self.exclude_list = exclude_list
    
    def plot_intensities(self, variant_name):
        GENOTYPE_MAPPING = {
            '00': {'name': 'homozygous A1', 'color': 'blue'},
            '01': {'name': 'missing', 'color': 'lightgrey'},
            '10': {'name': 'heterozygous', 'color': 'limegreen'},
            '11': {'name': 'homozygous A2', 'color': 'red'},
        }
        variant_index = self.variants.get_index(variant_name)
        genotypes = self.genotypes.get_genotypes(variant_index)
        xy = self.intensities.get_intensities(variant_index)
#         xy = self.
        fig, ax = plt.subplots(figsize=(10, 10))
        for k, v in GENOTYPE_MAPPING.iteritems():
            t = xy[genotypes == k]
            ax.scatter(t[:,0], t[:,1],
                       color=v['color'],
                       s=50,
                       alpha=0.75,
                       lw=0.333,
                       edgecolor='black',
                      )
        xymin = -0.1
        xymax = 1.1 * xy.max()
        ax.set_xlim([xymin,xymax])
        ax.set_ylim([xymin,xymax])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
#         ax.spines['bottom'].set_visible(False)
#         ax.spines['left'].set_visible(False)
#         fig.savefig('test.png')
        return ax

if __name__ == '__main__':        
    bfile_path = '/lustre/scratch115/realdata/mdt0/teams/barrett/users/dr9/gwas3_eva/coreex_gaibdc.1'
    int_path = '/lustre/scratch115/realdata/mdt0/teams/barrett/users/dr9/gwas3_eva/coreex_gaibdc.synced.int'
    e = EvokerLite(bfile_path, int_path=int_path)
    e.plot_intensities('rs1039823')

    