import numpy as np
import matplotlib.pyplot as plt
from math import ceil
from struct import unpack
import seaborn as sns
import re
import os
import shutil


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
    
#     def get_counts(self, variant_index):
#         A1 = get_A2(variant_index)
#         A1 = get_A2(variant_index)
        
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
        batches = []
        with open(self.file_path) as f:
            for line in f:
                tokens = line.split()
                sample, batch = tokens[0], tokens[-1]
                samples.append(sample)
                batches.append(batch)
        self.__samples = samples
        self.__batches = np.array(batches)
        
    def get_samples(self):
        return self.__samples
    
    def get_batches(self):
        return self.__batches
    
    def get_n_samples(self):
        return len(self.__samples)
    
    def get_index(self, sample_id):
        return self.__samples.index(sample_id)


class Genotypes(object):
    
    def __init__(self, file_path, n_samples, n_variants):
        self.file_path = file_path
        self.n_samples = n_samples
        self.n_variants = n_variants
        self.HEADER_MAGIC = (108, 27, 1)
#         self.check_file()
        
        
    def check_file(self):
        with open(self.file_path, 'rb') as f:
            _bytes = f.read(3)
            magic = unpack('3B', _bytes)
            assert magic == self.HEADER_MAGIC

    def get_offset(self, variant_index):
        offset = len(self.HEADER_MAGIC)
        offset += ceil(self.n_samples / 4.) * variant_index
#         offset = ceil(self.n_samples / 4.) * variant_index
        return int(offset)
    
    def get_genotypes(self, variant_index):
        offset = self.get_offset(variant_index)
        n_bytes = int(ceil(self.n_samples / 4))
        with open(self.file_path, 'rb') as f:
            f.seek(offset)
            _bytes = f.read(n_bytes)
        _bytes = unpack('B'*n_bytes, _bytes)
        genotypes = []
        for byte in _bytes:
            byte = '{:08b}'.format(byte)
            genotypes.extend([byte[i:i+2] for i in range(6, -2, -2)])
        return np.array(genotypes[:self.n_samples])


class BinaryIntensity(object):
    
    def __init__(self, file_path, variants, samples):
        self.file_path = file_path
        self.variants = variants
        self.samples = samples
        self.HEADER_MAGIC = (26, 49)

        self.open_file()
#         self.check_file()
        
    def open_file(self):
        self.f = open(self.file_path, 'rb')
        
    def check_file(self):
        self.f.seek(0)
        magic = self.f.read(2)
        magic = unpack('2B', magic)
        assert  magic == self.HEADER_MAGIC

    def get_offset(self, variant_index):
#         offset = len(self.HEADER_MAGIC)
#         offset += 8 * variant_index * self.samples.get_n_samples()
        offset = 8 * variant_index * self.samples.get_n_samples()
        return offset
    
    def get_n_bytes(self):
        return 8 * self.samples.get_n_samples()
    
    def get_intensities_for_variant(self, variant_index):
        offset = self.get_offset(variant_index)
        self.f.seek(offset)
        n_bytes = self.get_n_bytes()
        _bytes = self.f.read(n_bytes)
        AB = unpack('f'*int(n_bytes/4), _bytes) # 'f' is IEEE 754 binary32 or 4 bytes (=4*8=32)
        AB = np.array(AB, dtype='float')
        A,B = AB[0::2], AB[1::2]
        x = np.log2(A/B)    # Contrast (x-axis) = log2(A/B)
        y = np.log2(A*B)/2  # Strength (y-axis) = (log2(A*B))/2
        xy = np.column_stack((x,y))
        return xy

    def get_intensities_for_sample(self, sample_id):
        xy = []
        sample_offset = 8 * self.samples.get_index(sample_id)
        for variant_index in range(self.variants.get_n_variants()):
            offset = self.get_offset(variant_index)
            offset += sample_offset
            self.f.seek(offset)
            _bytes = self.f.read(8)
            xy.append(unpack('ff', _bytes))
        return xy
    
class EvokerLite(object):
    
    def __init__(self, bfile_path=None, bnt_path=None, int_path= None, fam_path=None, bim_path=None, bed_path=None, exclude_list_path=None):
        self.samples = Samples(fam_path or (bfile_path + '.fam'))
        self.variants = Variants(bim_path or (bfile_path + '.bim'))
        self.genotypes = Genotypes(bed_path or (bfile_path + '.bed'),
                                   self.samples.get_n_samples(),
                                   self.variants.get_n_variants(),
                                  )
        
        if bnt_path:
            self.intensities = BinaryIntensity(bnt_path, self.variants, self.samples)
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
    
    def plot_intensities(self, variant_name, batch=False):
        variant_index = self.variants.get_index(variant_name)
        A1 = self.variants.get_A1(variant_index)
        A2 = self.variants.get_A2(variant_index)
        GENOTYPE_MAPPING = (
            {'code':'00', 'name':'homozygous A1', 'color':'blue', 'label': '{A1}{A1}'.format(A1=A1,A2=A2)},
            {'code':'10', 'name':'heterozygous', 'color':'limegreen', 'label': '{A1}{A2}|{A2}{A1}'.format(A1=A1,A2=A2)},
            {'code':'11', 'name':'homozygous A2', 'color':'red', 'label': '{A2}{A2}'.format(A1=A1,A2=A2)},
            {'code':'01', 'name':'missing', 'color':'lightgrey', 'alpha':0.9, 'lw':0.75, 'label': 'Missing'},
        )
        genotypes = self.genotypes.get_genotypes(variant_index)
        xy = self.intensities.get_intensities_for_variant(variant_index)
        if batch:
            batches = self.samples.get_batches()
            batch_indices = batches == batch
            genotypes = genotypes[batch_indices]
            xy = xy[batch_indices]
            
        fig, ax = plt.subplots(figsize=(10, 10))
        for m in GENOTYPE_MAPPING:
            code = m['code']
            t = xy[genotypes == code]
            ax.scatter(t[:,0], t[:,1],
                       color=m['color'],
                       s=50,
                       alpha=m.get('alpha', 0.75),
                       lw=m.get('lw', 0.25),
                       edgecolor='black',
                       label=m['label']
                      )
#         xymin = min(-0.1, xy.min())
#         xymax = 1.1 * xy.max()
#         ax.set_xlim([xymin,xymax])
#         ax.set_ylim([xymin,xymax])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ax.legend()
        return fig

    def plot_save_all_batches(self, variant_name, parent_directory):
        directory = os.path.join(parent_directory, variant_name)
        try:
            shutil.rmtree(directory)
        except FileNotFoundError:
            print('FileNotFoundError')
        os.mkdir(directory)
        batches = list(set(self.samples.get_batches()))
        for batch in batches:
            fig = self.plot_intensities(variant_name, batch)
            filename = os.path.join(directory, batch + '.png')
            fig.savefig(filename)
            plt.close('all')


if __name__ == '__main__':
    bfile_path = '/lustre/scratch115/realdata/mdt0/teams/barrett/users/dr9/UKBioBankEvoker/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.chr22'
    bnt_path = '/lustre/scratch115/realdata/mdt0/teams/barrett/users/dr9/UKBioBankEvoker/V2_QCed.export.intensity.chr22.bin'
    bb = EvokerLite(bfile_path, bnt_path=bnt_path)
    bb.plot_save_all_batches('rs62224618', '/lustre/scratch115/realdata/mdt0/teams/barrett/users/dr9/UKBioBankEvoker/plots')