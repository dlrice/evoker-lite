from numpy import column_stack, reshape, array, log2
from struct import unpack
import re


class TextIntensity(object):

    def __init__(self, file_path, variants, samples):
        self.__file_path = file_path
        self.__variants = variants
        self.__samples = samples
        self.map_sample_ids()

    def map_sample_ids(self):
        with open(self.__file_path) as f:
            header = f.readline()
        tokens = header.split()
        header_samples = tokens[3::2]
        header_samples = [x[:-1] for x in header_samples]

        # Check if these variants exist in fam file
        # s = set(header_samples)
        # t = set(self.samples.get_samples())
        # intersection = s & t
        # if len(intersection) < min(len(s), len(t)):
        mapping = []
        for fam_sample in self.__samples.get_samples():
            for index, header_sample in enumerate(header_samples):
                if fam_sample in header_sample:
                    mapping.append(index)
                    break
            else:
                raise Exception(
                    'Sample from fam: {} not found!'.format(fam_sample))
        self.mapping = np.array(mapping, dtype=int)

    def get_intensities_for_variant(self, variant_index):
        variant_name = self.__variants.get_name(variant_index)
        p = re.compile(variant_name)
        with open(self.__file_path) as f:
            for line in f:
                if p.match(line):
                    data = line.split()
                    break
        xy = array(data[3:], dtype='float')
        xy = xy.reshape((-1, 2))
        xy = xy[self.mapping]
        return xy


class BinaryIntensity(object):

    def __init__(self, file_path, samples, ukbiobank=False):
        self.__samples = samples
        if ukbiobank:
            self.__HEADER_MAGIC = ()
        else:
            self.__HEADER_MAGIC = (26, 49)
        self.open_file(file_path)

    def open_file(self, file_path):
        self.__f = open(file_path, 'rb')

    def check_file(self):
        self.__f.seek(0)
        magic = self.__f.read(len(self.__HEADER_MAGIC))
        magic = unpack('2B', magic)
        assert magic == self.__HEADER_MAGIC

    def get_offset(self, variant_index):
        offset = len(self.__HEADER_MAGIC)
        offset += 8 * variant_index * self.__samples.get_n_samples()
        return offset

    def get_n_bytes(self):
        return 8 * self.__samples.get_n_samples()

    def get_intensities_for_variant(self, variant_index, transform=None):
        offset = self.get_offset(variant_index)
        self.__f.seek(offset)
        n_bytes = self.get_n_bytes()
        _bytes = self.__f.read(n_bytes)
        # 'f' is IEEE 754 binary32 or 4 bytes (=4*8=32)
        AB = unpack('f'*int(n_bytes/4), _bytes)
        AB = array(AB, dtype='float')
        A, B = AB[0::2], AB[1::2]
        if transform:
            # Contrast (x-axis) = log2(A/B)
            x = log2(A/B)
            # Strength (y-axis) = (log2(A*B))/2
            y = log2(A*B)/2
        else:
            x = A
            y = B
        xy = column_stack((x, y))
        return xy

    # def get_intensities_for_sample(self, sample_id):
    #     xy = []
    #     sample_offset = 8 * self.samples.get_index(sample_id)
    #     for variant_index in range(self.variants.get_n_variants()):
    #         offset = self.get_offset(variant_index)
    #         offset += sample_offset
    #         self.f.seek(offset)
    #         _bytes = self.f.read(8)
    #         xy.append(unpack('ff', _bytes))
    #     return xy
