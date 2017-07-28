class Variants(object):
    
    def __init__(self, file_path):
        self.load_variants(file_path)

    def load_variants(self, file_path):
        variants = []
        with open(file_path) as f:
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
   
    def get_variant(self, variant_index):
        return self.variants[variant_index]


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
    
    def get_chrom(self):
        return self.chrom

    def __str__(self):
        return '{}:{} {}'.format(self.chrom, self.pos, self.name)
