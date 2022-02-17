from numpy import array


class Samples(object):
    def __init__(self, file_path, ukbiobank=False):
        self.__ukbiobank = ukbiobank
        self.load_samples(file_path)

    def load_samples(self, file_path):
        samples = []
        if self.__ukbiobank:
            batches = []
            sex = []
        with open(file_path) as f:
            for line in f:
                tokens = line.split()
                samples.append(tokens[0])
                if self.__ukbiobank:
                    sex.append(tokens[4])
                    batches.append(tokens[5])
        self.__samples = samples
        if self.__ukbiobank:
            self.__sex = array(sex, dtype=int)
            self.__batches = array(batches)

    def get_samples(self):
        return self.__samples

    def get_n_samples(self):
        return len(self.__samples)

    def get_index(self, sample_id):
        return self.__samples.index(sample_id)

    def get_sex(self):
        return self.__sex

    def get_batches(self):
        if not self.__ukbiobank:
            raise Exception("Cannot access batches if not UKBiobank.")
        return self.__batches

    def get_n_batches(self):
        if not self.__ukbiobank:
            raise Exception("Cannot access batches if not UKBiobank.")
        return len(self.__batches)

    def __str__(self):
        L = [f"samples:{self.__samples}"]
        if self.__ukbiobank:
            L += [f"sex:{self.__sex}", f"batches:{self.__batches}"]
        return "".join(L)
