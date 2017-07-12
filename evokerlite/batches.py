class Batches:

    def __init__(self, file_path):
        self.load_batches(file_path)
        
    def load_batches(self, file_path):
        with open(file_path) as f:
            self.batches = {x.strip(): line_number for line_number, x in enumerate(f.readlines())}
           
    def get_index(self, batch_name):
        return self.batches[batch_name]
    
    def get_n_batches(self):
        return len(self.batches)