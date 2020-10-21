import numpy as np
import os

if __name__ == "__main__":

    files = [f for f in os.listdir(os.curdir) if os.path.isfile(f) and f.split('.')[1] == 'dat']
    print(files)
    pass