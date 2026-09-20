"""Native Phase 29 entry point; all experiments remain unperformed."""
import os
import sys
from pathlib import Path

os.environ.setdefault('OMP_NUM_THREADS', '2')
os.environ.setdefault('MKL_NUM_THREADS', '2')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
sys.path.insert(0, str(Path(__file__).resolve().parent / 'code'))
from phase29 import main

if __name__ == '__main__':
    main()
