#!/usr/bin/env python3

from caveman_wrapper.core import *
from caveman_wrapper.utils import *

def main():
    # Parse command line arguments
    parser = CavemanFlags.parser()
    args = parser.parse_args()
    
    # Initialise and run CaVEMan
    runner = CavemanRunner(**vars(args))
    runner.run_caveman()

if __name__ == "__main__":
    main()
