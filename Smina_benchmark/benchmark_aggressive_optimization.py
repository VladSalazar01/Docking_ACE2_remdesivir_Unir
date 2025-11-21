#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent))

from benchmark_smina_vinardo_v2 import *

class AggressiveConfig(Config):
    CONFIGS = {
        "local_refine": {
            "USE_LOCAL_ONLY": True,
            "EXHAUSTIVENESS": 64,
            "AUTOBOX_ADD": 4,
            "NUM_MODES": 10,
            "SCORING_FUNCTIONS": ["vinardo"]
        },
        "high_exhaustiveness": {
            "USE_LOCAL_ONLY": False,
            "EXHAUSTIVENESS": 128,
            "AUTOBOX_ADD": 6,
            "NUM_MODES": 30,
            "SCORING_FUNCTIONS": ["vinardo", "vina"]
        },
        "tight_box": {
            "USE_LOCAL_ONLY": False,
            "EXHAUSTIVENESS": 64,
            "AUTOBOX_ADD": 3,
            "NUM_MODES": 20,
            "SCORING_FUNCTIONS": ["vinardo"]
        },
        "minimize": {
            "USE_LOCAL_ONLY": False,
            "USE_MINIMIZE": True,
            "EXHAUSTIVENESS": 48,
            "AUTOBOX_ADD": 6,
            "NUM_MODES": 20,
            "SCORING_FUNCTIONS": ["vinardo"]
        }
    }


class AggressiveBenchmark(SminaBenchmark):
    def run_aggressive_tests(self):
        print("\n" + "=" * 70)
        print("AGGRESSIVE OPTIMIZATION TESTS")
        print("=" * 70)
        
        original_config = {
            'USE_LOCAL_ONLY': Config.USE_LOCAL_ONLY,
            'USE_MINIMIZE': Config.USE_MINIMIZE,
            'EXHAUSTIVENESS': Config.EXHAUSTIVENESS,
            'AUTOBOX_ADD': Config.AUTOBOX_ADD,
            'NUM_MODES': Config.NUM_MODES
        }
        
        for config_name, config_params in AggressiveConfig.CONFIGS.items():
            print("\n" + "*" * 70)
            print("TESTING: " + config_name.upper())
            print("*" * 70)
            print("Parameters: " + str(config_params))
            
            for key, value in config_params.items():
                if key != 'SCORING_FUNCTIONS':
                    setattr(Config, key, value)
            
            for scoring in config_params['SCORING_FUNCTIONS']:
                run_name = config_name + "_" + scoring
                self.run_redocking(run_name, scoring)
        
        for key, value in original_config.items():
            setattr(Config, key, value)


def main():
    print("""
================================================================
         AGGRESSIVE RMSD OPTIMIZATION
         Multiple strategies to improve docking quality
================================================================
""")
    
    benchmark = AggressiveBenchmark()
    
    if not benchmark.setup():
        return 1
    
    try:
        benchmark.run_aggressive_tests()
        benchmark.save_results()
        benchmark.print_summary()
        
        print("\n" + "=" * 70)
        print("OPTIMIZATION COMPLETE")
        print("=" * 70)
        
        return 0
    except Exception as e:
        print("\nError: " + str(e))
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
