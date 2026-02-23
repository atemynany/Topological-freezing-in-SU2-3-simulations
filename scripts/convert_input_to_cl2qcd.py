#!/usr/bin/env python3
"""
Convert cl2qcd_heatbath.txt input format to CL2QCD native input format.

Usage:
    python convert_input_to_cl2qcd.py input/cl2qcd_heatbath.txt output_dir/cl2qcd_input
"""

import sys
import os

def parse_input_file(filepath):
    """Parse the cl2qcd_heatbath.txt input file."""
    params = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue
            # Parse key-value pairs
            parts = line.split()
            if len(parts) >= 2:
                key = parts[0]
                value = parts[1]
                params[key] = value
    return params

def generate_cl2qcd_input(params, output_path):
    """Generate CL2QCD native input file."""
    
    # Map our parameter names to CL2QCD names
    cl2qcd_content = f"""# CL2QCD Input File (auto-generated)
# Generated from: {params.get('_source_file', 'unknown')}

# Hardware and lattice parameters
nSpace={params.get('nSpace', '8')}
nTime={params.get('nTime', '8')}
useGPU={params.get('useGPU', '0')}
precision={params.get('precision', '64')}

# Physics parameters
startCondition={params.get('startCondition', 'cold')}
beta={params.get('beta', '6.0')}

# Monte Carlo parameters
nThermalizationSteps={params.get('nThermalizationSteps', '0')}
nHeatbathSteps={params.get('nHeatbathSteps', '1000')}
nOverrelaxationSteps={params.get('nOverrelaxationSteps', '1')}

# Output parameters
onlineMeasureEvery={params.get('onlineMeasureEvery', '1')}
createCheckpointEvery={params.get('saveFrequency', '100')}
writeFrequency={params.get('saveFrequency', '100')}
overwriteTemporaryCheckpointEvery={params.get('saveFrequency', '100')}

# Config file naming
configPrefix={params.get('configPrefix', 'conf_')}
configPostfix=
configNumberDigits=5
"""
    
    with open(output_path, 'w') as f:
        f.write(cl2qcd_content)
    
    print(f"Generated CL2QCD input: {output_path}")
    return output_path

def generate_topcharge_input(params, output_path):
    """Generate input file for topological charge measurement."""
    
    nSpace = int(params.get('nSpace', '8'))
    nTime = int(params.get('nTime', '8'))
    save_freq = int(params.get('saveFrequency', '10'))
    n_therm = int(params.get('nThermalizationSteps', '100'))
    n_heatbath = int(params.get('nHeatbathSteps', '500'))
    
    # Calculate config numbers (after thermalization)
    start_conf = n_therm + save_freq
    end_conf = n_therm + n_heatbath
    
    content = f"""# Topological Charge Measurement Input (auto-generated)
# For converted CL2QCD configs

# Lattice
T                   {nTime}
L                   {nSpace}

# Config range
start_conf          {start_conf}
end_conf            {end_conf}
conf_step           {save_freq}

# Smearing
smear_steps         {params.get('smearSteps', '20')}
smear_alpha         {params.get('smearAlpha', '0.3')}

# Directories
config_dir          {params.get('convertedDir', 'output/configs_su3')}
output_file         {params.get('resultsDir', 'output/results')}/topcharge_cl2qcd.dat
"""
    
    with open(output_path, 'w') as f:
        f.write(content)
    
    print(f"Generated topcharge input: {output_path}")
    return output_path

def main():
    if len(sys.argv) < 2:
        print("Usage: python convert_input_to_cl2qcd.py <input_file> [output_dir]")
        print("Example: python convert_input_to_cl2qcd.py input/cl2qcd_heatbath.txt")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "."
    
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found")
        sys.exit(1)
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Parse input
    params = parse_input_file(input_file)
    params['_source_file'] = input_file
    
    # Generate CL2QCD input
    cl2qcd_input = os.path.join(output_dir, "cl2qcd_native_input")
    generate_cl2qcd_input(params, cl2qcd_input)
    
    # Generate topcharge measurement input
    topcharge_input = os.path.join(output_dir, "topcharge_input.txt")
    generate_topcharge_input(params, topcharge_input)
    
    # Print summary
    print("\n" + "="*60)
    print("Conversion Summary:")
    print("="*60)
    print(f"  Lattice: {params.get('nTime')} x {params.get('nSpace')}^3")
    print(f"  Beta: {params.get('beta')}")
    print(f"  Thermalization: {params.get('nThermalizationSteps')} sweeps")
    print(f"  Production: {params.get('nHeatbathSteps')} sweeps")
    print(f"  Save frequency: {params.get('saveFrequency')}")
    print(f"  GPU: {'Yes' if params.get('useGPU') == '1' else 'No'}")
    print("="*60)
    
    return params

if __name__ == "__main__":
    main()
