#!/bin/bash

# Set environment variables
ENV_NAME="spectra_codec_env"
YAML_FILE="environment.yml"

# Check if conda is installed
if ! command -v conda &> /dev/null; then
  echo "Error: conda is not installed. Please install Miniconda or Anaconda."
  exit 1
fi

# Check if the YAML file exists
if [ ! -f "$YAML_FILE" ]; then
  echo "Error: environment.yml not found. Please create it."
  exit 1
fi

# Create the conda environment
conda env create -f "$YAML_FILE" -n "$ENV_NAME" &> conda_output.txt

# Check the exit code of conda
if [ $? -ne 0 ]; then
  echo "Error creating conda environment. Check conda_output.txt for details."
  exit 1
fi

# Check if the environment was created successfully
if conda info --envs | grep "$ENV_NAME"; then
  echo "Conda environment '$ENV_NAME' created successfully."
else
  echo "Error: Conda environment '$ENV_NAME' creation failed. Check conda_output.txt for details."
  exit 1
fi

echo "Environment file: $YAML_FILE"