import os
import sys

# spectra_codec.py lives two directories up from this examples folder
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))

from huggingface_hub import list_repo_files, snapshot_download
from spectra_codec import SpectraCodec
import pandas as pd
import pymzml
import numpy as np
from scipy.interpolate import interp1d


def multiprocessing_encode_file(filename_dict):
    from spectra_codec import SpectraCodec

    message = filename_dict['message']
    original_filename = filename_dict['original_filename']
    output_filename = filename_dict['output_filename']
    encoder = SpectraCodec()
    encoder.encode_message_to_file(message, original_filename, output_filename, method='hilbert')
    encoder = SpectraCodec()
    decoded_message = encoder.decode_message_from_file(output_filename)#, parser='pymzml', method='hilbert')
    # do an assert message == decoded_message and catch any exceptions
    if message != decoded_message:
        print(f"File {output_filename} encoding/decoding failed")
        raise ValueError(f"Decoded message does not match original: {message} != {decoded_message}")
    
    
def set_huggingface_env(token_file, cache_dir):
    """Set Hugging Face environment variables."""
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)
    os.environ['HF_HOME'] = cache_dir
    try:
        with open(token_file, 'r') as f:
            token = f.read().strip()
        if token:
            os.environ['HUGGINGFACE_HUB_TOKEN'] = token
            os.environ['HF_TOKEN'] = token
            print("Environment variables set.")
        else:
            print("Token is empty. Environment variables not set.")
    except FileNotFoundError:
        print(f"Error: Token file not found at {token_file}")
    except Exception as e:
        print(f"An error occurred while reading the token file: {e}")

def check_repo(repo_id):
    """Check repository type and list files."""
    for repo_type in ["dataset", "model"]:
        try:
            files = list_repo_files(repo_id, repo_type=repo_type)
            # print(f"Repository type: {repo_type}")
            # print("Files:")
            # for file in files:
                # print(f"  - {file}")
            return repo_type, files
        except:
            continue
    print(f"Repository '{repo_id}' not found or not accessible")
    return None, None

def download_repo(repo_id, repo_type="dataset", cache_dir="./hf_cache"):
    """Download repository files to local cache."""
    # see if environment variable is set
    if 'HF_HOME' in os.environ:
        cache_dir = os.environ['HF_HOME']
    else:
        print("Using default cache directory:", cache_dir)
    try:
        local_dir = snapshot_download(
            repo_id=repo_id,
            repo_type=repo_type,
            cache_dir=cache_dir
        )

        return local_dir
    except Exception as e:
        print(f"Error downloading: {e}")
        return None

def get_paths_downloaded_files(local_dir,file_extension='.mzML'):
    """Processes all files in the downloaded directory."""
    all_files = []
    for root, _, files in os.walk(local_dir):
        for file in files:
            filepath = os.path.join(root, file)
            if file.endswith(file_extension):
                # print(f"Found file: {filepath}")
                all_files.append(filepath)
    return all_files


def get_eic(parameters,min_intensity=200000):

    """Extract the EIC for a target m/z from a mzML file.
    
    Parameters:
        parameters (dict): A dictionary containing the following
        keys:
            - 'filename': Path to the mzML file.
            - 'target_mz': The target m/z value to extract the EIC for.
            - 'mz_tolerance': The m/z tolerance for peak detection (in Da).
    Returns:
        list: A list of tuples containing the scan time, intensity, m/z value,
              and filename for each detected peak in the EIC.
    """
    mzml_file = parameters['filename']
    target_mz = parameters['target_mz']
    target_tolerance = parameters['mz_tolerance']
    run = pymzml.run.Reader(mzml_file, use_index=True, MS1_only=True)
    eic = []
    for spectrum in run:
        mz = np.array(spectrum.mz)
        intensity = np.array(spectrum.i)
        idx = np.where((mz > target_mz-target_tolerance) & (mz < target_mz+target_tolerance) & (intensity > min_intensity))
        if len(idx[0]) > 0:
            mz_vals = mz[idx]
            intensity_vals = intensity[idx]
            max_idx = np.argmax(intensity_vals)
            mz_at_max = mz_vals[max_idx]
            intensity_at_max = intensity_vals[max_idx]
            scan_time = spectrum.scan_time
            eic.append((scan_time, intensity_at_max, mz_at_max, mzml_file))
    return eic



def get_ms1_data(mzml_file):

    """
    Extracts MS1 data from a mzML file.
    Parameters:
        mzml_file (str): Path to the mzML file.
    Returns:
        dict: A dictionary containing the following keys:
            - 'mz': A numpy array of m/z values.
            - 'intensity': A numpy array of intensity values.
            - 'rt': A numpy array of retention times corresponding to the m/z and intensity values
    """
    run = pymzml.run.Reader(mzml_file, use_index=True, MS1_only=True)
    mz = []
    intensity = []
    rt = []
    for spectrum in run:
        mz_one_spec = spectrum.mz
        intensity_one_spec = spectrum.i
        if len(mz_one_spec) > 0:
            scan_time = spectrum.scan_time[0]
            mz.extend(mz_one_spec)
            intensity.extend(intensity_one_spec)
            rt.extend([scan_time] * len(mz_one_spec))
    output = {
        'mz': mz,
        'intensity': intensity,
        'rt': rt
    }
    return output

def get_peak_height(data, metabolites_df, min_rt=0.6, ppm_tolerance=5):
    """
    Extracts peak heights from MS1 data based on m/z values of metabolites.
    Parameters:
        data (dict): A dictionary containing 'mz', 'intensity', and 'rt' keys
        metabolites_df (pd.DataFrame): A DataFrame containing metabolite information with a 'pos_mz' column.
        min_rt (float): Minimum retention time to consider for peak extraction.
        ppm_tolerance (float): PPM tolerance for m/z matching.
    Returns:
        pd.DataFrame: A DataFrame containing detected peaks information
    """
    # Used for manuscript figures only; requires https://github.com/biorack/metatlas
    try:
        from metatlas.io.feature_tools import (
            group_consecutive,
            map_mzgroups_to_data,
            filter_raw_data_using_atlas,
        )
    except ImportError as e:
        raise ImportError(
            "get_peak_height requires the metatlas package "
            "(https://github.com/biorack/metatlas); it is only needed to "
            "reproduce the manuscript figures.") from e

    # Create DataFrame from data and filter by minimum RT
    df = pd.DataFrame({'rt': data['rt'], 'mz': data['mz'], 'i': data['intensity']})
    df = df[df['rt'] > min_rt]
    
    # Create atlas with necessary columns for the filtering process
    atlas = pd.DataFrame({
        'mz': metabolites_df['pos_mz'].values,
        'rt_min': 0,  # Set to 0 to only filter by m/z, not RT
        'rt_max': float('inf'),  # Set to infinity to only filter by m/z, not RT
        'label': metabolites_df.index.values,
        'extra_time': 0,
        'ppm_tolerance': ppm_tolerance
    })
    
    # Group consecutive m/z values that are within ppm_tolerance
    atlas['group_index'] = group_consecutive(atlas['mz'].values, 
                                           stepsize=ppm_tolerance, 
                                           do_ppm=True)
    
    # Map mz groups to data
    df['group_index'] = map_mzgroups_to_data(atlas['mz'].values,
                                           atlas['group_index'].values,
                                           df['mz'].values)
    
    # Filter raw data using the atlas
    filtered_df = filter_raw_data_using_atlas(atlas, df)
    filtered_df = filtered_df[filtered_df['in_feature'] == True]
    
    # For each metabolite/group, find the peak with the highest intensity
    result = []
    for label, group_data in filtered_df.groupby('label'):
        if len(group_data) > 0:
            idx = group_data['i'].idxmax()
            peak = group_data.loc[idx]
            area = group_data['i'].sum()
            # Get the original metabolite index
            metabolite_idx = metabolites_df.index[metabolites_df.index == label][0]
            
            result.append({
                'mz': peak['mz'],
                'intensity': peak['i'],
                'area': area,
                'rt': peak['rt'],
                'ppm_diff': abs(peak['mz'] - peak['mz_atlas']) / peak['mz_atlas'] * 1e6,
                'target_mz': peak['mz_atlas'],
                'target_index': metabolite_idx
            })
    
    # Create DataFrame from results
    feature_df = pd.DataFrame(result) if result else pd.DataFrame(
        columns=['mz', 'intensity', 'rt', 'ppm_diff', 'target_mz', 'target_index']
    )
    
    # Merge with metabolites_df to get the metabolite names
    if len(feature_df) > 0:
        feature_df = feature_df.merge(metabolites_df, left_on='target_index', right_index=True, how='left')
    
    # Sort by intensity and remove duplicates
    feature_df.sort_values(by='intensity', ascending=False, inplace=True)
    feature_df.drop_duplicates('target_index', inplace=True)
    feature_df.reset_index(drop=True, inplace=True)
    
    return feature_df


def get_compound_data_from_files(my_file, metabolites_df):
    """
    Get data from files in the files_dict.
    """
    data = get_ms1_data(my_file)
    feature_df = get_peak_height(data, metabolites_df, min_rt=0.6)
    feature_df['filename'] = my_file
    return feature_df