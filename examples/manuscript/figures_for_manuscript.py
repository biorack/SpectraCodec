import matplotlib.pyplot as plt
import numpy as np
import string
import json
import os

import random
from matplotlib.ticker import LogLocator, ScalarFormatter
import matplotlib.gridspec as gridspec
from spectra_codec import SpectraCodec

# Set up the figure with publication-quality settings
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['xtick.major.width'] = 1
plt.rcParams['ytick.major.width'] = 1
plt.rcParams['xtick.minor.width'] = 0.5
plt.rcParams['ytick.minor.width'] = 0.5

def make_spectrum_figures_showing_encoding(original_filename, output_filename):
    # get the first spectrum o the original file
    from pymzml import run
    original_run = run.Reader(original_filename)
    for s in original_run:
        original_spectrum = s
        break  # Get the first spectrum only

    output_run = run.Reader(output_filename)
    for s in output_run:
        encoded_spectrum = s
        break  # Get the first spectrum only

    fig,ax = plt.subplots(1,2, figsize=(8,4), sharey=True, sharex=True)
    ax[0].vlines(original_spectrum.mz,0, original_spectrum.i, color='blue', label='Original Spectrum')
    ax[0].set_xlabel('m/z', fontsize=18)
    ax[0].set_ylabel('Intensity', fontsize=18)
    ax[0].set_title('a', fontsize=24, loc='left', pad=20,x=0.015,y=0.8)
    ax[1].vlines(encoded_spectrum.mz, 0, encoded_spectrum.i, color='red', label='Encoded Spectrum')
    ax[1].set_xlabel('m/z', fontsize=18)
    ax[1].set_title('b', fontsize=24, loc='left', pad=20,x=0.015,y=0.8)

    # ax[1].set_ylabel('Intensity')

    for a in ax:
        a.set_xlim(0,100)
        a.set_ylim(1,1e8)
        a.set_yscale('log')
        # increase tick label size
        a.tick_params(axis='both', which='major', labelsize=14)
    plt.tight_layout()
    return fig
    # encoded_spectrum.peaks.plot(ax=ax[1], title='Encoded Spectrum')
    

def make_scaling_figure(raw_data_folder, converted_data_folder, metadata_df, original_manuscript,min_length=1,max_length=7.5,num_points=10):
    # Create figure with custom gridspec for better control
    fig = plt.figure(figsize=(7, 8))
    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1], hspace=0.3)

    # Create subplots
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax3 = fig.add_subplot(gs[2])

    # Function to generate random unicode text
    def generate_random_unicode_text(length):
        # Mix of ASCII and some common unicode characters
        chars = string.ascii_letters + string.digits + string.punctuation
        unicode_chars = 'αβγδεζηθικλμνξοπρστυφχψω∑∏∫√≈≠≤≥±∞'
        all_chars = chars + unicode_chars
        return ''.join(random.choice(all_chars) for _ in range(length))

    # Generate message lengths (log scale from 10 to 10000)
    message_lengths = np.logspace(min_length,max_length,num_points, dtype=int)

    # Placeholder lists for results
    mz_ranges = []
    unique_mz_coords = []
    unique_intensity_coords = []
    hilbert_orders = []

    # Simulate encoding for each message length
    for length in message_lengths:
        # Generate random message
        message = generate_random_unicode_text(length)
        encoder = SpectraCodec()  # Assuming SpectraCodec is defined elsewhere
        message = encoder.string_to_binary_matrix(message)
        encoder.determine_hilbert_curve_order(len(message))
        coords = np.column_stack((encoder.x_indices, encoder.y_indices))
        indices = np.argwhere(np.asarray(encoder.encoded_message_vector)!=0).flatten()
        # print(coords[indices,:])
        mz, intensity = encoder.scale_coords(coords, indices)
        
        mz_range = mz.max() - mz.min()
        mz_ranges.append(mz_range)
        
        unique_mz = len(set(coords[:,0]))
        unique_intensity = len(set(coords[:,1]))    
        unique_mz_coords.append(unique_mz)
        unique_intensity_coords.append(unique_intensity)
        
        
        
        hilbert_orders.append(encoder.order)

    # make the calcs for a real metadata message
    my_file = '20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_51_S40-D89_C_Rg80to1200-CE102040-soil-S1_Run227.mzML'
    original_filename = os.path.join(raw_data_folder, my_file)
    output_filename = os.path.join(converted_data_folder, my_file)
    message = metadata_df[my_file].to_dict()
    message['manuscript'] = original_manuscript
    message = json.dumps(message)
    real_message_length = len(message)
    encoder = SpectraCodec()  # Assuming SpectraCodec is defined elsewhere
    message = encoder.string_to_binary_matrix(message)
    encoder.determine_hilbert_curve_order(len(message))
    coords = np.column_stack((encoder.x_indices, encoder.y_indices))
    indices = np.argwhere(np.asarray(encoder.encoded_message_vector)!=0).flatten()
    # print(coords[indices,:])
    mz, intensity = encoder.scale_coords(coords, indices)
    real_mz_range = mz.max() - mz.min()
    real_num_unique_mz = len(set(coords[:,0]))
    real_order = encoder.order



    # Plot 1: m/z range vs message length
    ax1.semilogx(message_lengths, mz_ranges, 'o-', color='#1f77b4', 
            markersize=6, linewidth=2, markeredgecolor='white', markeredgewidth=0.5)
    ax1.axhline(real_mz_range, color='gray', linestyle='--', linewidth=1.5, 
                label=f'Real m/z Range: {real_mz_range:.2f}')
    ax1.axvline(real_message_length, color='gray', linestyle='--', linewidth=1.5, 
                label=f'Real Message Length: {real_message_length} characters')
    # ax1.legend(loc='upper left',
    ax1.set_ylabel('m/z Range (Da)', fontsize=14, fontweight='normal')
    ax1.grid(True, which='both', alpha=0.3, linestyle='-', linewidth=0.5)

    # lower the text 
    ax1.set_title('a', loc='left', fontsize=24, fontweight='normal', pad=10, y=0.8,x=0.015)

    # # Add trend line
    # z = np.polyfit(np.log(message_lengths), np.log(mz_ranges), 1)
    # p = np.poly1d(z)
    # ax1.loglog(message_lengths, np.exp(p(np.log(message_lengths))), 
    #            '--', color='red', alpha=0.7, linewidth=1.5, 
    #            label=f'Slope: {z[0]:.2f}')
    # ax1.legend(loc='upper left', frameon=True, fancybox=False, edgecolor='black')

    # Plot 2: Number of unique coordinates
    ax2.semilogx(message_lengths, unique_mz_coords, 'o-', color='#ff7f0e', 
                markersize=6, linewidth=2, label='Either m/z or intensity dimension',
                markeredgecolor='white', markeredgewidth=0.5)
    ax2.axhline(real_num_unique_mz, color='gray', linestyle='--', linewidth=1.5,
                label=f'Real Unique m/z Coordinates: {real_num_unique_mz}')
    ax2.axvline(real_message_length, color='gray', linestyle='--', linewidth=1.5, 
                label=f'Real Message Length: {real_message_length} characters')
    # ax2.semilogx(message_lengths, unique_intensity_coords, 's-', color='#2ca02c', 
    #              markersize=6, linewidth=2, label='Intensity dimension',
    #              markeredgecolor='white', markeredgewidth=0.5)
    ax2.set_ylabel('Unique Coordinates', fontsize=14, fontweight='normal')
    ax2.grid(True, which='both', alpha=0.3, linestyle='-', linewidth=0.5)
    # ax2.legend(loc='upper left', frameon=True, fancybox=False, edgecolor='black')
    ax2.set_title('b', loc='left', fontsize=24, fontweight='normal', pad=10, y=0.8,x=0.015)

    # Plot 3: Hilbert order required
    ax3.semilogx(message_lengths, hilbert_orders, 'o-', color='#d62728', 
                markersize=6, linewidth=2, markeredgecolor='white', 
                markeredgewidth=0.5, drawstyle='steps-post')
    ax3.axhline(real_order, color='gray', linestyle='--', linewidth=1.5, 
                label=f'Real Hilbert Order: {real_order}')
    ax3.axvline(real_message_length, color='gray', linestyle='--', linewidth=1.5, 
                label=f'Real Message Length: {real_message_length} characters')
    ax3.set_xlabel('Message Length (characters)', fontsize=14, fontweight='normal')
    ax3.set_ylabel('Hilbert Order', fontsize=14, fontweight='normal')
    ax3.grid(True, which='both', alpha=0.3, linestyle='-', linewidth=0.5)
    ax3.set_title('c', loc='left', fontsize=24, fontweight='normal', pad=10, y=0.8,x=0.015)

    # Set integer y-ticks for Hilbert order
    ax3.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

    # Format x-axis for all subplots
    for ax in [ax1, ax2, ax3]:
        ax.set_xlim([message_lengths[0] * 0.8, message_lengths[-1] * 1.2])
        ax.xaxis.set_major_locator(LogLocator(base=10))
        ax.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
        # make tick labels bigger
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
        # Remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Make remaining spines thicker
        ax.spines['left'].set_linewidth(1.2)
        ax.spines['bottom'].set_linewidth(1.2)

    # Only show x-label on bottom plot
    ax1.set_xticklabels([])
    ax2.set_xticklabels([])

    # Add main title
    # fig.suptitle('SpectraCodec Scaling Properties', fontsize=14, fontweight='bold', y=0.98)

    # Adjust layout
    plt.tight_layout()

    # Save figure
    # plt.savefig('spectracodec_scaling.pdf', dpi=300, bbox_inches='tight')
    # plt.savefig('spectracodec_scaling.png', dpi=300, bbox_inches='tight')

    plt.show()
    return fig
# from spectra_codec import SpectraCodec
# import numpy as np
# import matplotlib.pyplot as plt
# def make_comparison_figure():
#     import numpy as np
#     import matplotlib.pyplot as plt

#     # Parameters
#     min_mz = 5
#     mz_scale = 0.0001
#     min_intensity = 100
#     intensity_scale = 0.001

#     # Message lengths from 1 to 500k
#     message_lengths = np.logspace(0, np.log10(500000), 100).astype(int)
#     orders = [5, 10, 13, 15]
#     colors = ['green', 'orange', 'purple', 'blue']

#     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

#     # Panel A: Required m/z Range
#     seq_mz_range = (message_lengths - 1) * mz_scale
#     ax1.plot(message_lengths, seq_mz_range, 'r-', label='Sequential-ASCII', linewidth=2)

#     # Hilbert curve encoding
#     order = 13
#     max_chars = 4**order // 128
#     valid_lengths = message_lengths[message_lengths <= max_chars]
#     hilbert_extent = np.minimum(np.sqrt(valid_lengths * 128), 2**order - 1)
#     hilb_mz_max = hilbert_extent * mz_scale

#     ax1.plot(valid_lengths, hilb_mz_max, '-', color='purple', 
#             label=f'Hilbert Order {order}', linewidth=2)
#     ax1.set_xscale('log')
#     ax1.set_ylabel('m/z range (Da)', fontsize=14)
#     ax1.set_xlabel('Message length (characters)', fontsize=14)
#     ax1.set_ylim(0, 5)
#     # ax1.grid(True, alpha=0.3)
#     ax1.text(0.02, 0.98, 'a', transform=ax1.transAxes, fontsize=24, fontweight='normal', 
#             verticalalignment='top')

#     # Panel B: Encoding Capacity vs. Space Constraints
#     max_dimensions = np.logspace(0, 5, 50)
#     seq_capacity = max_dimensions

#     ax2.loglog(max_dimensions, seq_capacity, 'r-', linewidth=3, label='Sequential-ASCII')
#     line_data_x = []
#     line_data_y = []
#     for i, order in enumerate(orders):
#         hilbert_capacity = 4**order // 128
#         one_dimensional_capacity = 2**order - 1

#         ax2.loglog([one_dimensional_capacity], [hilbert_capacity], 'o', color=colors[i], markersize=14,
#                 label=f'Hilbert Order {order}')
#         line_data_x.append(one_dimensional_capacity)
#         line_data_y.append(hilbert_capacity)

#         if order == 13:
#             # Add a horizontal line for the maximum capacity
#             ax2.axhline(y=hilbert_capacity, color=colors[i], linestyle='--', linewidth=1.5)
#             ax2.axvline(x=one_dimensional_capacity, color=colors[i], linestyle='--', linewidth=1.5)
#             ax2.text(1.05, hilbert_capacity * 1.025,
#                     f'{hilbert_capacity} Characters', color='black', fontsize=14, ha='left', va='bottom')
            
#     ax2.plot(line_data_x, line_data_y, '-', color='black', linewidth=2,alpha=0.5)

#     ax2.set_xlabel('Maximum coordinates in one dimension', fontsize=14)
#     ax2.set_ylabel('Character capacity', fontsize=14)
#     # ax2.grid(True, alpha=0.3)
#     ax2.text(0.02, 0.98, 'b', transform=ax2.transAxes, fontsize=24, fontweight='normal', 
#             verticalalignment='top')

#     # put legend outside the plot
#     ax2.legend(loc='lower right', fontsize=14)

#     # increase font size for both axes
#     ax1.tick_params(axis='both', which='major', labelsize=14)
#     ax2.tick_params(axis='both', which='major', labelsize=14)

#     # remove top and right spines
#     for ax in (ax1, ax2):
#         ax.spines['top'].set_visible(False)
#         ax.spines['right'].set_visible(False)

#     plt.tight_layout()
#     plt.show()

# # def make_abcd_figure():
# #     # message = "the quick brown fox jumped over the lazy dog"
# #     message = [chr(i) for i in range(97,123)]
# #     message = ''.join(message)
# #     encoder = SpectraCodec(order=13)
# #     mz,intensity = encoder.sequential_encode_message(message)

# #     mz = np.asarray(mz)
# #     intensity = np.asarray(intensity) - encoder.min_intensity
# #     import matplotlib.pyplot as plt
# #     fig,(ax1,ax2) = plt.subplots(ncols=2,nrows=1,figsize=(12,4))
# #     intensity = intensity - intensity.min() + 1
# #     ax1.vlines(mz, 0, intensity, color='blue', lw=2)
# #     ax1.set_xlabel('m/z', fontsize=14)
# #     ax1.set_ylabel('Scaled Intensity', fontsize=14)
# #     # ax1.set_title('Linear Encoding of Message')

# #     for i,c in enumerate(message):
# #         # print a character for each point
# #         ax1.text(mz[i],intensity[i]*1.02,message[i],ha='center', fontsize=14)
# #     ax1.text(0.02, 0.98, 'a', transform=ax1.transAxes, fontsize=24, fontweight='normal', 
# #                 verticalalignment='top')
# #     ax1.set_ylim(0,intensity.max()*1.1)
# #     encoder.convert_text_to_vector(message)
# #     # encoder.convert_text_to_vector('abc')
# #     coords = np.column_stack((encoder.x_indices, encoder.y_indices))
# #     indices = np.argwhere(np.asarray(encoder.encoded_ascii_vector)!=0).flatten()
# #     # # print(coords[indices,:])
# #     mz, intensity = encoder.scale_coords(coords, indices)
# #     intensity = intensity - encoder.min_intensity
# #     ax2.vlines(mz,0,intensity)
# #     ax2.set_xlabel('m/z', fontsize=14)
# #     ax2.set_ylabel('Scaled Intensity', fontsize=14)
# #     for i,c in enumerate(message):
# #         # print a character for each point
# #         ax2.text(mz[i],intensity[i]*1.02,message[i],ha='center', fontsize=14)
# #     ax2.text(0.02, 0.98, 'b', transform=ax2.transAxes, fontsize=24, fontweight='normal', 
# #                 verticalalignment='top')
# #     ax2.set_ylim(0,intensity.max()*1.1)
# #         # increase font size for both axes
# #     ax1.tick_params(axis='both', which='major', labelsize=14)
# #     ax2.tick_params(axis='both', which='major', labelsize=14)

# #     # remove top and right spines
# #     for ax in (ax1, ax2):
# #         ax.spines['top'].set_visible(False)
# #         ax.spines['right'].set_visible(False)

# #     plt.tight_layout()
# #     plt.show()