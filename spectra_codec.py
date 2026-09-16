    
import zlib
import base64
import numpy as np
from lxml import etree
import requests
import pandas as pd
from pyteomics import mzml
import pymzml
import json
import hashlib
import struct
import re
import os

__version__ = "1.0.0"

class SpectraCodec:
    def __init__(self, order=9, encoding_chars=7, min_mz=5, mz_scale=0.0001, 
                 max_mz=12, min_intensity=100, intensity_scale=0.001, max_intensity=1000,
                   repetitions=1,max_order=14):
        """
        Initialize the spectrum_text_encoder class.

        
        Parameters:
        order (int): The order of the Hilbert curve.
        encoding_chars (int): The number of bit positions to use for encoding each character.
        min_mz (float): Minimum m/z value.
        mz_scale (float): Scaling factor for m/z values.
        min_intensity (float): Minimum intensity value.
        intensity_scale (float): Scaling factor for intensity values.
        repetitions (int): Number of times to repeat each character for redundancy.
        """
        self.order = order
        self.encoding_chars = encoding_chars
        self.min_mz = min_mz
        self.max_mz = max_mz
        self.max_intensity = max_intensity
        self.mz_scale = mz_scale
        self.min_intensity = min_intensity
        self.intensity_scale = intensity_scale
        self.repetitions = repetitions
        self.max_order = max_order
        coords = self.generate_full_hilbert_curve(order)
    
        self.base64_alphabet = (
            'ABCDEFGHIJKLMNOPQRSTUVWXYZ' +
            'abcdefghijklmnopqrstuvwxyz' +
            '0123456789+/='
        )
        self.char_to_index = {char: idx for idx, char in enumerate(self.base64_alphabet)}
        
    def determine_hilbert_curve_order(self, message_length):
        """
        Determine the order of the Hilbert curve based on the message length.
        The order is calculated such that the Hilbert curve can accommodate the message length.
        """
        for order in range(self.max_order+1):
            if (message_length+2) < (2**order * 2**order):
                break
        if order != self.order:
            self.order = order
            _ = self.generate_full_hilbert_curve(self.order)
        # print('The order of the Hilbert curve is:', self.order)

    def generate_full_hilbert_curve(self,order):
        """Generates a Hilbert curve of the specified order, starting at (0, 0)."""
        if order == 0:
            return np.array([[0], [0]])

        n = 2**(order -1)
        x, y = self.generate_full_hilbert_curve(order - 1)

        # Quadrant 1
        x1 = y
        y1 = x

        # Quadrant 2
        x2 = x
        y2 = y + n

        # Quadrant 3
        x3 = x + n
        y3 = y + n

        # Quadrant 4
        x4 = n - 1 - y
        y4 = n - 1 - x
        x4 = x4 + n

        # Combine quadrants and return
        # self.x_indices = np.concatenate((x1, x2, x3, x4), axis=1)
        # self.y_indices = np.concatenate((y1, y2, y3, y4), axis=1)
        coords = np.concatenate((np.array([x1, y1]), np.array([x2, y2]), np.array([x3, y3]), np.array([x4, y4])), axis=1)
        self.x_indices = coords[0,:]
        self.y_indices = coords[1,:]
        return coords
    
    def convert_coords_to_chars(self, coords):
        """
        Convert coordinates to characters.

        Parameters:
        coords (numpy.ndarray): The coordinates array.

        Returns:
        str: The decoded text.
        """
        df_full = pd.DataFrame()
        df_full['x'] = self.x_indices
        df_full['y'] = self.y_indices
        df_full.index.name = 'hilbert_index'
        df_full.reset_index(inplace=True,drop=False)
        df_spectra = pd.DataFrame()
        df_spectra['x'] = coords[:,0]
        df_spectra['y'] = coords[:,1]    
        df = df_full.merge(df_spectra, on=['x','y'], how='inner')
        z_decoded = df['hilbert_index'].values
        num_chars = int(np.ceil((z_decoded.max()+1) / self.encoding_chars))
        # num_chars has to be a multiple of 4 so make it that way
        if num_chars % 4 != 0:
            num_chars += 4 - (num_chars % 4)
        length_binary_array = num_chars * self.encoding_chars
        # number of characters * self.encoding_chars has to be the multiple of 7
        if length_binary_array % self.encoding_chars != 0:
            length_binary_array += self.encoding_chars - (length_binary_array % self.encoding_chars)

        decoded_text = np.zeros(length_binary_array)
        decoded_text[z_decoded] = 1
        decoded_text = self.binary_array_to_string(decoded_text)
        return decoded_text
    
    def hilbert_xy_to_index(self, x, y):
        #find the closest x and y in the curve 
        idx = np.argwhere((self.x_indices == x) & (self.y_indices == y)).flatten()
        return idx[0]
    
    def scale_coords(self, coords, idx):
        """
        Scale coordinates to m/z and intensity values.

        Parameters:
        coords (numpy.ndarray): The coordinates array.
        idx (int): The index of the coordinate to scale.

        Returns:
        tuple: The scaled m/z and intensity values.
        """
        mz = coords[idx, 0] * self.mz_scale
        mz = mz + self.min_mz
        intensity = (coords[idx, 1] * self.intensity_scale) + self.min_intensity
        
        return mz, intensity
    
    def make_spectra_to_coords(self, mz, intensity):
        """
        Convert m/z and intensity values to coordinates.

        Parameters:
        mz (numpy.ndarray): The m/z values.
        intensity (numpy.ndarray): The intensity values.

        Returns:
        numpy.ndarray: The coordinates array.
        """

        mz = (mz - self.min_mz) / self.mz_scale
        intensity = (intensity - self.min_intensity) / self.intensity_scale
        mz_rounded = mz.round().astype(int)
        intensity_rounded = intensity.round().astype(int)
        coords = np.column_stack((mz_rounded, intensity_rounded))        
        return coords
    
    def decode_message_from_file(self, input_mzML,parser='pymzml',method='hilbert'):
        """ 
        Decode a message from an mzML file using either pymzml or pyteomics.
        Parameters:
        input_mzML (str): The path to the mzML file.
        parser (str): The parser to use, either 'pymzml' or 'pyteomics'. Default is 'pymzml'.
        method (str): The decoding method to use, either 'hilbert' or 'sequential'. Default is 'hilbert'.
        """
        if parser=='pymzml':
            run = pymzml.run.Reader(input_mzML)
            # run is a generator.  We only want the first spectrum
            for spectrum in run:
                mz = np.array(spectrum.mz)
                intensity = np.array(spectrum.i)

                self.handle_decoding_details(mz, intensity, method=method)
                break

        elif parser=='pyteomics':
            with mzml.read(input_mzML) as reader:
                for spectrum in reader:
                    mz = np.array(spectrum['m/z array'])
                    intensity = np.array(spectrum['intensity array'])

                    self.handle_decoding_details(mz, intensity, method=method)
                    break

        else:
            print("Error: Invalid parser.")
            return None
        return self.decoded_message

    def handle_decoding_details(self,mz, intensity, method='hilbert'):
        # mz and intensity are parallel arrays (one peak per index), so they
        # must be filtered with a single joint mask to keep the pairs aligned
        keep = (mz >= self.min_mz) & (mz <= self.max_mz) & \
               (intensity >= self.min_intensity) & (intensity <= self.max_intensity)
        mz = mz[keep]
        intensity = intensity[keep]
        if len(mz) == 0 or len(intensity) == 0:
            print("Error: No valid message was found.")
            return None
        if method=='hilbert':
            coords = self.make_spectra_to_coords(mz, intensity)
            max_coords = coords.max()
            # for i in range(0,self.max_order+1):
                # my_1d_coords = 2**i
                # if my_1d_coords>max_coords:
                    # num_possible_points = my_1d_coords * my_1d_coords#2**i * 2**i
                    # self.determine_hilbert_curve_order(num_possible_points)
                    # print("Hilbert curve order:", self.order)
                    # break
            self.determine_hilbert_curve_order(max_coords*max_coords)
            message = self.convert_coords_to_chars(coords)
        elif method=='sequential':
            message = self.sequential_decode_message(mz, intensity)
        else:
            print("Error: Invalid method.")
            return None

        self.decoded_message = message

    def encode_message_to_file(self, original_message, input_mzML, output_mzML,method='hilbert',
                               signing_key_path='auto', key_password=None, key_id=None):
        """
        Encode a message to the first spectrum in an mzml file.

        Parameters:
        message (str): The message to encode.
        input_mzML (str): The input mzML file path.
        output_mzML (str): The output mzML file path.
        method (str): The encoding method to use. Default is 'hilbert'. Options are 'hilbert' or 'sequential'.
        signing_key_path: 'auto' (default) signs with $SPECTRACODEC_SIGNING_KEY or
            ~/secrets/spectracodec_private_key if either exists, and encodes
            unsigned otherwise; a path forces signing with that Ed25519 PEM key;
            None disables signing. Signing adds a "provenance" block with
            payload_signature and spectra_signature (see sign_message).

        Returns the exact message string that was encoded.
        """

        if method not in ['hilbert', 'sequential']:
            print("Error: Invalid encoding method. Use 'hilbert' or 'sequential'.")
            return None

        if signing_key_path == 'auto':
            signing_key_path = resolve_signing_key_path()
            if signing_key_path is None:
                print("Note: no signing key found; encoding UNSIGNED. "
                      f"Set ${PRIVATE_KEY_ENV_VAR} or create {DEFAULT_PRIVATE_KEY_PATH} to sign.")
        if signing_key_path is not None:
            signed = sign_message(original_message, input_mzML,
                                  private_key_path=signing_key_path,
                                  key_password=key_password, key_id=key_id)
            original_message = json.dumps(signed, separators=(",", ":"))
        elif not isinstance(original_message, str):
            original_message = json.dumps(original_message, separators=(",", ":"))
        # Convert the message to a one-hot encoded vector
        if method == 'hilbert':
            # determine the order of the Hilbert curve based on the message length
            # Convert the message to a one-hot encoded vector
            message = self.string_to_binary_matrix(original_message)
            self.determine_hilbert_curve_order(len(message))
            # Scale the coordinates to m/z and intensity values
            coords = np.column_stack((self.x_indices, self.y_indices))
            indices = np.argwhere(np.asarray(self.encoded_message_vector)!=0).flatten()
            # print(coords[indices,:])
            mz, intensity = self.scale_coords(coords, indices)
            self.mz = mz
            self.intensity = intensity
        elif method == 'sequential':
            # Sequential encoding
            mz, intensity = self.sequential_encode_message(message)


        # Add the message to the spectrum
        self.add_values_to_first_spectrum(input_mzML, output_mzML, mz, intensity)
        self.rebuild_mzmlfile_index(output_mzML)
        self.update_checksum(output_mzML)
        decoded_message = self.decode_message_from_file(output_mzML, parser='pymzml', method=method)
        # raise error if decoded message is not the same as the original message
        assert decoded_message == original_message, "Decoded message does not match the original message."
        return original_message
    
    def string_to_binary_matrix(self, message):
        """
        Convert a string to a binary matrix representation using base64 encoding.
        Each character in the base64 string is represented by 7 bits.
        """
        message = message.encode('utf-8')
        message = zlib.compress(message)
        message = base64.b64encode(message).decode('ascii')
        
        # Store the original base64 string length for proper reconstruction
        original_length = len(message)
        # Convert each character to its index in the base64 alphabet
        indices = [self.char_to_index[char] for char in message]
        
        # Create a 2D array where each row represents the 7-bit binary sequence
        binary_matrix = np.array([
            [(idx >> i) & 1 for i in range(self.encoding_chars)]  # Gets bits from LSB to MSB
            for idx in indices
        ], dtype=np.uint8)
        
        # Flatten the array to a 1d array
        binary_array = binary_matrix.flatten()
        self.encoded_message_vector = binary_array
        return binary_array

    def binary_array_to_string(self, binary_array):
        """
        Convert a binary array back to the original string.
        """
        # Reshape the 1D array back to 2D where each row is 7 bits
        num_chars = len(binary_array) // self.encoding_chars
        binary_matrix = binary_array.reshape(num_chars, self.encoding_chars)
        
        # Convert each 7-bit row back to its index value
        indices = []
        for row in binary_matrix:
            idx = sum(int(bit) << i for i, bit in enumerate(row))
            indices.append(idx)
        
        # Convert indices back to base64 characters
        base64_chars = [self.base64_alphabet[idx] for idx in indices]
        base64_string = ''.join(base64_chars)
        compressed_bytes = base64.b64decode(base64_string)
        decompressed_bytes = zlib.decompress(compressed_bytes)
        original_string = decompressed_bytes.decode('utf-8')
        return original_string
        
    def make_spectrum_dict(self, mzml_file):
        with open(mzml_file, 'r') as f:
            mzml_text = f.read()

        pat = r'<spectrum index="\d+" id="([^"]+)"'

        # find all occurrences of the pattern in the text
        matches = re.finditer(pat, mzml_text)
        # create a list of dictionaries to hold the scan number, controllerType and byte location
        spectrum_dict = {}
        # loop through the matches and extract the scan number, controllerType and byte location
        for match in matches:
            # get the scan number, controllerType and byte location
            # scan_number = match.group(1)
            scan_id = match.group(1)
            byte_location = match.start()
            # create a dictionary to hold the scan number, controllerType and byte location
            spectrum_dict[scan_id] = int(byte_location)
            
        return spectrum_dict

    def rebuild_mzmlfile_index(self,mzml_file):
        # huge_tree allows text nodes >10MB (large encoded messages)
        tree = etree.parse(mzml_file, parser=etree.XMLParser(huge_tree=True))
        root = tree.getroot()
        ns = {'mzml': 'http://psi.hupo.org/ms/mzml'}

        # Find the indexList element
        index_list = root.find('.//mzml:indexList', namespaces=ns)
        if index_list is not None:
            spectrum_dict = self.make_spectrum_dict(mzml_file)
            # update the text for each offset element
            for index in index_list.findall('.//mzml:index', namespaces=ns):
                for offset in index.findall('.//mzml:offset', namespaces=ns):
                
                    # get the idRef attribute
                    id_ref = offset.get('idRef')
                    # get the controller type and scan number from the idRef attribute
                    if id_ref in spectrum_dict:
                        byte_location = spectrum_dict[id_ref]
                        # update the text for the offset element
                        offset.text = str(byte_location)
                    else:
                        continue
        tree.write(mzml_file, pretty_print=True, xml_declaration=True, encoding="UTF-8")
            
            
    def update_checksum(self,mzml_file):
        """
        Update the checksum of the mzML file.
        """
        # Calculate the SHA-1 checksum
        sha1 = hashlib.sha1()
        with open(mzml_file, 'rb') as f:
            while chunk := f.read(8192):
                sha1.update(chunk)
        checksum = sha1.hexdigest()

        # Parse the mzML file
        tree = etree.parse(mzml_file, parser=etree.XMLParser(huge_tree=True))
        root = tree.getroot()
        ns = {'mzml': 'http://psi.hupo.org/ms/mzml'}

        # Find the fileChecksum element
        file_checksum = root.find('.//mzml:fileChecksum', namespaces=ns)
        if file_checksum is not None:
            file_checksum.text = checksum

        # Save the modified mzML file
        tree.write(mzml_file, pretty_print=True, xml_declaration=True, encoding="UTF-8")

    def add_values_to_first_spectrum(self,input_mzML, output_mzML, new_mz_values, new_intensity_values):
        # Parse the mzML file
        tree = etree.parse(input_mzML, parser=etree.XMLParser(huge_tree=True))
        root = tree.getroot()
        ns = {'mzml': 'http://psi.hupo.org/ms/mzml'}

        # Find the first spectrum
        first_spectrum = root.find('.//mzml:spectrum', namespaces=ns)

        if first_spectrum is not None:
            # Find the m/z array and intensity array binary data elements
            mz_array_binary = first_spectrum.xpath('.//mzml:binaryDataArray[mzml:cvParam[@name="m/z array"]]', namespaces=ns)
            intensity_array_binary = first_spectrum.xpath('.//mzml:binaryDataArray[mzml:cvParam[@name="intensity array"]]', namespaces=ns)

            if mz_array_binary and intensity_array_binary:
                mz_array_binary = mz_array_binary[0]
                intensity_array_binary = intensity_array_binary[0]

                # Determine the precision for m/z and intensity arrays
                mz_precision = mz_array_binary.xpath('.//mzml:cvParam[@name="64-bit float"]', namespaces=ns)
                if mz_precision:
                    mz_format = 'd'  # 64-bit float
                else:
                    mz_format = 'f'  # 32-bit float

                intensity_precision = intensity_array_binary.xpath('.//mzml:cvParam[@name="64-bit float"]', namespaces=ns)
                if intensity_precision:
                    intensity_format = 'd'  # 64-bit float
                else:
                    intensity_format = 'f'  # 32-bit float
                    
                # Check for compression
                is_mz_compressed = len(mz_array_binary.xpath('.//mzml:cvParam[@name="zlib compression"]', namespaces=ns)) > 0
                is_intensity_compressed = len(intensity_array_binary.xpath('.//mzml:cvParam[@name="zlib compression"]', namespaces=ns)) > 0

                # Find the binary elements
                mz_binary = mz_array_binary.find('.//mzml:binary', namespaces=ns)
                intensity_binary = intensity_array_binary.find('.//mzml:binary', namespaces=ns)

                if mz_binary is not None and intensity_binary is not None:
                    # Decode the existing data
                    existing_mz_data = base64.b64decode(mz_binary.text)
                    existing_intensity_data = base64.b64decode(intensity_binary.text)
                    
                    # Decompress if necessary
                    if is_mz_compressed:
                        existing_mz_data = zlib.decompress(existing_mz_data)
                    if is_intensity_compressed:
                        existing_intensity_data = zlib.decompress(existing_intensity_data)

                    # Unpack the binary data
                    if mz_format == 'f':
                        existing_mz_values = struct.unpack('<' + 'f' * (len(existing_mz_data) // 4), existing_mz_data)
                    else:
                        existing_mz_values = struct.unpack('<' + 'd' * (len(existing_mz_data) // 8), existing_mz_data)

                    if intensity_format == 'f':
                        existing_intensity_values = struct.unpack('<' + 'f' * (len(existing_intensity_data) // 4), existing_intensity_data)
                    else:
                        existing_intensity_values = struct.unpack('<' + 'd' * (len(existing_intensity_data) // 8), existing_intensity_data)
                    
                    # Combine the existing data with the new data
                    combined_mz_values = list(existing_mz_values) + list(new_mz_values)
                    combined_intensity_values = list(existing_intensity_values) + list(new_intensity_values)

                    # Sort by m/z values
                    combined_mz_values = np.array(combined_mz_values)
                    combined_intensity_values = np.array(combined_intensity_values)
                    idx = np.argsort(combined_mz_values)
                    combined_mz_values = combined_mz_values[idx].tolist()
                    combined_intensity_values = combined_intensity_values[idx].tolist()
                    
                    # Update the defaultArrayLength attribute
                    first_spectrum.set('defaultArrayLength', str(len(combined_mz_values)))

                    # Encode the data with or without compression
                    if is_mz_compressed:
                        packed_mz_data = struct.pack('<' + mz_format * len(combined_mz_values), *combined_mz_values)
                        compressed_mz_data = zlib.compress(packed_mz_data)
                        combined_mz_encoded = base64.b64encode(compressed_mz_data).decode('ascii')
                    else:
                        packed_mz_data = struct.pack('<' + mz_format * len(combined_mz_values), *combined_mz_values)
                        combined_mz_encoded = base64.b64encode(packed_mz_data).decode('ascii')
                    
                    if is_intensity_compressed:
                        packed_intensity_data = struct.pack('<' + intensity_format * len(combined_intensity_values), *combined_intensity_values)
                        compressed_intensity_data = zlib.compress(packed_intensity_data)
                        combined_intensity_encoded = base64.b64encode(compressed_intensity_data).decode('ascii')
                    else:
                        packed_intensity_data = struct.pack('<' + intensity_format * len(combined_intensity_values), *combined_intensity_values)
                        combined_intensity_encoded = base64.b64encode(packed_intensity_data).decode('ascii')

                    # Update the binary elements with the combined encoded data
                    mz_binary.text = combined_mz_encoded
                    intensity_binary.text = combined_intensity_encoded
                    
                    # CRITICAL FIX: Update the encodedLength attributes
                    if 'encodedLength' in mz_array_binary.attrib:
                        mz_array_binary.set('encodedLength', str(len(combined_mz_encoded)))
                    if 'encodedLength' in intensity_array_binary.attrib:
                        intensity_array_binary.set('encodedLength', str(len(combined_intensity_encoded)))
                    
                    # Also check for arrayLength attributes
                    if 'arrayLength' in mz_array_binary.attrib:
                        mz_array_binary.set('arrayLength', str(len(combined_mz_values)))
                    if 'arrayLength' in intensity_array_binary.attrib:
                        intensity_array_binary.set('arrayLength', str(len(combined_intensity_values)))
        
       
        
        # Save the modified mzML file do not proceed until write is complete

        # Save the modified mzML file
        tree.write(output_mzML, pretty_print=True, xml_declaration=True, encoding="UTF-8")


# ----------------------------------------------------------------------------
# Signing / verification (tamper evidence)
#
# Two Ed25519 signatures are embedded in every signed message, in a
# top-level "provenance" block:
#   payload_signature  - over the canonical JSON of the whole message minus
#                        "provenance" (i.e. unique_file_id + payload):
#                        proves the embedded metadata/documents are authentic.
#   spectra_signature  - over unique_file_id + the spectral digest (scd-1):
#                        proves the acquired spectra are untampered.
# The spectral digest covers every spectrum EXCEPT the first (the carrier
# that holds the encoded message), so it is identical before and after
# encoding - that is what makes signing-then-embedding possible.
# ----------------------------------------------------------------------------

SIGNING_ALGORITHM = "Ed25519"
SPECTRAL_DIGEST_SPEC = "scd-1"
DEFAULT_PRIVATE_KEY_PATH = "~/secrets/spectracodec_private_key"
PRIVATE_KEY_ENV_VAR = "SPECTRACODEC_SIGNING_KEY"
KEY_PASSWORD_ENV_VAR = "SPECTRACODEC_KEY_PASSWORD"
_PAYLOAD_SIG_CONTEXT = b"spectracodec-payload-v1:"
_SPECTRA_SIG_CONTEXT = b"spectracodec-spectra-v1:"


def canonical_json_bytes(obj):
    """Deterministic JSON serialization used for everything that gets signed."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":")).encode("utf-8")


def compute_spectral_digest(mzml_file, spec=SPECTRAL_DIGEST_SPEC):
    """
    Canonical digest of the acquired spectral data (scd-1).

    Chained SHA-256 over every spectrum EXCEPT the first (the SpectraCodec
    carrier), in file order. Per spectrum it hashes: position among the
    digested spectra, native spectrum ID, peak count, and the decoded m/z and
    intensity arrays as little-endian float64 - so the digest is invariant to
    re-compression / base64 details, and to anything done to spectrum 1.
    """
    if spec != SPECTRAL_DIGEST_SPEC:
        raise ValueError(f"Unknown spectral digest spec: {spec}")
    h = hashlib.sha256()
    h.update(spec.encode("utf-8") + b"\x00")
    run = pymzml.run.Reader(mzml_file)
    n = 0
    for i, spectrum in enumerate(run):
        if i == 0:  # carrier spectrum: excluded by design
            continue
        mz = np.asarray(spectrum.mz, dtype="<f8")
        intensity = np.asarray(spectrum.i, dtype="<f8")
        h.update(struct.pack(">Q", n))
        h.update(str(spectrum.ID).encode("utf-8") + b"\x00")
        h.update(struct.pack(">Q", len(mz)))
        h.update(mz.tobytes())
        h.update(intensity.tobytes())
        n += 1
    h.update(b"end" + struct.pack(">Q", n))
    return h.hexdigest()


def key_fingerprint(verification_key):
    """SHA256 fingerprint (base64, unpadded) of the raw 32-byte public key."""
    from cryptography.hazmat.primitives import serialization
    raw = verification_key.public_bytes(
        serialization.Encoding.Raw, serialization.PublicFormat.Raw)
    return "SHA256:" + base64.b64encode(hashlib.sha256(raw).digest()).decode().rstrip("=")


def generate_signing_keypair(private_key_path, verification_key_path, password=None):
    """Mint a new Ed25519 keypair (PEM). Returns the fingerprint."""
    from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey
    from cryptography.hazmat.primitives import serialization
    for p in (private_key_path, verification_key_path):
        if os.path.exists(os.path.expanduser(p)):
            raise FileExistsError(f"refusing to overwrite existing key: {p}")
    priv = Ed25519PrivateKey.generate()
    enc = (serialization.BestAvailableEncryption(password.encode())
           if password else serialization.NoEncryption())
    fd = os.open(os.path.expanduser(private_key_path),
                 os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    with os.fdopen(fd, "wb") as f:
        f.write(priv.private_bytes(serialization.Encoding.PEM,
                                   serialization.PrivateFormat.PKCS8, enc))
    with open(os.path.expanduser(verification_key_path), "wb") as f:
        f.write(priv.public_key().public_bytes(
            serialization.Encoding.PEM, serialization.PublicFormat.SubjectPublicKeyInfo))
    return key_fingerprint(priv.public_key())


def resolve_signing_key_path():
    """Locate a signing key: $SPECTRACODEC_SIGNING_KEY, else ~/secrets/spectracodec_private_key."""
    env = os.environ.get(PRIVATE_KEY_ENV_VAR)
    if env and os.path.exists(os.path.expanduser(env)):
        return os.path.expanduser(env)
    default = os.path.expanduser(DEFAULT_PRIVATE_KEY_PATH)
    if os.path.exists(default):
        return default
    return None


def load_private_key(path, password=None):
    from cryptography.hazmat.primitives import serialization
    data = open(os.path.expanduser(path), "rb").read()
    if password is None:
        password = os.environ.get(KEY_PASSWORD_ENV_VAR)
    pw = password.encode() if isinstance(password, str) else password
    try:
        return serialization.load_pem_private_key(data, password=pw)
    except TypeError:
        if pw is not None:  # password given but key is unencrypted
            return serialization.load_pem_private_key(data, password=None)
        raise ValueError(
            f"Private key {path} is encrypted: pass key_password= or set "
            f"${KEY_PASSWORD_ENV_VAR}.")


def load_verification_key(path_or_pem):
    """Load a public key from a PEM file path or an in-message PEM string."""
    from cryptography.hazmat.primitives import serialization
    if isinstance(path_or_pem, str) and "BEGIN PUBLIC KEY" in path_or_pem:
        return serialization.load_pem_public_key(path_or_pem.encode("utf-8"))
    return serialization.load_pem_public_key(
        open(os.path.expanduser(path_or_pem), "rb").read())


def sign_message(message, input_mzML, private_key_path=None, key_password=None,
                 key_id=None):
    """
    Add a "provenance" block with the two signatures to a message dict.

    The spectral digest is computed from input_mzML (the file about to be
    encoded); because scd-1 excludes the carrier spectrum, the digest of the
    encoded output will match. Returns a new dict with provenance placed
    right after unique_file_id.
    """
    from cryptography.hazmat.primitives import serialization
    if isinstance(message, str):
        message = json.loads(message)
    if not isinstance(message, dict):
        raise ValueError("Signing requires the message to be a JSON object.")
    if private_key_path is None:
        private_key_path = resolve_signing_key_path()
        if private_key_path is None:
            raise FileNotFoundError(
                f"No signing key found (${PRIVATE_KEY_ENV_VAR} or "
                f"{DEFAULT_PRIVATE_KEY_PATH}).")
    priv = load_private_key(private_key_path, key_password)
    pub = priv.public_key()
    fingerprint = key_fingerprint(pub)

    unsigned = {k: v for k, v in message.items() if k != "provenance"}
    digest = {"spec": SPECTRAL_DIGEST_SPEC, "algorithm": "sha256",
              "value": compute_spectral_digest(input_mzML)}
    payload_sig = priv.sign(_PAYLOAD_SIG_CONTEXT + canonical_json_bytes(unsigned))
    spectra_sig = priv.sign(_SPECTRA_SIG_CONTEXT + canonical_json_bytes(
        {"unique_file_id": unsigned.get("unique_file_id"),
         "spectral_digest": digest}))

    from datetime import datetime, timezone
    provenance = {
        "signing_scheme_version": "1.0",
        "algorithm": SIGNING_ALGORITHM,
        "key_id": key_id if key_id is not None else fingerprint,
        "verification_key": pub.public_bytes(
            serialization.Encoding.PEM,
            serialization.PublicFormat.SubjectPublicKeyInfo).decode("utf-8"),
        "verification_key_fingerprint": fingerprint,
        "signed_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "spectral_digest": digest,
        "payload_signature": base64.b64encode(payload_sig).decode(),
        "spectra_signature": base64.b64encode(spectra_sig).decode(),
    }
    signed = {}
    if "unique_file_id" in unsigned:
        signed["unique_file_id"] = unsigned["unique_file_id"]
    signed["provenance"] = provenance
    for k, v in unsigned.items():
        signed.setdefault(k, v)
    return signed


def verify_signed_file(mzml_file, verification_key_path=None):
    """
    Decode an encoded mzML and check both signatures.

    verification_key_path: PEM file of the key you TRUST (fetched from the
    lab's published location). If omitted, the key embedded in the message is
    used - that only proves internal consistency, not origin, and the report
    says so via trusted_key_source == "embedded_untrusted".

    Returns a report dict; report["valid"] is the overall verdict.
    """
    from cryptography.exceptions import InvalidSignature
    raw = SpectraCodec().decode_message_from_file(mzml_file)
    try:
        message = json.loads(raw)
    except (TypeError, ValueError):
        return {"signed": False, "valid": False,
                "error": "decoded message is not JSON", "message": raw}
    if not isinstance(message, dict) or "provenance" not in message:
        return {"signed": False, "valid": False,
                "error": "message carries no provenance block (unsigned or legacy)",
                "message": message}

    prov = message["provenance"]
    if verification_key_path is not None:
        pub = load_verification_key(verification_key_path)
        trusted_source = "user_supplied"
    else:
        pub = load_verification_key(prov["verification_key"])
        trusted_source = "embedded_untrusted"
    fingerprint = key_fingerprint(pub)
    checks = {}
    checks["verification_key_matches_supplied"] = (
        fingerprint == prov.get("verification_key_fingerprint")
        if verification_key_path is not None else None)

    unsigned = {k: v for k, v in message.items() if k != "provenance"}
    try:
        pub.verify(base64.b64decode(prov["payload_signature"]),
                   _PAYLOAD_SIG_CONTEXT + canonical_json_bytes(unsigned))
        checks["payload_signature"] = True
    except InvalidSignature:
        checks["payload_signature"] = False

    computed = compute_spectral_digest(mzml_file, prov["spectral_digest"]["spec"])
    checks["spectral_digest"] = (computed == prov["spectral_digest"]["value"])
    try:
        pub.verify(base64.b64decode(prov["spectra_signature"]),
                   _SPECTRA_SIG_CONTEXT + canonical_json_bytes(
                       {"unique_file_id": unsigned.get("unique_file_id"),
                        "spectral_digest": prov["spectral_digest"]}))
        checks["spectra_signature"] = True
    except InvalidSignature:
        checks["spectra_signature"] = False

    return {
        "signed": True,
        "valid": all(v for v in checks.values() if v is not None),
        "checks": checks,
        "trusted_key_source": trusted_source,
        "verification_key_fingerprint": fingerprint,
        "key_id": prov.get("key_id"),
        "signed_at": prov.get("signed_at"),
        "spectral_digest_expected": prov["spectral_digest"]["value"],
        "spectral_digest_computed": computed,
        "unique_file_id": message.get("unique_file_id"),
        "message": message,
    }


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="SpectraCodec signing tools")
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("verify", help="verify the signatures in an encoded mzML")
    v.add_argument("mzml")
    v.add_argument("--key", default=None,
                   help="trusted verification key (PEM); omit to use the embedded key")
    k = sub.add_parser("keygen", help="mint a new Ed25519 signing keypair")
    k.add_argument("private_key_path")
    k.add_argument("verification_key_path")
    args = ap.parse_args()
    if args.cmd == "keygen":
        fp = generate_signing_keypair(args.private_key_path, args.verification_key_path)
        print(f"wrote {args.private_key_path} and {args.verification_key_path}")
        print(f"fingerprint: {fp}")
    elif args.cmd == "verify":
        report = verify_signed_file(args.mzml, verification_key_path=args.key)
        for field in ("signed", "valid", "trusted_key_source", "key_id",
                      "signed_at", "unique_file_id", "error"):
            if field in report:
                print(f"{field}: {report[field]}")
        for name, ok in report.get("checks", {}).items():
            print(f"  {name}: {'PASS' if ok else 'n/a' if ok is None else 'FAIL'}")
        if report.get("trusted_key_source") == "embedded_untrusted":
            print("WARNING: verified against the key embedded in the file itself; "
                  "this proves internal consistency, NOT origin. Pass --key with "
                  "the published verification key to prove origin.")
        raise SystemExit(0 if report.get("valid") else 1)



# import numpy as np
# import matplotlib.pyplot as plt
# import base64
# import re

# import struct
# import zlib
# from lxml import etree
# import requests
# import pandas as pd
# from pyteomics import mzml
# import pymzml
# import json
# import hashlib

# class SpectraCodec:
#     def __init__(self, order=10, ascii_chars=128, min_mz=5, mz_scale=0.0001, 
#                  max_mz=12, min_intensity=100, intensity_scale=0.001, max_intensity=1000,
#                    repetitions=1):
#         """
#         Initialize the spectrum_text_encoder class.
#         Critical to note that you really do not want to change the order of the 
#         Hilbert curve, as this will change the encoding and decoding of the messages.
        
#         Parameters:
#         order (int): The order of the Hilbert curve.
#         ascii_chars (int): The number of ASCII characters to consider.
#         min_mz (float): Minimum m/z value.
#         mz_scale (float): Scaling factor for m/z values.
#         min_intensity (float): Minimum intensity value.
#         intensity_scale (float): Scaling factor for intensity values.
#         repetitions (int): Number of times to repeat each character for redundancy.
#         """
#         self.order = order
#         self.ascii_chars = ascii_chars
#         self.min_mz = min_mz
#         self.max_mz = max_mz
#         self.max_intensity = max_intensity
#         self.mz_scale = mz_scale
#         self.min_intensity = min_intensity
#         self.intensity_scale = intensity_scale
#         self.repetitions = repetitions
#         coords = self.generate_full_hilbert_curve(order)


#     def determine_hilbert_curve_order(self, message_length):
#         """
#         Determine the order of the Hilbert curve based on the message length.
#         The order is calculated such that the Hilbert curve can accommodate the message length.
#         """

#         order  = int(np.ceil(np.emath.logn(4,message_length*128)))
#         self.order = order
#         print(f"Determined Hilbert curve order: {self.order} for message length: {message_length}")
#         _ = self.generate_full_hilbert_curve(self.order)

        

#     def sequential_encode_message(self,message):
#         byte_locations = [ord(char) for char in message if char.isascii() and char.isprintable()]
#         mz_locations = [self.min_mz + (i * self.mz_scale) for i,loc in enumerate(byte_locations)]
#         intensity_locations = [self.min_intensity + loc for loc in byte_locations]
#         return mz_locations, intensity_locations

#     def sequential_decode_message(self, mz, intensity):
#         """
#         Decode a message from m/z and intensity values.

#         Parameters:
#         mz (numpy.ndarray): The m/z values.
#         intensity (numpy.ndarray): The intensity values.

#         Returns:
#         str: The decoded message.
#         """
#         # Convert m/z and intensity to byte locations
#         byte_locations = [(mz[i] - self.min_mz) / self.mz_scale for i in range(len(mz))]
#         byte_locations = np.asarray(byte_locations)
#         byte_locations = np.round(byte_locations).astype(int)
#         loc_diff = np.diff(byte_locations)
#         # find the first difference that is greater than 1
#         first_diff = np.where(loc_diff > 1)[0]
#         # trim intensity to the first difference
#         if len(first_diff) > 0:
#             intensity = intensity[:first_diff[0]+1]
#             byte_locations = byte_locations[:first_diff[0]+1]
#         intensity_to_ascii = [(intensity[i] - self.min_intensity) for i in range(len(intensity))]
#         ascii_to_char = [chr(int(loc)) for loc in intensity_to_ascii]
#         decoded_message = ''.join(ascii_to_char)
#         return decoded_message

#     def generate_full_hilbert_curve(self,order):
#         """Generates a Hilbert curve of the specified order, starting at (0, 0)."""
#         if order == 0:
#             return np.array([[0], [0]])

#         n = 2**(order -1)
#         x, y = self.generate_full_hilbert_curve(order - 1)

#         # Quadrant 1
#         x1 = y
#         y1 = x

#         # Quadrant 2
#         x2 = x
#         y2 = y + n

#         # Quadrant 3
#         x3 = x + n
#         y3 = y + n

#         # Quadrant 4
#         x4 = n - 1 - y
#         y4 = n - 1 - x
#         x4 = x4 + n

#         # Combine quadrants and return
#         # self.x_indices = np.concatenate((x1, x2, x3, x4), axis=1)
#         # self.y_indices = np.concatenate((y1, y2, y3, y4), axis=1)
#         coords = np.concatenate((np.array([x1, y1]), np.array([x2, y2]), np.array([x3, y3]), np.array([x4, y4])), axis=1)
#         self.x_indices = coords[0,:]
#         self.y_indices = coords[1,:]
#         return coords
    
#     def convert_text_to_vector(self, text):
#         """
#         Convert text to a one-hot encoded vector.

#         Parameters:
#         text (str): The input text to be converted.

#         Returns:
#         numpy.ndarray: The one-hot encoded vector.
#         """
#         my_text = text.encode('ascii', 'replace').decode('ascii')
#         my_list = [ord(c) for c in my_text]
#         if self.repetitions > 1:
#             my_redundant_list = self.interleaved_redundancy(my_list)
#             z = np.zeros((len(my_redundant_list), self.ascii_chars))
#         else:
#             my_redundant_list = my_list
#             z = np.zeros((len(my_list), self.ascii_chars))
#         # Track space runs in the redundant list
#         in_space_run = False
#         for i, c in enumerate(my_redundant_list):
#             if c == 32:  # ASCII space
#                 if not in_space_run:
#                     # Only encode the first space in a run
#                     z[i, c] = 1
#                     in_space_run = True
#             else:
#                 # Always encode non-space characters
#                 z[i, c] = 1
#                 in_space_run = False
#         z = z.reshape(-1)
#         # z = list(z)
#         # tile z so that it fills the entire Hilbert curve
#         self.encoded_ascii_vector = z
#         if len(self.encoded_ascii_vector) > len(self.x_indices):
#             print("Error: Text is too long to encode.")
#             return None
    
#         return z
    
#     def interleaved_redundancy(self,data):
#         interleaved = []
#         for i in range(len(data)):
#             for _ in range(self.repetitions):
#                 interleaved.append(data[i])
#         return interleaved

#     def deinterleave_robust(self,data):
#         expected_length = len(data) // self.repetitions * self.repetitions #Expected length without spurious signals
#         if len(data) != expected_length:
#             print("Warning: Spurious signals detected. Attempting correction...")

#         deinterleaved = []
#         for i in range(0, expected_length, self.repetitions):  #Only process expected length
#             counts = {}
#             for j in range(i, i + self.repetitions):
#                 val = data[j]
#                 counts[val] = counts.get(val, 0) + 1
#             deinterleaved.append(max(counts, key=counts.get))
#         return deinterleaved
    


#     def scale_coords(self, coords, idx):
#         """
#         Scale coordinates to m/z and intensity values.

#         Parameters:
#         coords (numpy.ndarray): The coordinates array.
#         idx (int): The index of the coordinate to scale.

#         Returns:
#         tuple: The scaled m/z and intensity values.
#         """
#         mz = coords[idx, 0] * self.mz_scale
#         mz = mz + self.min_mz
#         intensity = (coords[idx, 1] * self.intensity_scale) + self.min_intensity

#         # mz_spike = []
#         # intensity_spike = []
#         # for i in range(self.len_spike):
#         #     mz_spike.append(self.min_mz - 0.001 - (i * 0.001))
#         #     intensity_spike.append(self.min_intensity+1)

#         # mz = np.concatenate([mz_spike, mz])
#         # intensity = np.concatenate([intensity_spike, intensity])
        
#         return mz, intensity

#     def make_spectra_to_coords(self, mz, intensity):
#         """
#         Convert m/z and intensity values to coordinates.

#         Parameters:
#         mz (numpy.ndarray): The m/z values.
#         intensity (numpy.ndarray): The intensity values.

#         Returns:
#         numpy.ndarray: The coordinates array.
#         """
#         # mz = mz[self.len_spike:]
#         # intensity = intensity[self.len_spike:]
#         mz = (mz - self.min_mz) / self.mz_scale
#         intensity = (intensity - self.min_intensity) / self.intensity_scale
#         mz_rounded = mz.round().astype(int)
#         intensity_rounded = intensity.round().astype(int)
#         mz_diff = abs(mz - mz.round(4))
#         intensity_diff = abs(intensity - intensity.round(4))
#         diff_count = sum((mz_diff<(10*self.mz_scale)) & (intensity_diff<(10*self.intensity_scale)))
#         # set the hilbert curve order based on the number of differences
#         if diff_count > 0:
#             self.determine_hilbert_curve_order(diff_count)
        
#         # idx = np.argwhere(((mz - mz_rounded)<0.1) & ((intensity - intensity_rounded)<0.1) & (mz_rounded >= 0) & (intensity_rounded >= 0)).flatten()
#         # print(mz.shape,intensity.shape,mz_rounded.shape,intensity_rounded.shape,mz_diff.shape,intensity_diff.shape)
#         # coords = np.column_stack((mz_rounded[idx], intensity_rounded[idx], mz_diff[idx], intensity_diff[idx]))
        
#         coords = np.column_stack((mz_rounded, intensity_rounded))
#         diffs = np.column_stack((mz_diff, intensity_diff))
        
#         return coords,diffs

#     def convert_coords_to_chars(self, coords):
#         """
#         Convert coordinates to characters.

#         Parameters:
#         coords (numpy.ndarray): The coordinates array.

#         Returns:
#         str: The decoded text.
#         """
#         df_full = pd.DataFrame()
#         df_full['x'] = self.x_indices
#         df_full['y'] = self.y_indices
#         df_full.index.name = 'hilbert_index'
#         df_full.reset_index(inplace=True,drop=False)
#         df_spectra = pd.DataFrame()
#         df_spectra['x'] = coords[:,0]
#         df_spectra['y'] = coords[:,1]    
#         df = df_full.merge(df_spectra, on=['x','y'], how='inner')
#         z_decoded = df['hilbert_index'].values
    
        
#         num_chars = int(np.ceil(z_decoded.max() / self.ascii_chars))

#         decoded_text = np.zeros(num_chars * self.ascii_chars)
#         decoded_text[z_decoded] = 1
#         decoded_text = decoded_text.reshape(num_chars, self.ascii_chars)
#         # decoded_text = self.deinterleave_robust(list(decoded_text))
#         # print(np.argwhere(decoded_text))
#         decoded_text = [chr(np.argmax(c)) for c in decoded_text]
#         decoded_text = self.deinterleave_robust(decoded_text)
#         decoded_text = "".join(decoded_text).replace('\x00', '')

#         return decoded_text
    
#     def hilbert_xy_to_index(self, x, y):
#         #find the closest x and y in the curve 
#         idx = np.argwhere((self.x_indices == x) & (self.y_indices == y)).flatten()
#         return idx[0]

#     def filename_to_dict(self, filename):
#         """
#         Convert a filename to a dictionary.

#         Parameters:
#         filename (str): The filename to convert.

#         Returns:
#         dict: The dictionary.
#         """
#         FIELDS = [
#             {"code": "DATE", "label": "Date"},
#             {"code": "NORTHENLABINITIALS", "label": "Initials of LCMS user"},
#             {"code": "COLLABINITIALS", "label": "Initials of sample sumbitter"},
#             {"code": "PROJ", "label": "Project"},
#             {"code": "EXP", "label": "Experiment"},
#             {"code": "SAMPSET", "label": "Sample set"},
#             {"code": "SYSTEM", "label": "System"},
#             {"code": "COLUMN-method", "label": "Chromatography and optional method parameters"},
#             {"code": "SERIAL", "label": "Column serial number"},
#             {"code": "POL", "label": "Polarity"},
#             {"code": "ACQ", "label": "Acquisition type"},
#             {"code": "SAMPLE#", "label": "Sample number"},
#             {"code": "SAMPLEGROUP", "label": "Sample group"},
#             {"code": "REP", "label": "Replicate number"},
#             {"code": "OPTIONAL", "label": "Additonal parameters"},
#             {"code": "SEQ", "label": "Sequence injection #"},
#         ]
#         filename = filename.split("_")
#         d = {}
#         for i, field in enumerate(FIELDS):
#             d[field["label"]] = filename[i]
#         return d
    
# def fetch_spectrum_from_gnps2(self,usi, mz_min=None, mz_max=None, max_intensity=None, 
#                             grid=True, fragment_mz_tolerance=0.1):
#     """
#     Fetch spectrum data from GNPS2 with additional parameters
#     """
#     base_url = "https://metabolomics-usi.gnps2.org/json/"
    
#     params = {
#         "usi1": usi,
#     }
    
#     # Filter out None values
#     params = {k: v for k, v in params.items() if v is not None}
    
#     response = requests.get(base_url, params=params)
    
#     if response.status_code == 200:
#         return response.json()
#     else:
#         print(f"Error: {response.status_code}")
#         print(response.text)
#         return None

#     def decode_message_from_gnps2_spectrum(self,usi):
#         """
#         Decode a message from a GNPS2 spectrum.
#         """
#         data = self.fetch_spectrum_from_gnps2(usi)
#         if data is None:
#             return None
        
#         mz = np.array(data['peaks'])[:,0]
#         intensity = np.array(data['peaks'])[:,1]
        
#         coords, diffs = self.make_spectra_to_coords(mz, intensity)
#         message = self.convert_coords_to_chars(coords)
#         return message

#     def decode_message_from_file(self, input_mzML,parser='pymzml',method='hilbert'):
#         """ 
#         Decode a message from an mzML file using either pymzml or pyteomics.
#         Parameters:
#         input_mzML (str): The path to the mzML file.
#         parser (str): The parser to use, either 'pymzml' or 'pyteomics'. Default is 'pymzml'.
#         method (str): The decoding method to use, either 'hilbert' or 'sequential'. Default is 'hilbert'.
#         """
#         if parser=='pymzml':
#             run = pymzml.run.Reader(input_mzML)
#             # run is a generator.  We only want the first spectrum
#             for spectrum in run:
#                 mz = np.array(spectrum.mz)
#                 intensity = np.array(spectrum.i)

#                 self.handle_decoding_details(mz, intensity, method=method)
#                 break

#         elif parser=='pyteomics':
#             with mzml.read(input_mzML) as reader:
#                 for spectrum in reader:
#                     mz = np.array(spectrum['m/z array'])
#                     intensity = np.array(spectrum['intensity array'])

#                     self.handle_decoding_details(mz, intensity, method=method)
#                     break

#         else:
#             print("Error: Invalid parser.")
#             return None

#     def handle_decoding_details(self,mz, intensity, method='hilbert'):
#         mz = mz[mz >= self.min_mz]
#         mz = mz[mz <= self.max_mz]
#         intensity = intensity[intensity >= self.min_intensity]
#         intensity = intensity[intensity <= self.max_intensity]    
#         if len(mz) == 0 or len(intensity) == 0:
#             print("Error: No valid message was found.")
#             return None
#         if method=='hilbert':
#             coords, diffs = self.make_spectra_to_coords(mz, intensity)
#             message = self.convert_coords_to_chars(coords)
#         elif method=='sequential':
#             message = self.sequential_decode_message(mz, intensity)
#         else:
#             print("Error: Invalid method.")
#             return None
#         # check if the message is a JSON string
#         try:
#             message_dict = json.loads(message)
#         except json.JSONDecodeError:
#             message_dict = message

#         # If the message is a dictionary, decode the manuscript
#         if isinstance(message_dict, dict):
#             # Decode from base64 and then decompress to get the original text
#             if 'manuscript_compressed_b64' in message_dict:
#                 compressed_bytes = base64.b64decode(message_dict['manuscript_compressed_b64'])
#                 original_manuscript_bytes = zlib.decompress(compressed_bytes)
#                 message_dict['manuscript'] = original_manuscript_bytes.decode('utf-8')

#         self.decoded_message = message_dict


#     def make_spectrum_dict(self, mzml_file):
#         with open(mzml_file, 'r') as f:
#             mzml_text = f.read()

#         pat = r'<spectrum index="\d+" id="([^"]+)"'

#         # find all occurrences of the pattern in the text
#         matches = re.finditer(pat, mzml_text)
#         # create a list of dictionaries to hold the scan number, controllerType and byte location
#         spectrum_dict = {}
#         # loop through the matches and extract the scan number, controllerType and byte location
#         for match in matches:
#             # get the scan number, controllerType and byte location
#             # scan_number = match.group(1)
#             scan_id = match.group(1)
#             byte_location = match.start()
#             # create a dictionary to hold the scan number, controllerType and byte location
#             spectrum_dict[scan_id] = int(byte_location)
            
#         return spectrum_dict

#     def rebuild_mzmlfile_index(self,mzml_file):
#         tree = etree.parse(mzml_file)
#         root = tree.getroot()
#         ns = {'mzml': 'http://psi.hupo.org/ms/mzml'}

#         # Find the indexList element
#         index_list = root.find('.//mzml:indexList', namespaces=ns)
#         if index_list is not None:
#             spectrum_dict = self.make_spectrum_dict(mzml_file)
#             # update the text for each offset element
#             for index in index_list.findall('.//mzml:index', namespaces=ns):
#                 for offset in index.findall('.//mzml:offset', namespaces=ns):
                
#                     # get the idRef attribute
#                     id_ref = offset.get('idRef')
#                     # get the controller type and scan number from the idRef attribute
#                     if id_ref in spectrum_dict:
#                         byte_location = spectrum_dict[id_ref]
#                         # update the text for the offset element
#                         offset.text = str(byte_location)
#                     else:
#                         continue
#         tree.write(mzml_file, pretty_print=True, xml_declaration=True, encoding="UTF-8")
            
            
#     def update_checksum(self,mzml_file):
#         """
#         Update the checksum of the mzML file.
#         """
#         # Calculate the SHA-1 checksum
#         sha1 = hashlib.sha1()
#         with open(mzml_file, 'rb') as f:
#             while chunk := f.read(8192):
#                 sha1.update(chunk)
#         checksum = sha1.hexdigest()

#         # Parse the mzML file
#         tree = etree.parse(mzml_file)
#         root = tree.getroot()
#         ns = {'mzml': 'http://psi.hupo.org/ms/mzml'}

#         # Find the fileChecksum element
#         file_checksum = root.find('.//mzml:fileChecksum', namespaces=ns)
#         if file_checksum is not None:
#             file_checksum.text = checksum

#         # Save the modified mzML file
#         tree.write(mzml_file, pretty_print=True, xml_declaration=True, encoding="UTF-8")

#     def add_values_to_first_spectrum(self,input_mzML, output_mzML, new_mz_values, new_intensity_values):
#         # Parse the mzML file
#         tree = etree.parse(input_mzML)
#         root = tree.getroot()
#         ns = {'mzml': 'http://psi.hupo.org/ms/mzml'}

#         # Find the first spectrum
#         first_spectrum = root.find('.//mzml:spectrum', namespaces=ns)

#         if first_spectrum is not None:
#             # Find the m/z array and intensity array binary data elements
#             mz_array_binary = first_spectrum.xpath('.//mzml:binaryDataArray[mzml:cvParam[@name="m/z array"]]', namespaces=ns)
#             intensity_array_binary = first_spectrum.xpath('.//mzml:binaryDataArray[mzml:cvParam[@name="intensity array"]]', namespaces=ns)

#             if mz_array_binary and intensity_array_binary:
#                 mz_array_binary = mz_array_binary[0]
#                 intensity_array_binary = intensity_array_binary[0]

#                 # Determine the precision for m/z and intensity arrays
#                 mz_precision = mz_array_binary.xpath('.//mzml:cvParam[@name="64-bit float"]', namespaces=ns)
#                 if mz_precision:
#                     mz_format = 'd'  # 64-bit float
#                 else:
#                     mz_format = 'f'  # 32-bit float

#                 intensity_precision = intensity_array_binary.xpath('.//mzml:cvParam[@name="64-bit float"]', namespaces=ns)
#                 if intensity_precision:
#                     intensity_format = 'd'  # 64-bit float
#                 else:
#                     intensity_format = 'f'  # 32-bit float
                    
#                 # Check for compression
#                 is_mz_compressed = len(mz_array_binary.xpath('.//mzml:cvParam[@name="zlib compression"]', namespaces=ns)) > 0
#                 is_intensity_compressed = len(intensity_array_binary.xpath('.//mzml:cvParam[@name="zlib compression"]', namespaces=ns)) > 0

#                 # Find the binary elements
#                 mz_binary = mz_array_binary.find('.//mzml:binary', namespaces=ns)
#                 intensity_binary = intensity_array_binary.find('.//mzml:binary', namespaces=ns)

#                 if mz_binary is not None and intensity_binary is not None:
#                     # Decode the existing data
#                     existing_mz_data = base64.b64decode(mz_binary.text)
#                     existing_intensity_data = base64.b64decode(intensity_binary.text)
                    
#                     # Decompress if necessary
#                     if is_mz_compressed:
#                         existing_mz_data = zlib.decompress(existing_mz_data)
#                     if is_intensity_compressed:
#                         existing_intensity_data = zlib.decompress(existing_intensity_data)

#                     # Unpack the binary data
#                     if mz_format == 'f':
#                         existing_mz_values = struct.unpack('<' + 'f' * (len(existing_mz_data) // 4), existing_mz_data)
#                     else:
#                         existing_mz_values = struct.unpack('<' + 'd' * (len(existing_mz_data) // 8), existing_mz_data)

#                     if intensity_format == 'f':
#                         existing_intensity_values = struct.unpack('<' + 'f' * (len(existing_intensity_data) // 4), existing_intensity_data)
#                     else:
#                         existing_intensity_values = struct.unpack('<' + 'd' * (len(existing_intensity_data) // 8), existing_intensity_data)
                    
#                     # Combine the existing data with the new data
#                     combined_mz_values = list(existing_mz_values) + list(new_mz_values)
#                     combined_intensity_values = list(existing_intensity_values) + list(new_intensity_values)

#                     # Sort by m/z values
#                     combined_mz_values = np.array(combined_mz_values)
#                     combined_intensity_values = np.array(combined_intensity_values)
#                     idx = np.argsort(combined_mz_values)
#                     combined_mz_values = combined_mz_values[idx].tolist()
#                     combined_intensity_values = combined_intensity_values[idx].tolist()
                    
#                     # Update the defaultArrayLength attribute
#                     first_spectrum.set('defaultArrayLength', str(len(combined_mz_values)))

#                     # Encode the data with or without compression
#                     if is_mz_compressed:
#                         packed_mz_data = struct.pack('<' + mz_format * len(combined_mz_values), *combined_mz_values)
#                         compressed_mz_data = zlib.compress(packed_mz_data)
#                         combined_mz_encoded = base64.b64encode(compressed_mz_data).decode('ascii')
#                     else:
#                         packed_mz_data = struct.pack('<' + mz_format * len(combined_mz_values), *combined_mz_values)
#                         combined_mz_encoded = base64.b64encode(packed_mz_data).decode('ascii')
                    
#                     if is_intensity_compressed:
#                         packed_intensity_data = struct.pack('<' + intensity_format * len(combined_intensity_values), *combined_intensity_values)
#                         compressed_intensity_data = zlib.compress(packed_intensity_data)
#                         combined_intensity_encoded = base64.b64encode(compressed_intensity_data).decode('ascii')
#                     else:
#                         packed_intensity_data = struct.pack('<' + intensity_format * len(combined_intensity_values), *combined_intensity_values)
#                         combined_intensity_encoded = base64.b64encode(packed_intensity_data).decode('ascii')

#                     # Update the binary elements with the combined encoded data
#                     mz_binary.text = combined_mz_encoded
#                     intensity_binary.text = combined_intensity_encoded
                    
#                     # CRITICAL FIX: Update the encodedLength attributes
#                     if 'encodedLength' in mz_array_binary.attrib:
#                         mz_array_binary.set('encodedLength', str(len(combined_mz_encoded)))
#                     if 'encodedLength' in intensity_array_binary.attrib:
#                         intensity_array_binary.set('encodedLength', str(len(combined_intensity_encoded)))
                    
#                     # Also check for arrayLength attributes
#                     if 'arrayLength' in mz_array_binary.attrib:
#                         mz_array_binary.set('arrayLength', str(len(combined_mz_values)))
#                     if 'arrayLength' in intensity_array_binary.attrib:
#                         intensity_array_binary.set('arrayLength', str(len(combined_intensity_values)))
        
       
        
#         # Save the modified mzML file do not proceed until write is complete
        
#         # Save the modified mzML file
#         tree.write(output_mzML, pretty_print=True, xml_declaration=True, encoding="UTF-8")


#     def encode_message_to_file(self, message, input_mzML, output_mzML,method='hilbert'):
#         """
#         Encode a message to the first spectrum in an mzml file.

#         Parameters:
#         message (str): The message to encode.
#         input_mzML (str): The input mzML file path.
#         output_mzML (str): The output mzML file path.
#         method (str): The encoding method to use. Default is 'hilbert'. Options are 'hilbert' or 'sequential'.
        
#         """
#         if method not in ['hilbert', 'sequential']:
#             print("Error: Invalid encoding method. Use 'hilbert' or 'sequential'.")
#             return None
#         # Convert the message to a one-hot encoded vector
#         if method == 'hilbert':
#             # determine the order of the Hilbert curve based on the message length
#             self.determine_hilbert_curve_order(len(message))
#             # Convert the message to a one-hot encoded vector
#             self.convert_text_to_vector(message)
#             # Scale the coordinates to m/z and intensity values
#             coords = np.column_stack((self.x_indices, self.y_indices))
#             indices = np.argwhere(np.asarray(self.encoded_ascii_vector)!=0).flatten()
#             # print(coords[indices,:])
#             mz, intensity = self.scale_coords(coords, indices)
#         elif method == 'sequential':
#             # Sequential encoding
#             mz, intensity = self.sequential_encode_message(message)


#         # Add the message to the spectrum
#         self.add_values_to_first_spectrum(input_mzML, output_mzML, mz, intensity)
#         self.rebuild_mzmlfile_index(output_mzML)
#         self.update_checksum(output_mzML)
#         return None
    


