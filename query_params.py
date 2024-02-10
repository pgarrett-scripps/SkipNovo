import dataclasses
import urllib
from typing import List, Tuple, Dict

import requests

import constants

from msms_compression import SpectrumCompressorF32LzstringUri, SpectrumCompressorUrl

from msms_compression import BaseCompressor, url_encoder, SpectrumCompressorI32, brotli_compressor

lossy_compressor = BaseCompressor(SpectrumCompressorI32(2, 1), brotli_compressor, url_encoder)


class InvalidQueryParam(Exception):
    pass


@dataclasses.dataclass
class QueryParams:
    seed_size: int
    step_size: int
    max_steps: int
    ppm: float
    precursor_mz: float
    precursor_charge: int
    compression_algorithm: str
    spectra: List[Tuple[float, float]]
    min_intensity: float
    top_n: int
    peak_picker: bool
    peak_picker_min_intensity: float
    peak_picker_mass_tolerance: float
    min_charge: int
    max_charge: int
    keep_n_rows: int

    n_term_static_mod: float
    c_term_static_mod: float
    num_static_mods: int
    static_mods: Dict[str, float]

    n_term_var_mod: float
    c_term_var_mod: float
    num_var_mods: int
    var_mods: Dict[str, float]


def serialize_spectra(spectra: List[Tuple[float, float]], compression_algorithm: str) -> str:
    if not spectra:
        mzs = []
        ints = []
    else:
        mzs, ints = zip(*spectra)

    if compression_algorithm == 'lzstring':
        return SpectrumCompressorF32LzstringUri.compress(mzs, ints)
    elif compression_algorithm == 'brotli':
        return SpectrumCompressorUrl.compress(mzs, ints)
    elif compression_algorithm == 'lossy':
        return lossy_compressor.compress(mzs, ints)
    elif compression_algorithm == 'key':
        # Assuming 'spectra' is a list with one element that is the key
        data = {"mzs": list(mzs), "intensities": list(ints)}

        response = requests.post(f'{constants.COMP_API}/store', json=data)
        if response.status_code == 200:
            response_data = response.json()
            return response_data.get("key")
        else:
            raise Exception("Failed to store spectra on the server")

    else:
        raise ValueError(f'Invalid compression algorithm: {compression_algorithm}')


def deserialize_spectra(s: str, compression_algorithm: str) -> List[Tuple[float, float]]:
    if compression_algorithm == 'lzstring':
        mzs, ints = SpectrumCompressorF32LzstringUri.decompress(s)
        return list(zip(mzs, ints))
    elif compression_algorithm == 'brotli':
        mzs, ints = SpectrumCompressorUrl.decompress(s)
        return list(zip(mzs, ints))
    elif compression_algorithm == 'lossy':
        mzs, ints = lossy_compressor.decompress(s)
        return list(zip(mzs, ints))
    elif compression_algorithm == 'key':
        # Assuming 'spectra' is a list with one element that is the key
        response = requests.get(f'{constants.COMP_API}/retrieve/{s}')
        if response.status_code == 200:
            response_data = response.json()
            return list(zip(response_data.get("mzs"), response_data.get("intensities")))
        else:
            raise Exception("Failed to retrieve spectra from the server")
    else:
        raise ValueError(f'Invalid compression algorithm: {compression_algorithm}')


def validate_float(value):
    try:
        return float(value)
    except ValueError:
        raise ValueError(f"Must be a float")


def validate_int(value):
    try:
        return int(value)
    except ValueError:
        raise ValueError(f"Must be an int")


def validate_bool(value):
    if isinstance(value, bool):
        return value
    if value in ['True', 'true', 1]:
        return True
    if value in ['False', 'false', 0]:
        return False
    raise ValueError(f"Must be a boolean")


def parse_query_params(query_params: dict) -> QueryParams:
    seed_size = validate_int(query_params.get('seed_size', 5))
    step_size = validate_int(query_params.get('step_size', 2))
    max_steps = validate_int(query_params.get('max_steps', 10))
    ppm = validate_float(query_params.get('ppm', 20.0))
    precursor_mz = validate_float(query_params.get('precursor_mz', 690.3821))
    precursor_charge = validate_int(query_params.get('precursor_charge', 2))

    compression_algorithm = query_params.get('compression_algorithm', 'key')
    spectra = query_params.get('spectra', serialize_spectra(constants.DEFAULT_SPECTRA, compression_algorithm))
    spectra = deserialize_spectra(spectra, compression_algorithm)

    for s in spectra:
        if s[0] < 0 or s[1] < 0:
            InvalidQueryParam("Spectra must be non-negative numbers.")

    min_intensity = validate_float(query_params.get('min_intensity', 0.0))
    top_n = validate_int(query_params.get('top_n', 100))
    peak_picker = validate_bool(query_params.get('peak_picker', True))
    peak_picker_min_intensity = validate_float(query_params.get('peak_picker_min_intensity', 1.0))
    peak_picker_mass_tolerance = validate_float(query_params.get('peak_picker_mass_tolerance', 0.002))
    min_charge = validate_int(query_params.get('min_charge', 1))
    max_charge = validate_int(query_params.get('max_charge', 1))

    keep_n_rows = validate_int(query_params.get('keep_n_rows', 100))

    n_term_static_mod = validate_float(query_params.get('n_term_static_mod', 0.0))
    c_term_static_mod = validate_float(query_params.get('c_term_static_mod', 0.0))
    num_static_mods = validate_int(query_params.get('num_static_mods', 1))
    query_static_mods_str = query_params.get('static_mods', 'C:57.02146')
    static_mods = [(s.split(':')[0], float(s.split(':')[1])) for s in query_static_mods_str.split(';') if s]

    n_term_var_mod = validate_float(query_params.get('n_term_var_mod', 0.0))
    c_term_var_mod = validate_float(query_params.get('c_term_var_mod', 0.0))
    num_var_mods = validate_int(query_params.get('num_var_mods', 1))
    query_var_mods_str = query_params.get('var_mods', 'M:15.99491')
    var_mods = [(s.split(':')[0], float(s.split(':')[1])) for s in query_var_mods_str.split(';') if s]

    return QueryParams(
        seed_size=seed_size,
        step_size=step_size,
        max_steps=max_steps,
        ppm=ppm,
        precursor_mz=precursor_mz,
        precursor_charge=precursor_charge,
        compression_algorithm=compression_algorithm,
        spectra=spectra,
        min_intensity=min_intensity,
        top_n=top_n,
        peak_picker=peak_picker,
        peak_picker_min_intensity=peak_picker_min_intensity,
        peak_picker_mass_tolerance=peak_picker_mass_tolerance,
        min_charge=min_charge,
        max_charge=max_charge,
        keep_n_rows=keep_n_rows,
        n_term_static_mod=n_term_static_mod,
        c_term_static_mod=c_term_static_mod,
        num_static_mods=num_static_mods,
        static_mods=static_mods,
        n_term_var_mod=n_term_var_mod,
        c_term_var_mod=c_term_var_mod,
        num_var_mods=num_var_mods,
        var_mods=var_mods,
    )


def generate_app_url(qp: QueryParams) -> str:

    # flip the dictionary
    static_mods_rev = {}
    for aa, mod in qp.static_mods.items():
        static_mods_rev.setdefault(mod, []).append(aa)
    static_mod_str = ';'.join([f'{"".join(aa)}:{mod}' for mod, aa in static_mods_rev.items()])

    var_mods_rev = {}
    for aa, mod in qp.var_mods.items():
        var_mods_rev.setdefault(mod, []).append(aa)
    var_mod_str = ';'.join([f'{"".join(aa)}:{mod}' for mod, aa in var_mods_rev.items()])


    params = {
        'seed_size': qp.seed_size,
        'step_size': qp.step_size,
        'max_steps': qp.max_steps,
        'ppm': qp.ppm,
        'precursor_mz': qp.precursor_mz,
        'precursor_charge': qp.precursor_charge,
        'spectra': serialize_spectra(qp.spectra, qp.compression_algorithm),
        'min_intensity': qp.min_intensity,
        'top_n': qp.top_n,
        'peak_picker': qp.peak_picker,
        'peak_picker_min_intensity': qp.peak_picker_min_intensity,
        'peak_picker_mass_tolerance': qp.peak_picker_mass_tolerance,
        'min_charge': qp.min_charge,
        'max_charge': qp.max_charge,
        'keep_n_rows': qp.keep_n_rows,
        'n_term_static_mod': qp.n_term_static_mod,
        'c_term_static_mod': qp.c_term_static_mod,
        'num_static_mods': qp.num_static_mods,
        'static_mods': static_mod_str,
        'n_term_var_mod': qp.n_term_var_mod,
        'c_term_var_mod': qp.c_term_var_mod,
        'num_var_mods': qp.num_var_mods,
        'var_mods': var_mod_str,

    }

    # Serialize parameters to query string, ensuring values are properly encoded
    query_string = '&'.join(f'{key}={urllib.parse.quote(str(value))}' for key, value in params.items())
    return f'{constants.BASE_URL}?{query_string}'