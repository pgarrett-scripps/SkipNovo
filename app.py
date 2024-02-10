import ms_deisotope

import streamlit as st
from peptacular.fragment import fragment
from peptacular.mass import calculate_mz
from peptacular.sequence import split_sequence, product_sequence, sort_sequence, calculate_sequence_length, \
    sequence_counter, permutate_sequence, combinate_with_replacement_sequence

import constants
from constants import AMINO_ACIDS
from intervaltree import Interval, IntervalTree
import pandas as pd
from msms_compression import SpectrumCompressorUrl

from query_params import parse_query_params, InvalidQueryParam, generate_app_url, QueryParams

st.set_page_config(page_title="SkipNovo", page_icon=":rock:", layout="wide")

try:
    qp = parse_query_params(st.query_params)
except InvalidQueryParam as e:
    st.error(str(e))
    st.stop()

with st.sidebar:
    st.header('SkipNovo')
    st.markdown('A tool for denovo peptide sequencing.')

    c1, c2, c3 = st.columns(3)

    seed_size = c1.number_input(
        label='N-mer Size',
        value=qp.seed_size
    )

    step_size = c2.number_input(
        label='Step Size',
        value=qp.step_size
    )

    max_steps = c3.number_input(
        label='Max Steps',
        value=qp.max_steps
    )

    keep_n_rows = st.number_input(
        label='Keep N Rows',
        value=qp.keep_n_rows
    )

    ppm = st.number_input(
        label='PPM',
        value=qp.ppm
    )

    c1, c2 = st.columns(2)
    min_intensity = c1.number_input(
        label='Minimum Intensity',
        value=qp.min_intensity
    )

    top_n = c2.number_input(
        label='Top N',
        value=qp.top_n
    )

    peak_picker = st.checkbox(
        label='Peak Picker',
        value=qp.peak_picker
    )

    peak_picker_min_intensity = qp.peak_picker_min_intensity
    peak_picker_mass_tolerance = qp.peak_picker_mass_tolerance
    if peak_picker:
        c1, c2 = st.columns(2)
        peak_picker_min_intensity = c1.number_input(
            label='Peak Picker Min Intensity',
            value=qp.peak_picker_min_intensity,
            min_value=0.0,
            max_value=1e9)
        peak_picker_mass_tolerance = c2.number_input(
            label='Peak Picker Mass Tolerance',
            value=qp.peak_picker_mass_tolerance,
            min_value=0.0,
            max_value=1e9, )

    c1, c2 = st.columns(2)

    min_charge = c1.number_input(
        label='Min Charge',
        value=qp.min_charge
    )

    max_charge = c2.number_input(
        label='Max Charge',
        value=qp.max_charge
    )

    c1, c2 = st.columns(2)
    precursor_mz = c1.number_input(
        label='Precursor MZ',
        value=qp.precursor_mz
    )

    precursor_charge = c2.number_input(
        label='Precursor Charge',
        value=qp.precursor_charge
    )

    spectra = st.text_area(
        label='Spectra',
        value='\n'.join(f'{s[0]} {s[1]}' for s in qp.spectra)
    )

    spectra = [tuple(map(float, s.split())) for s in spectra.split('\n')]

    with st.expander('Static Modifications'):

        c1, c2 = st.columns(2)
        n_term_static_mod = c1.number_input(label='N-term mod',
                                            value=qp.n_term_static_mod,
                                            help='Apply a static modification to the N-terminus')

        c_term_static_mod = c2.number_input(label='C-term mod',
                                            value=qp.c_term_static_mod,
                                            help='Apply a static modification to the C-terminus')

        num_static_mods = st.number_input(label='Number of unique static modifications',
                                          min_value=0,
                                          max_value=10,
                                          value=qp.num_static_mods,
                                          step=1,
                                          help='Add another modification row')

        # columns to lay out the inputs
        grid = st.columns([3, 2])


        def add_static_modification(r):

            aas = list(qp.static_mods[r][0]) if r < len(qp.static_mods) else []
            mod = qp.static_mods[r][1] if r < len(qp.static_mods) else 0.0

            with grid[0]:
                st.multiselect(label='Amino acids',
                               key=f'static_mod_residue{r}',
                               options=list(AMINO_ACIDS),
                               help='Select amino acids for which to apply the static modification',
                               default=aas)
            with grid[1]:
                st.number_input(label='Modification Mass (Da)',
                                step=0.00001,
                                key=f'static_mod_mass{r}',
                                help='The mass of the modification (in daltons)',
                                value=mod,
                                format='%.5f')


        # Loop to create rows of input widgets
        for r in range(num_static_mods):
            add_static_modification(r)

        static_mods = {}
        for r in range(num_static_mods):
            mod = "{:.5f}".format(st.session_state[f'static_mod_mass{r}'])
            for residue in st.session_state[f'static_mod_residue{r}']:
                static_mods[residue] = mod

    with st.expander('Variable Modifications'):

        c1, c2 = st.columns(2)
        n_term_var_mod = c1.number_input(label='N-term var mod',
                                         value=qp.n_term_var_mod,
                                         help='Apply a static modification to the N-terminus')

        c_term_var_mod = c2.number_input(label='C-term var mod',
                                         value=qp.c_term_var_mod,
                                         help='Apply a static modification to the C-terminus')

        num_var_mods = st.number_input(label='Number of unique var modifications',
                                       min_value=0,
                                       max_value=10,
                                       value=qp.num_var_mods,
                                       step=1,
                                       help='Add another modification row')

        # columns to lay out the inputs
        grid = st.columns([3, 2])


        def add_var_modification(r):

            aas = list(qp.var_mods[r][0]) if r < len(qp.var_mods) else []
            mod = qp.var_mods[r][1] if r < len(qp.var_mods) else 0.0

            with grid[0]:
                st.multiselect(label='Amino acids',
                               key=f'var_mod_residue{r}',
                               options=list(AMINO_ACIDS),
                               help='Select amino acids for which to apply the var modification',
                               default=aas)
            with grid[1]:
                st.number_input(label='Modification Mass (Da)',
                                step=0.00001,
                                key=f'var_mod_mass{r}',
                                help='The mass of the modification (in daltons)',
                                value=mod,
                                format='%.5f')


        # Loop to create rows of input widgets
        for r in range(num_var_mods):
            add_var_modification(r)

        var_mods = {}
        for r in range(num_var_mods):
            mod = "{:.5f}".format(st.session_state[f'var_mod_mass{r}'])
            for residue in st.session_state[f'var_mod_residue{r}']:
                var_mods[residue] = mod

qp_new = QueryParams(
    seed_size=seed_size,
    step_size=step_size,
    max_steps=max_steps,
    ppm=ppm,
    precursor_mz=precursor_mz,
    precursor_charge=precursor_charge,
    compression_algorithm='key',
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
    var_mods=var_mods
)

url = generate_app_url(qp_new)
url_chars = len(url)
st.write(f'##### [Analysis URL]({url}) (copy me and send to your friends!)')
st.write(f'Url Length: {url_chars} characters, {round(url_chars / 1024, 2)} KB')

if not st.button('Run'):
    st.stop()

precursor_mass = precursor_mz * precursor_charge - precursor_charge * 1.007276466812
precursor_mh = precursor_mass + 1.007276466812

st.write(f'Precursor MH: {precursor_mh}')

with st.spinner('filtering spectra...'):
    mzs = [s[0] for s in spectra]
    ints = [s[1] for s in spectra]

    if peak_picker:
        peaks = [(mz, intensity) for mz, intensity in zip(mzs, ints)]
        deconvoluted_peaks, _ = ms_deisotope.deconvolute_peaks(peaks, averagine=ms_deisotope.peptide,
                                                               scorer=ms_deisotope.MSDeconVFitter(
                                                                   peak_picker_min_intensity,
                                                                   peak_picker_mass_tolerance))
        mzs = [float(peak.mz) for peak in deconvoluted_peaks]
        ints = [float(peak.intensity) for peak in deconvoluted_peaks]

    spectra = [s for s in spectra if s[1] > min_intensity]

    # take top n
    spectra = sorted(spectra, key=lambda x: x[1], reverse=True)[:top_n]

amino_acids = list(AMINO_ACIDS)

# apply static mods
for residue, mod in static_mods.items():
    for i, aa in enumerate(amino_acids):
        if aa == residue:
            amino_acids[i] = f'{aa}({mod})'

# apply var mods
new_amino_acids = []
for residue, mod in var_mods.items():
    for i, aa in enumerate(amino_acids):
        if aa == residue:
            new_amino_acids.append(f'{aa}({mod})')

amino_acids += new_amino_acids

aa_sequence = ''.join(amino_acids)

st.write(f'### Amino Acid Sequence: {aa_sequence}')

with st.spinner('Generating combinations...'):
    sequences = set()
    for i in range(1, seed_size + 1):  # Start from 1 since peptides have at least one amino acid
        for combination in combinate_with_replacement_sequence(aa_sequence, i):
            sequences.add(''.join(combination))

    # get all modified/unmodified K and R residues
    tryptic_residues = {aa for aa in amino_acids if aa[0] in {'K', 'R'}}

    y_sequences = tryptic_residues.copy()

    for tryptic_residue in tryptic_residues:
        y_sequences.update({s + tryptic_residue for s in sequences})

    y_seqs = {}
    for sequence in y_sequences:
        charge_mzs = []
        for charge in range(min_charge, max_charge + 1):
            charge_mzs.append(calculate_mz(sequence=sequence, charge=charge, ion_type='y'))
        y_seqs[sequence] = charge_mzs

with st.spinner('Creating interval tree...'):
    spectra_tree = IntervalTree()
    for seq_mz, intensity in spectra:
        spectra_tree.add(
            Interval(seq_mz - seq_mz * ppm / 1_000_000, seq_mz + seq_mz * ppm / 1_000_000, (seq_mz, intensity)))

    wide_spectra_tree = IntervalTree()
    offset = precursor_mh * ppm / 1_000_000
    for seq_mz, intensity in spectra:
        wide_spectra_tree.add(
            Interval(seq_mz - offset, seq_mz + offset, (seq_mz, intensity)))

# add precursor_mh to spectra_tree
# spectra_tree.add(Interval(precursor_mh - precursor_mh * ppm / 1_000_000, precursor_mh + precursor_mh * ppm / 1_000_000, (precursor_mh, 1_000_000)))

with st.spinner('Finding matches...'):
    def find_matches(seqs):
        matches = []
        for seq in seqs:
            for charge, seq_mz in enumerate(seqs[seq], start=min_charge):
                for match in spectra_tree[seq_mz]:
                    matches.append((seq, match.data[0], match.data[1], seq_mz, charge))
        return matches


    y_matches = find_matches(y_seqs)


def make_comb_df(matches):
    df = pd.DataFrame(matches, columns=['sequence', 'mz', 'intensity', 'theoretical_mz', 'charge'])
    df['error'] = df['mz'] - df['theoretical_mz']
    df['ppm'] = df['error'] / df['mz'] * 1_000_000
    return df


y_df = make_comb_df(y_matches)


# Function to check if counter2 is a subsequence of counter1
def is_counter_subsequence(counter1, counter2):
    """Check if counter2 is a subsequence of counter1."""
    return all(counter1[elem] >= count for elem, count in counter2.items())


def get_subsequence_counts(df):
    sequences = df['sequence'].unique()

    # sort the sequences by length
    sequences = sorted(sequences, key=calculate_sequence_length, reverse=True)

    sequence_counters = {sequence: sequence_counter(sequence) for sequence in sequences}
    subsequence_counts = {sequence: 0 for sequence in sequences}

    for sequence, seq_counter in sequence_counters.items():
        for other_sequence, other_seq_counter in sequence_counters.items():
            if calculate_sequence_length(sequence) < calculate_sequence_length(other_sequence):
                break
            if sequence != other_sequence and is_counter_subsequence(seq_counter, other_seq_counter):
                subsequence_counts[sequence] += 1

    return subsequence_counts


with st.spinner('Finding subsequences...'):
    y_subsequence_counts = get_subsequence_counts(y_df)


def update_df(df, subsequence_counts):
    df['subset_count'] = df['sequence'].map(subsequence_counts)
    df = df[df['subset_count'] > 0]

    # keep only longest peptides use calculate_sequence_length to sort by length
    df['length'] = [calculate_sequence_length(seq) for seq in df['sequence']]
    max_len = df['length'].max()
    df = df[df['length'] == max_len]
    return df


y_df = update_df(y_df, y_subsequence_counts)

# keep top 100 rows based on subset_count
y_df = y_df.sort_values('subset_count', ascending=False).head(100)

with st.expander('combinations'):
    st.dataframe(y_df)

st.write(f'Found {len(y_df)} y sequences')


def generate_perms(sequences):
    permutations = set()
    for sequence in sequences:
        # Generate all permutations for each sequence and add them to the set
        for perm in permutate_sequence(sequence, calculate_sequence_length(sequence)):
            permutations.add(''.join(perm))
    return permutations


y_perms = generate_perms(y_df['sequence'].unique())

# ensure that y_perms ends in either k or R
y_perms = {perm for perm in y_perms if perm[-1] in {'K', 'R'}}

st.write(f'Found {len(y_perms)} y permutations')


def fragment_perms(sequences, ion_type):
    fragment_ions = {}
    for perm in sequences:
        fragment_ions[perm] = fragment(sequence=perm, ion_types=(ion_type,),
                                       charges=[c for c in range(min_charge, max_charge + 1)])
    return fragment_ions


y_fragment_ions = fragment_perms(y_perms, 'y')


def get_matches(fragment_ions):
    scores = {}
    int_scores = {}
    for perm, ions in fragment_ions.items():
        matches = 0
        intensity = 0
        for ion in ions:
            for match in spectra_tree[ion]:
                matches += 1
                intensity += match.data[1]
        scores[perm] = matches
        int_scores[perm] = intensity

    return scores, int_scores


def get_complement_matches(fragment_ions):
    scores = {}
    int_scores = {}
    for perm, ions in fragment_ions.items():
        matches = 0
        intensity = 0
        for ion in ions:
            ion = precursor_mh - ion + 1.007276
            for match in wide_spectra_tree[ion]:
                matches += 1
                intensity += match.data[1]
        scores[perm] = matches
        int_scores[perm] = intensity

    return scores, int_scores


y_scores, y_int_scores = get_matches(y_fragment_ions)

y_complete_scores, y_complete_int_scores = get_complement_matches(y_fragment_ions)


def make_perm_df(scores, int_scores):
    df = pd.DataFrame(scores.items(), columns=['sequence', 'matches'])
    df['intensity'] = df['sequence'].map(int_scores)
    # remove rows which have a subset count less than 1
    df = df[df['matches'] > 0]

    return df


def update_complete_df(df, complete_scores, complete_int_scores):
    df['complement_matches'] = df['sequence'].map(complete_scores)
    df['complement_intensity'] = df['sequence'].map(complete_int_scores)
    return df


y_df = make_perm_df(y_scores, y_int_scores)

y_df = update_complete_df(y_df, y_complete_scores, y_complete_int_scores)

y_df['comb_seq'] = [sort_sequence(seq) for seq in y_df['sequence']]

# sort by matches and intensity
y_df = y_df.sort_values(by=['matches', 'intensity'], ascending=False)

# filter by top 100
y_df = y_df.head(100)

mzs = [s[0] for s in spectra]
ints = [s[1] for s in spectra]
compresseed_spectra = SpectrumCompressorUrl.compress(mzs, ints)


def make_clickable(sequence, ion, final=False):
    sequence_mz = calculate_mz(sequence=sequence, ion_type='y', charge=1)
    missing_mz = precursor_mh - sequence_mz

    fragment_types = [f'{c}{i}' for c in range(min_charge, max_charge + 1) for i in ['b', 'y']]
    fragment_types = ';'.join(fragment_types)

    missing_mz_str = ''
    if abs(missing_mz) > abs(precursor_mh - sequence_mz) * ppm / 1e6:
        missing_mz_str = f'X({missing_mz})'

    if final is True:
        missing_mz_str = ''

    if ion == 'b':
        sequence = sequence + missing_mz_str
    else:
        sequence = missing_mz_str + sequence

    return constants.SPECTRA_VIEWER_URL + f'?sequence={sequence}&fragment_types={fragment_types}&spectra={compresseed_spectra}&compression_algorithm=brotli&mass_tolerance_type=ppm&mass_tolerance={ppm}&aa_masses=X0.0'


y_df['link'] = [make_clickable(seq, 'y') for seq in y_df['sequence']]

st.write('SEEDS')
st.dataframe(
    y_df,
    column_config={
        "link": st.column_config.LinkColumn(
            display_text="ions"),
    },
    hide_index=True,
    use_container_width=True
)

valid_rows = []

step_sequences = set()
for i in range(1, step_size + 1):
    step_sequences.update(set(product_sequence(aa_sequence, i)))


def fast_score(old_sequence, new_step):
    y_scores, y_int_scores = 0, 0
    y_complete_scores, y_complete_int_scores = 0, 0

    new_ions = []
    new_sequence = old_sequence
    for a in split_sequence(new_step):
        new_sequence = a + new_sequence
        for c in range(min_charge, max_charge + 1):
            y_ion = calculate_mz(sequence=new_sequence, ion_type='y', charge=c)

            for match in spectra_tree[y_ion]:
                y_scores += 1
                y_int_scores += match.data[1]

            # TODO: Fix this issue... no worky im dumb im dumb
            inverse_ion = (precursor_mass + 1.007276 * c) / c - y_ion + 1.007276

            for match in wide_spectra_tree[inverse_ion]:
                y_complete_scores += 1
                y_complete_int_scores += match.data[1]

    return y_scores, y_int_scores, y_complete_scores, y_complete_int_scores, new_sequence


with st.spinner('Calculating valid sequences...'):
    with st.expander('Steps'):
        for s in range(max_steps):

            data = []
            for seq, m, i, cm, ci in y_df[
                ['sequence', 'matches', 'intensity', 'complement_matches', 'complement_intensity']].values:
                for new_step in step_sequences:
                    y_scores, y_int_scores, y_complete_scores, y_complete_int_scores, new_sequence = fast_score(seq,
                                                                                                                new_step)
                    data.append([new_sequence, m + y_scores, i + y_int_scores, cm + y_complete_scores,
                                 ci + y_complete_int_scores])

            y_df = pd.DataFrame(data, columns=['sequence', 'matches', 'intensity', 'complement_matches',
                                               'complement_intensity'])

            # new_sequences = [step_sequence + seq for seq in y_df['sequence'] for step_sequence in step_sequences]
            # y_fragment_ions = fragment_perms(new_sequences, 'y')
            # y_scores, y_int_scores = get_matches(y_fragment_ions)
            # y_complete_scores, y_complete_int_scores = get_complement_matches(y_fragment_ions)
            # y_df = make_perm_df(y_scores, y_int_scores)
            # y_df = update_complete_df(y_df, y_complete_scores, y_complete_int_scores)

            y_df['mz'] = [calculate_mz(sequence=seq, ion_type='y', charge=1) for seq in y_df['sequence']]
            y_df['comb_seq'] = [sort_sequence(seq) for seq in y_df['sequence']]
            y_df['len'] = [calculate_sequence_length(seq) for seq in y_df['sequence']]

            # score should be matches * intensity + (complement_matches * complement_intensity)*0.5
            y_df['score'] = (y_df['matches'] * y_df['intensity'] + (
                    y_df['complement_matches'] * y_df['complement_intensity'])) / y_df['len'] * 0.5

            valid_row_flags = (y_df['mz'] <= precursor_mh + precursor_mh * ppm / 1_000_000) & (
                    y_df['mz'] >= precursor_mh - precursor_mh * ppm / 1_000_000)

            valid_df = y_df[valid_row_flags]
            valid_rows.append(valid_df)
            y_df = y_df[~valid_row_flags]

            y_df = y_df[y_df['mz'] <= precursor_mh + precursor_mh * ppm / 1_000_000]
            y_df = y_df.sort_values(by=['score'], ascending=False)
            y_df = y_df.groupby('mz').first().reset_index()
            y_df = y_df.sort_values(by=['score'], ascending=False)

            y_df = y_df.head(100)

            st.write(f'Step {s}')

            y_df['link'] = [make_clickable(seq, 'y') for seq in y_df['sequence']]

            st.dataframe(
                y_df,
                column_config={
                    "link": st.column_config.LinkColumn(
                        display_text="ions"),
                },
                hide_index=True,
                use_container_width=True
            )

            if len(y_df) == 0:
                break

st.write('done')

st.write('valid sequences')

valid_df = pd.concat(valid_rows)
valid_df = valid_df.drop_duplicates(subset=['sequence'])
valid_df = valid_df.sort_values(by=['matches', 'intensity'], ascending=False)
valid_df['link'] = [make_clickable(seq, 'y') for seq in valid_df['sequence']]

# refragment and rescore

final_data = []
for sequence, mz in valid_df[['sequence', 'mz']].values:
    frags = fragment(sequence=sequence, ion_types=['y', 'b'], charges=[c for c in range(min_charge, max_charge + 1)])

    matches = 0
    intensity = 0
    for ion in frags:
        for match in spectra_tree[ion]:
            matches += 1
            intensity += match.data[1]

    final_data.append({
        'sequence': sequence,
        'matches': matches,
        'intensity': intensity,
        'link': make_clickable(sequence, 'y', True)
    })

final_df = pd.DataFrame(final_data)

if len(final_df) == 0:
    st.write('No valid sequences')
    st.stop()

final_df['score'] = final_df['matches'] * final_df['intensity']
final_df = final_df.sort_values(by=['score'], ascending=False)

st.dataframe(
    final_df,
    column_config={
        "link": st.column_config.LinkColumn(
            display_text="ions"),
    },
    hide_index=True,
    use_container_width=True
)
