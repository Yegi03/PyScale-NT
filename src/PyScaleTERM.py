#!/usr/bin/env python3
"""
PyScale - Terminal Version
Nucleotide propensity-scale profiling with Lyapunov exponent estimation.

Authors: Joshua Davies, Yeganeh Abdollahinejad, Darren R. Flower, Amit K. Chattopadhyay
License: MIT
"""

import numpy as np
import os
import sys

# Configuration
arguments = sys.argv[1:]
base_directory = os.path.join(os.path.dirname(__file__), '..')

if '-h' in arguments:
    with open(os.path.join(base_directory, 'docs', 'help.txt'), 'r') as f:
        print(f.read())
    sys.exit()

# Parse command-line options
all_output_override = '-ao' in arguments
complement_strand = '-c' in arguments
flip_read = '-r' in arguments
periodicity_assumption = '-p' in arguments
window_averaging = '-wa' in arguments

try:
    input_location = arguments[arguments.index('-f')+1]
except ValueError:
    input_location = os.path.join(base_directory, 'data', 'examples')

try:
    output_location = arguments[arguments.index('-o')+1]
except ValueError:
    output_location = os.path.join(base_directory, 'data', 'output')

try:
    output_type = tuple(arguments[arguments.index('-t')+1])
except ValueError:
    output_type = 'a', 'i'

try:
    scale_location = arguments[arguments.index('-s')+1]
except:
    scale_location = os.path.join(base_directory, 'data', 'scales')

try:
    step_size = int(arguments[arguments.index('-ss')+1])
except:
    step_size = 2

try:
    window_size = int(arguments[arguments.index('-ws')+1])
except:
    window_size = 7

valid_extensions = '.txt', '.dat'

# Lyapunov parameters
try:
    lyap_m = int(arguments[arguments.index('--lyap_m')+1])
except:
    lyap_m = 2

try:
    lyap_tau = int(arguments[arguments.index('--lyap_tau')+1])
except:
    lyap_tau = 1

try:
    lyap_kmax = int(arguments[arguments.index('--lyap_kmax')+1])
except:
    lyap_kmax = 30

try:
    lyap_theiler = int(arguments[arguments.index('--lyap_theiler')+1])
except:
    lyap_theiler = 0

try:
    lyap_eps = float(arguments[arguments.index('--lyap_eps')+1])
except:
    lyap_eps = 0.0

try:
    lyap_fit_start = int(arguments[arguments.index('--lyap_fit_start')+1])
except:
    lyap_fit_start = 1

try:
    lyap_fit_end = int(arguments[arguments.index('--lyap_fit_end')+1])
except:
    lyap_fit_end = 10

try:
    lyap_knn = int(arguments[arguments.index('--lyap_knn')+1])
except:
    lyap_knn = 15

lyap_auto_fit = '--lyap_autofit' in arguments

try:
    lyap_r2min = float(arguments[arguments.index('--lyap_r2min')+1])
except:
    lyap_r2min = 0.98

try:
    lyap_minlen = int(arguments[arguments.index('--lyap_minlen')+1])
except:
    lyap_minlen = 6

# Utility functions
clear = lambda: os.system('clear')
divider = lambda length: '\n' + '─' * length + '\n'
files_per = lambda scales: {True: len(range(window_size, 1, -(step_size))) + 1, False: 1}[window_averaging] * len(scales)
setting_display = lambda: '\n'.join([var + ': ' + str(state) for var, state in (('Periodicity', periodicity_assumption), ('Flip', flip_read), ('Complement', complement_strand))])


class DAT_PIECE():

    def __init__(self, title, raw_array):
        self.title = title
        self.original_raw = np.array(list(''.join(raw_array)))

    def create_complement(self):
        complementary_bases = {'A': 'T', 'T': 'A', 'U': 'A', 'G': 'C', 'C': 'G'}
        self.complement_raw = np.array([complementary_bases[original_base] for original_base in self.original_raw])

    def convert(self, scale):
        self.original_convert = np.array([scale[base] for base in self.original_raw])
        try:
            self.complement_convert = np.array([scale[base] for base in self.complement_raw])
        except AttributeError:
            pass

    def window_average(self, width, periodicity, reverse, complement):
        N = len(self.original_raw)
        half_width = width // 2

        if complement:
            working_array = self.complement_convert
        else:
            working_array = self.original_convert
        if reverse:
            working_array = working_array[::-1]
        if periodicity:
            sample_range = range(N)
        else:
            sample_range = range(half_width, N - half_width)

        output_array = np.array([np.mean(np.take(working_array, range(i - half_width, i + 1 + half_width), mode='wrap')) for i in sample_range])

        if not periodicity:
            return np.concatenate((np.zeros(half_width), output_array, np.zeros(half_width)))
        return output_array


class File():

    def __init__(self, path):
        self.file = path
        self.file_name = os.path.basename(path)
        with open(self.file, 'r') as f:
            self.raw_lines = f.read().splitlines()

    def parse_lines(self, current_contents):
        for i, line in enumerate(current_contents[1:]):
            i += 1
            if line.startswith('>'):
                return [DAT_PIECE(current_contents[0][1:], current_contents[1:i])] + self.parse_lines(current_contents[i:])
        return [DAT_PIECE(current_contents[0][1:], current_contents[1:])]

    def parse_scales(self, current_contents):
        for i, line in enumerate(current_contents[1:]):
            i += 1
            if line.startswith('>'):
                new_scale = {current_contents[0][1:]: dict(zip([i.split()[0] for i in current_contents[1:i]], [float(i.split()[1]) for i in current_contents[1:i]]))}
                new_scale.update(self.parse_scales(current_contents[i:]))
                return new_scale
        new_scale = {current_contents[0][1:]: dict(zip([i.split()[0] for i in current_contents[1:]], [float(i.split()[1]) for i in current_contents[1:]]))}
        return new_scale

    def output(self, data_array):
        self.progress_bar = DualProgressBar(len(data_array), files_per(scale_data), 55, ['File ' + self.file_name + ' & data piece ' + i.title for i in data_array])
        for data_piece in data_array:
            if complement_strand:
                data_piece.create_complement()
            file_name = ''.join([character for character, state in (('P', periodicity_assumption), ('F', flip_read), ('C', complement_strand)) if state == True])
            if file_name == '':
                file_name = 'vanilla'
            for scale_title in scale_data:
                scale = scale_data[scale_title]
                data_piece.convert(scale)
                output_arrays = [data_piece.window_average(i, periodicity_assumption, flip_read, complement_strand) for i in range(window_size, 2, -(step_size))]
                file_path = os.path.join(output_location, data_piece.title, scale_title)
                self.write_file(output_arrays[0], file_name, file_path, window_size, periodicity_assumption)
                if window_averaging:
                    for i, array in enumerate(output_arrays[1:]):
                        width = window_size - ((i + 1) * step_size)
                        self.write_file(array, file_name, file_path, width, periodicity_assumption)
                    self.write_file(average_arrays(output_arrays), file_name, file_path, 'avg', periodicity_assumption)

    def write_file(self, output_array, file_name, file_path, width, periodicity):
        self.progress_bar.step()
        file_name = file_name + 'W' + str(width) + '.txt'
        file_path = os.path.join(file_path, file_name)
        if width != 'avg':
            lyap_val, lyap_meta = return_lyapunov(
                output_array, periodicity, width,
                m=lyap_m, tau=lyap_tau, k_max=lyap_kmax,
                theiler=lyap_theiler, eps=lyap_eps,
                fit_start=lyap_fit_start, fit_end=lyap_fit_end,
                k_neighbors=lyap_knn, auto_fit=lyap_auto_fit,
                r2_min=lyap_r2min, min_fit_len=lyap_minlen
            )
            header_lines = ['# lyapunov = ' + str(lyap_val) + '\n',
                            '# ' + str(lyap_meta) + '\n']
        else:
            header_lines = ['# lyapunov = unavailable for averaged data sets\n']
        if not os.path.exists(os.path.dirname(file_path)):
            os.makedirs(os.path.dirname(file_path))
        with open(file_path, 'w') as f:
            output_lines = header_lines + [str(x) + ' ' + str(y) + '\n' for x, y in enumerate(output_array)]
            f.writelines(output_lines)


class DualProgressBar():

    def __init__(self, master_max, sub_max, length, headers):
        self.__dict__.update(locals())
        self.sub_progress = 0
        self.master_progress = 0

    def step(self):
        self.sub_progress += 1
        self.master_progress += self.sub_progress // self.sub_max
        if self.master_progress != self.master_max:
            self.sub_progress %= self.sub_max
        self.update()

    def update(self):
        clear()
        try:
            print('Processing ' + self.headers[self.master_progress])
        except IndexError:
            print('Done!')
        print(divider(55) + 'Current settings: ' + divider(55) + setting_display() + divider(55))
        for progress, maximum in (self.sub_progress, self.sub_max), (self.master_progress, self.master_max):
            completed = int((progress / maximum) * self.length)
            remaining = self.length - completed
            print('\n' + str(progress) + '/' + str(maximum) + divider(55) + '▮' * completed + '-' * remaining + divider(55))


def average_arrays(array_tuple):
    zipped = list(zip(*array_tuple))
    output = [np.mean(i) for i in zipped]
    return output


def average_dictionaries(dictionary_tuple):
    output = {key: np.mean([dictionary[key] for dictionary in dictionary_tuple]) for key in dictionary_tuple[0]}
    return output


def merge_dictionaries(dictionary_tuple):
    [dictionary_tuple[0].update(i) for i in dictionary_tuple[1:]]
    return dictionary_tuple[0]


def return_lyapunov(array, periodicity, width,
                    m=2, tau=1, k_max=30, theiler=0, eps=0.0,
                    fit_start=1, fit_end=10,
                    k_neighbors=15, auto_fit=True, r2_min=0.98, min_fit_len=6):
    """
    Estimate the largest Lyapunov exponent using the Rosenstein algorithm
    with K-nearest neighbor averaging and optional automatic fit-region selection.
    """
    x = np.asarray(array, dtype=float)

    if (not periodicity) and (width != 'avg'):
        half_width = int(width) // 2
        if len(x) > 2 * half_width:
            x = x[half_width:-half_width]

    N = len(x)
    if N < 20:
        return (float('nan'), "unavailable: series too short")

    if theiler is None or theiler == 0:
        theiler = m * tau

    eps_cutoff = None if (eps is None or eps <= 0.0) else float(eps)
    k_neighbors = max(1, min(k_neighbors, 30))

    if m < 1 or tau < 1:
        return (float('nan'), "unavailable: m and tau must be >=1")

    M = N - (m - 1) * tau
    if M <= 2:
        return (float('nan'), "unavailable: series too short for embedding")

    if M <= k_max + 1:
        k_max = M - 2
    if k_max < 2:
        return (float('nan'), "unavailable: k_max too large for series length")

    # Delay embedding
    idx = np.arange(M)[:, None] + np.arange(m)[None, :] * tau
    X = x[idx]

    # K-nearest neighbor search with Theiler exclusion
    nn = np.full((M, k_neighbors), -1, dtype=int)

    for i in range(M):
        d = np.linalg.norm(X - X[i], axis=1)
        lo = max(0, i - theiler)
        hi = min(M, i + theiler + 1)
        d[lo:hi] = np.inf

        finite_mask = np.isfinite(d)
        if not finite_mask.any():
            continue

        sorted_idxs = np.argsort(d)

        count = 0
        for j in sorted_idxs:
            if count >= k_neighbors:
                break
            if not np.isfinite(d[j]):
                break
            if eps_cutoff is not None and d[j] > eps_cutoff:
                continue
            nn[i, count] = j
            count += 1

    valid = np.where(nn[:, 0] >= 0)[0]
    if len(valid) < max(5, M // 20):
        return (float('nan'), "unavailable: too few neighbor pairs")

    # Mean log divergence curve
    small = 1e-12
    y = np.full(k_max + 1, np.nan)

    for k in range(k_max + 1):
        log_divs = []
        for i in valid:
            neighbors = nn[i, :]
            neighbors = neighbors[neighbors >= 0]

            if len(neighbors) == 0:
                continue

            dists = []
            for j in neighbors:
                if i + k < M and j + k < M:
                    dist = np.linalg.norm(X[i + k] - X[j + k])
                    dists.append(max(dist, small))

            if len(dists) > 0:
                mean_dist = np.mean(dists)
                log_divs.append(np.log(mean_dist))

        if len(log_divs) < 5:
            break
        y[k] = float(np.mean(log_divs))

    # Fit slope
    if auto_fit:
        best_fit = None
        best_score = -1

        for a in range(k_max + 1):
            for b in range(a + min_fit_len, k_max + 1):
                segment = y[a:b + 1]
                if not np.isfinite(segment).all():
                    continue
                if len(segment) < min_fit_len:
                    continue

                t = np.arange(a, b + 1)
                slope, intercept = np.polyfit(t, segment, 1)

                if slope <= 0:
                    continue

                y_pred = slope * t + intercept
                ss_res = np.sum((segment - y_pred) ** 2)
                ss_tot = np.sum((segment - np.mean(segment)) ** 2)
                if ss_tot < 1e-12:
                    continue
                r2 = 1 - ss_res / ss_tot

                if r2 < r2_min:
                    continue

                seg_len = b - a + 1
                if seg_len > best_score or (seg_len == best_score and best_fit is not None and r2 > best_fit[2]):
                    best_fit = (a, b, r2, slope)
                    best_score = seg_len

        if best_fit is None:
            return (float('nan'), f"unavailable: no linear region (r2_min={r2_min}, min_len={min_fit_len})")

        fit_start, fit_end, r2, slope = best_fit
        meta = f"method=rosenstein_knn m={m} tau={tau} theiler={theiler} k={k_neighbors} k_max={k_max} eps={eps_cutoff} fit=[{fit_start},{fit_end}] r2={r2:.3f}"
        return (float(slope), meta)

    else:
        if fit_start < 0 or fit_end > k_max or fit_start >= fit_end:
            return (float('nan'), "unavailable: invalid fit range")

        segment = y[fit_start:fit_end + 1]
        if not np.isfinite(segment).all():
            return (float('nan'), "unavailable: insufficient divergence data in fit range")

        t = np.arange(fit_start, fit_end + 1)
        slope, intercept = np.polyfit(t, segment, 1)

        y_pred = slope * t + intercept
        ss_res = np.sum((segment - y_pred) ** 2)
        ss_tot = np.sum((segment - np.mean(segment)) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0

        meta = f"method=rosenstein_knn m={m} tau={tau} theiler={theiler} k={k_neighbors} k_max={k_max} eps={eps_cutoff} fit=[{fit_start},{fit_end}] r2={r2:.3f}"
        return (float(slope), meta)


def option_menu(titles, container, prompt):
    print('\na Use all files')
    print(''.join(['\n' + str(i) + ' ' + title + '\n' for i, title in enumerate(titles)]))
    while True:
        try:
            user_selection = input(prompt).split()
            if user_selection == []:
                raise ValueError
            if 'a' in user_selection:
                return container
            if all([int(i) in range(len(container)) for i in user_selection]):
                return np.take(container, user_selection)
        except ValueError:
            print('\nInvalid entry! ensure your entries are present in the menu')


# Main execution
scale_files = [File(os.path.join(scale_location, scale_file)) for scale_file in os.listdir(scale_location)]
scale_data = merge_dictionaries([scale_file.parse_scales(scale_file.raw_lines) for scale_file in scale_files])
if 'a' in output_type and not ('i' in output_type):
    scale_data = {'Averaged': average_dictionaries(tuple(scale_data.values()))}
elif 'a' in output_type:
    scale_data['Averaged'] = average_dictionaries(tuple(scale_data.values()))

data_files = [File(os.path.join(input_location, data_file)) for data_file in os.listdir(input_location)]
data_files = option_menu([data_file.file_name for data_file in data_files], data_files, 'Multiple data files detected, enter one or more options separated by a space: ')
if not all_output_override:
    [data_file.output(data_file.parse_lines(data_file.raw_lines)) for data_file in data_files]
else:
    window_averaging = True
    for periodicity_assumption in (True, False):
        for flip_read in (True, False):
            for complement_strand in (True, False):
                [data_file.output(data_file.parse_lines(data_file.raw_lines)) for data_file in data_files]
