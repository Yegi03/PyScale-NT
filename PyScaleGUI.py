#!/usr/bin/env python3
"""
PyScale - GUI Version
Nucleotide propensity-scale profiling with Lyapunov exponent estimation.

Authors: Joshua Davies, Yeganeh Abdollahinejad, Darren R. Flower, Amit K. Chattopadhyay
License: MIT
"""

import numpy as np
import os
import tkinter.filedialog as filedialog
import tkinter.messagebox as messagebox
from tkinter import *
from tkinter import ttk
from shutil import rmtree

root = Tk()
root.resizable(width=False, height=False)

plot_calls = 0

# GUI Layout
option_toolbar = Frame(root, bd=1, relief='solid')
option_toolbar.grid(column=0, columnspan=2, sticky=(N, E, W), padx=5, pady=5)

file_browser = Frame(root, bd=1, relief='solid')
file_browser.grid(column=0, row=1, sticky=(N, S, W), padx=5, pady=5)

file_tree = ttk.Treeview(file_browser, height=15)
file_tree.grid(column=0, row=4, padx=5, pady=5)

progress_updater = Frame(root, bd=1, relief='solid')
progress_updater.grid(column=1, row=1, sticky=(N, S, W), padx=5, pady=5)
progress_output = Text(progress_updater, height=30, width=224, state=DISABLED)
progress_output.grid(column=0, row=0)

# Default option values
all_output_override = False
complement_strand = False
flip_read = False
input_location = 'resources/data_files'
output_location = 'resources/output'
output_type = 'ai'
periodicity_assumption = False
scale_location = 'resources/scales'
step_size = 2
valid_extensions = '.txt', '.dat'
window_averaging = False
window_size = 7

# Lyapunov parameters
compute_lyapunov = True
lyap_m = 2
lyap_tau = 1
lyap_kmax = 30
lyap_theiler = 0
lyap_eps = 0.0
lyap_fit_start = 1
lyap_fit_end = 10
lyap_knn = 15
lyap_auto_fit = True
lyap_r2min = 0.98
lyap_minlen = 6

# Utility functions
clear = lambda: progress_output.delete(1.0, END)
divider = lambda master, length: Frame(master, height=3, bd=1, width=length, relief='ridge')
textdivider = lambda length: '\n' + '─' * length + '\n'
files_per = lambda scales: {True: len(range(window_size, 1, -(step_size))) + 1, False: 1}[window_averaging] * len(scales)
setting_display = lambda: '\n'.join([var + ': ' + str(state) for var, state in (
    ('Periodicity', periodicity_assumption), ('Flip', flip_read), ('Complement', complement_strand),
    ('Lyapunov', compute_lyapunov), ('lyap_m', lyap_m), ('lyap_tau', lyap_tau), ('lyap_kmax', lyap_kmax),
    ('lyap_fit', f'{lyap_fit_start}-{lyap_fit_end}'), ('lyap_eps', lyap_eps), ('lyap_theiler', lyap_theiler),
    ('lyap_knn', lyap_knn), ('lyap_auto_fit', lyap_auto_fit), ('lyap_r2min', lyap_r2min))])


class Options():

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, {
                str: StringVar(),
                bool: BooleanVar(),
                int: IntVar(),
                float: DoubleVar()
            }[type(value)])
            getattr(self, key).set(value)

    def globalise(self):
        current_values = {key: self.__dict__[key].get() for key in self.__dict__}
        globals().update(current_values)


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

    def output(self, data_array, scale_data):
        self.progress_bar = DualProgressBar(len(data_array), files_per(scale_data), 160, ['File ' + self.file_name + ' & data piece ' + i.title for i in data_array])
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
        if width != 'avg' and compute_lyapunov:
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
        elif width == 'avg':
            header_lines = ['# lyapunov = unavailable for averaged data sets\n']
        else:
            header_lines = ['# lyapunov = disabled\n']
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
        progress_output.configure(state='normal')
        clear()
        try:
            progress_output.insert('end', 'Processing :' + self.headers[self.master_progress] + textdivider(160))
        except IndexError:
            progress_output.insert('end', 'Done!' + textdivider(160))
        progress_output.insert('end', 'Current settings: \n' + setting_display() + textdivider(160))
        for progress, maximum in (self.sub_progress, self.sub_max), (self.master_progress, self.master_max):
            completed = int((progress / maximum) * self.length)
            remaining = self.length - completed
            progress_output.insert('end', '\n' + str(progress).zfill(2) + '/' + str(maximum).zfill(2) + textdivider(160) + '▮' * completed + '-' * remaining + textdivider(160))
        progress_output.configure(state='disabled')
        root.update()


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


def purge_directory():
    confirm = messagebox.askquestion('Purge Directory', 'Are you sure you want to delete all files and folders present in the output directory?')
    if confirm == "yes":
        for i in os.listdir(output_location):
            rmtree(os.path.join(output_location, i))
        messagebox.showinfo('Cleared', 'Directory purged')
        scan_output(0, output_location)


def produce_output():
    options.globalise()

    scale_files = [File(os.path.join(scale_location, scale_file)) for scale_file in os.listdir(scale_location)]
    scale_data = merge_dictionaries([scale_file.parse_scales(scale_file.raw_lines) for scale_file in scale_files])
    if 'a' in output_type and not ('i' in output_type):
        scale_data = {'Averaged': average_dictionaries(tuple(scale_data.values()))}
    elif 'a' in output_type:
        scale_data['Averaged'] = average_dictionaries(tuple(scale_data.values()))

    data_files = [File(os.path.join(input_location, data_file)) for data_file in os.listdir(input_location)]
    if not all_output_override:
        [data_file.output(data_file.parse_lines(data_file.raw_lines), scale_data) for data_file in data_files]
    else:
        options.window_averaging.set(True)
        for periodicity_assumption in (True, False):
            options.periodicity_assumption.set(periodicity_assumption)
            for flip_read in (True, False):
                options.flip_read.set(flip_read)
                for complement_strand in (True, False):
                    options.complement_strand.set(complement_strand)
                    options.globalise()
                    [data_file.output(data_file.parse_lines(data_file.raw_lines), scale_data) for data_file in data_files]

    scan_output(0, output_location)


def scan_output(depth, current_folder):
    if depth == 3:
        return
    if depth == 0:
        file_tree.delete(*file_tree.get_children())
    for i in sorted(os.listdir(current_folder)):
        if i.startswith('.'):
            continue
        new_folder = os.path.join(current_folder, i)
        parent = current_folder
        if depth == 0:
            parent = ''
        try:
            file_tree.insert(parent, index='end', iid=new_folder, text=i[:20] + '..')
        except:
            pass
        scan_output(depth + 1, new_folder)


def plot_selected():
    global plot_calls
    selected_file = file_tree.focus()
    format_string = ["--", "-.", "-x", "-D", "-h", "-p", "-o"]
    try:
        if not selected_file.endswith('.txt'):
            raise FileNotFoundError
        data = np.loadtxt(selected_file)
        plt.plot(data[:, 0], data[:, 1], format_string[plot_calls % len('markers')], label=selected_file)
        plt.xlabel('Sequence Position')
        plt.ylabel('Propensity Scale')
        plt.show(block=False)
        plot_calls += 1
    except (FileNotFoundError, NameError) as current_error:
        if type(current_error) == FileNotFoundError:
            messagebox.showerror('Error', 'Ensure you selected a valid data file!')
        else:
            messagebox.showerror('Error', 'Plotting disabled due to no matplotlib install')


# Initialize options
variable_names = 'all_output_override', 'complement_strand', 'flip_read', 'input_location', 'output_location', 'output_type', 'periodicity_assumption', \
                 'scale_location', 'step_size', 'window_averaging', 'window_size', \
                 'compute_lyapunov', 'lyap_m', 'lyap_tau', 'lyap_kmax', 'lyap_theiler', 'lyap_eps', 'lyap_fit_start', 'lyap_fit_end', \
                 'lyap_knn', 'lyap_auto_fit', 'lyap_r2min', 'lyap_minlen'

options = Options(**{name: eval(name) for name in variable_names})
variable_zip = zip(variable_names, [i.replace('_', ' ').capitalize() for i in variable_names])

# Toolbar population
Label(option_toolbar, text='Generic options:').grid(column=0, columnspan=5, row=0, padx=5, pady=5, sticky=W)
divider(option_toolbar, 695).grid(column=0, columnspan=5, row=1, padx=5, pady=5)

checkbuttons = [Checkbutton(option_toolbar, text=variable_title, variable=getattr(options, variable_name)) for variable_name, variable_title in variable_zip if type(eval(variable_name)) == bool]
for i, button in enumerate(checkbuttons):
    button.grid(column=i, row=2, padx=5, pady=5)

for column, title in enumerate(('Step size:      ', 'Window size:')):
    Label(option_toolbar, text=title).grid(column=6 + column, row=0, padx=5, pady=5, sticky=W)
    divider(option_toolbar, 100).grid(column=6 + column, row=1, padx=5, pady=5)

optionmenus = [OptionMenu(option_toolbar, getattr(options, variable_name), *option_range) for variable_name, option_range in (('step_size', range(2, 24, 2)), ('window_size', range(3, 25, 2)))]
for i, menu in enumerate(optionmenus):
    menu.grid(column=6 + i, row=2, padx=5, pady=5, sticky=(E, W))

for column, title in enumerate(('Input file location:', 'Output file location:', 'Scale file location:')):
    Label(option_toolbar, text=title).grid(column=8 + (column * 2), columnspan=2, row=0, padx=5, pady=5, sticky=W)
    divider(option_toolbar, 215).grid(column=8 + (column * 2), columnspan=2, row=1, padx=5, pady=5)

for i, variable_name in enumerate([variable_name for variable_name in variable_names if type(eval(variable_name)) == str and os.path.exists(eval(variable_name))]):
    Entry(option_toolbar, textvariable=getattr(options, variable_name)).grid(column=8 + (i * 2), row=2, padx=5, pady=5)
    Button(option_toolbar, text='⤶', command=lambda variable_name=variable_name: getattr(options, variable_name).set(filedialog.askdirectory())).grid(column=9 + (i * 2), row=2, padx=5, pady=5)

# Lyapunov options UI
Label(option_toolbar, text='Lyapunov (Rosenstein) options:').grid(column=0, columnspan=8, row=3, padx=5, pady=5, sticky=W)
divider(option_toolbar, 695).grid(column=0, columnspan=8, row=4, padx=5, pady=5)

Label(option_toolbar, text='m').grid(column=0, row=5, padx=5, pady=5, sticky=W)
OptionMenu(option_toolbar, options.lyap_m, *range(1, 11)).grid(column=1, row=5, padx=5, pady=5, sticky=(E, W))

Label(option_toolbar, text='tau').grid(column=2, row=5, padx=5, pady=5, sticky=W)
OptionMenu(option_toolbar, options.lyap_tau, *range(1, 11)).grid(column=3, row=5, padx=5, pady=5, sticky=(E, W))

Label(option_toolbar, text='kmax').grid(column=4, row=5, padx=5, pady=5, sticky=W)
OptionMenu(option_toolbar, options.lyap_kmax, *range(5, 101, 5)).grid(column=5, row=5, padx=5, pady=5, sticky=(E, W))

Label(option_toolbar, text='fit start').grid(column=6, row=5, padx=5, pady=5, sticky=W)
OptionMenu(option_toolbar, options.lyap_fit_start, *range(0, 31)).grid(column=7, row=5, padx=5, pady=5, sticky=(E, W))

Label(option_toolbar, text='fit end').grid(column=8, row=5, padx=5, pady=5, sticky=W)
OptionMenu(option_toolbar, options.lyap_fit_end, *range(1, 51)).grid(column=9, row=5, padx=5, pady=5, sticky=(E, W))

Label(option_toolbar, text='eps').grid(column=10, row=5, padx=5, pady=5, sticky=W)
Entry(option_toolbar, textvariable=options.lyap_eps, width=8).grid(column=11, row=5, padx=5, pady=5, sticky=(E, W))

Label(option_toolbar, text='theiler (0=auto)').grid(column=12, row=5, padx=5, pady=5, sticky=W)
OptionMenu(option_toolbar, options.lyap_theiler, *([0] + list(range(1, 51)))).grid(column=13, row=5, padx=5, pady=5, sticky=(E, W))

Label(option_toolbar, text='K neighbors').grid(column=0, row=6, padx=5, pady=5, sticky=W)
OptionMenu(option_toolbar, options.lyap_knn, *range(1, 31)).grid(column=1, row=6, padx=5, pady=5, sticky=(E, W))

Checkbutton(option_toolbar, text='Auto fit region', variable=options.lyap_auto_fit).grid(column=2, columnspan=2, row=6, padx=5, pady=5, sticky=W)

Label(option_toolbar, text='R² min').grid(column=4, row=6, padx=5, pady=5, sticky=W)
Entry(option_toolbar, textvariable=options.lyap_r2min, width=6).grid(column=5, row=6, padx=5, pady=5, sticky=(E, W))

Label(option_toolbar, text='min fit len').grid(column=6, row=6, padx=5, pady=5, sticky=W)
OptionMenu(option_toolbar, options.lyap_minlen, *range(3, 21)).grid(column=7, row=6, padx=5, pady=5, sticky=(E, W))

# Option overrides
optionmenus[0].configure(state=DISABLED)
getattr(options, 'all_output_override').trace('w', lambda name, index, mode: [i.configure(state={True: DISABLED, False: NORMAL}[options.all_output_override.get()]) for i in checkbuttons[1:]])
getattr(options, 'all_output_override').trace('w', lambda name, index, mode: [i.select() if options.all_output_override.get() else i.deselect() for i in checkbuttons[1:]])
getattr(options, 'window_averaging').trace('w', lambda name, index, mode: optionmenus[0].configure(state={True: NORMAL, False: DISABLED}[options.window_averaging.get()]))

output_button = Button(option_toolbar, text='Output with these settings', command=produce_output, height=4).grid(column=14, row=0, rowspan=3, padx=5, pady=5)

# File browser buttons
refresh_directory = Button(file_browser, text='Search directory', width=22, command=lambda: scan_output(0, output_location))
refresh_directory.grid(column=0, row=0, padx=5, pady=5)

purge_directory = Button(file_browser, text='Purge directory', width=22, command=purge_directory)
purge_directory.grid(column=0, row=1, padx=5, pady=5)

plot_selected = Button(file_browser, text='Plot selected file', command=plot_selected, width=22)
plot_selected.grid(column=0, row=2, padx=5, pady=5)

divider(file_browser, 200).grid(column=0, row=3, padx=5, pady=5)

# Check for matplotlib
try:
    import matplotlib.pyplot as plt
    plt.ion()
    font = {'family': 'Helvetica', 'weight': 'regular', 'size': 18}
    plt.rc('font', **font)
except ImportError:
    messagebox.showerror('Error!', 'Matplotlib not found, plotting functionality disabled')

root.mainloop()
