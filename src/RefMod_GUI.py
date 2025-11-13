# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 11:20:47 2025

@author: alaguillog
"""

import FreeSimpleGUI as sg
from PIL import Image
from pathlib import Path
import configparser
import subprocess
import threading
import io
import os
import signal
import json
import sys

ROOT_DIR = Path(__file__).parent.parent
SETTINGS_FILE = ROOT_DIR/"config/RefMod_GUI_settings.json" # User preferences file
DOCS_PATH = "docs/docs.pdf" # Documentation

def load_settings():
    """Load saved user settings"""
    if os.path.exists(SETTINGS_FILE):
        try:
            with open(SETTINGS_FILE, "r") as f:
                return json.load(f)
        except Exception:
            pass
    return {}

def save_settings(values):
    """Save current user inputs to JSON"""
    try:
        data = {k: v for k, v in values.items() if isinstance(v, (str, int, bool))}
        with open(SETTINGS_FILE, "w") as f:
            json.dump(data, f, indent=2)
    except Exception:
        pass

def load_ini(file_path):
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(file_path)
    return config

def ini_to_dict(config):
    """Build dict from configparser for showing in the GUI"""
    data = {}
    for section in config.sections():
        for key, value in config.items(section):
            data[f"{section}.{key}"] = value
    return data

def dict_to_ini(data, file_path):
    """Build ini from dict for saving to file"""
    config = configparser.ConfigParser()
    for composite_key, value in data.items():
        section, key = composite_key.split(".", 1)
        if section not in config:
            config[section] = {}
        config[section][key] = value
    with open(file_path, "w") as f:
        config.write(f)
        
def custom_popup(title, message, image_file, window_size=(500, 450)):
    with Image.open(image_file) as img:
        img_ratio = img.width / img.height
        target_width = window_size[0] - 40  # leave a little padding
        target_height = int(target_width / img_ratio)
        img_resized = img.resize((target_width, target_height))
        bio = io.BytesIO()
        img_resized.save(bio, format="PNG")
        image_data = bio.getvalue()
        
    layout = [
        [sg.Column(
            [[sg.Image(data=image_data, key="-REFMOD_LOGO-")]],
            justification="center",
            element_justification="center",
            expand_x=True
        )],
        
        [sg.Multiline(
            message,
            size=(None, None),
            expand_x=True,
            expand_y=True,
            border_width=0,
            background_color=sg.theme_background_color(),
            text_color=sg.theme_text_color(),
            no_scrollbar=True,
            disabled=True,
            pad=(10, 10)
        )],
        
        [sg.Column(
            [[sg.Button("Close", size=(10, 1))]],
            justification="center",
            element_justification="center",
            expand_x=True
        )],
    ]

    window = sg.Window(
        title,
        layout,
        modal=True,
        size=window_size,
        element_justification="center",
        finalize=True,
        icon=ROOT_DIR/'assets/RefMod_icon.ico'
    )

    while True:
        event, _ = window.read()
        if event in (sg.WINDOW_CLOSED, "Close"):
            break
    window.close()

def run_script(values, window):
    """Build command and run RefMod"""
    cmd = [sys.executable, ROOT_DIR/"src/RefMod.py"]

    # Required arguments
    cmd += ["-i", values["-INFILE-"]]
    cmd += ["-r", values["-RAWFILE-"]]
    cmd += ["-d", values["-DMFILE-"]]

    # Optional arguments
    if values["-DIA-"]:
        cmd += ["-a", values["-DIA-"]]
    # if values["-SCANRANGE-"]:
    #     cmd += ["-s", values["-SCANRANGE-"]]
    scan_start = values.get("-SCAN_START-", "")
    scan_end = values.get("-SCAN_END-", "")
    scan_range = ""
    if scan_start and scan_end:
        scan_range = f"{scan_start},{scan_end}"
    elif scan_start:
        scan_range = scan_start
    elif scan_end:
        scan_range = scan_end
    if scan_range:
        cmd += ["-s", scan_range]
    if values["-OUTDIR-"]:
        cmd += ["-o", values["-OUTDIR-"]]
    if values["-CONFIG_TO_RUN-"]:
        cmd += ["-c", values["-CONFIG_TO_RUN-"]]
    if values["-WORKERS-"]:
        cmd += ["-w", str(values["-WORKERS-"])]
    if values["-VERBOSE-"]:
        cmd += ["-v"]

    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        universal_newlines=True,
        creationflags=subprocess.CREATE_NEW_PROCESS_GROUP
    )

    window.write_event_value('-PROCESS-', process)
    
    error_lines = []

    progress = 0
    progress_current = 1
    progress_total = 0
    if os.path.isfile(values["-INFILE-"]): progress_total = 1
    elif os.path.isdir(values["-INFILE-"]): progress_total = len([f for f in os.listdir(values["-INFILE-"]) if any(f.lower().endswith(s) for s in [".tsv"])]) # , ".txt"])])
    window["-PROGRESS_BAR-"].update(0, progress_total)
    
    for line in iter(process.stdout.readline, ''):
        if '%|' in line:
            # tqdm line
            window.write_event_value('-UPDATE-', line)
        else:
            if " - INFO - Reading MSFragger file (" in line: # update progress label
                window["-PROGRESS_LABEL-"].update("Searching file " + str(progress_current) + " out of " + str(progress_total) + " ...")
                progress += 1
                progress_current += 1
            elif " - INFO - Done." in line: # update progress bar
                window["-PROGRESS_BAR-"].update(progress)
            elif " - INFO - ERROR: " in line:
                error_lines.append(line)
            #elif "INFO - end script" in line:
            # normal line
            window.write_event_value('-APPEND-', line)

    process.wait()
    window.write_event_value('-ERRORS-', error_lines)
    window.write_event_value('-DONE-', process.returncode)


# Load saved settings
settings = load_settings()
sg.theme('TealMono')

# GUI layout
italic = (sg.DEFAULT_FONT[0], sg.DEFAULT_FONT[1], "italic")
bold = (sg.DEFAULT_FONT[0], sg.DEFAULT_FONT[1], "bold")
amino_acids = [
    ("Alanine", "A"),
    ("Arginine", "R"),
    ("Asparagine", "N"),
    ("Aspartic Acid", "D"),
    ("Cysteine", "C"),
    ("Glutamic Acid", "E"),
    ("Glutamine", "Q"),
    ("Glycine", "G"),
    ("Histidine", "H"),
    ("Isoleucine", "I"),
    ("Leucine", "L"),
    ("Lysine", "K"),
    ("Methionine", "M"),
    ("Phenylalanine", "F"),
    ("Proline", "P"),
    ("Serine", "S"),
    ("Threonine", "T"),
    ("Selenocysteine", "U"),
    ("Tryptophan", "W"),
    ("Tyrosine", "Y"),
    ("Valine", "V"),
    ("Pyrrolysine", "O"),
    ("Ambiguous E/Q", "Z")
] # TODO: Handle adding custom amino acids

default_masses = {
    "A": 71.037114, "R": 156.101111, "N": 114.042927, "D": 115.026943,
    "C": 103.009185, "E": 129.042593, "Q": 128.058578, "G": 57.021464,
    "H": 137.058912, "I": 113.084064, "L": 113.084064, "K": 128.094963,
    "M": 131.040485, "F": 147.068414, "P": 97.052764, "S": 87.032028,
    "T": 101.047679, "U": 150.953630, "W": 186.079313, "Y": 163.063329,
    "V": 99.068414, "O": 132.089878, "Z": 129.042594
}

default_mods = {
    "A": 0, "R": 0, "N": 0, "D": 0, "C": 57.021464, "E": 0, "Q": 0, "G": 0,
    "H": 0, "I": 0, "L": 0, "K": 0, "M": 0, "F": 0, "P": 0, "S": 0,
    "T": 0, "U": 0, "W": 0, "Y": 0, "V": 0, "O": 0, "Z": 0
}

aa_rows = []
for name, code in amino_acids:
    mass = default_masses.get(code.upper(), "")
    mods = default_mods.get(code.upper(), "")
    aa_rows.append([
        sg.Text(f"{name} ({code})", size=(25,1)),
        sg.Input(default_text=mass, key=f"-{code}_MASS-", size=(20,1), disabled=True, justification="right"),
        sg.Input(default_text=mods, key=f"-{code}_FM-", size=(20,1), justification="right")
    ])

iniedit_layout = [
    [sg.Text("Configuration File:"), sg.Input(settings.get("-CONFIG-", ROOT_DIR/"config/RefMod.ini"), key="-CONFIG-"),
     sg.FileBrowse(file_types = (('INI Files', '*.ini;*.INI'),), initial_folder = "."),
     sg.Button("Load Config"),
     sg.Button("Save Config", key="-SAVE_CONFIG-"),
     sg.Button("Reset To Default", key="-RESET_CONFIG-")],
    [sg.Text("Hover over the name of each parameter to show a brief description.", font=italic)],
    [sg.Column([
        [sg.Text("\nSEARCH PARAMETERS", font=bold)],
        [sg.Text("Batch Size", size=(25,1), tooltip=" Size (number of PSMs) of each task that will be submitted to a CPU core. "),
         sg.Input(default_text=settings.get("-BATCH_SIZE-", 1000), key="-BATCH_SIZE-", size=(10,1), enable_events=True)],
        [sg.Text("Fragment Tolerance", size=(25,1), tooltip=" Fragment mass tolerance, in parts-per-million. "),
         sg.Input(default_text=settings.get("-FRAGMENT_TOLERANCE-", 20.0), key="-FRAGMENT_TOLERANCE-", size=(10,1), enable_events=True),
         sg.Text("ppm")],
        [sg.Text("Theoretical Δmass Tolerance", size=(25,1), tooltip=" Tolerance for matching of theoretical and experimental Δmasses, in Dalton. \n This is an absolute value. "),
         sg.Input(default_text=settings.get("-DELTAMASS_TOLERANCE-", 3.0), key="-DELTAMASS_TOLERANCE-", size=(10,1), enable_events=True),
         sg.Text("Da")],
        [sg.Text("Score Mode", size=(25,1), tooltip=" The method for hyperscore calculation. "),
         sg.Radio("MOD-Hyperscore", "MODE_GROUP", key="-MODE_A-", default=True, tooltip=" Equivalent to MSFragger hyperscore. "),
         sg.Radio("HYB-Hyperscore", "MODE_GROUP", key="-MODE_B-", tooltip=" Attempts to match all non-modified fragment ions \n regardless of the position of the modification. ")],
        [sg.Text("Y-series Matching", size=(25,1), tooltip=" How to use the y-series for fragment matching. "),
         sg.Radio("Exclude y\u00b9", "Y_GROUP", key="-Y_A-", default=True, tooltip=" Exclude the y\u00b9 ion. Equivalent to MSFragger. "),
         sg.Radio("Use Full Series", "Y_GROUP", key="-Y_B-", tooltip=" Include the full y-series up to y\u207f. ")],
        [sg.Text("Δmass Preference", size=(25,1), tooltip=" Which Δmass candidate should be reported \n if both have the same score. "),
         sg.Radio("Experimental", "PREF_GROUP", key="-PREF_A-", default=True, tooltip=" Δmass determined by MSFragger. "),
         sg.Radio("Theoretical", "PREF_GROUP", key="-PREF_B-", tooltip=" Highest-scoring Δmass from the curated list. ")],
        
        [sg.Text("\nSPECTRUM PROCESSING", font=bold)],
        [sg.Text("Top N", size=(25,1), tooltip=" Maximum number of peaks (sorted by intensity) \n to keep from each spectrum. "),
         sg.Input(default_text=settings.get("-TOP_N-", 150), key="-TOP_N-", size=(10,1), enable_events=True)],
        [sg.Text("Minimum Intensity Ratio", size=(25,1), tooltip=" Remove peaks less intense than this multiple \n of the base peak intensity. "),
         sg.Input(default_text=settings.get("-MIN_RATIO-", 0.01), key="-MIN_RATIO-", size=(10,1), enable_events=True)],
        [sg.Text("Bin Top N", size=(25,1), tooltip=" Bin spectra according to the average aminoacid mass \n and keep the Top N peaks in each bin. "),
         sg.Checkbox("", default=settings.get("-BIN_TOP_N-", False), key="-BIN_TOP_N-")],
        [sg.Text("Minimum Fragment m/z", size=(25,1), tooltip=" Remove peaks with m/z lower than or equal to this value. "),
         sg.Input(default_text=settings.get("-MIN_FRAG_MZ-", 0), key="-MIN_FRAG_MZ-", size=(10,1), enable_events=True)],
        [sg.Text("Maximum Fragment m/z", size=(25,1), tooltip=" Remove peaks with m/z greater than or equal to this value. \n A value of 0 ignores this parameter. "),
         sg.Input(default_text=settings.get("-MAX_FRAG_MZ-", 0), key="-MAX_FRAG_MZ-", size=(10,1), enable_events=True)],
        [sg.Text("Deisotope", size=(25,1), tooltip=" Remove non-monoisotopic peaks up to 3 Carbon-13, \n with a tolerance of 0.005 Th. \n This is an experimental parameter. "),
         sg.Checkbox("", default=settings.get("-DEISO-", False), key="-DEISO-")],
        
        [sg.Text("\nFDR PARAMETERS", font=bold)],
        [sg.Text("Protein Column", size=(25,1), tooltip=" Name of the column containing protein names. "),
         sg.Input(default_text=settings.get("-PROTEIN-", "protein"), key="-PROTEIN-", size=(40,1))],
        [sg.Text("Decoy Prefix", size=(25,1), tooltip=" The prefix that marks decoy protein IDs. "),
         sg.Input(default_text=settings.get("-DECOY-", "DECOY"), key="-DECOY-", size=(40,1))],
        [sg.Text("Filter Targets", size=(25,1), tooltip=" Remove Decoys from output. "),
         sg.Checkbox("", default=settings.get("-FILTER_TARGET-", False), key="-FILTER_TARGET-")],
        [sg.Text("Filter FDR", size=(25,1), tooltip=" Remove PSMs above this FDR threshold. \n A value of 0 ignores this parameter. "),
         sg.Spin([round(x * 0.01, 2) for x in range(0, 101)], initial_value=settings.get("-FILTER_FDR-", 0), key="-FILTER_FDR-", size=(10,1))],
        #[sg.Column([], scrollable=True, vertical_scroll_only=True, size=(600,300), key="-INI_INPUTS-")],
        
        [sg.Text("\nAMINO ACIDS", font=bold)],
        [sg.Text("Adding custom amino acids through the GUI is currently unsupported. It can still be done by editing the INI file directly.", font=italic)],
        [sg.Text("", size=(25,1)), sg.Text("Amino Acid Mass", size=(18,1)), sg.Text("Fixed Modifications", size=(20,1))],
        *aa_rows,
        
        [sg.Text("\nOTHER MASSES", font=bold)],
        [sg.Text("Proton", size=(25,1)), sg.Input(default_text=1.007276, key="-PROTON_MASS-", size=(20,1), disabled=True, justification="right")],
        [sg.Text("Hydrogen", size=(25,1)), sg.Input(default_text=1.007825, key="-HYDROGEN_MASS-", size=(20,1), disabled=True, justification="right")],
        [sg.Text("Oxygen", size=(25,1)), sg.Input(default_text=15.994915, key="-OXYGEN_MASS-", size=(20,1), disabled=True, justification="right")],
        
        [sg.Text("\nLOGGING", font=bold)],
        [sg.Text("Create Log", size=(25,1), tooltip=" Create a log file. "),
         sg.Checkbox("", default=settings.get("-CREATE_LOG-", True), key="-CREATE_LOG-")],
        [sg.Text("Create INI", size=(25,1), tooltip=" Create a copy of the configuration file in the input directory. "),
         sg.Checkbox("", default=settings.get("-CREATE_INI-", True), key="-CREATE_INI-")],
        [sg.Text("Create Summary", size=(25,1), tooltip=" Create summary files containing metadata about each search. "),
         sg.Checkbox("", default=settings.get("-CREATE_SUMMARY-", True), key="-CREATE_SUMMARY-")],
        
        [sg.Text("\nDEBUG", font=bold)],
        [sg.Text("Debug Scores", size=(25,1)), sg.Checkbox("", default=settings.get("-DEBUG_SCORES-", False), key="-DEBUG_SCORES-")],
    ], scrollable=True, vertical_scroll_only=True, size=(760,640))]
]

run_layout = [
    [sg.Text("MSFragger File(s)", size=(25,1), justification='right'),
     sg.Input(settings.get("-INFILE-", ""), key="-INFILE-", size=(60,1), enable_events=True),
     sg.FileBrowse(button_text = "Load File", target="-INFILE-", file_types = (('Tab-separated Text Files', '*.tsv;*.txt'),), key="-BROWSE_FILE-"),
     sg.FolderBrowse(button_text = "Load Folder", target="-INFILE-", key="-BROWSE_FOLDER-")],
    [sg.Text("MS Data File(s)", size=(25,1), justification='right'),
     sg.Input(settings.get("-RAWFILE-", ""), key="-RAWFILE-", size=(60,1), enable_events=True),
     sg.FileBrowse(button_text = "Load File", file_types = (('MS Data Files', '*.mzML;*.MGF;*.mzml;*.mgf'),), key="-BROWSE_RAW-"),
     sg.FolderBrowse(button_text = "Load Folder", target="-RAWFILE-", key="-BROWSE_RAW_FOLDER-")],
    [sg.Text("Δmass File", size=(25,1), justification='right'),
     sg.Input(settings.get("-DMFILE-", ROOT_DIR/"data/dm_list.tsv"), key="-DMFILE-", size=(60,1), enable_events=True),
     sg.FileBrowse(button_text = "Load File", file_types = (('Tab-separated Text Files', '*.tsv;*.txt'),), key="-BROWSE_DM-")],
    # [sg.Text("_chN Files", size=(25,1), justification='right'),
    #  sg.Input(settings.get("-DIA-", ""), key="-DIA-", size=(60,1), enable_events=True)],
    [sg.Input("", key="-DIA-", visible=False)],
    [sg.Text("Scan Range", size=(25,1), justification='right'),
     sg.Spin([i for i in range(0, 1000000)], initial_value=int(settings.get("-SCAN_START-", 0)), key="-SCAN_START-", enable_events=True, size=(8,1)),
     sg.Text("-", pad=(0,0)),
     sg.Spin([i for i in range(0, 1000000)], initial_value=int(settings.get("-SCAN_END-", 0)), key="-SCAN_END-", enable_events=True, size=(8,1)),
     sg.Text("A value of 0 ignores these parameters.", font=italic)],
    [sg.Text("Output Directory", size=(25,1), justification='right'),
     sg.Input(settings.get("-OUTDIR-", ""), key="-OUTDIR-", size=(60,1), enable_events=True),
     sg.FolderBrowse(button_text = "Select Folder", key="-BROWSE_OUTPUT-")],
    [sg.Text("Configuration File", size=(25,1), justification='right'),
     sg.Input(settings.get("-CONFIG_TO_RUN-", ROOT_DIR/"config/RefMod.ini"), key="-CONFIG_TO_RUN-", size=(60,1), enable_events=True),
     sg.FileBrowse(button_text = "Load File", key="-BROWSE_CONFIG-")],
     #sg.Button("Load Config", key="-LOAD_CONFIG-")],
    [sg.Text("Number of workers", size=(25,1), justification='right'),
     sg.Spin([i for i in range(0, os.cpu_count()+1)], initial_value=os.cpu_count(), key="-WORKERS-", size=(8,1), enable_events=True)],
    [sg.Text("", size=(25,1)), sg.Checkbox("Verbose", default=settings.get("-VERBOSE-", False), key="-VERBOSE-")],
    [sg.Column([[sg.Multiline(size=(90, 25), key='-OUTPUT-', autoscroll=True, write_only=True, font=('Courier', 10))]], element_justification='center', expand_x=True)],
    [sg.Column([[sg.ProgressBar(100, orientation='h', size=(45, 20), bar_color=('green', 'white'), key='-PROGRESS_BAR-')]], pad=((30, 5), (10)), element_justification='left', expand_x=False),
     sg.Text("", justification="left", key="-PROGRESS_LABEL-")], # TODO get max value from input file list and update these values
    [sg.Column([[
        sg.Button("Run", bind_return_key=True),
        sg.Button("Stop", disabled=True, button_color=('white','red')),
        sg.Button("Exit")]], element_justification='center', expand_x=True)] # TODO add progress bar for 1 out of n files, etc
]
layout = [
    [sg.Text("RefMod v1.0"+" "*67, font=(sg.DEFAULT_FONT[0], sg.DEFAULT_FONT[1]*2, "bold")), sg.Button("About"), sg.Button("Help")], # TODO get version from script
    [sg.TabGroup([
        [sg.Tab('Parameters', iniedit_layout, key='-INI_TAB-'), sg.Tab('Run RefMod', run_layout, key='-RUN_TAB-')]
    ])]
]
window = sg.Window("RefMod GUI", layout, icon=ROOT_DIR/'assets/RefMod_icon.ico')

buffer = []
process = None

# Event loop
while True:
    event, values = window.read()
    if event in (sg.WINDOW_CLOSED, "Exit"):
        if process and process.poll() is None:
            try:
                process.terminate()
            except Exception:
                pass
        # Save settings before closing
        #save_settings(values)
        break

    elif event == "Run":
        if process and process.poll() is None:
            sg.popup_error("A process is already running.")
            continue

        # Save settings when running
        save_settings(values)

        window["-OUTPUT-"].update("")
        # window["-PROGRESS_BAR-"].update(0)
        buffer.clear()
        window["Run"].update(disabled=True)
        window["Stop"].update(disabled=False)
        window["Exit"].update(disabled=True)
        window["-PROGRESS_BAR-"].update(bar_color=("green", "white"))
        
        disable_keys = [
            "-INI_TAB-", "-BROWSE_FILE-", "-BROWSE_FOLDER-", "-BROWSE_RAW-", "-BROWSE_RAW_FOLDER-",
            "-BROWSE_DM-", "-BROWSE_OUTPUT-", "-BROWSE_CONFIG-",
            "-INFILE-", "-RAWFILE-", "-DMFILE-", "-DIA-", "-SCAN_START-",
            "-SCAN_END-", "-OUTDIR-", "-CONFIG_TO_RUN-", "-WORKERS-"
        ]
        for key in disable_keys:
            try:
                window[key].update(disabled=True, text_color="gray")
            except TypeError:
                window[key].update(disabled=True)

        threading.Thread(target=run_script, args=(values, window), daemon=True).start() # TODO: always save INI showing in GUI?
        
    elif event == "Load Config":
        ini_path = values["-CONFIG-"]
        if ini_path:
            config = configparser.ConfigParser(inline_comment_prefixes='#')
            config.read(ini_path)
            # SEARCH
            window["-BATCH_SIZE-"].update(int(config._sections['Search']['batch_size']))
            window["-FRAGMENT_TOLERANCE-"].update(float(config._sections['Search']['f_tol']))
            window["-DELTAMASS_TOLERANCE-"].update(float(config._sections['Search']['dm_tol']))
            radio = int(config._sections['Search']['score_mode'])
            window["-MODE_A-"].update(value=(radio == 0))
            window["-MODE_B-"].update(value=(radio == 1))
            radio = int(config._sections['Search']['full_y'])
            window["-Y_A-"].update(value=(radio == 1))
            window["-Y_B-"].update(value=(radio == 0))
            radio = int(config._sections['Search']['preference'])
            window["-PREF_A-"].update(value=(radio == 0))
            window["-PREF_B-"].update(value=(radio == 1))
            # SPECTRUM PROCESSING
            window["-TOP_N-"].update(int(config._sections['Spectrum Processing']['top_n']))
            window["-MIN_RATIO-"].update(float(config._sections['Spectrum Processing']['min_ratio']))
            window["-BIN_TOP_N-"].update(bool(int((config._sections['Spectrum Processing']['bin_top_n']))))
            window["-MIN_FRAG_MZ-"].update(float(config._sections['Spectrum Processing']['min_fragment_mz']))
            window["-MAX_FRAG_MZ-"].update(float(config._sections['Spectrum Processing']['max_fragment_mz']))
            window["-DEISO-"].update(bool(int((config._sections['Spectrum Processing']['deisotope']))))
            # FDR
            window["-PROTEIN-"].update(str(config._sections['FDR']['prot_column']))
            window["-DECOY-"].update(str(config._sections['FDR']['decoy_prefix']))
            window["-FILTER_TARGET-"].update(bool(int(config._sections['FDR']['filter_target'])))
            window["-FILTER_FDR-"].update(float(config._sections['FDR']['filter_fdr']))
            # AMINO ACIDS
            AAs = dict(config._sections['Aminoacids'])
            MODs = dict(config._sections['Fixed Modifications'])
            for name, code in amino_acids:
                if code.lower() in AAs:
                    window[f"-{code}_MASS-"].update(AAs[code.lower()])
                if code.lower() in MODs:
                    window[f"-{code}_FM-"].update(MODs[code.lower()])
            # MASSES
            window["-PROTON_MASS-"].update(str(config._sections['Masses']['m_proton']))
            window["-HYDROGEN_MASS-"].update(str(config._sections['Masses']['m_hydrogen']))
            window["-OXYGEN_MASS-"].update(str(config._sections['Masses']['m_oxygen']))
            # LOGGING
            window["-CREATE_LOG-"].update(bool(int((config._sections['Logging']['create_log']))))
            window["-CREATE_INI-"].update(bool(int((config._sections['Logging']['create_ini']))))
            window["-CREATE_SUMMARY-"].update(bool(int((config._sections['Logging']['create_summary']))))
            # DEBUG
            window["-DEBUG_SCORES-"].update(bool(int((config._sections['Debug']['debug_scores']))))
            
    elif event == "-SAVE_CONFIG-":
        save_path = sg.popup_get_file("Save INI file as...", save_as=True,
                                      file_types=(("INI Files", "*.ini;*.INI"),),
                                      default_extension=".ini",
                                      initial_folder=".") # TODO current path in -CONFIG-
        if save_path:
                config_data = {
                    "Search.batch_size": str(values["-BATCH_SIZE-"]),
                    "Search.f_tol": str(values["-FRAGMENT_TOLERANCE-"]),
                    "Search.dm_tol": str(values["-DELTAMASS_TOLERANCE-"]),
                    "Search.score_mode": "0" if values["-MODE_A-"] else "1",
                    "Search.full_y": "1" if values["-Y_A-"] else "0",
                    "Search.preference": "0" if values["-PREF_A-"] else "1",
                    "Spectrum Processing.top_n": str(values["-TOP_N-"]),
                    "Spectrum Processing.min_ratio": str(values["-MIN_RATIO-"]),
                    "Spectrum Processing.bin_top_n": "1" if values["-BIN_TOP_N-"] else "0",
                    "Spectrum Processing.min_fragment_mz": str(values["-MIN_FRAG_MZ-"]),
                    "Spectrum Processing.max_fragment_mz": str(values["-MAX_FRAG_MZ-"]),
                    "Spectrum Processing.deisotope": "1" if values["-DEISO-"] else "0",
                    "FDR.prot_column": values["-PROTEIN-"],
                    "FDR.decoy_prefix": values["-DECOY-"],
                    "FDR.filter_target": "1" if values["-FILTER_TARGET-"] else "0",
                    "FDR.filter_fdr": str(values["-FILTER_FDR-"]),
                    "Logging.create_log": "1" if values["-CREATE_LOG-"] else "0",
                    "Logging.create_ini": "1" if values["-CREATE_INI-"] else "0",
                    "Logging.create_summary": "1" if values["-CREATE_SUMMARY-"] else "0",
                    "Debug.debug_scores": "1" if values["-DEBUG_SCORES-"] else "0",
                    "Masses.m_proton": values["-PROTON_MASS-"],
                    "Masses.m_hydrogen": values["-HYDROGEN_MASS-"],
                    "Masses.m_oxygen": values["-OXYGEN_MASS-"]
                }
                for name, code in amino_acids:
                    config_data[f"Aminoacids.{code.upper()}"] = str(values[f"-{code}_MASS-"])
                    config_data[f"Fixed Modifications.{code.upper()}"] = str(values[f"-{code}_FM-"])
                
                dict_to_ini(config_data, save_path)
                #sg.popup("Configuration saved successfully!", title="Save Config")
                window["-CONFIG-"].update(save_path)
                window["-CONFIG_TO_RUN-"].update(save_path)
        
    elif event == "-RESET_CONFIG-":
        default_masses = {
            "A": 71.037114, "R": 156.101111, "N": 114.042927, "D": 115.026943,
            "C": 103.009185, "E": 129.042593, "Q": 128.058578, "G": 57.021464,
            "H": 137.058912, "I": 113.084064, "L": 113.084064, "K": 128.094963,
            "M": 131.040485, "F": 147.068414, "P": 97.052764, "S": 87.032028,
            "T": 101.047679, "U": 150.953630, "W": 186.079313, "Y": 163.063329,
            "V": 99.068414, "O": 132.089878, "Z": 129.042594
        }
        default_mods = {
            "A": 0, "R": 0, "N": 0, "D": 0, "C": 57.021464, "E": 0, "Q": 0, "G": 0,
            "H": 0, "I": 0, "L": 0, "K": 0, "M": 0, "F": 0, "P": 0, "S": 0,
            "T": 0, "U": 0, "W": 0, "Y": 0, "V": 0, "O": 0, "Z": 0
        }
        # window["-CONFIG-"].update(settings.get("-CONFIG-", ""))
        window["-BATCH_SIZE-"].update(1000)
        window["-FRAGMENT_TOLERANCE-"].update(20.0)
        window["-DELTAMASS_TOLERANCE-"].update(3.0)
        window["-MODE_A-"].update(True)
        window["-MODE_B-"].update(False)
        window["-Y_A-"].update(True)
        window["-Y_B-"].update(False)
        window["-PREF_A-"].update(True)
        window["-PREF_B-"].update(False)
        window["-TOP_N-"].update(150)
        window["-MIN_RATIO-"].update(0.01)
        window["-BIN_TOP_N-"].update(False)
        window["-MIN_FRAG_MZ-"].update(0)
        window["-MAX_FRAG_MZ-"].update(0)
        window["-DEISO-"].update(False)
        window["-PROTEIN-"].update("protein")
        window["-DECOY-"].update("DECOY")
        window["-FILTER_TARGET-"].update(False)
        window["-FILTER_FDR-"].update(0)
        window["-CREATE_LOG-"].update(True)
        window["-CREATE_INI-"].update(True)
        window["-CREATE_SUMMARY-"].update(True)
        window["-DEBUG_SCORES-"].update(False)
        for name, code in amino_acids:
                mass = default_masses.get(code.upper(), "")
                mods = default_mods.get(code.upper(), "")
                window[f"-{code}_MASS-"].update(mass)
                window[f"-{code}_FM-"].update(mods)
        window["-PROTON_MASS-"].update(1.007276)
        window["-HYDROGEN_MASS-"].update(1.007825)
        window["-OXYGEN_MASS-"].update(15.994915)
            
    elif event == "Help":
        try:
            if sys.platform.startswith('darwin'):  # macOS
                subprocess.call(('open', ROOT_DIR / DOCS_PATH))
            elif os.name == 'nt':  # Windows
                os.startfile(ROOT_DIR / DOCS_PATH)
            elif os.name == 'posix':  # Linux
                subprocess.call(('xdg-open', ROOT_DIR / DOCS_PATH))
        except Exception as e:
            sg.popup_error(f"Could not open help file:\n{e}")
            
    elif event == "About":
        about_text = "\nRefMod is an implementation of the ReCom concept (Laguillo-Gómez et al., 2023) designed to be run as a post-processing step after an open MSFragger search. It is compatible with both DDA and DIA data. When using RefMod, DIA data can be searched in a “pseudo-DDA” workflow, using a curated list of theoretical Δmass values to correct errors caused by the uncertainty in precursor masses contained within the same fragmentation window. \n\nRefMod has been developed at the Cardiovascular Proteomics Lab / Proteomics Unit at CNIC (Spanish National Centre for Cardiovascular Research). "
        version_text = "\n\nVersion v1.0"
        custom_popup(title="About RefMod", message=about_text+version_text, image_file=ROOT_DIR/"assets/RefMod_logo_text.png")

    elif event == "-PROCESS-":
        process = values[event]

    elif event == "-APPEND-":
        buffer.append(values[event])
        window["-OUTPUT-"].update(''.join(buffer))

    elif event == "-UPDATE-":
        if buffer:
            buffer[-1] = values[event]
        else:
            buffer.append(values[event])
        window["-OUTPUT-"].update(''.join(buffer))
    
    elif event == "-ERRORS-":
        error_lines = values[event]
        
    elif event in ["-BATCH_SIZE-", "-TOP_N-", "-SCAN_START-", "-SCAN_END-"]:
        # Keep only digits
        val = values[event]
        if not val.isdigit():
            new_val = "".join(c for c in val if c.isdigit())
            window[event].update(new_val)
            
    elif event in ["-FRAGMENT_TOLERANCE-", "-DELTAMASS_TOLERANCE-", "-MIN_RATIO-", "-FILTER_FDR-"]:
        # Keep only digits and a single decimal point
        val = values[event]
        new_val = ""
        decimal_found = False
        for c in val:
            if c.isdigit():
                new_val += c
            elif c == '.' and not decimal_found:
                new_val += c
                decimal_found = True
        if new_val != val:
            window[event].update(new_val)

    # elif event == "-PROGRESS-":
    #     window["-PROGRESS_BAR-"].update(values[event])
    #     window.refresh()

    elif event == "Stop":
        if process and process.poll() is None:
            try:
                if os.name == 'nt':
                    process.send_signal(signal.CTRL_BREAK_EVENT)
                else:
                    os.killpg(os.getpgid(process.pid), signal.SIGTERM)
            except Exception:
                try:
                    process.terminate()
                except Exception:
                    pass
            window["-OUTPUT-"].print("\nRefMod stopped by user.\n", text_color='red')
        window["Stop"].update(disabled=True)
        window["Run"].update(disabled=False)
        window["Exit"].update(disabled=False)
        enable_keys = [
            "-INI_TAB-", "-BROWSE_FILE-", "-BROWSE_FOLDER-", "-BROWSE_RAW-", "-BROWSE_RAW_FOLDER-",
            "-BROWSE_DM-", "-BROWSE_OUTPUT-", "-BROWSE_CONFIG-",
            "-INFILE-", "-RAWFILE-", "-DMFILE-", "-DIA-", "-SCAN_START-",
            "-SCAN_END-", "-OUTDIR-", "-CONFIG_TO_RUN-", "-WORKERS-"
        ]
        for key in enable_keys:
            try:
                window[key].update(disabled=False, text_color="black")
            except TypeError:
                window[key].update(disabled=False)
            
    elif event == "-DONE-":
        code = values[event]
        if error_lines:
            window["-PROGRESS_BAR-"].update(bar_color=("orange", "white"))
            window["-OUTPUT-"].print("\nWarnings detected during execution:\n", text_color='orange')
            for err in error_lines:
                    window["-OUTPUT-"].print(err, text_color='orange')
        if code == 0:
            window["-OUTPUT-"].print("\nRefMod finished successfully!\n", text_color='green')
        else:
            window["-OUTPUT-"].print("\nRefMod finished with an error or was stopped.\n", text_color='red')
        
        current_text = window["-PROGRESS_LABEL-"].get()
        updated_text = current_text.replace("Searching file ", "Searched ")
        updated_text = updated_text.replace("...", "files")
        window["-PROGRESS_LABEL-"].update(updated_text)

        window["Run"].update(disabled=False)
        window["Stop"].update(disabled=True)
        window["Exit"].update(disabled=False)
        enable_keys = [
            "-INI_TAB-", "-BROWSE_FILE-", "-BROWSE_FOLDER-", "-BROWSE_RAW-", "-BROWSE_RAW_FOLDER-",
            "-BROWSE_DM-", "-BROWSE_OUTPUT-", "-BROWSE_CONFIG-",
            "-INFILE-", "-RAWFILE-", "-DMFILE-", "-DIA-", "-SCAN_START-",
            "-SCAN_END-", "-OUTDIR-", "-CONFIG_TO_RUN-", "-WORKERS-"
        ]
        for key in enable_keys:
            try:
                window[key].update(disabled=False, text_color="black")
            except TypeError:
                window[key].update(disabled=False)
            
        process = None
        
window.close()
